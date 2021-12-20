# **************************************************************************
# *
# * Authors: Yunior C. Fonseca Reyna    (cfonseca@cnb.csic.es)
# *
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
import getpass
import os
import shutil
import time

from pkg_resources import parse_version

import pyworkflow.utils as pwutils
from pwem.constants import SCIPION_SYM_NAME
from pwem.constants import (SYM_CYCLIC, SYM_TETRAHEDRAL,
                            SYM_OCTAHEDRAL, SYM_I222, SYM_I222r)
from pwem.convert import Ccp4Header

from . import Plugin
from .constants import (CS_SYM_NAME, SYM_DIHEDRAL_Y, CRYOSPARC_USER,
                        CRYO_PROJECTS_DIR, V2_14_0, V2_13_0, CRYOSPARC_HOME,
                        CRYOSPARC_USE_SSD, V_UNKNOWN, V3_0_0, V2_15_0,
                        CRYOSPARC_MASTER)

VERSION = 'version'

STATUS_FAILED = "failed"
STATUS_ABORTED = "aborted"
STATUS_COMPLETED = "completed"
STATUS_KILLED = "killed"
STATUS_RUNNING = "running"
STATUS_QUEUED = "queued"
STATUS_LAUNCHED = "launched"
STATUS_STARTED = "started"
STATUS_BUILDING = "building"

STOP_STATUSES = [STATUS_ABORTED, STATUS_COMPLETED, STATUS_FAILED, STATUS_KILLED]
ACTIVE_STATUSES = [STATUS_QUEUED, STATUS_RUNNING, STATUS_STARTED,
                   STATUS_LAUNCHED, STATUS_BUILDING]

# Module variables
_csVersion = None  # Lazy variable: never use it directly. Use getCryosparcVersion instead


def getCryosparcDir(*paths):
    """
    Get the root directory where cryoSPARC code and dependencies are installed.
    """
    return Plugin.getHome(*paths)


def getCryosparcProgram(mode="cli"):
    """
    Get the cryosparc program to launch any command
    """
    csDir = getCryosparcDir()

    # TODO Find a better way to do that
    if csDir is not None:
        if os.path.exists(os.path.join(csDir, CRYOSPARC_MASTER, "bin")):
            # Case of CS v3.X.X is instaled
            return os.path.join(csDir, CRYOSPARC_MASTER, "bin",
                                'cryosparcm %s' % mode)
        else:
            # Case of CS v2.X.X is instaled
            return os.path.join(csDir, 'cryosparc2_master', "bin",
                                'cryosparcm %s' % mode)

    return None


def cryosparcExists():
    """
    Determine if scipion can find cryosparc
    :returns True if found, False otherwise
    """
    csDir = getCryosparcDir()
    return csDir is not None and os.path.exists(csDir)


def isCryosparcRunning():
    """
    Determine if cryosparc services are running
    :returns True if running, false otherwise
    """
    import subprocess
    status = -1
    if getCryosparcProgram() is not None:
        test_conection_cmd = (getCryosparcProgram() +
                              ' %stest_connection()%s ' % ("'", "'"))
        test_conection = subprocess.getstatusoutput(test_conection_cmd)
        status = test_conection[0]

    return status == 0


def cryosparcValidate():
    """
    Validates some cryo properties that must be satisfy
    """
    if not cryosparcExists():
        return ["cryoSPARC software not found at %s. Please, fill %s variable "
                "in scipion's config file." % (getCryosparcDir(),
                                               CRYOSPARC_HOME)]

    if not isCryosparcRunning():
        return ['Failed to connect to cryoSPARC. Please, make sure cryoSPARC '
                'is running.\nRunning: *%s* might fix this.'
                % getCryosparcProgram("start")]

    cryosparcVersion = parse_version(getCryosparcVersion())
    supportedVersions = Plugin.getSupportedVersions()
    minSupportedVersion = parse_version(supportedVersions[0])
    maxSupportedVersion = parse_version(supportedVersions[-1])

    # If version lower than first one
    if minSupportedVersion > cryosparcVersion:
        return ['The installed Cryosparc version is not '
                'compatible with the plugin. This can cause a '
                'malfunction of the protocol. Please install one of '
                'these versions: ' + str(supportedVersions).replace('\'', '')]

    elif maxSupportedVersion < cryosparcVersion:
        print(pwutils.yellowStr("cryoSPARC %s is newer than those we've tested %s. Instead of blocking the "
                                "execution, we are allowing this to run assuming compatibility is not broken."
                                "If it fails, please consider:\n A - upgrade the plugin, there might be an update.\n "
                                "B - downgrade cryosparc version.\n C - Contact plugin maintainers"
                                " at https://github.com/scipion-em/scipion-em-cryosparc2"
                                % (cryosparcVersion, str(supportedVersions).replace('\'', ''))))

    return []


def gpusValidate(gpuList, checkSingleGPU=False):
    """
    Validate a gpu list
    """
    # Case in which the protocol allow a single GPU
    if checkSingleGPU and len(gpuList) > 1:
        return ['This protocol can only be run on a single GPU.']
    return []


def getCryosparcEnvInformation(envVar=VERSION):
    """
    Get the cryosparc environment information
    """
    import ast
    system_info = getSystemInfo()
    dictionary = ast.literal_eval(system_info[1])
    envVariable = str(dictionary[envVar])
    return envVariable


def getCryosparcVersion():
    """ Gets cryosparc version 1st, from a variable if populated,
     2nd from the version txt file, if fails, asks CS using getCryosparcEnvInformation"""
    global _csVersion
    if _csVersion is None:
        try:
            _csVersion = _getCryosparcVersionFromFile()
        except Exception:
            try:
                _csVersion = getCryosparcEnvInformation(VERSION)
            except Exception:
                print("Couldn't get Cryosparc's version. Please review your config (%s)" % Plugin.getUrl())
                _csVersion = V_UNKNOWN
    return _csVersion.rstrip('\n')


def _getCryosparcVersionFromFile():
    versionFile = getCryosparcDir(CRYOSPARC_MASTER, "version")
    # read the version file
    with open(versionFile, "r") as fh:
        return fh.readline()


def getCryosparcUser():
    """
    Get the full name of the initial admin account
    """
    return os.path.basename(os.environ.get(CRYOSPARC_USER, ""))


def getCryosparcProjectsDir():
    """
    Get the path on the worker node to a writable directory
    """
    # Make a join in case is relative it will prepend getHome.
    cryoProject_Dir = os.path.join(Plugin.getHome(),
                                   Plugin.getVar(CRYO_PROJECTS_DIR))

    if not os.path.exists(cryoProject_Dir):
        os.mkdir(cryoProject_Dir)

    return cryoProject_Dir


def getProjectName(scipionProjectName):
    """ returns the name of the cryosparc project based on
    scipion project name and  a hash based on the user name"""

    username = getpass.getuser()
    return "%s-%s" % (scipionProjectName, username)


def getProjectPath(projectDir):
    """
    Gets all projects of given path .
    projectDir: Folder path to get sub folders.
    returns: Set with all sub folders.
    """
    folderPaths = os.listdir(projectDir)
    return folderPaths


def getJobLog(projectDirName, projectName, job):
    """
    Return the job log
    """
    return os.path.join(getCryosparcProjectsDir(), projectDirName, projectName,
                        job, 'job.log')


def createEmptyProject(projectDir, projectTitle):
    """
    create_empty_project(owner_user_id, project_container_dir, title=None,
                            desc=None)
    """

    create_empty_project_cmd = (getCryosparcProgram() +
                                ' %screate_empty_project("%s", "%s", "%s")%s '
                                % ("'", str(getCryosparcUser()),
                                   str(projectDir), str(projectTitle), "'"))

    return runCmd(create_empty_project_cmd, printCmd=False)


def createProjectDir(project_container_dir):
    """
    Given a "root" directory, create a project (PXXX) dir if it doesn't already
     exist
    :param project_container_dir: the "root" directory in which to create the
                                  project (PXXX) directory
    :returns: str - the final path of the new project dir with shell variables
              still in the returned path (the path should be expanded every
              time it is used)
    """
    create_project_dir_cmd = (getCryosparcProgram() +
                              ' %scheck_or_create_project_container_dir("%s")%s '
                              % ("'", project_container_dir, "'"))
    return runCmd(create_project_dir_cmd, printCmd=False)


def createEmptyWorkSpace(projectName, workspaceTitle, workspaceComment):
    """
    create_empty_workspace(project_uid, created_by_user_id,
                           created_by_job_uid=None,
                           title=None, desc=None)
    returns the new uid of the workspace that was created
    """
    create_work_space_cmd = (getCryosparcProgram() +
                             ' %screate_empty_workspace("%s", "%s", "%s", "%s", "%s")%s '
                             % ("'", projectName, str(getCryosparcUser()),
                                "None", str(workspaceTitle),
                                str(workspaceComment), "'"))
    return runCmd(create_work_space_cmd, printCmd=False)


def doImportParticlesStar(protocol):
    """
    do_import_particles_star(puid, wuid, uuid, abs_star_path,
                             abs_blob_path=None, psize_A=None)
    returns the new uid of the job that was created
    """
    print(pwutils.yellowStr("Importing particles..."), flush=True)
    className = "import_particles"
    params = {"particle_meta_path": str(os.path.join(os.getcwd(),
                                                     protocol._getFileName('input_particles'))),
              "particle_blob_path": str(os.path.join(os.getcwd(),
                                                     protocol._getTmpPath())),
              "psize_A": str(protocol._getInputParticles().getSamplingRate())
              }

    import_particles = enqueueJob(className, protocol.projectName, protocol.workSpaceName,
                                  str(params).replace('\'', '"'), '{}', protocol.lane)

    waitForCryosparc(protocol.projectName.get(), import_particles.get(),
                     "An error occurred importing particles. "
                     "Please, go to cryosPARC software for more "
                     "details.")

    return import_particles


def doImportVolumes(protocol, refVolumePath, refVolume, volType, msg):
    """
    :return:
    """
    print(pwutils.yellowStr(msg), flush=True)
    className = "import_volumes"
    params = {"volume_blob_path": str(refVolumePath),
              "volume_out_name": str(volType),
              "volume_psize": str(refVolume.getSamplingRate())}

    importedVolume = enqueueJob(className, protocol.projectName,
                                protocol.workSpaceName,
                                str(params).replace('\'', '"'), '{}',
                                protocol.lane)

    waitForCryosparc(protocol.projectName.get(), importedVolume.get(),
                     "An error occurred importing the volume. "
                     "Please, go to cryosPARC software for more "
                     "details."
                     )

    return importedVolume


def doJob(jobType, projectName, workSpaceName, params, input_group_connect):
    """
    do_job(job_type, puid='P1', wuid='W1', uuid='devuser', params={},
           input_group_connects={})
    """
    do_job_cmd = (getCryosparcProgram() +
                  ' %sdo_job("%s","%s","%s", "%s", %s, %s)%s' %
                  ("'", jobType, projectName, workSpaceName, getCryosparcUser(),
                   params, input_group_connect, "'"))

    return runCmd(do_job_cmd)


def enqueueJob(jobType, projectName, workSpaceName, params, input_group_connect,
               lane, gpusToUse=False, group_connect=None, result_connect=None):
    """
    make_job(job_type, project_uid, workspace_uid, user_id,
             created_by_job_uid=None, params={}, input_group_connects={})
    """
    from pyworkflow.object import String

    cryosparcVersion = getCryosparcVersion()
    # Create a compatible job to versions < v2.14.X
    make_job_cmd = (getCryosparcProgram() +
                    ' %smake_job("%s","%s","%s", "%s", "None", %s, %s)%s' %
                    ("'", jobType, projectName, workSpaceName, getCryosparcUser(),
                     params, input_group_connect, "'"))

    # Create a compatible job to versions >= v2.14.X
    if parse_version(cryosparcVersion) >= parse_version(V2_14_0):
        make_job_cmd = (getCryosparcProgram() +
                        ' %smake_job("%s","%s","%s", "%s", "None", "None", %s, %s)%s' %
                        ("'", jobType, projectName, workSpaceName,
                         getCryosparcUser(),
                         params, input_group_connect, "'"))

    # Create a compatible job to versions >= v3.0.X
    if parse_version(cryosparcVersion) >= parse_version(V3_0_0):
        make_job_cmd = (getCryosparcProgram() +
                        ' %smake_job("%s","%s","%s", "%s", "None", "None", %s, %s, "False", 0)%s' %
                        ("'", jobType, projectName, workSpaceName,
                         getCryosparcUser(),
                         params, input_group_connect, "'"))

    exitCode, cmdOutput = runCmd(make_job_cmd)

    # Extract the jobId
    jobId = String(cmdOutput.split()[-1])

    if group_connect is not None:
        for key, valuesList in group_connect.items():
            for value in valuesList:
                job_connect_group = (getCryosparcProgram() +
                                     ' %sjob_connect_group("%s", "%s", "%s")%s' %
                                     ("'", projectName, value, (str(jobId) + "." + key), "'"))
                runCmd(job_connect_group, printCmd=False)

    if result_connect is not None:
        for key, value in result_connect.items():
            job_connect_group = (getCryosparcProgram() +
                                 ' %sjob_connect_result("%s", "%s", "%s")%s' %
                                 ("'", projectName, value, (str(jobId) + "." + key), "'"))
            runCmd(job_connect_group, printCmd=True)

    print(pwutils.greenStr("Got %s for JobId" % jobId), flush=True)

    # Queue the job
    if parse_version(cryosparcVersion) < parse_version(V2_13_0):
        enqueue_job_cmd = (getCryosparcProgram() +
                           ' %senqueue_job("%s","%s","%s")%s' %
                           ("'", projectName, jobId,
                            lane, "'"))

    elif parse_version(cryosparcVersion) <= parse_version(V2_15_0):
        hostname = getCryosparcEnvInformation('master_hostname')
        if gpusToUse:
            gpusToUse = str(gpusToUse)
        enqueue_job_cmd = (getCryosparcProgram() +
                           ' %senqueue_job("%s","%s","%s", "%s", %s)%s' %
                           ("'", projectName, jobId,
                            lane, hostname, gpusToUse, "'"))

    elif parse_version(cryosparcVersion) >= parse_version(V3_0_0):
        hostname = getCryosparcEnvInformation('master_hostname')
        if gpusToUse:
            gpusToUse = str(gpusToUse)
        no_check_inputs_ready = False
        enqueue_job_cmd = (getCryosparcProgram() +
                           ' %senqueue_job("%s","%s","%s", "%s", %s, "%s")%s' %
                           ("'", projectName, jobId,
                            lane, hostname, gpusToUse,
                            no_check_inputs_ready, "'"))

    runCmd(enqueue_job_cmd)

    return jobId


def runCmd(cmd, printCmd=True):
    """ Runs a command and check its exit code. If different than 0 it raises an exception
    :parameter cmd command to run
    :parameter printCmd (default True) prints the command"""
    import subprocess
    if printCmd:
        print(pwutils.greenStr("Running: %s" % cmd), flush=True)
    exitCode, cmdOutput = subprocess.getstatusoutput(cmd)

    if exitCode != 0:
        raise Exception("%s failed --> Exit code %s, message %s" % (cmd, exitCode, cmdOutput))

    return exitCode, cmdOutput


def waitForCryosparc(projectName, jobId, failureMessage):
    """ Waits for cryosparc to finish or fail a job
    :parameter projectName: Cryosparc project name
    :parameter jobId: cryosparc job id
    :parameter failureMessage: Message for the exception thrown in case job fails
    :returns job Status
    :raises Exception when parsing cryosparc's output looks wrong"""

    # While is needed here, cause waitJob has a timeout of 5 secs.
    while True:
        status = getJobStatus(projectName, jobId)
        if status not in STOP_STATUSES:
            waitJob(projectName, jobId)
        else:
            break

    if status != STATUS_COMPLETED:
        raise Exception(failureMessage)

    return status


def getJobStatus(projectName, job):
    """
    Return the job status
    """
    get_job_status_cmd = (getCryosparcProgram() +
                          ' %sget_job_status("%s", "%s")%s'
                          % ("'", projectName, job, "'"))

    status = runCmd(get_job_status_cmd, printCmd=False)
    return status[-1]


def waitJob(projectName, job):
    """
    Wait while the job not finished
    """
    wait_job_cmd = (getCryosparcProgram() +
                    ' %swait_job_complete("%s", "%s")%s'
                    % ("'", projectName, job, "'"))
    runCmd(wait_job_cmd, printCmd=False)


def get_job_streamlog(projectName, job, fileName):
    get_job_streamlog_cmd = (getCryosparcProgram() +
                             ' %sget_job_streamlog("%s", "%s")%s%s'
                             % ("'", projectName, job, "'", ">" + fileName))

    runCmd(get_job_streamlog_cmd, printCmd=False)


def killJob(projectName, job):
    """
     Kill a Job (if running)
    :param projectName: the uid of the project that contains the job to kill
    :param job: the uid of the job to kill
    """
    kill_job_cmd = (getCryosparcProgram() +
                    ' %skill_job("%s", "%s")%s'
                    % ("'", projectName, job, "'"))
    runCmd(kill_job_cmd, printCmd=True)


def clearJob(projectName, job):
    """
         Clear a Job (if queued) to get it back to building state (do not clear
         params or inputs)
        :param projectName: the uid of the project that contains the job to clear
        :param job: the uid of the job to clear
        ** IMPORTANT: This method can be launch only if the job is queued
        """
    clear_job_cmd = (getCryosparcProgram() +
                     ' %sclear_job("%s", "%s")%s'
                     % ("'", projectName, job, "'"))
    runCmd(clear_job_cmd, printCmd=False)


def clearIntermediateResults(projectName, job, wait=3):
    """
     Clear the intermediate result from a specific Job
    :param projectName: the uid of the project that contains the job to clear
    :param job: the uid of the job to clear
    """
    print(pwutils.yellowStr("Removing intermediate results..."), flush=True)
    clear_int_results_cmd = (getCryosparcProgram() +
                             ' %sclear_intermediate_results("%s", "%s")%s'
                             % ("'", projectName, job, "'"))
    runCmd(clear_int_results_cmd, printCmd=False)
    # wait a delay in order to delete intermediate results correctly
    time.sleep(wait)


def getSystemInfo():
    """
    Returns system-related information related to the cryosparc app
    :returns: dict -- dictionary listing information about cryosparc environment
    {
        'master_hostname' : master_hostname,
        'port_webapp' : os.environ['CRYOSPARC_HTTP_PORT'],
        'port_mongo' : os.environ['CRYOSPARC_MONGO_PORT'],
        'port_command_core' : os.environ['CRYOSPARC_COMMAND_CORE_PORT'],
        'port_command_vis' : os.environ['CRYOSPARC_COMMAND_VIS_PORT'],
        'port_command_proxy' : os.environ['CRYOSPARC_COMMAND_PROXY_PORT'],
        'port_command_rtp' : os.environ['CRYOSPARC_COMMAND_RTP_PORT'],
        'port_rtp_webapp' : os.environ['CRYOSPARC_HTTP_RTP_PORT'],
        'version' : get_running_version(),
    }
    """
    system_info_cmd = (getCryosparcProgram() + " 'get_system_info()'")
    return runCmd(system_info_cmd, printCmd=False)


def addComputeSectionParams(form, allowMultipleGPUs=True):
    """
    Add the compute settings section
    """
    from pyworkflow.protocol.params import (BooleanParam, StringParam, NonEmpty,
                                            GPU_LIST)
    computeSSD = os.getenv(CRYOSPARC_USE_SSD)
    if computeSSD is None:
        computeSSD = False

    form.addParam('compute_use_ssd', BooleanParam, default=computeSSD,
                  label='Cache particle images on SSD',
                  help='Whether or not to copy particle images to the local '
                       'SSD before running. The cache is persistent, so after '
                       'caching once, particles should be available for '
                       'subsequent jobs that require the same data. Not '
                       'using an SSD can dramatically slow down processing.')

    # This is here because getCryosparcEnvInformation is failing in some machines
    try:
        versionAllowGPUs = parse_version(getCryosparcVersion()) >= parse_version(V2_13_0)
    # Code is failing to get CS info, either stop or some error
    except Exception:
        # ... we assume is a modern version
        versionAllowGPUs = True

    if versionAllowGPUs:
        if allowMultipleGPUs:
            form.addHidden(GPU_LIST, StringParam, default='0',
                           label='Choose GPU IDs:', validators=[NonEmpty],
                           help='This argument is necessary. By default, the '
                                'protocol will attempt to launch on GPU 0. You can '
                                'override the default allocation by providing a '
                                'list of which GPUs (0,1,2,3, etc) to use. '
                                'GPU are separated by ",". For example: "0,1,5"')
        else:
            form.addHidden(GPU_LIST, StringParam, default='0',
                           label='Choose GPU ID:', validators=[NonEmpty],
                           help='This argument is necessary. By default, the '
                                'protocol will attempt to launch on GPU 0. You can '
                                'override the default allocation by providing a '
                                'single GPU (0, 1, 2 or 3, etc) to use.')

    form.addParam('compute_lane', StringParam, default='default',
                  label='Lane name:',
                  help='The scheduler lane name to add the protocol execution')


def addSymmetryParam(form, help=""):
    """
    Add the symmetry param with the conventions
    :param form:
    :return:
    """
    from pyworkflow.protocol.params import (EnumParam, IntParam, Positive)
    form.addParam('symmetryGroup', EnumParam,
                  choices=[CS_SYM_NAME[SYM_CYCLIC] +
                           " (" + SCIPION_SYM_NAME[SYM_CYCLIC] + ")",
                           CS_SYM_NAME[SYM_DIHEDRAL_Y] +
                           " (" + SCIPION_SYM_NAME[SYM_DIHEDRAL_Y] + ")",
                           CS_SYM_NAME[SYM_TETRAHEDRAL] +
                           " (" + SCIPION_SYM_NAME[SYM_TETRAHEDRAL] + ")",
                           CS_SYM_NAME[SYM_OCTAHEDRAL] +
                           " (" + SCIPION_SYM_NAME[SYM_OCTAHEDRAL] + ")",
                           CS_SYM_NAME[SYM_I222] +
                           " (" + SCIPION_SYM_NAME[SYM_I222] + ")",
                           CS_SYM_NAME[SYM_I222r] +
                           " (" + SCIPION_SYM_NAME[SYM_I222r] + ")"],
                  default=SYM_CYCLIC,
                  label="Symmetry",
                  help="Symmetry as defined by cryosparc. Please note that "
                       "Dihedral symmetry in cryosparc is defined with respect"
                       "to y axis (Dyn).\n"
                       "If no symmetry is present, use C1.\n" +
                       help
                  )
    form.addParam('symmetryOrder', IntParam, default=1,
                  condition='symmetryGroup==%d or symmetryGroup==%d' %
                            (SYM_DIHEDRAL_Y - 1, SYM_CYCLIC),
                  label='Symmetry Order',
                  validators=[Positive],
                  help='Order of symmetry.')


def getSymmetry(symmetryGroup, symmetryOrder):
    """
    Get the symmetry(string) taking into account the symmetry convention
    """
    symmetry = {
        0: CS_SYM_NAME[SYM_CYCLIC][0] + str(symmetryOrder),  # Cn
        1: CS_SYM_NAME[SYM_DIHEDRAL_Y][0] + str(symmetryOrder),  # Dn
        2: CS_SYM_NAME[SYM_TETRAHEDRAL],  # T
        3: CS_SYM_NAME[SYM_OCTAHEDRAL],  # O
        4: CS_SYM_NAME[SYM_I222],  # I1
        5: CS_SYM_NAME[SYM_I222r]  # I2
    }
    return symmetry.get(symmetryGroup, "C1")


def calculateNewSamplingRate(newDims, previousSR, previousDims):
    """
    :param newDims:
    :param previousSR:
    :param previousDims:
    :return:
    """
    pX = previousDims[0]
    nX = newDims[0]
    return previousSR * pX / nX


def fixVolume(paths):
    """

    :param paths: accept a string or a list of strings
    :return:
    """
    if isinstance(paths, str):
        paths = [paths]
    for path in paths:
        ccp4header = Ccp4Header(path, readHeader=True)
        ccp4header.setISPG(1)
        ccp4header.writeHeader()


def copyFiles(src, dst, files=None):
    """
    Copy a list of files from src to dst. If files is None, all files of src are
    copied to dst
    :param src: source folder path
    :param dst: destiny folder path
    :param files: a list of files to be copied
    :return:
    """
    try:
        if files is None:
            shutil.copytree(src, dst)
        else:
            if isinstance(files, str):
                files = [files]
            for file in files:
                shutil.copy(os.path.join(src, file),
                            os.path.join(dst, file))
    except Exception as ex:
        print("Unable to execute the copy: Files or directory does not exist: ",
              ex)