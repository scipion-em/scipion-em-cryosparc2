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
import ast
import getpass
import logging
import os
import re
import shutil
import time

from pkg_resources import parse_version

import pyworkflow.utils as pwutils
from pwem.constants import SCIPION_SYM_NAME
from pwem.convert import Ccp4Header
from pyworkflow.protocol import IntParam

from . import Plugin
from .constants import *

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

# logging variable
logger = logging.getLogger(__name__)


class NestedDict:
    def __init__(self, depth=1):
        self.data = {}
        self.depth = depth

    def insert(self, keys, value):
        current = self.data
        for key in keys[:self.depth - 1]:
            current = current.setdefault(key, {})
        current[keys[self.depth - 1]] = value

    def search(self, keys):
        current = self.data
        for key in keys[:self.depth]:
            if key in current:
                current = current[key]
            else:
                return None
        return current


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
        logger.info(pwutils.yellowStr("cryoSPARC %s is newer than those we've tested %s. Instead of blocking the "
                                "execution, we are allowing this to run assuming compatibility is not broken."
                                "If it fails, please consider:\n A - upgrade the plugin, there might be an update.\n "
                                "B - downgrade cryosparc version.\n C - Contact plugin maintainers"
                                " at https://github.com/scipion-em/scipion-em-cryosparc2"
                                % (cryosparcVersion, str(supportedVersions).replace('\'', ''))))

    if cryosparcVersion >= parse_version(V4_1_0):
        if not os.environ.get(CRYOSPARC_USER):
            return ["You need to define the cryoSPARC user variable "
                    "(CRYOSPARC_USER) in the Scipion config file. Note that the "
                    "cryoSPARC username is the email address."]
        elif not userExist(os.environ.get(CRYOSPARC_USER)):
            return ["The user defined in the Scipion config file does not exist within CS."]

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
            _csVersion = _getCryosparcVersionFromFile().split('+')[0]
        except Exception:
            try:
                _csVersion = getCryosparcEnvInformation(VERSION).split('+')[0]
            except Exception as e:
                logger.error("Couldn't get Cryosparc's version. Please review your config (%s)" % Plugin.getUrl(), exc_info=e)
                _csVersion = V_UNKNOWN
    return _csVersion.rstrip('\n')


def _getCryosparcVersionFromFile():
    versionFile = getCryosparcDir(CRYOSPARC_MASTER, CRYOSPARC_VERSION_FILE)
    # read the version file
    with open(versionFile, "r") as fh:
        return fh.readline()


def _getLicenceFromFile():
    configFile = getCryosparcDir(CRYOSPARC_MASTER, CRYOSPARC_CONFIG_FILE)
    with open(configFile, 'r') as f:
        configContent = f.read().strip().split("\n")
        for variable in configContent:
            if CRYOSPARC_LICENSE_ID_VARIABLE in variable:
                return variable.split("=")[1].replace("\"", "")
        return None


def getCryosparcUser(userId=True):
    """
    Get the user
    """
    user = os.environ.get(CRYOSPARC_USER, "admin")
    cryosparcVersion = getCryosparcVersion()
    if parse_version(cryosparcVersion) >= parse_version(V4_1_0):
        if userId:
            user = getUserId(user)

    return user


def getCryosparcProjectsList():
    """
    Get list of all projects available
    :return: all projects available in the database
    """
    projects_list_cmd = (getCryosparcProgram() + ' %slist_projects()%s ' % ("'", "'"))
    cmd = runCmd(projects_list_cmd, printCmd=False)[1]
    projectList = ast.literal_eval(cmd)
    return projectList


def getCryosparcWorkSpaces(projectId):
    """ List all workspaces inside a given project (or all projects if not
    specified)

    :param projectId: target project UID, e.g. "P1", defaults to None
    :type projectId: str, optional
    :return: list of workpaces in a project or all projects if not specified
    :rtype: list
    """
    workspace_list_cmd = (getCryosparcProgram() + ' %slist_workspaces("%s")%s ' % ("'", str(projectId), "'"))
    cmd = runCmd(workspace_list_cmd, printCmd=False)[1]
    workspacesList = ast.literal_eval(cmd)
    return workspacesList


def isCryosparcStandalone():
    """
    Get the cryoSPARC installation mode. If True, we have a standalone installation
    else a cluster installation is considered. If the environment variable
    CRYOSPARC_STANDALONE_INSTALLATION isn't present, then we assume that we have
    a standalone installation and then, this method returns True
    """
    return os.environ.get(CRYOSPARC_STANDALONE_INSTALLATION, 'True') == 'True'


def getCryosparcDefaultLane():
    """
    Get the cryoSPARC default lane
    """
    return os.environ.get(CRYOSPARC_DEFAULT_LANE, None)


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


def getCryosparcProjectId(projectDir):
    """
    Get the project Id form project.json file.
    :param projectDir: project directory path
    """
    import json
    projectJsonFilePath = os.path.join(projectDir.get(), 'project.json')

    with open(projectJsonFilePath, 'r') as file:
        prjson = json.load(file)

    pId = prjson['uid']
    return pId


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


def getProjectInformation(project_uid, info='project_dir'):
    """
    Get information about a single project
    :param project_uid: the id of the project
    :return: the information related to the project that's stored in the database
    """
    import ast
    from cryosparc.tools import CryoSPARC
    getProject_cmd = (getCryosparcProgram() +
                                ' %sget_project("%s")%s '
                                % ("'", str(project_uid), "'"))

    project_info = runCmd(getProject_cmd, printCmd=False)
    dictionary = ast.literal_eval(project_info[1])
    return str(dictionary[info])


def getUserToken(email):
    get_user_cmd = (getCryosparcProgram() +
                                ' %sGetUser("%s")%s '
                                % ("'", str(email),"'"))

    return runCmd(get_user_cmd, printCmd=False)


def updateProjectDirectory(project_uid, new_project_dir):
    """
       Safely updates the project directory of a project given a directory. Checks
       if the directory exists, is readable, and writeable.
       :param project_uid: uid of the project to update
       :param new_project_dir_container: the new directory
       """
    updateProjectDirectory_cmd = (getCryosparcProgram() +
                      ' %supdate_project_directory("%s", "%s")%s '
                      % ("'", str(project_uid), str(new_project_dir), "'"))

    runCmd(updateProjectDirectory_cmd, printCmd=False)


def getOutputPreffix(projectName):
    cryosparcVersion = getCryosparcVersion()
    preffix = "cryosparc_" + projectName+"_" if parse_version(cryosparcVersion) < parse_version(V4_0_0) else ""
    return preffix


def createProjectContainerDir(project_container_dir):
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
                             % ("'", projectName, str(getCryosparcUser(userId=False)),
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
                                                     protocol._getPath())),
              "psize_A": str(protocol._getInputParticles().getSamplingRate())
              }

    import_particles = enqueueJob(className, protocol.projectName, protocol.workSpaceName,
                                  str(params).replace('\'', '"'), '{}', protocol.lane)

    waitForCryosparc(protocol.projectName.get(), import_particles.get(),
                     "An error occurred importing particles. "
                     "Please, go to cryoSPARC software for more "
                     "details.")

    return import_particles


def doImportVolumes(protocol, refVolumePath, refVolume, volType, msg):
    """
    :return:
    """
    logger.info(pwutils.yellowStr(msg))
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
                     "Please, go to cryoSPARC software for more "
                     "details."
                     )

    return importedVolume


def doImportMicrographs(protocol):
    print(pwutils.yellowStr("Importing micrographs..."), flush=True)
    className = "import_micrographs"
    micrographs = protocol._getInputMicrographs()
    acquisition = micrographs.getAcquisition()
    micList = list(micrographs.getFiles())

    micFolder = os.path.join(protocol._getExtraPath('micrographs'))
    os.makedirs(micFolder, exist_ok=True)

    for micPath in micList:
        micName = os.path.basename(micPath)
        micLink = os.path.join(micFolder, micName)
        os.symlink(os.path.abspath(micPath), micLink)
    micExt = '*%s' % os.path.splitext(micList[0])[1]

    params = {"blob_paths": str(os.path.join(os.getcwd(), micFolder, micExt)),
              "psize_A": str(micrographs.getSamplingRate()),
              "accel_kv": str(acquisition.getVoltage()),
              "cs_mm": str(acquisition.getSphericalAberration()),
              "total_dose_e_per_A2": str(0.1),
              "output_constant_ctf": "True"
              }

    import_particles = enqueueJob(className, protocol.projectName, protocol.workSpaceName,
                                  str(params).replace('\'', '"'), '{}', protocol.lane)

    waitForCryosparc(protocol.projectName.get(), import_particles.get(),
                     "An error occurred importing particles. "
                     "Please, go to cryoSPARC software for more "
                     "details.")

    return import_particles


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
    standaloneInstallation = isCryosparcStandalone()

    # Create a compatible job to versions < v2.14.X                DEPRECATED
    # make_job_cmd = (getCryosparcProgram() +
    #                 ' %smake_job("%s","%s","%s", "%s", "None", %s, %s)%s' %
    #                 ("'", jobType, projectName, workSpaceName, getCryosparcUser(),
    #                  params, input_group_connect, "'"))

    # Create a compatible job to versions >= v2.14.X               DEPRECATED
    # if parse_version(cryosparcVersion) >= parse_version(V2_14_0):
    #     make_job_cmd = (getCryosparcProgram() +
    #                     ' %smake_job("%s","%s","%s", "%s", "None", "None", %s, %s)%s' %
    #                     ("'", jobType, projectName, workSpaceName,
    #                      getCryosparcUser(),
    #                      params, input_group_connect, "'"))

    # Create a compatible job to versions >= v3.0.X < v4_3_1
    if parse_version(V3_0_0) <= parse_version(cryosparcVersion) < parse_version(V4_3_1):
        make_job_cmd = (getCryosparcProgram() +
                        ' %smake_job("%s","%s","%s", "%s", "None", "None", %s, %s, "False", 0)%s' %
                        ("'", jobType, projectName, workSpaceName,
                         getCryosparcUser(),
                         params, input_group_connect, "'"))

    # Create a compatible job to versions >= v4_3_1
    elif parse_version(cryosparcVersion) >= parse_version(V4_3_1):
        make_job_cmd = (getCryosparcProgram() +
                        ' %smake_job("%s","%s","%s", "%s", "None", "None", "None", %s, %s, "False", 0)%s' %
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

    logger.info(pwutils.greenStr("Got %s for JobId" % jobId))

    # Queue the job  DEPRECATED
    # if parse_version(cryosparcVersion) < parse_version(V2_13_0):
    #     enqueue_job_cmd = (getCryosparcProgram() +
    #                        ' %senqueue_job("%s","%s","%s")%s' %
    #                        ("'", projectName, jobId,
    #                         lane, "'"))
    #
    # elif parse_version(cryosparcVersion) <= parse_version(V2_15_0):
    #     if standaloneInstallation:
    #         hostname = getCryosparcEnvInformation('master_hostname')
    #         if gpusToUse:
    #             gpusToUse = str(gpusToUse)
    #         enqueue_job_cmd = (getCryosparcProgram() +
    #                            ' %senqueue_job("%s","%s","%s", "%s", %s)%s' %
    #                            ("'", projectName, jobId,
    #                             lane, hostname, gpusToUse, "'"))
    #     else:
    #         enqueue_job_cmd = (getCryosparcProgram() +
    #                            ' %senqueue_job("%s","%s","%s")%s' %
    #                            ("'", projectName, jobId, lane, "'"))

    if parse_version(cryosparcVersion) <= parse_version(V3_3_2):
        if standaloneInstallation:
            hostname = getCryosparcEnvInformation('master_hostname')
            if gpusToUse:
                gpusToUse = str(gpusToUse)
            no_check_inputs_ready = False
            enqueue_job_cmd = (getCryosparcProgram() +
                               ' %senqueue_job("%s","%s","%s", "%s", %s, "%s")%s' %
                               ("'", projectName, jobId,
                                lane, hostname, gpusToUse,
                                no_check_inputs_ready, "'"))
        else:
            enqueue_job_cmd = (getCryosparcProgram() +
                               ' %senqueue_job("%s","%s","%s")%s' %
                               ("'", projectName, jobId,
                                lane, "'"))
    elif parse_version(cryosparcVersion) >= parse_version(V4_0_0):
        user = getCryosparcUser()
        if standaloneInstallation:
            hostname = getCryosparcEnvInformation('master_hostname')
            if gpusToUse:
                gpusToUse = str(gpusToUse)
            no_check_inputs_ready = False
            enqueue_job_cmd = (getCryosparcProgram() +
                               ' %senqueue_job("%s","%s","%s", "%s", "%s", %s, "%s")%s' %
                               ("'", projectName, jobId,
                                lane, user, hostname, gpusToUse,
                                no_check_inputs_ready, "'"))
        else:
            enqueue_job_cmd = (getCryosparcProgram() +
                               ' %senqueue_job("%s","%s","%s","%s")%s' %
                               ("'", projectName, jobId,
                                lane, user, "'"))
    runCmd(enqueue_job_cmd)

    return jobId


def customLatentTrajectory(latentsPoints, projectId, workspaceId, trainingJobId):
    """Output the trajectory as a new output in CryoSPARC.
       The resulting trajectory may be used as input to the 3D Flex Generator job
       to generate a volume series along the trajectory."""
    from cryosparc.tools import CryoSPARC

    credentials = _getCredentials()
    if not credentials[0]:
        logger.error("Error obtaining cryoSPARC's credentials: %s" % credentials[1])
        raise Exception("Error obtaining cryoSPARC's credentials: %s" % credentials[1])

    credentials = credentials[1]
    cs = CryoSPARC(license=credentials['license'],
                   host=credentials['host'],
                   base_port=int(credentials['base_port']),
                   email=credentials['email'],
                   password=credentials['password'])

    project = cs.find_project(projectId)
    particles = project.find_job(trainingJobId).load_output("particles")
    numComponents = int(len([x for x in particles.fields() if "components_mode" in x]) / 2)
    slot_spec = [{"dtype": "components", "prefix": f"components_mode_{k}", "required": True} for k in
                 range(numComponents)]
    job = project.create_external_job(workspaceId, "Custom Latents")
    job.connect("particles", trainingJobId, "particles", slots=slot_spec)

    if len(latentsPoints.shape) == 1:
        latentsPoints = latentsPoints[None, ...]

    latentsDSet = job.add_output(
        type="particle",
        name="latents",
        slots=slot_spec,
        title="Latents",
        alloc=len(latentsPoints),
    )

    for k in range(numComponents):
        latentsDSet[f"components_mode_{k}/component"] = k
        latentsDSet[f"components_mode_{k}/value"] = latentsPoints[:, k]

    # Save the output
    with job.run():
        job.save_output("latents", latentsDSet)

    return job.uid


def runFlexGeneratorJob(trainingJobId, customLatentsJobId, projectId, workspaceId, gpu=0, lane='default'):
    """Generate a volume series along the trajectory using a flex model."""
    className = "flex_generate"
    gpusToUse = [gpu]
    input_group_connect = {"flex_model": "%s.flex_model" % trainingJobId,
                           "latents": "%s.latents" % customLatentsJobId}
    params = {}

    run3DFlexGeneratorJob = enqueueJob(className,
                                       projectId,
                                       workspaceId,
                                       str(params).replace('\'', '"'),
                                       str(input_group_connect).replace('\'', '"'),
                                       lane, gpusToUse)

    waitForCryosparc(projectId,
                     run3DFlexGeneratorJob,
                     "An error occurred in the 3D Flex Training process. "
                     "Please, go to cryoSPARC software for more "
                     "details.")
    clearIntermediateResults(projectId,
                             run3DFlexGeneratorJob)

    return run3DFlexGeneratorJob


def generateFlexVolumes(latentsPoints, projectId, workspaceId, trainingJobId, gpu=0):
    """Load particle latent coordinates from a 3D Flex Training job and use the 3D Flex Generator job to
        generate a volume series along the trajectory.
        This method allows(FlexUtils plugin) visualizing specific regions or pathways through the latent conformational distribution
        of the particle."""
    try:
        latentTrajectoryJob = customLatentTrajectory(latentsPoints,
                                                     projectId,
                                                     workspaceId,
                                                     trainingJobId)
        flexGeneratorJob = runFlexGeneratorJob(trainingJobId,
                                               latentTrajectoryJob,
                                               projectId,
                                               workspaceId,
                                               gpu)

        return flexGeneratorJob
    except Exception as ex:
        raise Exception("Error generating the flex volume : %s" % ex)


def runCmd(cmd, printCmd=True):
    """ Runs a command and check its exit code. If different from 0 it raises an exception
    :parameter cmd command to run
    :parameter printCmd (default True) prints the command"""
    import subprocess
    if printCmd:
        logger.info(pwutils.greenStr("Running: %s" % cmd))
    else:
        logger.debug(pwutils.greenStr("Running: %s" % cmd))

    exitCode, cmdOutput = subprocess.getstatusoutput(cmd)

    if exitCode != 0:
        raise Exception("%s failed --> Exit code %s, message %s" % (cmd, exitCode, cmdOutput))

    return exitCode, cmdOutput.split('\n')[-1]


def waitForCryosparc(projectName, jobId, failureMessage, protocol=None):
    """ Waits for cryosparc to finish or fail a job
    :parameter projectName: Cryosparc project name
    :parameter jobId: cryosparc job id
    :parameter failureMessage: Message for the exception thrown in case job fails
    :returns job Status
    :raises Exception when parsing cryosparc's output looks wrong"""

    # While is needed here, cause waitJob has a timeout of 5 secs.
    while True:
        try:
            status = getJobStatus(projectName, jobId)
            if status not in STOP_STATUSES:
                waitJob(projectName, jobId)
                if protocol is not None:
                    jobStreamLog = getJobStreamlog(projectName, jobId)
                    jobStreamLogList = eval(jobStreamLog[1])
                    jobLogLastLine = protocol.getLogLine()
                    lenLog = len(jobStreamLogList)
                    if lenLog > jobLogLastLine:
                        protocol.setLogLine(lenLog)
                        for line in range(jobLogLastLine, lenLog):
                            logDict = jobStreamLogList[line]
                            if logDict['type'] == 'text' and 'text' in logDict and logDict['text']:
                                logger.info(logDict['text'])
                    else:
                        jobLogLastLine = len(jobStreamLogList) - 1
                        while jobLogLastLine:
                            logDict = jobStreamLogList[jobLogLastLine]
                            if logDict['type'] == 'text' and 'text' in logDict and logDict['text']:
                                logger.info(logDict['text'])
                                break
                            jobLogLastLine -= 1
            else:
                break
        except Exception as e:
            logger.error("Can't query cryoSPARC about the job %s. Maybe it needs a restart ? We'll wait 5 minutes" % jobId, exc_info=e)
            import time
            time.sleep(300)  # wait 5 minutes

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


def getJob(projectName, job):
    """
       Return the job
       """
    get_job_status_cmd = (getCryosparcProgram() +
                          ' %sget_job("%s", "%s")%s'
                          % ("'", projectName, job, "'"))

    job = runCmd(get_job_status_cmd, printCmd=False)
    return job


def getJobLog(projectName, job):
    """
       Get the full contents of the given job's standard output log
       """
    get_job_log_cmd = (getCryosparcProgram() +
                          ' %sget_job_log("%s", "%s")%s'
                          % ("'", projectName, job, "'"))

    logStr = runCmd(get_job_log_cmd, printCmd=False)
    return logStr


def getJobStreamlog(projectName, job):
    """
       Get a list of dictionaries representing the given job's event log
       """
    get_job_stream_log_cmd = (getCryosparcProgram() +
                          ' %sget_job_streamlog("%s", "%s")%s'
                          % ("'", projectName, job, "'"))

    logList = runCmd(get_job_stream_log_cmd, printCmd=False)
    return logList


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
    logger.info(pwutils.yellowStr("Removing intermediate results..."))
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


def userExist(email):
    """
    Return if an user exist into cryoSPARC
    """
    getUser_cmd = (getCryosparcProgram() + ' %sUserExists("%s")%s' % ("'", email, "'"))
    return runCmd(getUser_cmd, printCmd=False)[1] == 'True'


def getUserId(email):
    """Get the user Id taking into account the user email"""
    import ast
    getUser_cmd = (getCryosparcProgram() + ' %sGetUser("%s")%s' % ("'", email, "'"))
    user = runCmd(getUser_cmd, printCmd=False)
    return ast.literal_eval(user[1])['_id']


def _getCredentials():
    licence = _getLicenceFromFile()
    if licence is None:
        return False, 'Error obtaining cryoSPARC license'

    csIsRunning = isCryosparcRunning()
    if csIsRunning:
        hostName = getCryosparcEnvInformation('master_hostname')
        basePort = getCryosparcEnvInformation('port_app')

        email = Plugin.getUser()
        if email is None:
            return False, 'Error obtaining the cryoSPARC user'

        password = Plugin.getUserPassword()
        if password is None:
            return False, 'Error obtaining the %s password ' % email

        return True, {'license': licence,
                    'host': hostName,
                    'base_port': basePort,
                    'email': email,
                    'password': password}

    return False, 'Cryosparc is not running'


def getSchedulerLanes():
    """
     Returns a list of lanes that are registered with the master scheduler
     list of dicts -- information about each lane
     """
    _csLanes = ['default']
    _defaultLane = _csLanes[0]
    csValidate = cryosparcValidate()
    if not csValidate:
        try:
            lanes_info_cmd = (getCryosparcProgram() + " 'get_scheduler_lanes()'")
            _csLanes = runCmd(lanes_info_cmd, printCmd=False)[1]
            lanes_dict_list = eval(_csLanes)
            _csLanes = []
            for lanes in lanes_dict_list:
                _csLanes.append(lanes.get('name'))
            defaultLane = getCryosparcDefaultLane()
            _defaultLane = _csLanes[0] if defaultLane is None else defaultLane
            if _defaultLane not in _csLanes:
                logger.error("Couldn't get the lane %s to the cryoSPARC installation" % _defaultLane)
                _defaultLane = _csLanes[0]
        except Exception:
           logger.error("Couldn't get Cryosparc's lanes")
    return _csLanes, _defaultLane


def addComputeSectionParams(form, allowMultipleGPUs=True, needGPU=True):
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
        if isCryosparcStandalone():
            versionAllowGPUs = parse_version(getCryosparcVersion()) >= parse_version(V3_0_0)
        else:
            versionAllowGPUs = False
    # Code is failing to get CS info, either stop or some error
    except Exception:
        # ... we assume is a modern version
        versionAllowGPUs = True

    if needGPU and versionAllowGPUs:
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

    defaultLane = getCryosparcDefaultLane()
    if defaultLane is None:
        defaultLane = 'default'
    form.addParam('compute_lane', StringParam, default=defaultLane,
                  label='Lane name:', readOnly=True,
                  help='The scheduler lane name to add the protocol execution')

    from .protocols import ProtCryo2D
    if not isCryosparcStandalone() and isinstance(form._protocol, ProtCryo2D):
        form.addParam('compute_num_gpus', IntParam, default=1,
                      label='Number of GPUs to compute:',
                      help='Number of GPUs to compute:')


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
        logger.error("Unable to execute the copy: Files or directory does not exist: ", exc_info=ex)


def matchItemRow(item, row):
    """
    Matches an item with a row from a dataset by comparing its index and filename.

    :param item: The item object containing metadata such as index and filename.
    :param row: A dictionary-like object representing a row of metadata,
                expected to contain RELIONCOLUMNS.rlnImageName.
    :return: True if the item's index and filename match the extracted row data, False otherwise.
    """
    imageName = row.get(RELIONCOLUMNS.rlnImageName.value)
    if not imageName:
        return True  # If no image name is found in row, assume a match

    try:
        index, filePath = imageName.split('@')
        rowFileName = os.path.splitext(os.path.basename(filePath).split('_', 1)[-1])[0]  # Remove leading digits and underscore and extension
        itemFileName = os.path.splitext(os.path.basename(item.getFileName()))[0]  # Remove extension
        return int(index) == item.getIndex() and itemFileName == rowFileName
    except (ValueError, AttributeError):
        return False  # In case of unexpected format, assume no match

