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
import shutil
import time
import shlex
import subprocess
from packaging.version import parse as parse_version
from packaging.version import Version


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
STATUS_WAITING = "waiting"

STOP_STATUSES = [STATUS_ABORTED, STATUS_COMPLETED, STATUS_FAILED, STATUS_KILLED]
ACTIVE_STATUSES = [STATUS_QUEUED, STATUS_RUNNING, STATUS_STARTED,
                   STATUS_LAUNCHED, STATUS_BUILDING, STATUS_WAITING]

# Module variables
_csVersion = None  # Lazy variable: never use it directly. Use getCryosparcVersion instead

# logging variable
logger = logging.getLogger(__name__)


def _normalizeCryosparcVersion(version):
    version = str(version or V_UNKNOWN).strip()
    version = version.split('+')[0].split('-')[0].replace('_', '.')
    version = version.lstrip('vV')
    if not version:
        version = '0.0.0'
    return f'v{version}'


def _getCryosparcVersionForRouting():
    try:
        return _normalizeCryosparcVersion(_getCryosparcVersionFromFile())
    except Exception:
        if _csVersion is not None:
            return _normalizeCryosparcVersion(_csVersion)
        return _normalizeCryosparcVersion(V_UNKNOWN)


def _isCryosparcV5OrNewer():
    version = _getCryosparcVersionForRouting().lstrip('vV')
    return parse_version(version) >= parse_version('5.0.0')


def getVersionedEnumValue(index, legacyValues, v5Values=None, switchVersion=V5_0_0):
    """
    Return the correct CryoSPARC enum value depending on the installed version.
    """
    values = legacyValues

    if v5Values is not None:
        cryosparcVersion = parse_version(getCryosparcVersion())
        if cryosparcVersion >= parse_version(switchVersion):
            values = v5Values

    return values[index]


def _runCommandRaw(cmd, printCmd=True):
    if printCmd:
        logger.info(pwutils.greenStr("Running: %s" % cmd))
    else:
        logger.debug(pwutils.greenStr("Running: %s" % cmd))

    exitCode, cmdOutput = subprocess.getstatusoutput(cmd)

    if exitCode != 0:
        raise Exception("%s failed --> Exit code %s, message %s" % (cmd, exitCode, cmdOutput))

    return exitCode, cmdOutput


def _runCliExpression(expression, printCmd=True):
    cmd = "%s %s" % (getCryosparcProgram(), shlex.quote(expression))
    return _runCommandRaw(cmd, printCmd=printCmd)


def _runCliValue(expression, printCmd=True):
    _, output = _runCliExpression(expression, printCmd=printCmd)
    output = output.strip()

    if output == "":
        return None

    try:
        return ast.literal_eval(output)
    except Exception:
        return output


def _pythonLiteral(value, default=None):
    if value is None:
        return {} if default is None else default

    if isinstance(value, (dict, list, tuple, bool, int, float)):
        return value

    text = str(value).strip()
    if text == "":
        return {} if default is None else default

    return ast.literal_eval(text)


def _parseResultTarget(target):
    sourceJobUid, sourceOutputName, sourceResultName = str(target).split(".", 2)
    return sourceJobUid, sourceOutputName, sourceResultName


def _runCryosparcmSubcommand(*parts, printCmd=False):
    cmd = " ".join(
        [getCryosparcProgram("")] +
        [shlex.quote(str(p)) for p in parts if p is not None]
    )
    return _runCommandRaw(cmd, printCmd=printCmd)


def _writeTextFile(path, text):
    with open(path, "w", encoding="utf-8", errors="replace") as fh:
        fh.write("" if text is None else str(text))


def _normalizeJobStatus(status):
    return str(status or "").strip().strip("'\"").lower()


def _normalizeSourceOutputName(sourceOutputName):
    name = str(sourceOutputName)

    legacyToV5OutputMap = {
        "imported_volume_1.map": "imported_volume_1",
        "imported_volume.map": "imported_volume",

        "imported_volume_1.map_half_A": "imported_volume_1",
        "imported_volume_1.map_half_B": "imported_volume_1",
        "imported_volume.map_half_A": "imported_volume",
        "imported_volume.map_half_B": "imported_volume",

        "imported_mask_1.map": "imported_mask_1",
        "imported_mask.map": "imported_mask",
    }

    return legacyToV5OutputMap.get(name, name)


def _parseConnectionTarget(target):
    sourceJobUid, sourceOutputName = str(target).split(".", 1)

    if _isCryosparcV5OrNewer():
        sourceOutputName = _normalizeSourceOutputName(sourceOutputName)

    return sourceJobUid, sourceOutputName

def _tryV5StructuredEventLog(projectUid, jobUid):
    """
    Try a few likely v5 API method names first.
    If one works, return a Python object that can be repr()-written and
    parsed later by protocol_base.py using ast.literal_eval.
    """
    expressions = [
        "[(e.model_dump() if hasattr(e, 'model_dump') else "
        "(e.dict() if hasattr(e, 'dict') else e)) "
        "for e in api.jobs.get_event_log(%s, %s)]"
        % (repr(str(projectUid)), repr(str(jobUid))),

        "[(e.model_dump() if hasattr(e, 'model_dump') else "
        "(e.dict() if hasattr(e, 'dict') else e)) "
        "for e in api.jobs.get_event_logs(%s, %s)]"
        % (repr(str(projectUid)), repr(str(jobUid))),

        "[(e.model_dump() if hasattr(e, 'model_dump') else "
        "(e.dict() if hasattr(e, 'dict') else e)) "
        "for e in api.jobs.get_streamlog(%s, %s)]"
        % (repr(str(projectUid)), repr(str(jobUid))),

        "[(e.model_dump() if hasattr(e, 'model_dump') else "
        "(e.dict() if hasattr(e, 'dict') else e)) "
        "for e in api.jobs.get_job_streamlog(%s, %s)]"
        % (repr(str(projectUid)), repr(str(jobUid))),
    ]

    for expr in expressions:
        try:
            data = _runCliValue(expr, printCmd=False)
            if data is not None:
                return data
        except Exception:
            pass

    return None

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
    Get the cryosparc program to launch any command.
    mode="cli" returns ".../cryosparcm cli"
    mode="" returns ".../cryosparcm"
    """
    csDir = getCryosparcDir()

    if csDir is not None:
        command = 'cryosparcm' if not mode else 'cryosparcm %s' % mode

        if os.path.exists(os.path.join(csDir, CRYOSPARC_MASTER, "bin")):
            return os.path.join(csDir, CRYOSPARC_MASTER, "bin", command)
        else:
            return os.path.join(csDir, 'cryosparc2_master', "bin", command)

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
    if getCryosparcProgram() is None:
        return False

    try:
        if _isCryosparcV5OrNewer():
            status = _runCliValue("api.health()", printCmd=False)
            return str(status).strip("'\"") == "OK"

        testConnectionCmd = (
            getCryosparcProgram() +
            ' %stest_connection()%s ' % ("'", "'")
        )
        exitCode, _ = subprocess.getstatusoutput(testConnectionCmd)
        return exitCode == 0
    except Exception:
        return False


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
    Get the cryoSPARC environment information.
    """
    import ast

    systemInfo = getSystemInfo()
    dictionary = systemInfo[1]

    if isinstance(dictionary, str):
        dictionary = ast.literal_eval(dictionary)

    return str(dictionary[envVar])


def getCryosparcVersion():
    """
    Gets cryosparc version:
    1) from version file
    2) from v5 api.config.get_version()
    3) from legacy system info
    """
    global _csVersion

    if _csVersion is None:
        try:
            _csVersion = _normalizeCryosparcVersion(_getCryosparcVersionFromFile())
        except Exception:
            try:
                if _isCryosparcV5OrNewer():
                    _csVersion = _normalizeCryosparcVersion(
                        _runCliValue("api.config.get_version()", printCmd=False)
                    )
                else:
                    _csVersion = _normalizeCryosparcVersion(
                        getCryosparcEnvInformation(VERSION)
                    )
            except Exception as e:
                logger.error(
                    "Couldn't get Cryosparc's version. Please review your config (%s)"
                    % Plugin.getUrl(),
                    exc_info=e
                )
                _csVersion = _normalizeCryosparcVersion(V_UNKNOWN)

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
    Get list of all projects available.
    Returns a legacy-compatible list of dicts.
    """
    if _isCryosparcV5OrNewer():
        rawProjects = _runCliValue(
            "[(p.uid, p.title, api.projects.get_directory(p.uid)) for p in api.projects.find()]",
            printCmd=False
        ) or []

        return [
            {
                "uid": uid,
                "title": title,
                "project_dir": projectDir,
            }
            for uid, title, projectDir in rawProjects
        ]

    projectsListCmd = (getCryosparcProgram() + ' %slist_projects()%s ' % ("'", "'"))
    cmd = runCmd(projectsListCmd, printCmd=False)[1]
    projectList = ast.literal_eval(cmd)
    return projectList


def getCryosparcWorkSpaces(projectId):
    """
    List all workspaces inside a given project.
    Returns a legacy-compatible list of dicts.
    """
    if _isCryosparcV5OrNewer():
        rawWorkspaces = _runCliValue(
            "[(w.uid, getattr(w, 'title', ''), getattr(w, 'description', None)) "
            "for w in api.workspaces.find(project_uid=%s)]" % repr(str(projectId)),
            printCmd=False
        ) or []

        return [
            {
                "uid": uid,
                "title": title,
                "description": description,
            }
            for uid, title, description in rawWorkspaces
        ]

    workspaceListCmd = (
        getCryosparcProgram() +
        ' %slist_workspaces("%s")%s ' % ("'", str(projectId), "'")
    )
    cmd = runCmd(workspaceListCmd, printCmd=False)[1]
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


def getCryosparcPreprocessLane():
    """
    Get the cryoSPARC lane for lightweight preprocessing tasks.
    If not defined, fallback to the default (heavy processing) lane.
    """
    return os.environ.get(CRYOSPARC_PREPROCESS_LANE, getCryosparcDefaultLane())


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
    Legacy-compatible project creation.
    Keep return contract so protocol_base.py does not need to change.
    """
    if _isCryosparcV5OrNewer():
        projectUid = _runCliValue(
            "api.projects.create(title=%s, description=None, parent_dir=%s).uid"
            % (repr(str(projectTitle)), repr(str(projectDir))),
            printCmd=False
        )
        return 0, "Created project %s" % projectUid

    createEmptyProjectCmd = (
        getCryosparcProgram() +
        ' %screate_empty_project("%s", "%s", "%s")%s '
        % ("'", str(getCryosparcUser()), str(projectDir), str(projectTitle), "'")
    )

    return runCmd(createEmptyProjectCmd, printCmd=False)


def getProjectInformation(project_uid, info='project_dir'):
    """
    Get information about a single project.
    """
    if _isCryosparcV5OrNewer():
        if info == 'project_dir':
            return str(
                _runCliValue(
                    "api.projects.get_directory(%s)" % repr(str(project_uid)),
                    printCmd=False
                )
            )

        return str(
            _runCliValue(
                "getattr(api.projects.find_one(%s), %s)"
                % (repr(str(project_uid)), repr(str(info))),
                printCmd=False
            )
        )

    getProjectCmd = (
        getCryosparcProgram() +
        ' %sget_project("%s")%s '
        % ("'", str(project_uid), "'")
    )

    projectInfo = runCmd(getProjectCmd, printCmd=False)
    dictionary = ast.literal_eval(projectInfo[1])
    return str(dictionary[info])


def getUserToken(email):
    if _isCryosparcV5OrNewer():
        rawUser = _runCliValue(
            "api.users.find_one(%s)" % repr(str(email)),
            printCmd=False
        )
        return 0, str(rawUser)

    getUserCmd = (
        getCryosparcProgram() +
        ' %sGetUser("%s")%s ' % ("'", str(email), "'")
    )
    return runCmd(getUserCmd, printCmd=False)


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
    Given a root directory, ensure it exists and is usable for project creation.
    Keep legacy return contract: (exitCode, path)
    """
    os.makedirs(project_container_dir, exist_ok=True)

    if _isCryosparcV5OrNewer():
        _runCliValue(
            "api.projects.check_directory(path=%s)" % repr(str(project_container_dir)),
            printCmd=False
        )
        return 0, project_container_dir

    createProjectDirCmd = (
        getCryosparcProgram() +
        ' %scheck_or_create_project_container_dir("%s")%s '
        % ("'", project_container_dir, "'")
    )
    return runCmd(createProjectDirCmd, printCmd=False)


def createEmptyWorkSpace(projectName, workspaceTitle, workspaceComment):
    """
    Legacy-compatible workspace creation.
    Keep return contract so protocol_base.py does not need to change.
    """
    if _isCryosparcV5OrNewer():
        workspaceUid = _runCliValue(
            "api.workspaces.create(%s, title=%s, description=%s).uid"
            % (
                repr(str(projectName)),
                repr(str(workspaceTitle)),
                repr(str(workspaceComment)),
            ),
            printCmd=False
        )
        return 0, "Created workspace %s" % workspaceUid

    createWorkSpaceCmd = (
        getCryosparcProgram() +
        ' %screate_empty_workspace("%s", "%s", "%s", "%s", "%s")%s '
        % (
            "'",
            projectName,
            str(getCryosparcUser(userId=False)),
            "None",
            str(workspaceTitle),
            str(workspaceComment),
            "'"
        )
    )
    return runCmd(createWorkSpaceCmd, printCmd=False)


def _getProtocolPreprocessLane(protocol):
    preprocessLane = getattr(protocol, 'preprocessLane', None)
    if preprocessLane:
        return preprocessLane

    defaultPreprocessLane = getCryosparcPreprocessLane()
    if defaultPreprocessLane is not None:
        return str(defaultPreprocessLane)

    computeLane = getattr(protocol, 'lane', None)
    if computeLane:
        return computeLane

    if hasattr(protocol, 'getAttributeValue'):
        laneValue = protocol.getAttributeValue('compute_lane')
        if laneValue is not None:
            return str(laneValue)

    return None


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

    preprocessLane = _getProtocolPreprocessLane(protocol)
    import_particles = enqueueJob(className, protocol.projectName, protocol.workSpaceName,
                                  str(params).replace('\'', '"'), '{}', preprocessLane)

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

    preprocessLane = _getProtocolPreprocessLane(protocol)
    importedVolume = enqueueJob(className, protocol.projectName,
                                protocol.workSpaceName,
                                str(params).replace('\'', '"'), '{}',
                                preprocessLane)

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

    preprocessLane = _getProtocolPreprocessLane(protocol)
    import_particles = enqueueJob(className, protocol.projectName, protocol.workSpaceName,
                                  str(params).replace('\'', '"'), '{}', preprocessLane)

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
    Queue a CryoSPARC job.
    """
    from pyworkflow.object import String

    cryosparcVersion = getCryosparcVersion()
    standaloneInstallation = isCryosparcStandalone()

    if _isCryosparcV5OrNewer():
        paramsDict = _pythonLiteral(params, default={})
        inputGroupConnect = _pythonLiteral(input_group_connect, default={})
        gpus = [] if not gpusToUse else list(gpusToUse)

        projectUid = str(projectName)
        workspaceUid = str(workSpaceName)

        jobUid = _runCliValue(
            "api.jobs.create(%s, %s, %s, type=%s).uid"
            % (
                repr(projectUid),
                repr(workspaceUid),
                repr(paramsDict),
                repr(str(jobType)),
            ),
            printCmd=False
        )

        jobId = String(jobUid)

        if inputGroupConnect:
            for inputName, sourceTarget in inputGroupConnect.items():
                sourceJobUid, sourceOutputName = _parseConnectionTarget(sourceTarget)
                _runCliValue(
                    "api.jobs.connect(%s, %s, %s, source_output_name=%s, source_job_uid=%s)"
                    % (
                        repr(projectUid),
                        repr(str(jobId)),
                        repr(str(inputName)),
                        repr(str(sourceOutputName)),
                        repr(str(sourceJobUid)),
                    ),
                    printCmd=False
                )

        if group_connect is not None:
            for inputName, valuesList in group_connect.items():
                for sourceTarget in valuesList:
                    sourceJobUid, sourceOutputName = _parseConnectionTarget(sourceTarget)
                    _runCliValue(
                        "api.jobs.connect(%s, %s, %s, source_output_name=%s, source_job_uid=%s)"
                        % (
                            repr(projectUid),
                            repr(str(jobId)),
                            repr(str(inputName)),
                            repr(str(sourceOutputName)),
                            repr(str(sourceJobUid)),
                        ),
                        printCmd=False
                    )

        if result_connect is not None:
            for key, value in result_connect.items():
                inputName, inputSlot, resultName = _parseResultInputKey(key)
                sourceJobUid, sourceOutputName, sourceResultName = _parseResultTarget(value)

                _runCliValue(
                    "api.jobs.connect_result(%s, %s, %s, %s, %s, "
                    "source_output_name=%s, source_result_name=%s, source_job_uid=%s)"
                    % (
                        repr(projectUid),
                        repr(str(jobId)),
                        repr(str(inputName)),
                        inputSlot,
                        repr(str(resultName)),
                        repr(str(sourceOutputName)),
                        repr(str(sourceResultName)),
                        repr(str(sourceJobUid)),
                    ),
                    printCmd=False
                )

        hostname = None
        if standaloneInstallation:
            try:
                hostname = getCryosparcEnvInformation('master_hostname')
            except Exception:
                hostname = None

        _runCliValue(
            "api.jobs.enqueue(%s, %s, lane=%s, gpus=%s, no_check_inputs_ready=%s)"
            % (
                repr(projectUid),
                repr(str(jobId)),
                repr(str(lane)),
                repr(gpus),
                repr(False),
            ),
            printCmd=False
        )

        logger.info(pwutils.greenStr("Got %s for JobId" % jobId))
        return jobId

    # Legacy path kept as-is below
    if parse_version(V3_0_0) <= parse_version(cryosparcVersion) < parse_version(V4_3_1):
        make_job_cmd = (getCryosparcProgram() +
                        ' %smake_job("%s","%s","%s", "%s", "None", "None", %s, %s, "False", 0)%s' %
                        ("'", jobType, projectName, workSpaceName,
                         getCryosparcUser(),
                         params, input_group_connect, "'"))

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


def waitForCryosparc(projectName, jobName, errorMsg, protocol=None, sleepTime=15):
    """
    Wait until the cryoSPARC job reaches a stop status and stream job logs
    into the Scipion logger when a protocol instance is available.
    """
    projectName = str(projectName)
    jobName = str(jobName)

    if isinstance(protocol, (int, float)):
        sleepTime = protocol
        protocol = None

    status = None

    while True:
        try:
            status = getJobStatus(projectName, jobName)

            if status in STOP_STATUSES:
                break

            _logCryosparcJobEvents(projectName, jobName, protocol=protocol)

            if status in ACTIVE_STATUSES:
                waitJob(projectName, jobName)
            else:
                time.sleep(sleepTime)

        except Exception as ex:
            logger.error(
                "Can't query cryoSPARC about the job %s. Maybe it needs a restart? "
                "We'll wait 5 minutes" % jobName,
                exc_info=ex
            )
            time.sleep(300)

    _logCryosparcJobEvents(projectName, jobName, protocol=protocol)

    if status != STATUS_COMPLETED:
        raise Exception("%s Current status: %s" % (errorMsg, status))

    return status


def _getCryosparcJobLogEvents(projectName, jobName):
    """
    Return cryoSPARC job log events as a list.
    """
    if _isCryosparcV5OrNewer():
        events = _tryV5StructuredEventLog(projectName, jobName)
    else:
        events = getJobStreamlog(projectName, jobName)[1]

    if events is None:
        return []

    try:
        events = _pythonLiteral(events, default=[])
    except Exception:
        return []

    if isinstance(events, dict):
        events = events.get("events", events.get("logs", []))

    if not isinstance(events, list):
        return []

    return events


def _getCryosparcLogText(logEvent):
    """
    Extract a printable text message from a cryoSPARC log event.
    """
    if isinstance(logEvent, str):
        return logEvent.strip()

    if not isinstance(logEvent, dict):
        return None

    for key in ("text", "message", "msg", "description"):
        value = logEvent.get(key)
        if value:
            return str(value).strip()

    return None


def _logCryosparcJobEvents(projectName, jobName, protocol=None):
    """
    Log new cryoSPARC job events into Scipion logger.
    """
    events = _getCryosparcJobLogEvents(projectName, jobName)
    if not events:
        return

    if protocol is None:
        lastEvent = events[-1]
        text = _getCryosparcLogText(lastEvent)
        if text:
            logger.info(text)
        return

    jobLogLastLine = protocol.getLogLine()
    lenLog = len(events)

    if lenLog > jobLogLastLine:
        protocol.setLogLine(lenLog)
        for line in range(jobLogLastLine, lenLog):
            text = _getCryosparcLogText(events[line])
            if text:
                logger.info(text)
    else:
        jobLogLastLine = lenLog - 1
        while jobLogLastLine >= 0:
            text = _getCryosparcLogText(events[jobLogLastLine])
            if text:
                logger.info(text)
                break
            jobLogLastLine -= 1


def getJobStatus(projectName, jobId):
    """
    Get job status in a version-compatible way.
    """
    if _isCryosparcV5OrNewer():
        status = _runCliValue(
            "api.jobs.get_status(%s, %s)"
            % (repr(str(projectName)), repr(str(jobId))),
            printCmd=False
        )
        return _normalizeJobStatus(status)

    getStatusCmd = (
        getCryosparcProgram() +
        ' %sget_job_status("%s","%s")%s '
        % ("'", str(projectName), str(jobId), "'")
    )
    _, output = runCmd(getStatusCmd, printCmd=False)
    return _normalizeJobStatus(output)


def getJob(projectName, job):
    """
       Return the job
       """
    get_job_status_cmd = (getCryosparcProgram() +
                          ' %sget_job("%s", "%s")%s'
                          % ("'", projectName, job, "'"))

    job = runCmd(get_job_status_cmd, printCmd=False)
    return job


def get_job_log(projectName, jobId, outputFile):
    """
    Store stdout/stderr job log into outputFile.
    """
    if _isCryosparcV5OrNewer():
        _, rawLog = _runCryosparcmSubcommand(
            "job", "log", str(projectName), str(jobId),
            printCmd=False
        )
        _writeTextFile(outputFile, rawLog)
        return 0, outputFile

    legacyCmd = (
        getCryosparcProgram() +
        ' %sget_job_log("%s","%s")%s '
        % ("'", str(projectName), str(jobId), "'")
    )
    _, rawOutput = _runCommandRaw(legacyCmd, printCmd=False)
    _writeTextFile(outputFile, rawOutput)
    return 0, outputFile


def getJobStreamlog(projectName, job):
    """
       Get a list of dictionaries representing the given job's event log
       """
    get_job_stream_log_cmd = (getCryosparcProgram() +
                          ' %sget_job_streamlog("%s", "%s")%s'
                          % ("'", projectName, job, "'"))

    logList = runCmd(get_job_stream_log_cmd, printCmd=False)
    return logList


def waitJob(projectName, job, sleepTime=15):
    """
    Wait while the job is not finished.
    """
    projectName = str(projectName)
    job = str(job)

    if _isCryosparcV5OrNewer():
        while True:
            status = getJobStatus(projectName, job)
            if status in STOP_STATUSES:
                return status
            time.sleep(sleepTime)

    wait_job_cmd = (getCryosparcProgram() +
                    ' %swait_job_complete("%s", "%s")%s'
                    % ("'", projectName, job, "'"))
    runCmd(wait_job_cmd, printCmd=False)


def get_job_streamlog(projectName, jobId, outputFile):
    """
    Store the job event log into outputFile.

    Legacy path writes the old Python-literal representation.
    v5 path first tries to recover a structured API payload so downstream
    ast.literal_eval() continues to work. If that fails, it falls back to the
    documented `cryosparcm job events` command and writes plain text.
    """
    if _isCryosparcV5OrNewer():
        structured = _tryV5StructuredEventLog(projectName, jobId)

        if structured is not None:
            _writeTextFile(outputFile, repr(structured))
            return 0, outputFile

        _, rawEvents = _runCryosparcmSubcommand(
            "job", "events", str(projectName), str(jobId),
            printCmd=False
        )
        _writeTextFile(outputFile, rawEvents)
        return 0, outputFile

    legacyCmd = (
        getCryosparcProgram() +
        ' %sget_job_streamlog("%s","%s")%s '
        % ("'", str(projectName), str(jobId), "'")
    )
    _, rawOutput = _runCommandRaw(legacyCmd, printCmd=False)
    _writeTextFile(outputFile, rawOutput)
    return 0, outputFile


def killJob(projectName, jobId):
    """
    Kill a job.
    """
    if _isCryosparcV5OrNewer():
        return _runCryosparcmSubcommand(
            "job", "kill", str(projectName), str(jobId),
            printCmd=False
        )

    killJobCmd = (
        getCryosparcProgram() +
        ' %skill_job("%s","%s")%s '
        % ("'", str(projectName), str(jobId), "'")
    )
    return runCmd(killJobCmd, printCmd=False)


def clearJob(projectName, jobId):
    """
    Clear a job.
    """
    if _isCryosparcV5OrNewer():
        return _runCryosparcmSubcommand(
            "job", "clear", str(projectName), str(jobId),
            printCmd=False
        )

    clearJobCmd = (
        getCryosparcProgram() +
        ' %sclear_job("%s","%s")%s '
        % ("'", str(projectName), str(jobId), "'")
    )
    return runCmd(clearJobCmd, printCmd=False)


def clearIntermediateResults(projectName, jobName=None, workspaceName=None, alwaysKeepFinal=True):
    """
    Clear intermediate results in a version-compatible way.

    Supported use cases kept intentionally broad so existing call sites do not
    need to change:
      - clearIntermediateResults(projectUid, jobUid)
      - clearIntermediateResults(projectUid, workspaceName="W1")
      - clearIntermediateResults(projectUid)

    In v5 the exact public CLI/API entrypoint for this action is not clearly
    documented in the references we have, so we try several likely API forms
    and degrade gracefully if none are available.
    """
    projectUid = str(projectName) if projectName is not None else None
    jobUid = str(jobName) if jobName is not None else None
    workspaceUid = str(workspaceName) if workspaceName is not None else None

    if not projectUid:
        logger.warning("clearIntermediateResults called without project UID.")
        return 0, ""

    if _isCryosparcV5OrNewer():
        attempts = []

        # Most likely modern API shapes
        if jobUid is not None:
            attempts.extend([
                "api.jobs.clear_intermediate_results(%s, %s, always_keep_final=%s)"
                % (repr(projectUid), repr(jobUid), repr(bool(alwaysKeepFinal))),
                "api.jobs.clear_intermediate_results(project_uid=%s, job_uid=%s, always_keep_final=%s)"
                % (repr(projectUid), repr(jobUid), repr(bool(alwaysKeepFinal))),
                "api.projects.clear_intermediate_results(%s, job_uid=%s, always_keep_final=%s)"
                % (repr(projectUid), repr(jobUid), repr(bool(alwaysKeepFinal))),
                "api.projects.clear_intermediate_results(project_uid=%s, job_uid=%s, always_keep_final=%s)"
                % (repr(projectUid), repr(jobUid), repr(bool(alwaysKeepFinal))),
            ])

        if workspaceUid is not None:
            attempts.extend([
                "api.projects.clear_intermediate_results(%s, workspace_uid=%s, always_keep_final=%s)"
                % (repr(projectUid), repr(workspaceUid), repr(bool(alwaysKeepFinal))),
                "api.projects.clear_intermediate_results(project_uid=%s, workspace_uid=%s, always_keep_final=%s)"
                % (repr(projectUid), repr(workspaceUid), repr(bool(alwaysKeepFinal))),
            ])

        # Project-level clear
        attempts.extend([
            "api.projects.clear_intermediate_results(%s, always_keep_final=%s)"
            % (repr(projectUid), repr(bool(alwaysKeepFinal))),
            "api.projects.clear_intermediate_results(project_uid=%s, always_keep_final=%s)"
            % (repr(projectUid), repr(bool(alwaysKeepFinal))),
        ])

        # Legacy-like low-level entrypoints that may still be exposed in some v5 installs
        if jobUid is not None:
            attempts.extend([
                "clear_intermediate_results(project_uid=%s, job_uid=%s, always_keep_final=%s)"
                % (repr(projectUid), repr(jobUid), repr(bool(alwaysKeepFinal))),
                "api.clear_intermediate_results(project_uid=%s, job_uid=%s, always_keep_final=%s)"
                % (repr(projectUid), repr(jobUid), repr(bool(alwaysKeepFinal))),
            ])

        if workspaceUid is not None:
            attempts.extend([
                "clear_intermediate_results(project_uid=%s, workspace_uid=%s, always_keep_final=%s)"
                % (repr(projectUid), repr(workspaceUid), repr(bool(alwaysKeepFinal))),
                "api.clear_intermediate_results(project_uid=%s, workspace_uid=%s, always_keep_final=%s)"
                % (repr(projectUid), repr(workspaceUid), repr(bool(alwaysKeepFinal))),
            ])

        attempts.extend([
            "clear_intermediate_results(project_uid=%s, always_keep_final=%s)"
            % (repr(projectUid), repr(bool(alwaysKeepFinal))),
            "api.clear_intermediate_results(project_uid=%s, always_keep_final=%s)"
            % (repr(projectUid), repr(bool(alwaysKeepFinal))),
        ])

        lastError = None

        for expr in attempts:
            try:
                result = _runCliValue(expr, printCmd=False)
                logger.info(
                    pwutils.yellowStr(
                        "Intermediate results cleared with v5 expression: %s" % expr
                    )
                )
                return 0, "" if result is None else str(result)
            except Exception as ex:
                lastError = ex

        logger.warning(
            "Could not clear intermediate results for project=%s, job=%s, workspace=%s. "
            "Tried several v5 API forms and none succeeded. Continuing without failing the protocol. "
            "Last error: %s",
            projectUid, jobUid, workspaceUid, lastError
        )
        return 0, ""

    # Legacy <= 4.7 path
    if jobUid is not None:
        clearIntermediateExpr = (
                "clear_intermediate_results(project_uid=%s, job_uid=%s, always_keep_final=%s)"
                % (repr(projectUid), repr(jobUid), repr(bool(alwaysKeepFinal)))
        )
        return _runCliExpression(clearIntermediateExpr, printCmd=False)

    if workspaceUid is not None:
        clearIntermediateExpr = (
                "clear_intermediate_results(project_uid=%s, workspace_uid=%s, always_keep_final=%s)"
                % (repr(projectUid), repr(workspaceUid), repr(bool(alwaysKeepFinal)))
        )
        return _runCliExpression(clearIntermediateExpr, printCmd=False)

    clearIntermediateExpr = (
            "clear_intermediate_results(project_uid=%s, always_keep_final=%s)"
            % (repr(projectUid), repr(bool(alwaysKeepFinal)))
    )
    return _runCliExpression(clearIntermediateExpr, printCmd=False)


def getSystemInfo():
    """
    Get CryoSPARC system information.
    Returns the same contract as legacy code: (exitCode, repr(dictionary))
    """
    if _isCryosparcV5OrNewer():
        expr = (
            "(lambda s: "
            "s.model_dump() if hasattr(s, 'model_dump') else "
            "(s.dict() if hasattr(s, 'dict') else s.__dict__))"
            "(api.config.get_system_info())"
        )
        data = _runCliValue(expr, printCmd=False)
        return 0, repr(data)

    systemInfoCmd = (
        getCryosparcProgram() +
        ' %sget_system_info()%s ' % ("'", "'")
    )
    return runCmd(systemInfoCmd, printCmd=False)


def getCryosparcJobUrl(projectId, workspaceId=None, jobId=None):
    """
    Return the cryoSPARC browser URL for a project/workspace/job.
    """
    systemInfo = getSystemInfo()
    statusErrors = systemInfo[0]
    if statusErrors:
        return None

    systemInfo = systemInfo[1]
    if not isinstance(systemInfo, dict):
        systemInfo = _pythonLiteral(systemInfo, default={})
    if isinstance(systemInfo, str):
        systemInfo = _pythonLiteral(systemInfo, default={})
    if not isinstance(systemInfo, dict):
        return None

    masterHostname = systemInfo.get('master_hostname')
    portWebapp = systemInfo.get('port_webapp')
    portApp = systemInfo.get('port_app')
    version = systemInfo.get('version') or getCryosparcVersion()

    if not masterHostname:
        return None

    projectId = str(projectId)
    workspaceId = str(workspaceId) if workspaceId is not None else None
    jobId = str(jobId) if jobId is not None else None

    if parse_version(version) >= parse_version(V4_1_0):
        port = portApp or portWebapp
        if not port:
            return None

        if workspaceId:
            browseTarget = "%s-%s-J*" % (projectId, workspaceId)
        else:
            browseTarget = "%s-J*" % projectId

        url = "http://%s:%s/browse/%s" % (masterHostname, port, browseTarget)

        if jobId:
            url += "#job(%s-%s)" % (projectId, jobId)

        return url

    if not portWebapp:
        return None

    return "http://%s:%s/projects/%s/%s/%s" % (masterHostname, portWebapp, projectId, workspaceId, jobId)

def userExist(user):
    """
    Returns True if user exists.
    """
    if _isCryosparcV5OrNewer():
        cmd = "%s user exists --email %s" % (
            getCryosparcProgram(""),
            shlex.quote(str(user))
        )
        exitCode, _ = subprocess.getstatusoutput(cmd)
        return exitCode == 0

    userExistsCmd = (
        getCryosparcProgram() +
        ' %sUserExists("%s")%s ' % ("'", str(user), "'")
    )
    exitCode, _ = subprocess.getstatusoutput(userExistsCmd)
    return exitCode == 0


def getUserId(user):
    """
    Legacy compatibility helper.
    In v5, most user-facing APIs accept email directly, so keep email as-is.
    """
    if _isCryosparcV5OrNewer():
        return str(user)

    getUserIdCmd = (
        getCryosparcProgram() +
        ' %sget_id_by_email("%s")%s ' % ("'", str(user), "'")
    )
    return runCmd(getUserIdCmd, printCmd=False)[1]


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


def addPreprocessLaneParam(form):
    from pyworkflow.protocol.params import StringParam

    protocol = form._protocol

    defaultPreprocessLane = getCryosparcPreprocessLane()

    if not defaultPreprocessLane:
        defaultPreprocessLane = protocol.getAttributeValue('compute_lane')

    if not defaultPreprocessLane:
        defaultPreprocessLane = getCryosparcDefaultLane()

    if not defaultPreprocessLane:
        defaultPreprocessLane = 'default'

    defaultPreprocessLane = str(defaultPreprocessLane)

    form.addParam('preprocess_lane', StringParam,
                  default=defaultPreprocessLane,
                  label='Preprocessing lane name:',
                  readOnly=True,
                  help='Scheduler lane used for preprocessing imports '
                       '(particles, volumes, masks).')


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
    try:

        imageName = row.get(RELIONCOLUMNS.rlnImageName.value)
        if not imageName:
            return True  # If no image name is found in row, assume a match

        index, filePath = imageName.split('@')
        rowFileName = os.path.splitext(os.path.basename(filePath).split('_', 1)[-1])[0]  # Remove leading digits and underscore and extension
        itemFileName = os.path.splitext(os.path.basename(item.getFileName()))[0]  # Remove extension
        return int(index) == item.getIndex() and itemFileName == rowFileName
    except Exception:
        return False  # In case of unexpected format, assume no match

def parse_version(value: str) -> Version:
    value = str(value).strip()
    if value[:1] in {"V", "v"}:
        value = value[1:]
    return Version(value)

def _parseResultInputKey(key):
    """
    Parse result_connect keys like:
      volume.0.map_half_A

    Legacy code encoded the input slot inside the result name. In v5,
    connect_result expects the slot as a separate argument.
    """
    parts = str(key).split(".")

    inputName = parts[0]
    inputSlot = 0

    if len(parts) >= 3 and parts[1].isdigit():
        inputSlot = int(parts[1])
        resultName = ".".join(parts[2:])
    else:
        resultName = ".".join(parts[1:])

    return inputName, inputSlot, resultName
