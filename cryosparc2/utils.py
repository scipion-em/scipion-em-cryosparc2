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
import os
import commands
import pyworkflow.utils as pwutils
from pyworkflow.em import SCIPION_SYM_NAME
from pyworkflow.em.constants import (SYM_CYCLIC, SYM_TETRAHEDRAL,
                                     SYM_OCTAHEDRAL, SYM_I222, SYM_I222r)
from pyworkflow.protocol.params import EnumParam, IntParam, Positive
from cryosparc2.constants import CS_SYM_NAME, SYM_DIHEDRAL_Y

STATUS_FAILED = "failed"
STATUS_ABORTED = "aborted"
STATUS_COMPLETED = "completed"
STATUS_KILLED = "killed"
STATUS_RUNNING = "running"
STATUS_QUEUED = "queued"
STATUS_LAUNCHED = "launched"
STATUS_STARTED = "started"

STOP_STATUSES = [STATUS_ABORTED, STATUS_COMPLETED, STATUS_FAILED, STATUS_KILLED]
ACTIVE_STATUSES = [STATUS_QUEUED, STATUS_RUNNING, STATUS_STARTED,
                   STATUS_LAUNCHED]


def getCryosparcDir():
    """
    Get the root directory where cryoSPARC code and dependencies are installed.
    """
    return os.environ['CRYOSPARC_DIR']


def getCryosparcProgram():
    """
    Get the cryosparc program to launch any command
    """
    return os.path.join(getCryosparcDir(),
                        'cryosparc2_master/bin/cryosparcm cli')


def cryosparcExist():
    """
    Determine if cryosparc software exist
    """
    msg = []
    if not os.path.exists(getCryosparcDir()):
       msg.append(('The cryoSPARC software do not exist in %s. Please install it')
                  % str(os.environ['CRYOSPARC_DIR']))
    return msg


def isCryosparcRunning():
    """
    Determine if cryosparc services are running
    """
    msg = []
    test_conection_cmd = (getCryosparcProgram() +
                                ' %stest_connection()%s ' % ("'", "'"))
    test_conection = commands.getstatusoutput(test_conection_cmd)
    status = test_conection[0]

    if status != 0:
        msg = ['Failed to establish a new connection with cryoSPARC. Please, '
               'restart the cryoSPARC services.']

    return msg


def getCryosparcUser():
    """
    Get the full name of the initial admin account
    """
    return os.environ['CRYOSPARC_USER']


def getCryosparcSSD():
    """
    Get the path on the worker node to a writable directory residing on the
    local SSD
    """
    if os.environ.get('CRYOSSD_DIR') is None:
        cryoSSD_Dir = os.path.join(str(getCryosparcDir()), 'cryo_ssd')
        if not os.path.exists(cryoSSD_Dir):
            os.mkdir(cryoSSD_Dir)
    elif os.path.exists(os.environ['CRYOSSD_DIR']):
        cryoSSD_Dir = os.environ['CRYOSSD_DIR']
    else:
        cryoSSD_Dir = os.path.join(str(getCryosparcDir()), 'cryo_ssd')
        os.mkdir(cryoSSD_Dir)
    return cryoSSD_Dir


def getProjectPath(projectDir):
    """
    Gets all projects of given path .
    projectDir: Folder path to get subfolders.
    returns: Set with all subfolders.
    """
    folderPaths = os.listdir(projectDir)
    return folderPaths


def getJobLog(projectDirName, projectName, job):
    """
    Return the job log
    """
    return os.path.join(getCryosparcSSD(), projectDirName, projectName, job,
                        'job.log')


def createEmptyProject(projectDir, projectTitle):
    """
    create_empty_project(owner_user_id, project_container_dir, title=None,
                            desc=None)
    """

    create_empty_project_cmd = (getCryosparcProgram() +
                                ' %screate_empty_project("%s", "%s", "%s")%s '
                                % ("'", str(getCryosparcUser()),
                                   str(projectDir), str(projectTitle), "'"))

    return commands.getstatusoutput(create_empty_project_cmd)


def createProjectDir(project_container_dir):
    """
    Given a "root" directory, create a project (PXXX) dir if it doesn't already
     exist
    :param project_container_dir: the "root" directory in which to create the
                                  project (PXXX) directory
    :param projectName: the name of the project
    :returns: str - the final path of the new project dir with shell variables
              still in the returned path (the path should be expanded every
              time it is used)
    """
    create_project_dir_cmd = (getCryosparcProgram() +
                             ' %scheck_or_create_project_container_dir("%s")%s '
                             % ("'", project_container_dir, "'"))
    return commands.getstatusoutput(create_project_dir_cmd)


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
    return commands.getstatusoutput(create_work_space_cmd)


def doJob(jobType, projectName, workSpaceName, params, input_group_conect):
    """
    do_job(job_type, puid='P1', wuid='W1', uuid='devuser', params={},
           input_group_connects={})
    """
    do_job_cmd = (getCryosparcProgram() +
                  ' %sdo_job("%s","%s","%s", "%s", %s, %s)%s' %
                  ("'", jobType, projectName, workSpaceName, getCryosparcUser(),
                   params, input_group_conect, "'"))

    print(pwutils.greenStr(do_job_cmd))
    return commands.getstatusoutput(do_job_cmd)


def getJobStatus(projectName, job):
    """
    Return the job status
    """
    get_job_status_cmd = (getCryosparcProgram() +
                          ' %sget_job_status("%s", "%s")%s'
                          % ("'", projectName, job, "'"))

    status = commands.getstatusoutput(get_job_status_cmd)
    return status[-1].split()[-1]


def waitJob(projectName, job):
    """
    Wait while the job not finished
    """
    wait_job_cmd = (getCryosparcProgram() +
                    ' %swait_job_complete("%s", "%s")%s'
                    % ("'", projectName, job, "'"))
    commands.getstatusoutput(wait_job_cmd)


def get_job_streamlog(projectName, job, fileName):

    get_job_streamlog_cmd = (getCryosparcProgram() +
                             ' %sget_job_streamlog("%s", "%s")%s%s'
                             % ("'", projectName, job, "'", ">" + fileName))

    commands.getstatusoutput(get_job_streamlog_cmd)


def killJob(projectName, job):
    """
     Kill a Job (if running)
    :param projectName: the uid of the project that contains the job to delete
    :param job: the uid of the job to delete
    """
    kill_job_cmd = (getCryosparcProgram() +
                    ' %skill_job("%s", "%s")%s'
                    % ("'", projectName, job, "'"))
    commands.getstatusoutput(kill_job_cmd)


def getSystemInfo():
    """
    Get the cryoSPARC system information
    """
    system_info_cmd = (getCryosparcProgram() + ' %sget_system_info()%s') % ("'", "'")
    return commands.getstatusoutput(system_info_cmd)


def addSymmetryParam(form):
    """
    Add the symmetry param with the conventions
    :param form:
    :return:
    """
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
                       "If no symmetry is present, use C1. Enforcing symmetry "
                       "above C1 is not recommended for ab-initio reconstruction"
                  )
    form.addParam('symmetryOrder', IntParam, default=1,
                  condition='symmetryGroup==%d or symmetryGroup==%d' %
                            (SYM_DIHEDRAL_Y - 11, SYM_CYCLIC),
                  label='Symmetry Order',
                  validators=[Positive],
                  help='Order of cyclic symmetry.')


def getSymmetry(symmetryGroup, symmetryOrder):
    """
    Get the symmetry(string) taking into account the symmetry convention
    """
    if symmetryGroup == 1:
        symmetry = CS_SYM_NAME[SYM_DIHEDRAL_Y][0] + str(symmetryOrder)
    else:
        symmetry = CS_SYM_NAME[symmetryGroup][0]
        if symmetryGroup == SYM_CYCLIC:
            symmetry = symmetry + str(symmetryOrder)
        elif (symmetryGroup == SYM_I222 or
              symmetryGroup == SYM_I222r):
            symmetry += CS_SYM_NAME[symmetryGroup][1]
    return symmetry