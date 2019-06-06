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
from pyworkflow.utils import utils

def getCryosparcProgram():
    return os.path.join(os.environ['CRYOSPARC_DIR'],
                        'cryosparc2_master/bin/cryosparcm cli')


def getCryosparcUser():
    return os.environ['CRYOSPARC_USER']


def getCryosparcSSD():
    return os.environ['CRYOSSD_DIR']


def createEmptyProject():
    """ create_empty_project(owner_user_id, project_container_dir, title=None,
                            desc=None)
    """

    create_empty_project_cmd = (getCryosparcProgram() +
                                ' %screate_empty_project("%s+%s%s", "%s")%s '
                                % ("'", "'", str(getCryosparcUser()), "'",
                                str(getCryosparcSSD()), "'"))

    return commands.getstatusoutput(create_empty_project_cmd)


def createEmptyWorkSpace(projectName):
    """
    create_empty_workspace(project_uid, created_by_user_id,
                           created_by_job_uid=None,
                           title=None, desc=None)
    returns the new uid of the workspace that was created
    """

    create_work_space_cmd = (getCryosparcProgram() +
                             " \'create_empty_workspace(\"" +
                             projectName + "\", \"\'+" +
                             getCryosparcUser() + "\'\")\'")

    return commands.getstatusoutput(create_work_space_cmd)


def doJob(jobType, projectName, workSpaceName, params, input_group_conect):
    """
    do_job(job_type, puid='P1', wuid='W1', uuid='devuser', params={},
           input_group_connects={})
    """
    do_job_cmd = (getCryosparcProgram() +
                  ' %sdo_job("%s","%s","%s", "%s+%s%s", %s, %s)%s' %
                  ("'", jobType, projectName, workSpaceName, "'",
                   getCryosparcUser(), "'", params, input_group_conect, "'"))

    return commands.getstatusoutput(do_job_cmd)


def getJobStatus(projectName, job):
    """
    Return the job status
    """
    get_job_status_cmd = (getCryosparcProgram() +
                          ' %sget_job_status("%s", "%s")%s'
                          % ("'", projectName, job, "'"))

    return commands.getstatusoutput(get_job_status_cmd)[-1].split()[-1]


def waitJob(projectName, job):
    """
    Wait while the job not finished
    """
    wait_job_cmd = (getCryosparcProgram() +
                    ' %swait_job_complete("%s", "%s")%s'
                    % ("'", projectName, job, "'"))
    commands.getstatusoutput(wait_job_cmd)

