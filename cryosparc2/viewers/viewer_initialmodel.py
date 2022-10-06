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

import webbrowser

from pyworkflow.viewer import (ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO)
import pwem.viewers.showj as showj
from pwem.viewers import ObjectView
from pyworkflow.protocol.params import (LabelParam, EnumParam)

from ..constants import *
from ..utils import *
from ..protocols import ProtCryoSparcInitialModel, ProtCryoSparc3DClassification


class CryosPARCViewerInitialModel(ProtocolViewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj
    """
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [ProtCryoSparcInitialModel, ProtCryoSparc3DClassification]

    def _defineParams(self, form):
        self._env = os.environ.copy()
        form.addSection(label='Results')
        group = form.addGroup('3D Classes')

        group.addParam('show3DClasses', LabelParam,
                       label='Initial model output classes')

        group = form.addGroup('Volume')

        group.addParam('displayVol', EnumParam, choices=['dataViewer', 'cryoSPARC'],
                       default=DATA_VIEWER, display=EnumParam.DISPLAY_HLIST,
                       label='Display volume with',
                       help='*dataViewer*: display volumes as surface with Scipion.\n'
                            '*cryoSPARC: display volumes as surface with cryoSPARC')
        # '*slices*: display volumes as 2D slices along z axis.\n'

    def _getVisualizeDict(self):
        self._load()
        return {'show3DClasses': self._visualizeClasses,
                'displayVol': self._visualizeVolumes,
                }

    def _visualizeClasses(self, paramName=None):
        view = []

        obj = self.protocol.outputClasses
        fn = obj.getFileName()
        labels = 'enabled id _size _representative._filename '

        view.append(ObjectView(self._project, obj.strId(), fn,
                               viewParams={showj.MODE: showj.MODE_MD,
                                           showj.VISIBLE: labels,
                                           showj.RENDER: '_representative._filename'}))
        return view

    def _visualizeVolumes(self, paramName=None):
        if self.displayVol == DATA_VIEWER:
            return self._showVolumesDataViewer()
        elif self.displayVol == VOLUME_CRYOSPARC:
            return self._showCryoSPARVolume()

    def _showCryoSPARVolume(self, paramName=None):
        views = []
        system_info = getSystemInfo()
        status_errors = system_info[0]
        print(status_errors)
        if not status_errors:
            system_info = eval(system_info[1])
            master_hostname = system_info.get('master_hostname')
            port_webapp = system_info.get('port_webapp')

            projectId = self.protocol.projectName.get()
            workspaceId = self.protocol.workSpaceName.get()
            jobId = self.protocol.currenJob.get()
            url = os.path.join("http://",
                               master_hostname + ':' + port_webapp,
                               "projects", projectId, workspaceId,
                               jobId)

            browser = webbrowser.get()

            browser.open(url)

        return views

    def _showVolumesDataViewer(self):
        """ Visualize the volumes with the Scipion dataViewer. """
        view = []

        obj = self.protocol.outputVolumes
        fn = obj.getFileName()
        labels = 'id enabled comment _filename '
        objCommands = "'%s' '%s' '%s'" % (OBJCMD_CLASSAVG_PROJS,
                                          OBJCMD_PROJS,
                                          OBJCMD_INITVOL)

        view.append(ObjectView(self._project, obj.strId(), fn,
                               viewParams={showj.MODE: showj.MODE_MD,
                                           showj.VISIBLE: labels,
                                           showj.RENDER: '_filename',
                                           showj.OBJCMDS: objCommands}))
        return view

    def _load(self):
        """ Load the 3D classes and volumes for visualization mode. """
        self.protocol._defineFileNames()
