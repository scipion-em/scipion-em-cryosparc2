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
import sys
import webbrowser

import pwem.viewers.showj as showj
from pyworkflow.protocol.params import LabelParam
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER
from pwem.viewers import ClassesView, Classes3DView

from ..protocols import ProtCryo2D
from ..utils import *
from .. import Plugin

AX_X = 0
AX_Y = 1
AX_Z = 2


class CryosPARCViewer2D(ProtocolViewer):
    """
    Visualization tools for cryosPARC results.

    CryoSPARC is a backend and frontend software system that provides data
    processing and image analysis capabilities for single particle cryo-EM,
    along with a browser based user interface and command line tools.

    Please find the documentation at https://cryosparc.com
    """

    _environments = [DESKTOP_TKINTER]
    _targets = [ProtCryo2D]
    _label = 'viewer cryosPARC'

    def __init__(self, *args, **kwargs):
        ProtocolViewer.__init__(self, **kwargs)
        sys.path.append(Plugin.getVar(CRYOSPARC_HOME))

    def _defineParams(self, form):
        self._env = os.environ.copy()
        form.addSection(label='Visualization')
        group = form.addGroup('Particles')
        group.addParam('displayClass2D', LabelParam,
                       label='Display particle classes with Scipion')
        group.addParam('displayCryosPARC2D', LabelParam,
                       label='Display particle classes with cryosPARC GUI')

    def _getVisualizeDict(self):
        self._load()
        visualizeDict = {'displayClass2D': self._showScipionClasses,
                         'displayCryosPARC2D': self._showCryoSPARCClasses}

        # If the is some error during the load, just show that instead
        # of any viewer
        if self._errors:
            for k in visualizeDict.keys():
                visualizeDict[k] = self._showErrors

        return visualizeDict

    def _showErrors(self, param=None):
        views = []
        self.errorList(self._errors, views)
        return views

    def _getZoom(self):
        # Ensure that classes are shown at least at 128 px to
        # properly see the rlnClassDistribution label.
        dim = self.protocol.inputParticles.get().getDim()[0]
        if dim < 128:
            zoom = int(128*100/dim)
        else:
            zoom = 100
        return zoom

    def _showScipionClasses(self, paramName=None):
        views = []
        if getattr(self.protocol, 'outputClasses', None) is not None:
            fn = self.protocol.outputClasses.getFileName()
            v = self.createScipionView(fn)
            views.append(v)
        return views

    def _showCryoSPARCClasses(self, paramName=None):
        views = []
        system_info = getSystemInfo()
        status_errors = system_info[0]
        if not status_errors:
            system_info = eval(system_info[1])
            master_hostname = system_info.get('master_hostname')
            port_webapp = system_info.get('port_webapp')

            projectId = self.protocol.projectName.get()
            workspaceId = self.protocol.workSpaceName.get()
            jobId = self.protocol.currenJob.get()
            url = os.path.join("http://", master_hostname + ':' + port_webapp,
                               "projects", projectId, workspaceId, jobId)
            browser = webbrowser.get()

            browser.open(url)

        return views

    def _load(self):
        self._errors = []
        self.protocol._defineFileNames()

    def createScipionView(self, filename):
        labels = 'enabled id _size _representative._filename '
        labels += '_rlnclassDistribution _rlnAccuracyRotations _rlnAccuracyTranslations'
        viewParams = {showj.ORDER: labels,
                      showj.VISIBLE: labels,
                      showj.RENDER: '_representative._filename',
                      showj.SORT_BY: '_size desc',
                      showj.ZOOM: str(self._getZoom())
                      }
        inputParticlesId = self.protocol.inputParticles.get().strId()
        ViewClass = ClassesView if self.protocol.IS_2D else Classes3DView
        view = ViewClass(self._project,
                         self.protocol.strId(), filename, other=inputParticlesId,
                         env=self._env,
                         viewParams=viewParams)

        return view
