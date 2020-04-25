
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
from pwem.viewers import ObjectView

from ..protocols import (ProtCryoSparcSubtract,
                         ProtCryoSparcGlobalCtfRefinement,
                         ProtCryoSparcLocalCtfRefinement)
from ..constants import *
from ..utils import *
from .. import Plugin


class CryosPARCViewerPartSubtract(ProtocolViewer):
    """
    Visualization tools for cryosPARC results.

    CryoSPARC is a backend and frontend software system that provides data
    processing and image analysis capabilities for single particle cryo-EM,
    along with a browser based user interface and command line tools.

    Please find the documentation at https://cryosparc.com
    """

    _environments = [DESKTOP_TKINTER]
    _targets = [ProtCryoSparcSubtract, ProtCryoSparcGlobalCtfRefinement,
                ProtCryoSparcLocalCtfRefinement]
    _label = 'viewer cryosPARC'

    def __init__(self, *args, **kwargs):
        ProtocolViewer.__init__(self, **kwargs)
        sys.path.append(Plugin.getVar(CRYOSPARC_HOME))

    def _defineParams(self, form):
        self._env = os.environ.copy()
        form.addSection(label='Visualization')
        group = form.addGroup('Particles')
        group.addParam('displayDataViewer', LabelParam,
                       label='Display particle classes with Scipion')
        group.addParam('displayCryosPARC2D', LabelParam,
                       label='Display particle classes with cryosPARC GUI')

    def _getVisualizeDict(self):
        self._load()
        visualizeDict = {'displayDataViewer': self._showOutputParticles,
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

    def _showOutputParticles(self, paramName=None):
        views = []
        if getattr(self.protocol, 'outputParticles', None) is not None:
            fn = self.protocol.outputParticles.getFileName()
            v = self.createScipionPartView(fn)
            views.append(v)
        return views

    def createScipionPartView(self, filename, viewParams={}):
        inputParticlesId = self.protocol.inputParticles.get().strId()
        labels = 'enabled id _size _filename _transform._matrix'
        viewParams = {showj.ORDER: labels,
                      showj.VISIBLE: labels, showj.RENDER: '_filename',
                      'labels': 'id',
                      }
        return ObjectView(self._project,
                          self.protocol.strId(), filename,
                          other=inputParticlesId,
                          env=self._env, viewParams=viewParams)

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
        self.protocol._createFilenameTemplates()


