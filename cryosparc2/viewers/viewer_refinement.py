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

import pwem.viewers.showj as showj
from pwem import Domain
from pwem.objects import FSC, SetOfFSCs
from pyworkflow.gui.dialog import showInfo
from pyworkflow.protocol.constants import *
from pyworkflow.protocol.params import (LabelParam, FloatParam, EnumParam)
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO
from pwem.viewers import (ChimeraView, ChimeraClientView,
                          ObjectView, EmProtocolViewer, FscViewer)

from ..protocols import (ProtCryoSparcNonUniformRefine3D,
                         ProtCryoSparcRefine3D,
                         ProtCryoSparcLocalRefine, ProtCryoSparcHelicalRefine3D,
                         ProtCryoSparc3DHomogeneousRefine,
                         ProtCryoSparcNewNonUniformRefine3D,
                         ProtCryoSparcNaiveLocalRefine,
                         ProtCryoSparcHomogeneousReconstruct)
from ..constants import *
from ..utils import *


class CryosPARCViewer3DRefinement(EmProtocolViewer):
    """ Visualization of e2refine_easy results. """

    _targets = [ProtCryoSparcRefine3D, ProtCryoSparcNonUniformRefine3D,
                ProtCryoSparcLocalRefine, ProtCryoSparcHelicalRefine3D,
                ProtCryoSparc3DHomogeneousRefine, ProtCryoSparcNaiveLocalRefine,
                ProtCryoSparcNewNonUniformRefine3D, ProtCryoSparcHomogeneousReconstruct]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _label = 'viewer Refinement'

    def _defineParams(self, form):
        self._env = os.environ.copy()
        form.addSection(label='Results')

        group = form.addGroup('cryoSPARC')

        group.addParam('displayCS', LabelParam,
                       label='Display the processing with cryoSPARC')

        if self.protocol.isFinished():
            group = form.addGroup('Particles')

            group.addParam('showImagesAngularAssignment', LabelParam,
                           label='Particles angular assignment')

            group = form.addGroup('Volume')

            displayChoices = ['data viewer', 'chimera']
            label = 'Display volume with'
            help = '*data viewer: display volumes as surface with Scipion data viewer. \n ' \
                   '*chimera*: display volumes as surface with Chimera.'

            group.addParam('displayVol', EnumParam, choices=displayChoices,
                           default=DATA_VIEWER, display=EnumParam.DISPLAY_LIST,
                           label=label,
                           help=help)
            # '*slices*: display volumes as 2D slices along z axis.\n'

            if self.protocol.isFinished():
                group = form.addGroup('Resolution')

                self.choices = self.getChoices()

                group.addParam('resolutionPlotsFSC', EnumParam,
                               choices=list(self.choices),
                               default=0, display=EnumParam.DISPLAY_COMBO,
                               label='Display resolution plots (FSC)',
                               help='*unmasked*: display FSC of unmasked maps.\n'
                                    '*masked*: display FSC of masked maps.\n'
                                    '*masked tight*: display FSC of masked tight maps.')
                group.addParam('resolutionThresholdFSC', FloatParam, default=0.143,
                               expertLevel=LEVEL_ADVANCED,
                               label='Threshold ',
                               help='Threshold in resolution plots')


    def _getVisualizeDict(self):
        return {'showImagesAngularAssignment': self._showOutputParticles,
                'displayVol': self._showVolumes,
                'resolutionPlotsFSC': self._showFSC,
                'displayCS': self._showCryoSPARVolume
                }

    # =========================================================================
    # showImagesAngularAssignment
    # =========================================================================

    def _showOutputParticles(self, paramName=None):
        views = []

        if getattr(self.protocol, 'outputParticles', None) is not None:
            particles = self.protocol.outputParticles
            fn = particles.getFileName()
            labels = 'enabled id _filename _ctfModel._defocusU _ctfModel._defocusV _ctfModel._defocusAngle _transform._matrix'
            viewParams = {showj.ORDER: labels,
                          showj.VISIBLE: labels, showj.RENDER: '_filename',
                          'labels': 'id',
                          }
            v = self.createScipionPartView(fn, particles,
                                           viewParams=viewParams)
            views.append(v)
        return views

    def _showOutputVolume(self, paramName=None):
        views = []

        if getattr(self.protocol, 'outputVolume', None) is not None:
            volume = self.protocol.outputVolume
            fn = volume.getFileName()
            v = self.createScipionPartView(fn, volume)
            views.append(v)
        return views

    def createScipionPartView(self, filename, obj, viewParams={}):
        objId = obj.strId()
        return ObjectView(self._project,
                          self.protocol.strId(), filename,
                          other=objId,
                          env=self._env, viewParams=viewParams)


    # =========================================================================
    # ShowVolumes
    # =========================================================================
    def _showVolumes(self, paramName=None):
        if self.displayVol == VOLUME_CHIMERA:
            return self._showVolumesChimera()
        elif self.displayVol == DATA_VIEWER:
            return self._showOutputVolume()

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
            if jobId is None:
                jobId = ''
            url = os.path.join("http://",
                               master_hostname + ':' + port_webapp,
                               "projects", projectId, workspaceId,
                               jobId)

            browser = webbrowser.get()

            browser.open(url)

        return views

    def _showVolumesChimera(self):
        """ Create a chimera script to visualize selected volumes. """

        # Check if Chimera is installed
        view = []
        chimera = Domain.importFromPlugin('chimera')
        if chimera is not None:
            volumes = [self.protocol.outputVolume.getFileName()]
            if len(volumes) > 1:
                cmdFile = self.protocol._getExtraPath('chimera_volumes.cxc')
                f = open(cmdFile, 'w+')
                for vol in volumes:
                    # We assume that the chimera script will be generated
                    # at the same folder than eman volumes
                    if os.path.exists(vol):
                        localVol = os.path.relpath(vol,
                                                   self.protocol._getExtraPath())
                        f.write("open %s\n" % localVol)
                f.write('tile\n')
                f.close()
                view.append(ChimeraView(cmdFile))
            else:
                view.append(ChimeraClientView(volumes[0]))
        else:
            showInfo('Info', "Chimera plugin is not installed. Please, "
                             "install it to display the volume",
                     self.getTkRoot())

        return view

    def getChoices(self):
        choices = []
        output = self.protocol.outputFSC
        if isinstance(output, SetOfFSCs):
            self.setOfFSCs = self.protocol.outputFSC
            for fsc in self.setOfFSCs.iterItems():
                choices.append(fsc.getObjLabel())
        else:
            fscFile = "fsc.txt"
            fscFilePath = os.path.join(self.protocol._getExtraPath(), fscFile)
            inputParticles = self.protocol._getInputParticles()
            factor = inputParticles.getDim()[0] * inputParticles.getSamplingRate()
            self.setOfFSCs = self.protocol.getSetOfFCSsFromFile(fscFilePath, factor)
            self.protocol.deleteOutput(output)
            self.protocol._defineOutputs(outputFSC=self.setOfFSCs)
            for fsc in self.setOfFSCs.iterItems():
                choices.append(fsc.getObjLabel())
        choices.append('All')

        return choices

    # =========================================================================
    # plotFSC
    # =========================================================================
    def _showFSC(self, paramName=None):

        fscViewer = FscViewer(project=self.getProject(),
                              protocol=self.protocol)
        if self.resolutionPlotsFSC.get() == len(self.choices)-1:  # Case of all plot
            fscViewer.visualize(self.setOfFSCs)
        else:
            pos = 0
            for fsc in self.setOfFSCs.iterItems():
                if pos == self.resolutionPlotsFSC.get():
                    fscViewer.visualize(fsc)
                    break
                pos += 1
