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
from pwem.viewers import (ChimeraView, ObjectView, EmProtocolViewer, FscViewer)

from ..protocols import (ProtCryoSparcLocalRefine, ProtCryoSparcHelicalRefine3D,
                         ProtCryoSparc3DHomogeneousRefine,
                         ProtCryoSparcNewNonUniformRefine3D,
                         ProtCryoSparcHomogeneousReconstruct, ProtCryoSparc3DVariability,
                         ProtCryoSparc3DVariabilityDisplay, ProtCryoSparc3DFlexDataPrepare,
                         ProtCryoSparc3DFlexTraining, ProtCryoSparc3DFlexReconstruction, ProtCryoSparc3DFlexMeshPrepare)
from ..constants import *
from ..utils import *


class CryosPARCViewer3DRefinement(EmProtocolViewer):
    """ Visualization of e2refine_easy results. """

    _targets = [ProtCryoSparcLocalRefine, ProtCryoSparcHelicalRefine3D,
                ProtCryoSparc3DHomogeneousRefine, ProtCryoSparcNewNonUniformRefine3D,
                ProtCryoSparcHomogeneousReconstruct, ProtCryoSparc3DVariability,
                ProtCryoSparc3DVariabilityDisplay, ProtCryoSparc3DFlexDataPrepare,
                ProtCryoSparc3DFlexTraining, ProtCryoSparc3DFlexReconstruction]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _label = 'viewer Refinement/Flex'

    def _defineParams(self, form):
        self._env = os.environ.copy()
        form.addSection(label='Results')

        group = form.addGroup('cryoSPARC')

        group.addParam('displayCS', LabelParam,
                       label='Display the processing with cryoSPARC')

        isFlexProtocol = (isinstance(self.protocol, ProtCryoSparc3DFlexDataPrepare) or
                          isinstance(self.protocol, ProtCryoSparc3DFlexMeshPrepare) or
                          isinstance(self.protocol, ProtCryoSparc3DFlexTraining))

        if self.protocol.isFinished() and not isFlexProtocol:

            if (not isinstance(self.protocol, ProtCryoSparc3DVariabilityDisplay) or
                    not isinstance(self.protocol, ProtCryoSparc3DFlexReconstruction)):
                group = form.addGroup('Particles')

                group.addParam('showImagesAngularAssignment', LabelParam,
                               label='Particles angular assignment')

            is3dVariavilityProt = isinstance(self.protocol, ProtCryoSparc3DVariabilityDisplay)
            groupName = 'Volume' if not is3dVariavilityProt else 'Components'
            if is3dVariavilityProt and self.protocol.var_output_mode.get() == 0:
                groupName = 'Volume'

            group = form.addGroup(groupName)
            if not isinstance(self.protocol, ProtCryoSparc3DVariabilityDisplay):
                displayChoices = ['chimera', 'data viewer']
            else:
                if self.protocol.var_output_mode.get() != 0:
                    self.componetChoices = self.getComponetChoices()
                    group.addParam('component', EnumParam,
                                   choices=list(self.componetChoices),
                                   default=0, display=EnumParam.DISPLAY_COMBO,
                                   label='Select a component',
                                   help='Select a component to display the conformation')

                    group.addParam('maxFrameRate', IntParam,
                                   default=4,
                                   label='Playback rate',
                                   help='Specify a maximum playback rate in steps per second. By default, playback is as fast as possible, which can be fairly slow for large data. This option is used to slow playback when it is too fast.'
                                   )

                displayChoices = ['chimera']
            label = 'Display volume with' if not isinstance(self.protocol, ProtCryoSparc3DVariabilityDisplay) else 'Display component with'
            if is3dVariavilityProt and self.protocol.var_output_mode.get() == 0:
                label = 'Display volume with'
                displayChoices = ['data viewer']

            help = '*data viewer: display volumes as surface with Scipion data viewer. \n ' \
                   '*chimera*: display volumes as surface with Chimera.'

            group.addParam('displayVol', EnumParam, choices=displayChoices,
                           default=VOLUME_CHIMERA, display=EnumParam.DISPLAY_LIST,
                           label=label,
                           help=help)
            # '*slices*: display volumes as 2D slices along z axis.\n'

            if (self.protocol.isFinished() and parse_version(getCryosparcVersion()) != parse_version(V4_0_0) and
                    not isinstance(self.protocol, ProtCryoSparc3DVariabilityDisplay) and
                    not isinstance(self.protocol, ProtCryoSparc3DVariability)):
                group = form.addGroup('Resolution')

                self.choices = self.getChoices()

                group.addParam('resolutionPlotsFSC', EnumParam,
                               choices=list(self.choices),
                               default=0, display=EnumParam.DISPLAY_COMBO,
                               label='Display resolution plots (FSC)',
                               help=self._getFscChoiceHelp())
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

        if getattr(self.protocol, paramName, None) is not None:
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
            return self._showOutputVolume(paramName='outputVolume')

    def _showCryoSPARVolume(self, paramName=None):
        views = []
        url = getCryosparcJobUrl(self.protocol.projectName.get(), self.protocol.workSpaceName.get(),
                                 self.protocol.currenJob.get())
        if url:
            webbrowser.open(url)
        return views

    def _showVolumesChimera(self):
        """ Create a chimera script to visualize selected volumes. """

        # Check if Chimera is installed
        view = []
        chimera = Domain.importFromPlugin('chimera')
        if chimera is not None:
            if not isinstance(self.protocol, ProtCryoSparc3DVariabilityDisplay):
                volumes = [self.protocol.outputVolume.getFileName()]
                if len(volumes) > 1:
                    cmdFile = self.protocol._getExtraPath('chimera_volumes.cxc')
                    f = open(cmdFile, 'w+')
                    for vol in volumes:
                        # We assume that the chimera script will be generated
                        # at the same folder as eman volumes
                        if os.path.exists(vol):
                            localVol = os.path.relpath(vol,
                                                       self.protocol._getExtraPath())
                            f.write("open %s\n" % localVol)
                    f.write('tile\n')
                    f.close()
                    view.append(ChimeraView(cmdFile))
                else:
                    view.append(ChimeraView(volumes[0]))
            elif self.protocol.var_output_mode.get() != 0:
                component = self.componetChoices[self.component.get()]
                componentsPath = os.path.abspath(self.protocol._getExtraPath(component))
                cmdFile = self.protocol._getExtraPath('chimera_volumes.cxc')
                f = open(cmdFile, 'w+')
                f.write("open %s/*.mrc vseries true\n" % componentsPath)
                f.write("vseries play #1 loop true maxFrameRate %d\n" % self.maxFrameRate.get())
                f.close()
                view.append(ChimeraView(cmdFile))
            else:
                volume = self.protocol.outputVolumes
                fn = volume.getFileName()
                v = self.createScipionPartView(fn, volume)
                view.append(v)
        else:
            showInfo('Info', "Chimera plugin is not installed. Please, "
                             "install it to display the volume",
                     self.getTkRoot())

        return view

    def _normalizeFscViewerChoice(self, label):
        text = str(label or "").strip()
        key = text.lower()

        if "auto-tightened" in key and "correct" in key:
            return "Auto-tightened corrected"

        if "auto-tightened" in key and "resolution" in key:
            return "Auto-tightened resolution mask"

        if "resolution mask" in key:
            return "Resolution mask"

        if "input mask" in key and "correct" in key:
            return "Input mask corrected"

        if "input mask" in key:
            return "Input mask"

        if "correct" in key:
            return "Corrected"

        if "tight" in key:
            return "Tight"

        if "loose" in key:
            return "Loose"

        if "spherical" in key:
            return "Spherical"

        if "no mask" in key or "nomask" in key:
            return "No mask"

        return text

    def _getFscChoiceHelp(self):
        version = parse_version(getCryosparcVersion())

        if version >= parse_version(V5_0_0):
            return (
                "*No mask*: FSC without masking.\n"
                "*Spherical*: FSC with a soft spherical mask.\n"
                "*Resolution mask*: FSC using the v5 resolution mask.\n"
                "*Auto-tightened resolution mask*: FSC after v5 auto-tightening.\n"
                "*Auto-tightened corrected*: noise-substitution corrected FSC after auto-tightening.\n"
                "*Input mask*: FSC computed with the user-provided input mask when applicable.\n"
                "*Input mask corrected*: corrected FSC using the input mask when applicable."
            )

        return (
            "*No mask*: display FSC of unmasked maps.\n"
            "*Spherical*: display FSC of spherical masked maps.\n"
            "*Loose*: display FSC of loose masked maps.\n"
            "*Tight*: display FSC of tight masked maps.\n"
            "*Corrected*: display FSC corrected by noise substitution."
        )

    def getChoices(self):
        choices = []
        seen = set()

        output = self.protocol.outputFSC
        if isinstance(output, SetOfFSCs):
            self.setOfFSCs = self.protocol.outputFSC
        else:
            fscFile = "fsc.txt"
            fscFilePath = os.path.join(self.protocol._getExtraPath(), fscFile)
            inputParticles = self.protocol._getInputParticles()
            factor = inputParticles.getDim()[0] * inputParticles.getSamplingRate()
            self.setOfFSCs = self.protocol.getSetOfFCSsFromFile(fscFilePath, factor)
            self.protocol.deleteOutput(output)
            self.protocol._defineOutputs(outputFSC=self.setOfFSCs)

        for fsc in self.setOfFSCs.iterItems():
            label = self._normalizeFscViewerChoice(fsc.getObjLabel())
            if label not in seen:
                seen.add(label)
                choices.append(label)

        choices.append('All')
        return choices

    def getComponetChoices(self):
        if self.protocol.var_output_mode.get() == 2 and self.protocol.var_intermediate_output_frame_particles.get():
            return self.protocol._outputs[0:-1]
        else:
            return self.protocol._outputs


    # =========================================================================
    # plotFSC
    # =========================================================================
    def _showFSC(self, paramName=None):
        fscViewer = FscViewer(project=self.getProject(),
                              protocol=self.protocol)

        selectedIndex = self.resolutionPlotsFSC.get()
        if selectedIndex == len(self.choices) - 1:
            fscViewer.visualize(self.setOfFSCs)
            return

        selectedLabel = self.choices[selectedIndex]

        for fsc in self.setOfFSCs.iterItems():
            currentLabel = self._normalizeFscViewerChoice(fsc.getObjLabel())
            if currentLabel == selectedLabel:
                fscViewer.visualize(fsc)
                break
