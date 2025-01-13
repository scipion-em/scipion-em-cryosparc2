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
        system_info = getSystemInfo()
        status_errors = system_info[0]
        print(status_errors)
        if not status_errors:
            system_info = eval(system_info[1])
            master_hostname = system_info.get('master_hostname')
            port_webapp = system_info.get('port_webapp')
            port_app = system_info.get('port_app')
            version = system_info.get('version')
            projectId = self.protocol.projectName.get()
            workspaceId = self.protocol.workSpaceName.get()
            jobId = self.protocol.currenJob.get()
            if jobId is None:
                jobId = ''
            if parse_version(version) >= parse_version(V4_1_0):
                url = "http://%s:%s/browse/%s-%s-J*#job(%s-%s)" % (master_hostname, port_app,
                                                                   projectId, workspaceId, projectId, jobId)
            else:
                url = "http://%s:%s/projects/%s/%s/%s" % (master_hostname,
                                                          port_webapp, projectId,
                                                          workspaceId, jobId)
            browser = webbrowser.get()

            browser.open(url)

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
        if self.resolutionPlotsFSC.get() == len(self.choices)-1:  # Case of all plot
            fscViewer.visualize(self.setOfFSCs)
        else:
            pos = 0
            for fsc in self.setOfFSCs.iterItems():
                if pos == self.resolutionPlotsFSC.get():
                    fscViewer.visualize(fsc)
                    break
                pos += 1
