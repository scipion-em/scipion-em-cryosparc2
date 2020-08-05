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
from pyworkflow.protocol.constants import *
from pyworkflow.protocol.params import (LabelParam, FloatParam)
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO
from pwem.viewers import (ChimeraView, EmPlotter, ChimeraClientView,
                          ObjectView, EmProtocolViewer)

from ..protocols import (ProtCryoSparcNonUniformRefine3D,
                         ProtCryoSparcRefine3D,
                         ProtCryoSparcLocalRefine)
from ..constants import *
from ..utils import *


class CryosPARCViewer3DRefinement(EmProtocolViewer):
    """ Visualization of e2refine_easy results. """

    _targets = [ProtCryoSparcRefine3D, ProtCryoSparcNonUniformRefine3D,
                ProtCryoSparcLocalRefine]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _label = 'viewer Refinement'

    def _defineParams(self, form):
        self._env = os.environ.copy()
        form.addSection(label='Results')
        group = form.addGroup('Particles')

        group.addParam('showImagesAngularAssignment', LabelParam,
                       label='Particles angular assignment')

        group = form.addGroup('Volume')

        group.addParam('displayVol', EnumParam, choices=['chimera', 'cryoSPARC'],
                       default=VOLUME_CHIMERA, display=EnumParam.DISPLAY_LIST,
                       label='Display volume with',
                       help='*chimera*: display volumes as surface with Chimera.\n'
                            '*cryoSPARC: display volumes as surface with cryoSPARC')
        # '*slices*: display volumes as 2D slices along z axis.\n'

        group = form.addGroup('Resolution')

        group.addParam('resolutionPlotsFSC', EnumParam,
                       choices=['no mask', 'spherical', 'loose',
                                'tight', 'corrected', 'all'],
                       default=FSC_UNMASK, display=EnumParam.DISPLAY_COMBO,
                       label='Display resolution plots (FSC)',
                       help='*unmasked*: display FSC of unmasked maps.\n'
                            '*masked*: display FSC of masked maps.\n'
                            '*masked tight*: display FSC of masked tight maps.')
        group.addParam('resolutionThresholdFSC', FloatParam, default=0.143,
                       expertLevel=LEVEL_ADVANCED,
                       label='Threshold in resolution plots',
                       help='')

    def _getVisualizeDict(self):
        self._load()
        return {'showImagesAngularAssignment': self._showOutputParticles,
                'displayVol': self._showVolumes,
                'resolutionPlotsFSC': self._showFSC
                }

    # =========================================================================
    # showImagesAngularAssignment
    # =========================================================================

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

    # =========================================================================
    # ShowVolumes
    # =========================================================================
    def _showVolumes(self, paramName=None):
        if self.displayVol == VOLUME_CHIMERA:
            return self._showVolumesChimera()
        # elif self.displayVol == VOLUME_SLICES:
        #     return self._createVolumesSqlite()
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

    def _showVolumesChimera(self):
        """ Create a chimera script to visualize selected volumes. """
        volumes = self._getVolumeNames()

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
            view = ChimeraView(cmdFile)
        else:
            view = ChimeraClientView(volumes[0])

        return [view]

    def _createVolumesSqlite(self):
        """ Write an sqlite with the volumes selected for visualization. """

        path = self.protocol._getExtraPath('cryosparc_viewer_volume.sqlite')
        samplingRate = self.protocol.inputParticles.get().getSamplingRate()

        files = []
        volumes = self._getVolumeNames()
        for volFn in volumes:
            if os.path.exists(volFn.replace(':mrc', '')):
                files.append(volFn)
        self.createVolumesSqlite(files, path, samplingRate)
        return [ObjectView(self._project, self.protocol.strId(), path)]

    # =========================================================================
    # plotFSC
    # =========================================================================
    def _showFSC(self, paramName=None):
        threshold = self.resolutionThresholdFSC.get()
        gridsize = self._getGridSize(1)
        xplotter = EmPlotter(x=gridsize[0], y=gridsize[1],
                             windowTitle='Resolution FSC')

        plot_title = 'FSC'
        a = xplotter.createSubPlot(plot_title, 'Angstroms^-1', 'FSC',
                                   yformat=False)

        legends = []
        show = False
        fsc_path = self.protocol._getExtraPath('fsc.txt')
        if os.path.exists(fsc_path):
            show = True
            if self.resolutionPlotsFSC.get() == FSC_UNMASK:
                self._plotFSC(a, fsc_path, FSC_UNMASK)
                legends.append('No Mask')
                xplotter.showLegend(legends)
            elif self.resolutionPlotsFSC.get() == FSC_SPHERICALMASK:
                self._plotFSC(a, fsc_path, FSC_SPHERICALMASK)
                legends.append('Spherical')
                xplotter.showLegend(legends)
            elif self.resolutionPlotsFSC.get() == FSC_TIGHTMASK:
                self._plotFSC(a, fsc_path, FSC_TIGHTMASK)
                legends.append('Tight')
                xplotter.showLegend(legends)
            elif self.resolutionPlotsFSC.get() == FSC_LOOSEMASK:
                self._plotFSC(a, fsc_path, FSC_LOOSEMASK)
                legends.append('Loose')
                xplotter.showLegend(legends)
            elif self.resolutionPlotsFSC.get() == FSC_CORRECTEDMASK:
                self._plotFSC(a, fsc_path, FSC_CORRECTEDMASK+1)
                legends.append('Corrected')
                xplotter.showLegend(legends)
            elif self.resolutionPlotsFSC.get() == FSC_ALL:
                self._plotFSC(a, fsc_path, FSC_UNMASK)
                legends.append('No Mask')
                self._plotFSC(a, fsc_path, FSC_SPHERICALMASK)
                legends.append('Spherical')
                self._plotFSC(a, fsc_path, FSC_LOOSEMASK)
                legends.append('Loose')
                self._plotFSC(a, fsc_path, FSC_TIGHTMASK)
                legends.append('Tight')
                self._plotFSC(a, fsc_path, FSC_CORRECTEDMASK+1)
                legends.append('Corrected')
                xplotter.showLegend(legends)

        if show:
            if threshold < self.maxFrc:
                a.plot([self.minInv, self.maxInv], [threshold, threshold],
                       color='black', linestyle='--')
            a.grid(True)
        else:
            raise Exception("Set a valid iteration to show its FSC")

        return [xplotter]

    def _plotFSC(self, a, fscFn, col):
        resolution_inv = self._getColunmFromFilePar(fscFn, 0)
        frc = self._getColunmFromFilePar(fscFn, col+1)
        self.maxFrc = max(frc)
        self.minInv = min(resolution_inv)
        self.maxInv = max(resolution_inv)
        self.sampligRate = self.protocol._getInputParticles().getSamplingRate()
        factor = (2. * self.sampligRate * self.maxInv)
        resolution_inv = [x/factor for x in resolution_inv]
        self.minInv /= factor
        self.maxInv /= factor
        a.plot(resolution_inv, frc)
        a.xaxis.set_major_formatter(self._plotFormatter)
        a.set_ylim([-0.1, 1.1])

    # =========================================================================
    # Utils Functions
    # =========================================================================
    def _load(self):
        """ Load the 3D classes for visualization mode. """
        self.protocol._defineFileNames()
        from matplotlib.ticker import FuncFormatter
        self._plotFormatter = FuncFormatter(self._formatFreq)

    @staticmethod
    def _formatFreq(value, pos):
        """ Format function for Matplotlib formatter. """
        inv = 999.
        if value:
            inv = 1./value
        return "1/%0.2f" % inv

    def _getVolumeNames(self):
        vol = []
        vn = self.protocol.outputVolume.getFileName()
        vol.append(vn)
        return vol

    def _getGridSize(self, n=None):
        """ Figure out the layout of the plots given the number of
        references. """
        if n is None:
            n = len(self._refsList)

        if n == 1:
            gridsize = [1, 1]
        elif n == 2:
            gridsize = [2, 1]
        else:
            gridsize = [(n + 1) / 2, 2]

        return gridsize

    def _getColunmFromFilePar(self, fscFn, col):
        f1 = open(fscFn)
        f1.readline()
        value = []
        for l in f1:
            valList = l.split()
            val = float(valList[col])
            value.append(val)
        f1.close()
        return value
