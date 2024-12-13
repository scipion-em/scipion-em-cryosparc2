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

from ..protocols import ProtCryoSparc3DFlexMeshPrepare
from ..constants import *
from ..utils import *
import glob
from pwem.viewers.viewer_chimera import Chimera

class CryosPARCViewerShowMesh(EmProtocolViewer):
    """ Visualization of cryosparc flex mesh results. """

    _targets = [ProtCryoSparc3DFlexMeshPrepare]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _label = 'viewer Refinement/Flex'

    def _defineParams(self, form):
        self._env = os.environ.copy()
        form.addSection(label='Results')

        group = form.addGroup('cryoSPARC')

        group.addParam('displayCS', LabelParam,
                       label='Display the processing with cryoSPARC')

        group = form.addGroup('Mesh')
        group.addParam('showMesh', LabelParam,
                        label='Show tetrahedral mesh for 3DFlex',
                        help='Show the tetrahedral mesh used for 3DFlex')

    def _getVisualizeDict(self):
        return {'showMesh': self._showVolumesChimera,
                }

    # =========================================================================
    # showImagesAngularAssignment
    # =========================================================================


    def _showVolumesChimera(self, paramName=None):
        """ Create a chimera script to visualize the mesh """
        # get pdb file
        projectDir = self.protocol._getExtraPath()
        pattern = os.path.join(projectDir, '*_mesh_pdb.pdb')
        meshFn = glob.glob(pattern)
        # jobId = self.protocol.run3DFlexMeshPrepJob.get()
        # pdbMeshName = "%s_mesh_pdb.pdb" % jobId
        meshFn = os.path.abspath(meshFn[0])
        # get input volume (from dataPrepare)
        dataPrepare = self.protocol.dataPrepare.get()
        if dataPrepare is not None:
            vol = dataPrepare.refVolume.get()
            volumeFn = os.path.abspath(vol.getFileName())
            dim = vol.getDim()[0]
            sampling = vol.getSamplingRate()
            x, y, z = vol.getOrigin(force=True).getShifts()

        else:
            volumeFn = None
        # get mask
        mask = self.protocol.refMask.get()
        if mask is not None:
            maskFn = os.path.abspath(mask.getFileName())
        else:
            maskFn = None

        # Check if Chimera is installed
        view = []
        chimera = Domain.importFromPlugin('chimera')
        if chimera is not None:
            bildFileName = os.path.abspath(self.protocol._getExtraPath(
                "axis_output.bild"))
            Chimera.createCoordinateAxisFile(dim,
                                 bildFileName=bildFileName,
                                 sampling=sampling)
            cmdFile = self.protocol._getExtraPath('chimera_volumes.cxc')
            f = open(cmdFile, 'w+')
            id = 1
            f.write("open %s\n" % bildFileName)
            f.write("cofr 0,0,0\n")  # set center of coordinates
            if meshFn is not None:
                id += 1
                f.write("open %s\n" % meshFn)
                f.write(f"move {x:0.2f},{y:0.2f},{z:0.2f} models #{id}\n")
            if volumeFn is not None:
                f.write("open %s\n" % volumeFn)
                id += 1
                f.write("volume #%d  style surface voxelSize %f\n"
                        "volume #%d  origin %0.2f,%0.2f,%0.2f\n"
                        % (id, sampling, id, x, y, z))
            if maskFn is not None:
                f.write("open %s\n" % maskFn)
            f.close()
            view.append(ChimeraView(cmdFile))
        else:
            showInfo('Info', "Chimera plugin is not installed. Please, "
                             "install it to display the volume",
                     self.getTkRoot())

        return view

