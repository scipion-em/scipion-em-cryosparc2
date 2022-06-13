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
from pyworkflow.gui import TreeProvider
from pyworkflow.gui.dialog import ListDialog, showInfo
from pyworkflow.object import String
from .protocols import *
import pyworkflow.viewer as pwviewer

# Suggested number of images per class
from pyworkflow.wizard import Wizard
from .utils import getSchedulerLanes, cryosparcValidate

IMAGES_PER_CLASS = 200


# =============================================================================
# NUMBER OF CLASSES
# =============================================================================
class ProtCryo2DNumberOfClassesWizard(Wizard):
    _targets = [(ProtCryo2D, ['numberOfClasses'])]

    def _getNumberOfClasses(self, protocol):

        numberOfClasses = 64

        if protocol.inputParticles.hasValue():
            if protocol.inputParticles.get().getSize() > IMAGES_PER_CLASS:
                numberOfClasses = int(protocol.inputParticles.get().getSize()/IMAGES_PER_CLASS)

        return numberOfClasses

    def show(self, form, *args):
        form.setVar('numberOfClasses', self._getNumberOfClasses(form.protocol))


class ProtCryosparcLanesWizard(Wizard):
    _targets = [(ProtCryo2D, ['compute_lane']),
                (ProtCryoSparcRefine3D, ['compute_lane']),
                (ProtCryoSparcInitialModel, ['compute_lane']),
                (ProtCryoSparcNonUniformRefine3D, ['compute_lane']),
                (ProtCryoSparcSubtract, ['compute_lane']),
                (ProtCryoSparcNaiveLocalRefine, ['compute_lane']),
                (ProtCryoSparcLocalRefine, ['compute_lane']),
                (ProtCryoSparcGlobalCtfRefinement, ['compute_lane']),
                (ProtCryoSparcLocalCtfRefinement, ['compute_lane']),
                (ProtCryoSparcSharppening, ['compute_lane']),
                (ProtCryoSparc3DClassification, ['compute_lane']),
                (ProtCryoSparcHelicalRefine3D, ['compute_lane']),
                (ProtCryoSparc3DHomogeneousRefine, ['compute_lane']),
                (ProtCryoSparcNewNonUniformRefine3D, ['compute_lane']),
                (ProtCryoSparcSymmetryExpansion, ['compute_lane']),
                (ProtCryoSparcHomogeneousReconstruct, ['compute_lane']),
                (ProtCryoSparcNew3DClassification, ['compute_lane'])]

    def show(self, form, *args):
        protocol = form.protocol
        csValidate = cryosparcValidate()
        if not csValidate:
            d = LanesDialogView(form.root, protocol)
            dlg = d.show()
            if dlg.resultYes():
                form.setVar('compute_lane', str(dlg.values[0]))
        else:
            showInfo('Info', csValidate[0], form.root)


class LanesTreeProvider(TreeProvider):
    """ Model class that will retrieve the information from cryoSPARC lanes and
    prepare the columns/rows models required by the TreeDialog GUI.
    """
    COL_LANE = 'Lane name'

    def __init__(self, protocol):
        self.protocol = protocol
        TreeProvider.__init__(self)
        self.selectedDict = {}
        self.lanes = self._getComputeLanes()[0]

    def _getComputeLanes(self):
        return getSchedulerLanes()

    def getObjects(self):
        objects = []

        for lane in self.lanes:
            objects.append(String(lane))

        return objects

    def getColumns(self):
        cols = [(self.COL_LANE, 100)]
        return cols

    def getObjectInfo(self, obj):
        objId = obj.getObjValue()
        key = objId
        text = key

        return {
            'key': key,
            'text': text,
            'open': False,
            'value': [text],
            'selected': False,
            'parent': None,
        }


class LanesDialogView(pwviewer.View):
    """ This class implements a view using Tkinter ListDialog
    and the LanesTreeProvider.
    """
    def __init__(self, parent, protocol, **kwargs):
        self._tkParent = parent
        self._protocol = protocol
        self._provider = LanesTreeProvider(protocol)
        self.selectedLane = None

    def show(self):
        return ListDialog(self._tkParent, 'Lanes display', self._provider,
                          allowSelect=True, cancelButton=True,
                          selectOnDoubleClick=True)
