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
from .protocols import ProtCryo2D


# Suggested number of images per class
from pyworkflow.wizard import Wizard

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
