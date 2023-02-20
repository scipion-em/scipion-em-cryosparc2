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

from pyworkflow import BETA
from pyworkflow.protocol.params import (PointerParam, FloatParam, Positive,
                                        IntParam, BooleanParam, EnumParam)
from pwem.objects import Volume

from . import ProtCryosparcBase
from ..convert import *
from ..utils import *
from ..constants import *


class ProtCryoSparc3DFlexDataPrepare(ProtCryosparcBase):
    """
    Prepares particles for use in 3DFlex training and reconstruction
    """
    _label = '3D flex data prepare'
    _devStatus = BETA
    _className = "flex_prep"

    # --------------------------- DEFINE param functions ----------------------
    def _defineFileNames(self):
        """ Centralize how files are called within the protocol. """
        myDict = {
            'input_particles': self._getTmpPath('input_particles.star'),
            'out_particles': self._getPath() + '/output_particle.star',
            'stream_log': self._getPath() + '/stream.log'
        }
        self._updateFilenamesDict(myDict)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      label="Input particles",
                      important=True,
                      help='Particle stacks to use.')
        form.addParam('refVolume', PointerParam, pointerClass='Volume',
                      label='Initial volume',
                      important=True,
                      help='Initial volume raw data')

        # Data Preparation
        form.addParam('box_size_pix', IntParam, default=None,
                      allowsNull=True,
                      label="Crop box size (pix)",
                      help="Crop input particles and volume to this box size. "
                           "This is the box size that will be used for high "
                           "resolution reconstruction with 3D Flex. Particles "
                           "are cropped to this box size first, and then "
                           "downsampled to the Training Box Size for training "
                           "time. Default (None) means to keep the original "
                           "box size of the particles.")

        form.addParam('bin_size_pix', IntParam, default=128,
                      allowsNull=True,
                      label="Training box size (pix)",
                      help="Downsample cropped particles (via Fourier cropping)"
                           " to this box size for training the 3D Flex model. "
                           "This should be chosen to limit 3D Flex training to "
                           "a resolution below the gold-standard FSC resolution "
                           "of the consensus reconstruction, in order to ensure "
                           "high resolution reconstructions can be validated. "
                           "Box sizes over 256 may become prohibitively slow")

        form.addParam('alpha_min', IntParam, default=None,
                      allowsNull=True,
                      label="Min. scale to keep",
                      help="Only keep particles with scale factor above this "
                           "value. Useful for discarding particles that might "
                           "be junk.")

        form.addParam('keep_num_particles', IntParam, default=None,
                      allowsNull=True,
                      label="Num. particles to use",
                      help="Only keep the first X particles. The final number"
                           " of particles used during 3D Flex training "
                           "and reconstruction must be divisible by 1000. "
                           "If this is None (default) then the number of input "
                           "particles will be rounded down to the nearest 1000.")

        # --------------[Compute settings]---------------------------
        form.addSection(label="Compute settings")
        addComputeSectionParams(form, allowMultipleGPUs=False, needGPU=False)

    def _insertAllSteps(self):
        self._defineFileNames()
        self._defineParamsName()
        self._initializeCryosparcProject()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.processStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions ------------------------------
    def processStep(self):
        print(pwutils.yellowStr("3D Flex Data Preparation started..."))
        self.doRun3DFlexDataPrepare()

    def createOutputStep(self):
        pass

    # ------------------------- Utils methods ----------------------------------
    def _validate(self):
        validateMsgs = cryosparcValidate()
        if not validateMsgs:
            particles = self._getInputParticles()
            if not particles.hasCTF():
                validateMsgs.append(
                    "The Particles has not associated a CTF model")
                if not particles.hasAlignment3D():
                    validateMsgs.append("The Particles has not aligned")
        return validateMsgs

    def _defineParamsName(self):
        """ Define a list with all protocol parameters names"""
        self._paramsName = ['box_size_pix', 'bin_size_pix', 'alpha_min',
                            'keep_num_particles']
        self.lane = str(self.getAttributeValue('compute_lane'))

    def doRun3DFlexDataPrepare(self):
        input_group_connect = {"particles": self.particles.get(),
                               "volume": self.volume.get()}
        params = {}

        for paramName in self._paramsName:
            if self.getAttributeValue(paramName) is not None:
                params[str(paramName)] = str(self.getAttributeValue(paramName))

        run3DFlexDataPrepJob = enqueueJob(self._className, self.projectName.get(),
                                  self.workSpaceName.get(),
                                  str(params).replace('\'', '"'),
                                  str(input_group_connect).replace('\'', '"'),
                                  self.lane, False)

        self.run3DFlexDataPrepJob = String(run3DFlexDataPrepJob.get())
        self.currenJob.set(self.run3DFlexDataPrepJob.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(), self.run3DFlexDataPrepJob.get(),
                         "An error occurred in the 3D Flex Data Preparation process. "
                         "Please, go to cryoSPARC software for more "
                         "details.")
        clearIntermediateResults(self.projectName.get(), self.run3DFlexDataPrepJob.get())




