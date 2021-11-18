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
import os

from pwem.objects import Volume
from pwem.protocols import ProtAnalysis3D
import pyworkflow.utils as pwutils
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        LEVEL_ADVANCED, BooleanParam, IntParam,
                                        String)

from .protocol_base import ProtCryosparcBase
from ..utils import (addComputeSectionParams, cryosparcValidate, gpusValidate,
                     enqueueJob, waitForCryosparc, clearIntermediateResults,
                     fixVolume, copyFiles)


class ProtCryoSparcSharppening(ProtCryosparcBase, ProtAnalysis3D):
    """
    Wrapper protocol for the Cryosparc's to calculate the sharpened map.
    """
    _label = 'sharppening'
    _className = "sharpen"

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('refVolume', PointerParam, pointerClass='Volume',
                      important=True,
                      label="Input volume",
                      help='Provide a reference volume for sharpening')

        # -----------[Sharppening]------------------------
        form.addSection(label="Sharpening")
        form.addParam('sharp_bfactor', FloatParam, default=0.0,
                      label='B-Factor to apply',
                      help='Negative values sharpen.')

        form.addParam('sharp_use_FSC_full', BooleanParam, default=False,
                      label="Use full FSC",
                      help="False means to use the half-FSC, which is usually "
                           "an underestimate for sharpening.")

        form.addParam('sharp_falloff_order', IntParam, default=8,
                      expertLevel=LEVEL_ADVANCED,
                      label='Lowpass filter order',
                      help='Higher means faster falloff, 2 is usually best.')

        form.addParam('sharp_falloff_offset', IntParam, default=0,
                      expertLevel=LEVEL_ADVANCED,
                      label='Lowpass filter offset',
                      help='Offset for corner frequency from FSC resolution shell')

        form.addParam('sharp_generate_new_mask', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Generate new FSC mask",
                      help="Create a new mask for FSC and sharpening rather "
                           "than using the refinement FSC mask")

        form.addParam('sharp_mask_thresh', FloatParam, default=0.5,
                      expertLevel=LEVEL_ADVANCED,
                      label='Threshold',
                      help='Mask generation threshold')

        form.addParam('sharp_mask_near_A', IntParam, default=6,
                      expertLevel=LEVEL_ADVANCED,
                      label='Mask near (A)',
                      help='Mask dilation near (A)')

        form.addParam('sharp_mask_far_A', IntParam, default=12,
                      expertLevel=LEVEL_ADVANCED,
                      label='Mask far (A)',
                      help='Mask dilation far (A)')

        form.addParam('sharp_do_spherical_mask', BooleanParam, default=True,
                      label="Spherical mask final output",
                      help="Apply a spherical mask to the final map, to mask "
                           "out corners with over-sharpened density and to "
                           "reduce file size after compression.")

        form.addParam('sharp_do_expand_mask', BooleanParam, default=False,
                      label="Wide mask final output",
                      help="Also apply a wide mask to the final map, to reduce "
                           "file size after compression.")

        # --------------[Compute settings]---------------------------
        form.addSection(label="Compute settings")
        addComputeSectionParams(form, allowMultipleGPUs=False)

    # --------------------------- INSERT steps functions -----------------------

    def _insertAllSteps(self):
        self._defineParamsName()
        self._initializeCryosparcProject()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.processStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions ------------------------------

    def processStep(self):
        print(pwutils.yellowStr("Sharppening started..."), flush=True)
        self.doSharppening()

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to Relion file
        """
        self._initializeUtilsVariables()
        fnVolName = "cryosparc_%s_%s_map_sharp.mrc" % (self.projectName.get(),
                                                       self.runSharppening.get())
        # Copy the CS output sharpened volume to extra folder
        copyFiles(os.path.join(self.projectPath, self.projectName.get(),
                               self.runSharppening.get()),
                  self._getExtraPath(),
                  files=[fnVolName])

        fnVol = os.path.join(self._getExtraPath(), fnVolName)
        vol = Volume()
        fixVolume(fnVol)
        vol.setFileName(fnVol)
        vol.setSamplingRate(self._getInputVolume().getSamplingRate())
        self._defineOutputs(outputVolume=vol)

    def _defineParamsName(self):
        """ Define a list with all protocol parameters names"""
        self._paramsName = ['sharp_bfactor',
                            'sharp_use_FSC_full',
                            'sharp_falloff_order',
                            'sharp_falloff_offset',
                            'sharp_generate_new_mask',
                            'sharp_mask_thresh',
                            'sharp_mask_near_A',
                            'sharp_mask_far_A',
                            'sharp_do_spherical_mask',
                            'sharp_do_expand_mask']
        self.lane = str(self.getAttributeValue('compute_lane'))

    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        """ Should be overwritten in subclasses to
            return summary message for NORMAL EXECUTION.
        """
        validateMsgs = cryosparcValidate()
        if not validateMsgs:
            validateMsgs = gpusValidate(self.getGpuList(),
                                        checkSingleGPU=True)
            if not validateMsgs:
                if self.sharp_bfactor.get() >= 0.0:
                    validateMsgs.append('b-factor value must be negative')
        return validateMsgs

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputVolume'):
            summary.append("Output objects not ready yet.")
        else:
            summary.append("Input Refinement Protocol: %s" %
                           self.getObjectTag('inputRefinement'))
            summary.append("b-factor: %s" % self.sharp_bfactor.get())
            summary.append("------------------------------------------")

            summary.append("Output volume %s" %
                           self.getObjectTag('outputVolume'))
        return summary

    def doSharppening(self):

        input_group_connect = {"volume": self.volume.get()}

        input_result_connect = None
        if self._getInputVolume().hasHalfMaps():
            input_result_connect = {"volume.0.map_half_A": self.importVolumeHalfA.get(),
                                    "volume.0.map_half_B": self.importVolumeHalfB.get()}

        params = {}

        for paramName in self._paramsName:
            params[str(paramName)] = str(self.getAttributeValue(paramName))

        # Determinate the GPUs to use (in dependence of
        # the cryosparc version)
        try:
            gpusToUse = self.getGpuList()
        except Exception:
            gpusToUse = False

        runSharppeningJob = enqueueJob(self._className,
                                         self.projectName.get(),
                                         self.workSpaceName.get(),
                                         str(params).replace('\'', '"'),
                                         str(input_group_connect).replace('\'', '"'),
                                         self.lane, gpusToUse,
                                         result_connect=input_result_connect)

        self.runSharppening = String(runSharppeningJob.get())
        self.currenJob.set(runSharppeningJob.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(),
                         self.runSharppening.get(),
                         "An error occurred in the particles subtraction process. "
                         "Please, go to cryosPARC software for more "
                         "details.")
        clearIntermediateResults(self.projectName.get(), self.runSharppening.get())