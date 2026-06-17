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
                     fixVolume, copyFiles, getOutputPreffix)


class ProtCryoSparcSharppening(ProtCryosparcBase, ProtAnalysis3D):
    """
    Wrapper protocol for the Cryosparc's to calculate the sharpened map.
    """

    """
        CryoSPARC Sharpening Protocol — User Manual

        Overview

        The ProtCryoSparcSharppening protocol is a Scipion wrapper designed
        to execute cryoSPARC sharpening jobs inside an integrated cryo-EM
        workflow. Its primary purpose is to improve the visibility of
        high-resolution structural information contained within a reconstructed
        3D density map. In cryo-electron microscopy, reconstructed volumes
        frequently suffer from attenuation of high-frequency signal, which
        reduces the visibility of fine structural details. This protocol
        applies sharpening operations guided by B-factor correction and FSC
        information in order to restore part of this signal and improve map
        interpretability for structural analysis and atomic modeling.

        The protocol inherits from both ProtCryosparcBase and ProtAnalysis3D,
        allowing it to combine cryoSPARC execution management with Scipion’s
        generic 3D analysis framework. Internally, the class acts as a bridge
        between both environments by preparing the required parameters,
        launching the cryoSPARC sharpening job, monitoring execution status,
        and automatically importing the resulting sharpened map back into the
        Scipion project.

        Inputs and Workflow

        The protocol requires an input reference volume corresponding to the
        reconstructed cryo-EM map that will undergo sharpening. In most
        biological workflows, this volume originates from a previous
        refinement protocol and represents the final reconstructed density
        prior to post-processing. If half maps are available, the protocol
        also transfers them to cryoSPARC because FSC-based sharpening relies
        on independent half reconstructions to estimate signal quality and
        reduce overfitting risks.

        Execution is organized through Scipion’s step-based workflow
        architecture. The protocol initializes all parameter definitions,
        prepares the cryoSPARC project environment, converts and connects the
        input data, launches the sharpening process, and finally imports the
        generated sharpened map as an output volume. This modular structure
        allows integration into larger automated cryo-EM processing pipelines.

        B-Factor Sharpening

        The central parameter controlling the sharpening process is the
        sharpening B-factor (`sharp_bfactor`). Negative B-factor values are
        used to amplify high-frequency information and therefore sharpen the
        reconstructed density map. Positive values would instead blur the
        volume and are biologically inappropriate for sharpening workflows.
        For this reason, the protocol explicitly validates that the provided
        B-factor is negative before execution begins.

        From a biological perspective, B-factor selection strongly influences
        map interpretability. Moderate sharpening improves visualization of
        side chains, secondary structure elements, and ligand densities,
        whereas excessive sharpening may artificially amplify noise and
        generate misleading structural features. Conversely, insufficient
        sharpening may hide relevant high-resolution information required for
        accurate structural interpretation.

        FSC-Based Filtering and Resolution Control

        The protocol allows sharpening to be guided either by the full FSC
        curve or by the half-map FSC (`sharp_use_FSC_full`). Half-map FSC
        values are generally more conservative because they reduce the risk
        of overestimating resolution and introducing artificial details into
        the sharpened map. Full FSC-based sharpening may produce visually
        sharper reconstructions but can also increase the probability of
        noise amplification.

        Additional parameters such as `sharp_falloff_order` and
        `sharp_falloff_offset` control the behavior of the low-pass filtering
        stage near the estimated resolution limit. These settings determine
        how rapidly high-frequency components decay and therefore influence
        the balance between preserving structural detail and suppressing
        high-frequency noise. Biologically, this balance is essential because
        aggressive filtering may remove meaningful structural information,
        while insufficient filtering may retain artifacts and solvent noise.

        Masking Strategy

        The protocol includes advanced masking options specifically designed
        to improve sharpening robustness. Users may generate a new FSC mask
        optimized for sharpening rather than reusing the refinement mask.
        Parameters such as threshold values and mask dilation distances define
        how the molecular region is segmented and expanded around the density.

        In cryo-EM workflows, masking is one of the most important factors
        affecting post-processing quality because it determines which regions
        contribute to sharpening calculations. Proper masking restricts
        sharpening to biologically meaningful molecular density while
        suppressing empty solvent regions and background noise. Incorrect
        masking may distort flexible regions or artificially enhance noise.

        The protocol can additionally apply a spherical mask to the final map
        in order to suppress over-sharpened density near the corners of the
        reconstruction volume. An optional expanded mask may also be applied
        to reduce file size after compression while preserving the relevant
        structural region.

        GPU Execution and cryoSPARC Integration

        During execution, the protocol dynamically determines GPU allocation
        depending on whether the workflow is running locally or through a
        queue system. This behavior allows compatibility with both standalone
        workstations and distributed cryo-EM computing infrastructures.

        The sharpening job itself is launched through the cryoSPARC job
        scheduling interface using `enqueueJob`. Once submitted, the protocol
        continuously monitors execution status using `waitForCryosparc`,
        ensuring synchronization between Scipion and cryoSPARC. If execution
        fails, the protocol reports the corresponding error so that users can
        inspect additional details directly within cryoSPARC.

        Outputs and Biological Interpretation

        After successful execution, the sharpened map is imported back into
        the Scipion project as an output `Volume` object. The resulting MRC
        file preserves the original sampling rate and is corrected if needed
        before registration inside the workflow. The protocol summary reports
        the originating refinement protocol, the selected B-factor, and the
        generated sharpened output volume.

        From a biological perspective, sharpening should be understood as a
        visualization and interpretability enhancement step rather than a
        source of new structural information. Properly sharpened maps improve
        atomic model building, side-chain assignment, ligand identification,
        and structural comparison between conformational states. However,
        excessive sharpening may generate artificial densities and misleading
        high-resolution features, especially in flexible or poorly resolved
        regions. Careful parameter selection and visual inspection therefore
        remain essential for reliable cryo-EM structural interpretation.
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

        form.addParam('sharp_generate_new_mask', BooleanParam, default=True,
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
        self.info(pwutils.yellowStr("Sharpening started..."))
        self.doSharppening()

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to Relion file
        """
        self._initializeUtilsVariables()
        csOutputFolder = os.path.join(self.projectDir.get(),
                                      self.runSharppening.get())
        fnVolName = "%s%s_map_sharp.mrc" % (getOutputPreffix(self.projectName.get()),
                                            self.runSharppening.get())
        # Copy the CS output sharpened volume to extra folder
        copyFiles(csOutputFolder, self._getExtraPath(),
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
            if not self.useQueueForSteps() and not self.useQueue():  # not using queue system
                gpusToUse = self.getGpuList()
            else:  # using queue system
                gpusToUse = False
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
                         "Please, go to cryoSPARC software for more "
                         "details.", self)
        clearIntermediateResults(self.projectName.get(), self.runSharppening.get())