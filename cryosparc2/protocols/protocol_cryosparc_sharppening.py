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
    Performs post-processing sharpening of cryo-EM density maps using cryoSPARC-based
    refinement and filtering strategies. The protocol enhances high-resolution structural
    features that may be attenuated during reconstruction, helping users obtain maps with
    improved interpretability for visualization, atomic modeling, and structural analysis.

    AI Generated:

    CryoSPARC Sharpening (ProtCryoSparcSharppening) — User Manual
        Overview

        The CryoSPARC Sharpening protocol applies map sharpening procedures to reconstructed
        cryo-EM density maps in order to improve the visibility of structural details. In
        cryo-electron microscopy, reconstructed maps often contain dampened high-frequency
        information caused by experimental noise, particle heterogeneity, alignment
        uncertainty, or reconstruction regularization. Sharpening compensates for this loss
        and increases the contrast of fine structural features such as alpha helices,
        beta sheets, side chains, ligand densities, and nucleic acid backbones.

        For biological interpretation, sharpening is frequently one of the final stages
        before model building or publication-quality visualization. Proper sharpening can
        dramatically improve the readability of a map, while excessive sharpening may
        amplify noise and introduce misleading features. The protocol is therefore intended
        both as a visualization enhancement tool and as a preparation step for downstream
        structural analysis.

        Inputs and Biological Context

        The protocol requires a reconstructed three-dimensional volume as input. In most
        workflows, this volume originates from a consensus refinement, local refinement,
        or focused refinement procedure. The input map should already represent the best
        available reconstruction before sharpening is attempted.

        Biological users should ideally provide maps with reliable Fourier Shell Correlation
        information because sharpening depends strongly on the estimated signal quality at
        different spatial frequencies. When half maps are available, the protocol can use
        this information to estimate resolution-dependent attenuation more accurately and
        produce more reliable sharpening results.

        The quality of the sharpening outcome is closely tied to the biological homogeneity
        of the reconstruction. Flexible assemblies, partially occupied ligands, or mixed
        conformational states may respond differently to sharpening across regions of the
        map. Users should therefore evaluate sharpened results carefully and compare them
        against the original reconstruction.

        B-Factor Sharpening

        The central parameter in the protocol is the sharpening B-factor. Negative values
        enhance high-resolution frequencies and increase map contrast, while values that are
        too aggressive may amplify reconstruction noise or create fragmented densities.

        In practical cryo-EM workflows, moderate sharpening is often sufficient to reveal
        secondary structure elements and side-chain features. Stronger sharpening may become
        useful for very high-resolution maps where fine details are already present but
        partially attenuated. However, over-sharpening can distort weak densities, especially
        in flexible regions or membrane proteins with heterogeneous local resolution.

        Biological interpretation should always consider whether newly visible features are
        supported by the underlying signal rather than introduced by excessive enhancement.

        FSC-Based Filtering

        The protocol can use either full-map FSC information or half-map FSC estimates during
        sharpening. Half-map FSC is generally considered more conservative because it reflects
        independent reconstructions and provides a safer estimate of true signal. Full FSC
        may produce stronger enhancement but can sometimes overestimate resolution.

        For most biological analyses, conservative FSC usage is preferable when structural
        interpretation is uncertain or when maps contain heterogeneous regions. In highly
        refined datasets with stable reconstructions, full FSC approaches may provide more
        visually detailed results.

        Low-Pass Filtering and Frequency Falloff

        The sharpening workflow includes low-pass filtering controls that regulate how rapidly
        high-frequency information is attenuated beyond the estimated resolution limit.
        These parameters influence the visual smoothness of the final map and help suppress
        excessive noise amplification.

        Higher falloff values generate steeper transitions between preserved and attenuated
        frequencies, producing crisper visual features but potentially introducing ringing
        artifacts or noisy edges. Softer filtering generally yields more conservative and
        stable maps, especially for medium-resolution reconstructions.

        Biological users should balance visual sharpness against interpretability and avoid
        settings that generate fragmented or discontinuous densities unsupported by the
        experimental data.

        Mask Generation and Map Isolation

        The protocol can generate dedicated masks for sharpening and FSC calculations. Proper
        masking is biologically important because it isolates the macromolecular density from
        surrounding solvent noise. Accurate masking improves FSC estimation and reduces the
        risk of over-enhancing background regions.

        Threshold and mask expansion parameters determine how closely the mask follows the
        molecular envelope. Tight masks emphasize compact density regions, whereas broader
        masks better preserve peripheral domains or flexible extensions.

        For globular proteins, moderate mask expansion often provides stable results. Large
        assemblies with flexible domains, membrane-associated regions, or elongated complexes
        may require more permissive masks to avoid truncating biologically relevant density.

        Final Output Masking

        Optional spherical or expanded masking can be applied to the final sharpened map.
        These operations mainly improve visualization quality and reduce unnecessary solvent
        regions in the exported volume. Spherical masking is particularly useful for removing
        corner artifacts that may become prominent after sharpening.

        Expanded masking may provide cleaner outputs for visualization software and reduce
        storage requirements while preserving the biologically meaningful region of the map.

        Outputs and Interpretation

        The protocol produces a sharpened cryo-EM volume suitable for visualization,
        interpretation, and downstream atomic modeling. The output remains in the same
        coordinate system and sampling framework as the original reconstruction, allowing
        direct comparison with the input map and associated atomic models.

        Biologically, the sharpened map should be interpreted as an enhanced representation
        of the original reconstruction rather than a new reconstruction itself. Structural
        features that become visible after sharpening should always be validated against
        map continuity, local resolution estimates, and independent biological evidence.

        Practical Recommendations

        In routine cryo-EM practice, it is advisable to begin with moderate negative
        B-factors and visually inspect the resulting map before applying stronger sharpening.
        Excessively aggressive sharpening often produces disconnected densities and misleading
        structural features, particularly in flexible or low-resolution regions.

        Users working with heterogeneous assemblies should compare sharpened and unsharpened
        maps side by side to distinguish genuine biological features from amplified noise.
        Conservative masking and FSC usage generally provide the most reliable interpretation
        for challenging datasets.

        For high-resolution structures intended for atomic modeling, iterative adjustment of
        sharpening and masking parameters may substantially improve side-chain visibility and
        backbone continuity. Final interpretation should always integrate biological knowledge,
        local resolution information, and validation metrics.

        Final Perspective

        Map sharpening is one of the most influential post-processing steps in cryo-EM
        structural analysis because it directly affects how biological features are perceived
        and interpreted. Appropriate sharpening can transform a difficult reconstruction into
        a highly interpretable structural map, whereas excessive enhancement can obscure the
        true quality of the data. Careful parameter selection, conservative interpretation,
        and validation against experimental evidence are essential for obtaining reliable
        biological conclusions.
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