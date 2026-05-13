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
import emtable
import numpy

from pwem.objects import CTFModel
import pyworkflow.utils as pwutils
from pyworkflow import NEW
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        BooleanParam, IntParam,
                                        String)

from .protocol_base import ProtCryosparcBase
from .. import RELIONCOLUMNS
from ..convert import convertCs2Star
from ..utils import (addComputeSectionParams, cryosparcValidate,  enqueueJob, waitForCryosparc,
                     copyFiles)


class ProtCryoSparcPatchCTFEstimate(ProtCryosparcBase):
    """
    Patch-based CTF estimation automatically estimates defocus variation for tilted, bent,
    deformed samples and is accurate for all particle sizes and types including flexible and membrane proteins.

    AI Generated:

    Patch CTF Estimation (ProtCryoSparcPatchCTFEstimate) - User Manual
        Overview

        The Patch CTF Estimation protocol estimates the Contrast Transfer Function (CTF)
        parameters of cryo-EM micrographs using a patch-based strategy designed to remain
        robust across a wide range of biological specimens and imaging conditions. Its
        primary purpose is to characterize microscope-induced phase distortions and
        defocus variations so that downstream particle processing steps can accurately
        reconstruct high-resolution structures.

        In practical cryo-EM workflows, CTF estimation is one of the earliest and most
        critical preprocessing stages because the quality of these measurements directly
        influences particle alignment, classification, refinement, and map interpretability.
        Accurate defocus estimation becomes especially important for tilted datasets,
        flexible assemblies, membrane proteins, and heterogeneous samples where local
        optical conditions may vary substantially across the micrograph.

        Biological Context and Motivation

        Cryo-EM images are affected by microscope optics that modulate the recorded signal
        in a frequency-dependent manner. The CTF estimation process attempts to recover
        these optical parameters so that later computational corrections can compensate
        for information loss and phase inversions introduced during imaging.

        Traditional global estimation approaches may struggle when samples exhibit bending,
        local thickness variations, charging effects, uneven ice distribution, or tilted
        acquisition geometries. The patch-based strategy addresses this limitation by
        analyzing multiple local regions independently, improving robustness in difficult
        experimental conditions.

        This approach is particularly valuable for membrane proteins embedded in detergent
        or nanodiscs, filamentous systems, flexible molecular assemblies, and in situ
        datasets where local image conditions are often non-uniform. For many biological
        projects, reliable local defocus estimation significantly improves the quality of
        downstream refinements.

        Inputs and Experimental Considerations

        The protocol requires a set of micrographs as input. These micrographs should
        ideally originate from a well-calibrated acquisition workflow with consistent
        pixel size, acceleration voltage, and microscope metadata. Poor-quality images,
        severe contamination, crystalline ice, or extreme drift may reduce estimation
        accuracy regardless of the computational method used.

        The amplitude contrast parameter represents the fraction of electron scattering
        contributing to amplitude rather than phase contrast. Typical cryo-EM datasets
        commonly use values around 0.07 to 0.1. Although small inaccuracies in this
        parameter are often tolerated, biologically meaningful high-resolution analysis
        benefits from realistic values matching the imaging conditions.

        Resolution Search Range

        The protocol allows users to define minimum and maximum resolution limits for
        CTF estimation. These settings determine which spatial frequencies contribute
        to the fitting procedure.

        Lower-resolution limits help exclude very broad image features that are not
        informative for oscillatory CTF fitting, while upper-resolution limits define
        how far high-frequency information is considered. In routine biological workflows,
        default values are often appropriate, but challenging datasets may benefit from
        adjustment.

        For noisy datasets or thick ice conditions, restricting the highest resolution
        considered during fitting may improve stability. Conversely, exceptionally clean
        datasets with strong Thon rings can support higher-resolution fitting and more
        accurate optical characterization.

        Defocus Search Parameters

        The protocol includes configurable defocus search boundaries that determine the
        range explored during estimation. These values should approximately reflect the
        acquisition strategy used during microscopy.

        Wider search ranges improve robustness when acquisition conditions are uncertain,
        but they may increase runtime and occasionally introduce unstable fits. Narrower
        ranges are generally preferable when the microscope defocus settings are already
        known with confidence.

        Biological users should interpret defocus values in the context of their imaging
        goals. Lower defocus values generally preserve higher-resolution information but
        reduce image contrast, while higher defocus values improve visibility at the cost
        of high-frequency detail.

        Phase Shift Estimation

        The protocol supports phase shift estimation for datasets acquired using phase
        plates or other imaging modalities where additional phase modulation is present.
        Accurate phase-shift determination can substantially improve downstream refinement
        quality and map interpretability.

        The user may define the phase-shift search range and optionally restrict refinement
        exclusively to phase-shift optimization. This mode can be useful when defocus values
        are already reliable and only phase-related corrections need adjustment.

        For conventional cryo-EM datasets acquired without phase plates, default settings
        are generally sufficient. For Volta phase plate experiments or similar acquisition
        strategies, careful phase estimation becomes much more important.

        Outputs and Interpretation

        After execution, the protocol produces a set of estimated CTF models associated
        with the input micrographs. These outputs typically include defocus values,
        astigmatism parameters, phase shifts, and estimated resolution limits describing
        the quality of the fit.

        Biologically, these measurements help determine whether the dataset is suitable
        for high-resolution reconstruction. Large astigmatism, poor fit resolution, or
        inconsistent defocus behavior across micrographs may indicate acquisition problems
        that should be addressed before extensive downstream processing.

        The resulting CTF estimations are intended for subsequent particle extraction,
        refinement, and reconstruction workflows, where accurate optical correction is
        essential for preserving structural detail.

        Practical Recommendations

        In routine cryo-EM practice, users should visually inspect representative CTF
        fits and verify that estimated defocus values are physically reasonable relative
        to acquisition conditions. Extremely inconsistent measurements may indicate poor
        ice quality, contamination, incorrect metadata, or acquisition instability.

        For tilted or heterogeneous datasets, patch-based estimation is often preferable
        to simpler global methods because it better captures local variations across the
        micrograph. Membrane proteins, flexible complexes, and tomography-derived images
        particularly benefit from this local strategy.

        When processing high-quality single-particle datasets, the default parameters
        usually provide reliable results. More advanced tuning is typically reserved for
        difficult imaging conditions or specialized acquisition strategies.

        Final Perspective

        For most cryo-EM projects, accurate CTF estimation forms the optical foundation
        of the entire reconstruction workflow. Reliable characterization of defocus and
        phase behavior directly impacts the quality of particle alignment, classification,
        and final map reconstruction. Careful parameter selection, validation of estimated
        fits, and awareness of the biological sample conditions are essential for obtaining
        meaningful structural results.
    """
    _label = 'ctf_estimation'
    _className = "patch_ctf_estimation_multi"
    _devStatus = NEW

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputMicrographs', PointerParam, important=True,
                      label=pwutils.Message.LABEL_INPUT_MIC,
                      pointerClass='SetOfMicrographs')

        form.addParam('amp_contrast', FloatParam, default=0.1,
                      label='Amplitude Contrast',
                      help='Amplitude constrast to use. Typically 0.07 or 0.1 for cryo-EM data.')

        form.addParam('res_min_align', IntParam, default=25,
                      label='Minimum resolution (A)',
                      help='Minimum resolution (in A) to consider when estimating CTF.')

        form.addParam('res_max_align', IntParam, default=4,
                      label='Maximum resolution (A)',
                      help='Maximum resolution (in A) to consider when estimating CTF.')

        form.addParam('df_search_min', IntParam, default=1000,
                      label='Maximum resolution (A)',
                      help='Defocus range for gridsearch.')

        form.addParam('df_search_max', IntParam, default=40000,
                      label='Maximum resolution (A)',
                      help='Defocus range for gridsearch.')

        form.addParam('phase_shift_min', IntParam, default=0,
                      label='Min. search phase-shift (rad)',
                      help='Phase-shift range for gridsearch.')

        form.addParam('phase_shift_max', FloatParam, default=numpy.pi,
                      label='Min. search phase-shift (rad)',
                      help='Phase-shift range for gridsearch.')

        form.addParam('do_phase_shift_refine_only', BooleanParam, default=False,
                      label='Do phase refine only',
                      help='Whether to carry out refinement over phase shift only')

        """job.param_add('ctf_settings', "override_K_Y", base_value=None, title="Override knots Y", param_type="number",
                      hidden=False, advanced=True,
                      desc='Override automatically selected spline order for Y dimension (vertical)')
        job.param_add('ctf_settings', "override_K_X", base_value=None, title="Override knots X", param_type="number",
                      hidden=False, advanced=True,
                      desc='Override automatically selected spline order for X dimension (horizontal)')

        job.param_add_section('compute_settings', title='Compute settings', desc='')
        job.param_add('compute_settings', "compute_num_gpus", base_value=1, title="Number of GPUs to parallelize",
                      param_type="number", hidden=False, advanced=False,
                      desc='Number of GPUs over which to parallelize computation.')"""

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
        self.info(pwutils.yellowStr("Patch CTF estimate started..."))
        self.doPatchCTFEstimate()
        self.micrographs = String(str(self.runPatchCTF.get()) + '.exposures')

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to star file
        """
        self.info(pwutils.yellowStr("Create output started..."))
        self._initializeUtilsVariables()
        micSetPtr = self._getInputMicrographs()

        micList = {os.path.basename(mic.getFileName()): mic.clone() for mic in micSetPtr}

        # Copy the  CTF output to extra folder
        csOutputFolder = os.path.join(self.projectDir.get(),
                                      self.runPatchCTF.get())
        outputPath = os.path.join(self._getExtraPath(), self.runPatchCTF.get())
        copyFiles(csOutputFolder, outputPath)

        ctfEstimatedFileName = 'exposures_ctf_estimated.cs'
        csFile = os.path.join(outputPath, ctfEstimatedFileName)
        outputStarFn = self._getExtraPath('ctf.star')
        argsList = [csFile, outputStarFn]
        convertCs2Star(argsList)

        outputCtfSet = self._fillSetOfCTF(outputStarFn, micList)

        self._defineOutputs(outputCTF=outputCtfSet)
        self._defineSourceRelation(micSetPtr, outputCtfSet)

    def _fillSetOfCTF(self, outputCTFFn, micList):

        inputMics = self._getInputMicrographs()
        outputCtfSet = self._createSetOfCTF()
        outputCtfSet.setMicrographs(inputMics)
        mics = list(micList.values())

        ctf = CTFModel()
        mdFileName = '%s@%s' % ('micrograph', outputCTFFn)
        table = emtable.Table(fileName=outputCTFFn)

        for mic, row in enumerate(table.iterRows(mdFileName)):
            ctf.setDefocusU(row.get(RELIONCOLUMNS.rlnDefocusU.value))
            ctf.setDefocusV(row.get(RELIONCOLUMNS.rlnDefocusV.value))
            ctf.setPhaseShift(row.get(RELIONCOLUMNS.rlnPhaseShift.value))
            ctf.setResolution(row.get(RELIONCOLUMNS.rlnCtfMaxResolution.value))
            ctf.setDefocusAngle(row.get(RELIONCOLUMNS.rlnDefocusAngle.value))
            ctf.setMicrograph(mics[mic])
            outputCtfSet.append(ctf)

        return outputCtfSet

    def _defineParamsName(self):
        """ Define a list with all protocol parameters names"""

        self._paramsName = ['amp_contrast', 'res_min_align', 'res_max_align', 'df_search_min',
                            'df_search_max', 'phase_shift_min', 'phase_shift_max', 'do_phase_shift_refine_only']

        self.lane = str(self.getAttributeValue('compute_lane'))

    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        """ Should be overwritten in subclasses to
            return summary message for NORMAL EXECUTION.
        """
        validateMsgs = cryosparcValidate()
        return validateMsgs

    def _summary(self):
        summary = []
        return summary

    def doPatchCTFEstimate(self):
        input_group_connect = {"exposures": self.micrographs.get()}
        params = {'classic_mode': 'False'}
        try:
            if not self.useQueueForSteps() and not self.useQueue():  # not using queue system
                gpusToUse = self.getGpuList()
            else:  # using queue system
                gpusToUse = False
        except Exception:
            gpusToUse = False

        for paramName in self._paramsName:
            params[str(paramName)] = str(self.getAttributeValue(paramName))

        runPatchCTFJob = enqueueJob(self._className,
                                      self.projectName.get(),
                                      self.workSpaceName.get(),
                                      str(params).replace('\'', '"'),
                                      str(input_group_connect).replace('\'', '"'),
                                      self.lane, gpusToUse)

        self.runPatchCTF = String(runPatchCTFJob.get())
        self.currenJob.set(runPatchCTFJob.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(),
                         self.runPatchCTF.get(),
                         "An error occurred in the ctf estimation process. "
                         "Please, go to cryoSPARC software for more "
                         "details.", self)

