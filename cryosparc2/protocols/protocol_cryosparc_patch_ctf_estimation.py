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
    """
    """
        Patch CTF Estimation (ProtCryoSparcPatchCTFEstimate) — User Manual

        Overview

        The Patch CTF Estimation protocol performs patch-based Contrast Transfer
        Function (CTF) estimation for cryo-EM micrographs using cryoSPARC.
        The protocol is designed to estimate local defocus variations across
        micrographs, allowing accurate characterization of tilted, bent, or
        spatially deformed samples. This approach is particularly effective
        for heterogeneous datasets, membrane proteins, flexible complexes,
        and samples where global CTF estimation may not be sufficiently robust.

        In practical cryo-EM workflows, CTF estimation is one of the most
        important preprocessing steps because it directly influences particle
        quality, alignment accuracy, and the achievable reconstruction
        resolution. Reliable estimation of defocus and phase shift parameters
        is essential before particle extraction, classification, or refinement.

        Inputs and General Workflow

        The protocol requires a set of input micrographs that will be analyzed
        independently. During execution, cryoSPARC divides each micrograph
        into multiple local regions or patches and estimates the local CTF
        properties from the Fourier information contained in each region.
        This strategy improves robustness in datasets affected by specimen
        deformation, stage tilt, charging effects, or uneven ice thickness.

        The workflow begins by importing the input micrographs and defining
        the estimation parameters. Once the process starts, the protocol
        launches a cryoSPARC Patch CTF Estimation job and monitors execution
        until completion. After processing, the cryoSPARC outputs are converted
        into STAR-compatible metadata and transformed into Scipion CTF objects
        associated with the original micrographs.

        Resolution Range and Defocus Search

        The protocol allows the user to define the minimum and maximum
        resolution limits used during CTF estimation. These parameters determine
        the frequency range considered during fitting of the CTF model.

        The minimum resolution parameter controls the low-frequency boundary
        of the fitting process. Increasing this value may help reduce the
        influence of large-scale background variations or contamination.
        The maximum resolution parameter defines the highest spatial frequency
        included in the estimation and strongly influences the precision of
        the fitted CTF parameters.

        The protocol also provides a configurable defocus search range.
        Minimum and maximum defocus values define the interval explored
        during the grid search procedure. Broad search ranges are useful
        when the approximate defocus is unknown, while narrower ranges
        improve computational efficiency when acquisition conditions are
        already well characterized.

        Phase Shift Estimation

        The protocol supports phase-shift estimation for datasets collected
        using phase plates or imaging conditions where additional phase
        modulation is present. Users may define minimum and maximum phase-shift
        values to constrain the search interval.

        An optional refinement-only mode allows optimization exclusively
        over phase shift while preserving previously estimated parameters.
        This option may be useful when refining datasets acquired with
        stable defocus conditions but variable phase-plate behavior.

        Amplitude Contrast and Estimation Accuracy

        Amplitude contrast is another important parameter influencing the
        accuracy of the CTF model. Typical cryo-EM datasets commonly use
        values around 0.07 or 0.1 depending on the imaging conditions and
        specimen composition. Correct amplitude contrast estimation improves
        agreement between theoretical and experimental power spectra.

        Because the protocol operates locally using multiple patches,
        it is more tolerant to micrograph imperfections than traditional
        global estimation methods. This makes the protocol especially useful
        for challenging datasets where local variations significantly affect
        image quality.

        Execution and GPU Management

        During execution, the protocol automatically initializes the cryoSPARC
        project environment, converts input metadata, and launches the
        Patch CTF Estimation job using the selected computational resources.
        GPU allocation is handled dynamically depending on whether the workflow
        is executed locally or through a queue system.

        The protocol continuously monitors the cryoSPARC job status and waits
        until processing finishes successfully. If execution fails, an error
        message is generated instructing the user to inspect the cryoSPARC
        environment for additional diagnostic information.

        Outputs and Metadata Conversion

        After completion, the protocol copies the cryoSPARC outputs into the
        Scipion working directory and converts the generated .cs metadata files
        into STAR format. The resulting metadata are then used to populate
        Scipion CTF objects.

        Each output micrograph receives its associated CTF parameters,
        including defocus U, defocus V, astigmatism angle, phase shift,
        and estimated maximum resolution. The protocol preserves the original
        relationship between input micrographs and estimated CTF models,
        ensuring compatibility with downstream cryo-EM processing steps.

        Biological Interpretation

        From a biological perspective, accurate CTF estimation is fundamental
        for obtaining high-resolution reconstructions. Errors in defocus or
        phase-shift estimation can propagate through the workflow and reduce
        the quality of particle alignment and 3D refinement.

        Patch-based estimation is particularly beneficial for modern cryo-EM
        datasets containing tilted acquisitions, flexible membrane proteins,
        uneven ice distributions, or local beam-induced distortions. By
        modeling spatial variability across the micrograph, the protocol
        improves robustness and contributes to more reliable downstream
        structural interpretation.

        Practical Recommendations

        In routine workflows, the default parameter values are often sufficient
        for standard cryo-EM datasets. However, for challenging acquisitions,
        carefully adjusting the resolution limits and defocus search ranges
        can significantly improve estimation stability.

        Wide defocus ranges are recommended when acquisition conditions vary
        substantially between micrographs. Conversely, narrower ranges improve
        speed and stability for homogeneous datasets acquired under controlled
        imaging conditions.

        When processing phase-plate datasets, enabling phase-shift estimation
        is essential for obtaining physically meaningful CTF models. Users
        should also verify the resulting estimated resolutions and defocus
        distributions visually to detect potential outliers or problematic
        micrographs.

        Final Perspective

        Patch-based CTF estimation is not simply a technical preprocessing
        step but a critical component of reliable cryo-EM image analysis.
        Proper estimation of local defocus and phase-shift variations directly
        impacts particle quality assessment, alignment precision, and final
        map resolution. Careful parameter selection and validation of the
        resulting CTF models are therefore essential for robust structural
        interpretation and high-quality cryo-EM reconstructions.
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

