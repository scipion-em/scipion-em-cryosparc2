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
            gpusToUse = self.getGpuList()
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

