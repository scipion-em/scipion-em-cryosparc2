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

from pwem.objects import Coordinate, CTFModel
import pyworkflow.utils as pwutils
from pyworkflow import NEW
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        BooleanParam, IntParam,
                                        String)

from .protocol_base import ProtCryosparcBase
from .. import RELIONCOLUMNS
from ..convert import convertCs2Star
from ..utils import (addComputeSectionParams, cryosparcValidate,  enqueueJob, waitForCryosparc, clearIntermediateResults,
                     copyFiles)


class ProtCryoSparcBlobPicker(ProtCryosparcBase):
    """
    Automatically picks particles by searching for Gaussian signals.
    """
    _label = 'blob_picker'
    _className = "blob_picker_gpu"
    _devStatus = NEW

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputMicrographs', PointerParam, important=True,
                      label=pwutils.Message.LABEL_INPUT_MIC,
                      pointerClass='SetOfMicrographs')

        form.addParam('diameter', IntParam, default=None,
                      label='Minimum particle diameter (px)',
                      help='Minimum particle diameter (px)')

        form.addParam('diameter_max', IntParam, default=None,
                      label='Maximum particle diameter (px)',
                      help='Maximum particle diameter (px)')

        form.addParam('use_circle', BooleanParam, default=True,
                      label='Use circular blob',
                      help='Use three circular blobs, at minimum, average, and maximum particle diameters based on parameters above.')

        form.addParam('use_ellipse', BooleanParam, default=False,
                      label='Use elliptical blob',
                      help='Use an elliptical blob with minor diameter equal to minimum particle diamater param above, and major diameter equal to maximum particle diameter above. You may want to turn off circular blob with this on.')

        form.addParam('use_ring', BooleanParam, default=False,
                      label='Use ring blob',
                      help='Use a ring-shaped blob with inner diameter equal to minimum particle diameter param above, and outer diameter equal to maximum particle diameter above. You may want to turn off circular and elliptical blob with this on.')

        form.addParam('estimate_ctf', BooleanParam, default=False,
                      label='Estimate CTF before pick?',
                      help='Estimate CTF using cryoSPARC Patch CTF algorithm')

        # form.addParam('lowpass_res_template', IntParam, default=20,
        #               label='Lowpass filter to apply to templates (A)',
        #               help='Lowpass filter to apply to templates, (A)s')
        #
        # form.addParam('lowpass_res', IntParam, default=20,
        #               label='Lowpass filter to apply to micrographs (A)',
        #               help='Lowpass filter to apply to micrographs, (A)s')
        #
        # form.addParam('angular_spacing_deg', IntParam, default=5,
        #               label='Angular sampling (degrees)',
        #               help='Angular sampling of templates in degrees. Lower value will mean finer rotations.')

        form.addParam('min_distance', FloatParam, default=1.0,
                      label='Min. separation dist (diameters)',
                      help='Minimum distance between particles in units of particle diameter (min diameter for blob picker). The lower this value, the more and closer particles it picks.')

        form.addParam('num_process', IntParam, default=None,
                      allowsNull=True,
                      label='Number of mics to process',
                      help='Number of micrographs to process. None means all.')

        form.addParam('max_num_hits', IntParam, default=4000,
                      label='Maximum number of local maxima to consider',
                      help='Maximum number of local maxima (peaks) considered.')

        """
            job.param_add('template', "num_plot",             base_value=10,          title="Number of mics to plot",                                               param_type="number",    hidden=False,   advanced=False, desc='Number of micrographs to plot.')
            job.param_add('template', "recenter_templates",       base_value=True,        title="Recenter templates",                                                param_type="boolean",    hidden=False,   advanced=True, desc='Whether or not to recenter the input templates.')
        """

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
        if self.estimate_ctf.get():
            self.info(pwutils.yellowStr("Patch CTF estimate started..."))
            self.doPatchCTFEstimate()
            self.micrographs = String(str(self.runPatchCTF.get()) + '.exposures')

        self.info(pwutils.yellowStr("Blob picker started..."))
        self.doBlobPicker()

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to star file
        """
        self.info(pwutils.yellowStr("Create output started..."))
        self._initializeUtilsVariables()
        micSetPtr = self._getInputMicrographs()

        micList = {os.path.basename(mic.getFileName()): mic.clone() for mic in micSetPtr}

        csOutputFolder = os.path.join(self.projectDir.get(),
                                      self.runBlobPicker.get())
        # Copy the CS output coordinates to extra folder
        outputPath = os.path.join(self._getExtraPath(), self.runBlobPicker.get())
        copyFiles(csOutputFolder, outputPath)
        csPickedParticlesName = 'picked_particles.cs'

        csFile = os.path.join(outputPath, csPickedParticlesName)
        outputStarFn = self._getExtraPath('output_coordinates.star')
        argsList = [csFile, outputStarFn]
        convertCs2Star(argsList)

        outputCoords = self._fillSetOfCoordinates(micSetPtr, outputStarFn, micList)

        # Copy the  CTF output to extra folder
        if self.estimate_ctf.get():
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

        self._defineOutputs(outputCoordinates=outputCoords)
        self._defineSourceRelation(micSetPtr, outputCoords)

    def _fillSetOfCoordinates(self, micSetPtr, outputStarFn, micList):

        outputCoords = self._createSetOfCoordinates(micSetPtr)
        boxSixe = (self.diameter.get() + self.diameter_max.get()) / 2
        outputCoords.setBoxSize(int(boxSixe))

        coord = Coordinate()
        mdFileName = '%s@%s' % ('particles', outputStarFn)
        table = emtable.Table(fileName=outputStarFn)

        for row in table.iterRows(mdFileName):
            coord.setObjId(None)
            micName = os.path.basename(row.get(RELIONCOLUMNS.rlnMicrographName.value))
            splitMicName = micName.split('_')
            if len(splitMicName) > 1:
                micName = '_'.join(splitMicName[1:])
            else:
                micName = splitMicName[-1]
            coord.setMicrograph(micList[micName])
            x = row.get(RELIONCOLUMNS.rlnCoordinateX.value)
            y = row.get(RELIONCOLUMNS.rlnCoordinateY.value)
            dim = micList[micName].getDimensions()
            flipY = dim[1] - y
            coord.setPosition(x, flipY)
            # Add it to the set
            outputCoords.append(coord)

        return outputCoords

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
        self._paramsName = ['diameter', 'diameter_max', 'use_circle',
                            'use_ellipse', 'use_ring', 'min_distance',
                            'num_process', 'max_num_hits']

        self.lane = str(self.getAttributeValue('compute_lane'))

    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        """ Should be overwritten in subclasses to
            return summary message for NORMAL EXECUTION.
        """
        validateMsgs = cryosparcValidate()

        # if not validateMsgs:
        #     micrographs = self._getInputMicrographs()
        #     if micrographs is not None and not micrographs.hasCTF():
        #         validateMsgs.append("The micrographs has not associated a CTF model")

        return validateMsgs

    def _summary(self):
        summary = []
        return summary

    def doPatchCTFEstimate(self):
        input_group_connect = {"exposures": self.micrographs.get()}
        params = {'classic_mode': 'False'}
        className = 'patch_ctf_estimation_multi'
        try:
            gpusToUse = self.getGpuList()
        except Exception:
            gpusToUse = False

        runPatchCTFJob = enqueueJob(className,
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
                         "details.")

    def doBlobPicker(self):

        input_group_connect = {"micrographs": self.micrographs.get()}
        params = {}
        micSetPtr = self._getInputMicrographs()
        samplingRate = micSetPtr.getSamplingRate()
        for paramName in self._paramsName:
            if (paramName != 'diameter' and paramName != 'diameter_max' and
                    paramName != 'num_process'):
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif paramName == 'diameter' and self.diameter.get() is not None:
                params[str(paramName)] = str(int(self.diameter.get()*samplingRate))
            elif paramName == 'diameter_max' and self.diameter_max.get() is not None:
                params[str(paramName)] = str(int(self.diameter_max.get()*samplingRate))
            elif paramName == 'num_process' and self.num_process.get() is not None:
                params[str(paramName)] = str(self.num_process.get())

        # Determinate the GPUs to use (in dependence of
        # the cryosparc version)
        try:
            gpusToUse = self.getGpuList()
        except Exception:
            gpusToUse = False

        runBlobPickerJob = enqueueJob(self._className,
                                         self.projectName.get(),
                                         self.workSpaceName.get(),
                                         str(params).replace('\'', '"'),
                                         str(input_group_connect).replace('\'', '"'),
                                         self.lane, gpusToUse)

        self.runBlobPicker = String(runBlobPickerJob.get())
        self.currenJob.set(runBlobPickerJob.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(),
                         self.runBlobPicker.get(),
                         "An error occurred in the particles picking process. "
                         "Please, go to cryoSPARC software for more "
                         "details.", self)
        clearIntermediateResults(self.projectName.get(), self.runBlobPicker.get())