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
from pyworkflow.object import String
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        BooleanParam, IntParam)

from .protocol_base import ProtCryosparcBase
from .. import RELIONCOLUMNS
from ..convert import convertCs2Star
from ..utils import (addComputeSectionParams, cryosparcValidate,  enqueueJob, waitForCryosparc, clearIntermediateResults,
                     copyFiles)


class ProtCryoSparcBlobPicker(ProtCryosparcBase):
    """
    Automatically detects particle-like regions in cryo-EM micrographs by
    searching for intensity patterns that resemble Gaussian blobs. The protocol
    is designed as a rapid particle picking strategy for datasets where
    approximate particle size and shape are already known, allowing users to
    generate initial particle coordinates for downstream extraction,
    classification, and reconstruction workflows.

    AI Generated:

    Blob Picker (ProtCryoSparcBlobPicker) — User Manual
        Overview

        The Blob Picker protocol provides an automated method for identifying
        candidate particles directly from cryo-EM micrographs without requiring
        reference templates or previously reconstructed structures. Instead of
        relying on structural projections, the protocol searches for local image
        features that match expected particle dimensions and contrast patterns.
        This approach is especially useful during early stages of processing,
        rapid dataset screening, or exploratory analyses where no reliable
        templates are yet available.

        In practical cryo-EM workflows, blob-based picking is commonly used as
        an initial strategy to obtain a broad particle population before more
        specialized refinement or template-based approaches are introduced. It
        is particularly valuable for newly acquired datasets, heterogeneous
        samples, or projects where the biological target is poorly characterized.

        Inputs and Experimental Context

        The protocol requires a set of input micrographs representing the raw or
        motion-corrected cryo-EM images from which particles will be detected.
        Successful picking depends strongly on the quality of these micrographs.
        Images with severe contamination, crystalline ice, strong charging, or
        poor contrast may generate unreliable particle coordinates or excessive
        false positives.

        Users define approximate particle size limits through minimum and
        maximum particle diameters. These values establish the expected scale of
        detectable features and are among the most biologically important
        parameters in the protocol. Accurate diameter estimates improve the
        separation between real particles and background noise, while incorrect
        values may either miss true particles or include contaminants.

        In many biological projects, the minimum diameter should correspond to
        the smallest expected projection dimension, while the maximum diameter
        should reflect the largest visible orientation of the particle. Flexible
        assemblies or elongated complexes may require wider diameter ranges.

        Blob Shapes and Picking Strategies

        The protocol supports multiple blob geometries that adapt to different
        biological particle appearances. Circular blobs are typically the best
        starting point for globular proteins, compact assemblies, and symmetric
        complexes. Because many cryo-EM particles appear approximately circular
        in projection, this mode provides robust general-purpose performance.

        Elliptical blobs are more appropriate for elongated particles such as
        filaments, rod-like assemblies, or anisotropic macromolecular complexes.
        In these situations, introducing anisotropic shape information often
        improves detection sensitivity and reduces false positives generated by
        ice contamination or carbon edges.

        Ring-shaped blobs can be useful for hollow or toroidal structures whose
        strongest contrast appears near the particle perimeter rather than at
        the center. Examples include membrane channels, ring complexes, or
        cage-like assemblies with internal solvent cavities.

        In practice, users often begin with circular picking and only introduce
        elliptical or ring-based detection if the biological target exhibits
        obvious non-spherical morphology.

        Particle Separation and Density Control

        The protocol allows control over the minimum separation distance between
        detected particles. This parameter strongly influences particle density
        and helps avoid overlapping picks. Smaller separation distances increase
        the number of detected particles but may introduce duplicated or closely
        neighboring coordinates. Larger distances produce cleaner coordinate
        sets but may reject valid particles in crowded fields.

        From a biological perspective, optimal separation depends on specimen
        concentration and particle distribution. Sparse datasets tolerate larger
        separations, while densely packed samples often require more permissive
        values to avoid losing usable particles.

        The protocol also limits the total number of local maxima considered
        during picking. This helps stabilize performance in noisy micrographs
        and prevents excessive detection of contaminants or ice features. Very
        high values may increase false positives, whereas excessively low values
        may discard legitimate particles.

        Optional CTF Estimation

        An optional preprocessing stage estimates contrast transfer function
        parameters before particle picking. This is useful in workflows where
        micrographs lack prior CTF information or when users want a unified
        preprocessing pipeline directly within the same protocol.

        Accurate CTF estimation improves downstream processing quality and
        provides additional metadata for later refinement steps. However, it is
        important to understand that particle picking itself remains primarily
        driven by image features and expected particle geometry rather than by
        detailed structural interpretation.

        For high-quality datasets with previously validated CTF parameters,
        repeating estimation may be unnecessary. Conversely, newly collected or
        rapidly screened datasets often benefit from enabling this option.

        Outputs and Biological Interpretation

        After completion, the protocol produces a set of particle coordinates
        associated with the original micrographs. These coordinates define the
        particle centers that can later be used for extraction and downstream
        single-particle analysis.

        If CTF estimation was enabled, an additional set of CTF models is
        generated and linked to the corresponding micrographs. These outputs can
        be reused in later refinement and reconstruction stages.

        Biological users should carefully inspect the picked coordinates before
        continuing processing. Blob-based methods prioritize broad detection and
        therefore may include contaminants, carbon edges, ice features, or
        aggregated particles. Manual inspection of representative micrographs is
        strongly recommended before large-scale extraction.

        Practical Recommendations

        For most cryo-EM projects, a good starting strategy is to use circular
        blobs with conservative diameter estimates and moderate particle
        separation values. After visual inspection, parameters can be adjusted
        iteratively to improve particle coverage or reduce false positives.

        Highly elongated assemblies may benefit from elliptical blobs, while
        hollow complexes sometimes respond better to ring-shaped detection.
        When datasets contain substantial contamination, increasing particle
        separation and reducing the maximum number of candidate peaks often
        improves coordinate quality.

        Blob picking is especially effective as an initial dataset exploration
        method before training neural-network pickers or generating templates
        from 2D class averages. Many cryo-EM workflows begin with blob picking
        specifically because it avoids introducing structural bias during early
        processing stages.

        Final Perspective

        Automated blob picking serves as a fast and flexible entry point into
        cryo-EM particle selection. Although it lacks the structural specificity
        of template-based or deep learning approaches, its simplicity,
        robustness, and independence from prior models make it highly valuable
        for exploratory biological analysis. Careful adjustment of particle
        size, morphology, and separation parameters remains essential for
        obtaining biologically meaningful particle populations suitable for
        downstream reconstruction workflows.
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
            if not self.useQueueForSteps() and not self.useQueue():  # not using queue system
                gpusToUse = self.getGpuList()
            else:  # using queue system
                gpusToUse = False
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