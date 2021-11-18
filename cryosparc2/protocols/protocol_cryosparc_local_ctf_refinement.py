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

from pwem import ALIGN_PROJ
from pwem.protocols import ProtParticles
import pyworkflow.utils as pwutils
from pyworkflow import BETA
from pyworkflow.object import String
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        LEVEL_ADVANCED, IntParam, Positive)

from .protocol_base import ProtCryosparcBase
from ..convert import (defineArgs, convertCs2Star, createItemMatrix,
                       setCryosparcAttributes)
from ..utils import (addComputeSectionParams, cryosparcValidate, gpusValidate,
                     enqueueJob, waitForCryosparc, clearIntermediateResults,
                     copyFiles)
import pwem.emlib.metadata as md


class ProtCryoSparcLocalCtfRefinement(ProtCryosparcBase, ProtParticles):
    """
    Wrapper protocol for the Cryosparc's per-particle Local CTF refinement.
    Performs per-particle defocus estimation for each particle in a dataset,
    against a given 3D reference structure.
    """
    _label = 'local ctf refinement'
    _devStatus = BETA
    _className = "ctf_refine_local"

    def _initialize(self):
        self._createFilenameTemplates()

    def _createFilenameTemplates(self):
        """ Centralize how files are called. """
        myDict = {
            'input_particles': self._getTmpPath('input_particles.star'),
            'out_particles': self._getExtraPath('output_particle.star')
        }
        self._updateFilenamesDict(myDict)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      label="Input particles", important=True,
                      help='Provide a set of particles for local '
                           'CTF refinement.')
        form.addParam('refVolume', PointerParam, pointerClass='Volume',
                      important=True,
                      label="Input volume",
                      help='Provide a reference volume for local '
                           'CTF refinement.')
        form.addParam('refMask', PointerParam, pointerClass='VolumeMask',
                      label='Mask to be applied to this map',
                      important=True,
                      help="Provide a soft mask. if mask is present, use that, "
                           "otherwise use mask_refine if present, otherwise "
                           "fail")

        # -----------[Global CTF Refinement]------------------------
        form.addSection(label="Local CTF Refinement")
        form.addParam('crl_N', FloatParam, default=None,
                      allowsNull=True,
                      expertLevel=LEVEL_ADVANCED,
                      label='Refinement box size (Voxels)',
                      help='Size of reconstruction/image to use for refinement. '
                           'Blank means to use the particle box size '
                           '(upsampling input maps as needed)')
        form.addParam('crl_num_plots', IntParam, default=50,
                      validators=[Positive],
                      label='Num. groups to plot',
                      help='Number of exposure groups to make plots for. '
                           'After this many, stop plotting to save time.')

        form.addParam('crl_min_res_A', IntParam, default=20,
                      validators=[Positive],
                      label='Minimum Fit Res (A)',
                      help='The minimum resolution to use during refinement of '
                           'image aberrations.')

        form.addParam('crl_max_res_A', FloatParam, default=None,
                      label='Maximum Fit Res (A)',
                      expertLevel=LEVEL_ADVANCED,
                      allowsNull=True,
                      help='The maximum resolution to use during refinement of '
                           'image aberrations. If None, use input half-maps '
                           'to compute FSC and set max to FSC=0.5')

        form.addParam('crl_df_range', IntParam, default=2000,
                      label='Defocus Search Range (A +/-)',
                      help='Defocus search range in Angstroms, searching both '
                           'above and below the input defocus by this amount')

        # --------------[Compute settings]---------------------------
        form.addSection(label="Compute settings")
        addComputeSectionParams(form, allowMultipleGPUs=False)

    # --------------------------- INSERT steps functions -----------------------

    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._defineParamsName()
        self._initializeCryosparcProject()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.processStep)
        self._insertFunctionStep(self.createOutputStep)

    # -------------------------- UTILS functions ------------------------------

    def _defineParamsName(self):
        """ Define a list with all protocol parameters names"""
        self._paramsName = ['crl_N',
                            'crl_num_plots',
                            'crl_min_res_A',
                            'crl_max_res_A',
                            'crl_df_range',
                            'compute_use_ssd']
        self.lane = str(self.getAttributeValue('compute_lane'))

    def _getInputMask(self):
        if self.refMask.get() is not None:
            return self.refMask.get()
        else:
            inputProtocolMask = self._getInputPostProcessProtocol().refMask.get()
            if inputProtocolMask is not None:
                return inputProtocolMask

        return None

    # --------------------------- STEPS functions ------------------------------
    def processStep(self):
        print(pwutils.yellowStr("Local Ctf Refinement started..."), flush=True)
        self.doLocalCtfRefinement()

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to Relion file
        """
        self._initializeUtilsVariables()
        outputStarFn = self._getFileName('out_particles')
        csOutputFolder = os.path.join(self.projectPath, self.projectName.get(),
                                      self.runLocalCtfRefinement.get())
        csFileName = "particles.cs"

        # Copy the CS output particles to extra folder
        copyFiles(csOutputFolder, self._getExtraPath(), files=[csFileName])

        csFile = os.path.join(self._getExtraPath(), csFileName)

        argsList = [csFile, outputStarFn]

        parser = defineArgs()
        args = parser.parse_args(argsList)
        convertCs2Star(args)

        imgSet = self._getInputParticles()

        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(imgSet)
        self._fillDataFromIter(outImgSet)

        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(imgSet, outImgSet)

    def _fillDataFromIter(self, imgSet):
        outImgsFn = 'particles@' + self._getFileName('out_particles')
        imgSet.setAlignmentProj()
        imgSet.copyItems(self._getInputParticles(),
                         updateItemCallback=self._createItemMatrix,
                         itemDataIterator=md.iterRows(outImgsFn,
                                                      sortByLabel=md.RLN_IMAGE_ID))

    def _createItemMatrix(self, particle, row):
        createItemMatrix(particle, row, align=ALIGN_PROJ)
        setCryosparcAttributes(particle, row,
                               md.RLN_PARTICLE_RANDOM_SUBSET)

    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        """ Should be overwritten in subclasses to
            return summary message for NORMAL EXECUTION.
        """
        validateMsgs = cryosparcValidate()
        if not validateMsgs:
            validateMsgs = gpusValidate(self.getGpuList(), checkSingleGPU=True)
            if not validateMsgs:
                self._validateDim(self._getInputParticles(),
                                  self._getInputVolume(),
                                  validateMsgs, 'Input particles',
                                  'Input volume')
        return validateMsgs

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputParticles'):
            summary.append("Output Particles not ready yet.")
        else:
            summary.append("Input Particles: %s" %
                           self.getObjectTag('inputParticles'))
            summary.append("Reference Mask: %s" %
                           self.getObjectTag('refMask'))
            summary.append("--------------------------------------------------")
            summary.append("Output particles %s" %
                           self.getObjectTag('outputParticles'))
        return summary

    def doLocalCtfRefinement(self):
        """
         :return:
         """
        input_group_connect = {"particles": self.particles.get(),
                              "volume": self.volume.get(),
                              "mask": self.mask.get()}

        input_result_connect = None
        if self._getInputVolume().hasHalfMaps():
            input_result_connect = {"volume.0.map_half_A": self.importVolumeHalfA.get(),
                                    "volume.0.map_half_B": self.importVolumeHalfB.get()}

        params = {}

        for paramName in self._paramsName:
            if (paramName != 'crl_max_res_A' and paramName != 'crl_N'):
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif self.crl_max_res_A.get() is not None:
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif self.crl_N.get() is not None:
                params[str(paramName)] = str(self.getAttributeValue(paramName))

        # Determinate the GPUs to use (in dependence of
        # the cryosparc version)
        try:
            gpusToUse = self.getGpuList()
        except Exception:
            gpusToUse = False

        runLocalCtfRefinementJob = enqueueJob(self._className, self.projectName.get(),
                                                self.workSpaceName.get(),
                                                str(params).replace('\'', '"'),
                                                str(input_group_connect).replace('\'', '"'),
                                                self.lane, gpusToUse,
                                                result_connect=input_result_connect)

        self.runLocalCtfRefinement = String(runLocalCtfRefinementJob.get())
        self.currenJob.set(runLocalCtfRefinementJob.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(), self.runLocalCtfRefinement.get(),
                         "An error occurred in the particles subtraction process. "
                         "Please, go to cryosPARC software for more "
                         "details.")