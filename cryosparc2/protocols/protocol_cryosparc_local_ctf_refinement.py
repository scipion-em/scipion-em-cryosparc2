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
from pwem.protocols import ProtParticles
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        LEVEL_ADVANCED)

from . import ProtCryosparcBase
from ..convert import *
from ..utils import *
import pwem.emlib.metadata as md


class ProtCryoSparcLocalCtfRefinement(ProtCryosparcBase, ProtParticles):
    """
    Wrapper protocol for the Cryosparc's per-particle Local CTF refinement.
    Performs per-particle defocus estimation for each particle in a dataset,
    against a given 3D reference structure.
    """
    _label = 'local ctf refinement(BETA)'

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
                      help='Provide a set of particles for global '
                           'CTF refinement.')
        form.addParam('inputRefinement', PointerParam,
                      pointerClass='ProtCryoSparcRefine3D',
                      label="Select a Refinement protocol",
                      important=True,
                      help='Provide the refinement protocol that will be used. '
                           'If mask_refinement is present, use that, otherwise '
                           'use a selected mask. If mask does not present, '
                           'the protocol will fail')
        form.addParam('refMask', PointerParam, pointerClass='VolumeMask',
                      label='Mask to be applied to this map',
                      important=True,
                      allowsNull=True,
                      help="Provide a soft mask. if mask is present, use that, "
                           "otherwise use mask_refine if present, otherwise "
                           "fail")

        form.addParallelSection(threads=1, mpi=1)

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
        self._insertFunctionStep("convertInputStep")
        self._insertFunctionStep('processStep')
        self._insertFunctionStep('createOutputStep')

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

    def _getInputPostProcessProtocol(self):
        return self.inputRefinement.get()

    def _getInputVolume(self):
        return self._getInputPostProcessProtocol().outputVolume

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
        self.vol = self.importVolume.get() + '.imported_volume.map'
        self.mask = self.importMask.get() + '.imported_mask.mask'
        print(pwutils.yellowStr("Local Ctf Refinement started..."), flush=True)
        self.doLocalCtfRefinement()

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to Relion file
        """
        self._initializeUtilsVariables()
        outputStarFn = self._getFileName('out_particles')

        # Create the output folder
        os.system("cp -r " + self.projectPath + "/" + self.projectName.get() +
                  '/' + self.runLocalCtfRefinement.get() + " " + self._getExtraPath())

        csFileName = "particles.cs"
        csFile = os.path.join(self._getExtraPath(), self.runLocalCtfRefinement.get(),
                              csFileName)

        argsList = [csFile, outputStarFn]

        parser = defineArgs()
        args = parser.parse_args(argsList)
        convertCs2Star(args)

        imgSet = self._getInputParticles()

        outImgSet = self._createSetOfParticles()
        imgSet.setAlignmentProj()
        outImgSet.copyInfo(imgSet)
        self._fillDataFromIter(outImgSet)

        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(imgSet, outImgSet)

    def _fillDataFromIter(self, imgSet):
        outImgsFn = self._getFileName('out_particles')
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
        className = "ctf_refine_local"
        input_group_conect = {"particles": str(self.par),
                              "volume": str(self.vol),
                              "mask": str(self.mask)}
        # {'particles' : 'JXX.imported_particles' }
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

        self.runLocalCtfRefinement = enqueueJob(className, self.projectName.get(),
                                        self.workSpaceName.get(),
                                        str(params).replace('\'', '"'),
                                        str(input_group_conect).replace('\'',
                                                                        '"'),
                                        self.lane, gpusToUse)

        self.currenJob.set(self.runLocalCtfRefinement.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(), self.runLocalCtfRefinement.get(),
                         "An error occurred in the particles subtraction process. "
                         "Please, go to cryosPARC software for more "
                         "details.")