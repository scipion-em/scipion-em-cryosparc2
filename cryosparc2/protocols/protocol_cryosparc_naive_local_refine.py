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
from pwem.protocols import ProtOperateParticles

import pyworkflow.utils as pwutils
from pyworkflow.object import String
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        LEVEL_ADVANCED, IntParam, Positive,
                                        BooleanParam, EnumParam)
from pwem.objects import Volume

from .protocol_base import ProtCryosparcBase
from ..convert import (defineArgs, convertCs2Star, createItemMatrix,
                       setCryosparcAttributes)
from ..utils import (addComputeSectionParams, calculateNewSamplingRate,
                     cryosparcValidate, gpusValidate, enqueueJob,
                     waitForCryosparc, clearIntermediateResults, fixVolume,
                     copyFiles)
from ..constants import *


class ProtCryoSparcNaiveLocalRefine(ProtCryosparcBase, ProtOperateParticles):
    """ Signal subtraction protocol of cryoSPARC.
        Subtract projections of a masked volume from particles.
        """
    _label = 'naive local refinement(Legacy)'
    _className = "naive_local_refine"
    _fscColumns = 6

    def _initialize(self):
        self._defineFileNames()

    def _defineFileNames(self):
        """ Centralize how files are called. """
        myDict = {
            'input_particles': self._getTmpPath('input_particles.star'),
            'out_particles': self._getExtraPath('output_particle.star'),
            'stream_log': self._getPath() + '/stream.log'
        }
        self._updateFilenamesDict(myDict)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      label="Input particles", important=True,
                      help='Select the experimental particles.')
        form.addParam('refVolume', PointerParam, pointerClass='Volume',
                      label="Input map to be projected",
                      important=True,
                      help='Provide the input volume that will be used to '
                           'calculate projections, which will be subtracted '
                           'from the experimental particles. Make sure this '
                           'map was calculated by RELION from the same '
                           'particles as above, and preferably with those '
                           'orientations, as it is crucial that the absolute '
                           'greyscale is the same as in the experimental '
                           'particles.')
        form.addParam('refMask', PointerParam, pointerClass='VolumeMask',
                      label='Mask to be applied to this map',
                      important=True,
                      allowsNull=False,
                      help="Provide a soft mask where the protein density "
                           "you wish to subtract from the experimental "
                           "particles is white (1) and the rest of the "
                           "protein and the solvent is black (0). "
                           "That is: *the mask should INCLUDE the part of the "
                           "volume that you wish to SUBTRACT.*")

        # -----------[Local Refinement]------------------------
        form.addSection(label="Naive local refinement")

        form.addParam('local_align_extent_pix', IntParam, default=3,
                      validators=[Positive],
                      label='Local shift search extent (pix)',
                      help='The maximum extent of local shifts that will be '
                           'searched over, in pixels')

        form.addParam('local_align_extent_deg', IntParam, default=10,
                      label='Local rotation search extent (degrees)',
                      help='The maximum magnitude of the change in rotations '
                           'to search over, in degrees')

        form.addParam('local_align_max_align', FloatParam, default=0.5,
                      validators=[Positive],
                      label='Alignment resolution (degrees)',
                      help='Smallest search distance between angles, in degrees')

        form.addParam('local_align_grid_r', IntParam, default=9,
                      validators=[Positive],
                      label='Local shift search grid size',
                      help='The number of points on the search grid for local '
                           'shifts')

        form.addParam('local_align_grid_t', IntParam, default=9,
                      validators=[Positive],
                      label='Local rotation search grid size',
                      help='The number of points on the search grid for local '
                           'rotations')

        form.addParam('override_final_radwn', BooleanParam,
                      default=False,
                      label='Override final radwn')

        form.addParam('n_iterations', IntParam, default=1,
                      validators=[Positive],
                      label='Override number of iterations')

        # -----[Non Uniform Refinement]----------------------------------------

        form.addSection(label='Non-uniform refinement')
        form.addParam('NU-refine', BooleanParam, default=False,
                      label='Use Non-Uniform Refinement')

        # -----[Refinement]----------------------------------------

        form.addSection(label='Refinement')

        form.addParam('refine_num_final_iterations', IntParam, default=1,
                      label="Number of extra final passes",
                      help='Number of extra passes through the data to do '
                           'after the GS-FSC resolution has stopped improving')

        form.addParam('refine_res_init', IntParam, default=20,
                      validators=[Positive],
                      label="Initial lowpass resolution (A)",
                      help='Applied to input structure')

        form.addParam('refine_res_gsfsc_split', IntParam, default=20,
                      validators=[Positive],
                      label="GSFSC split resolution (A)",
                      help='Resolution beyond which two GS-FSC halves are '
                           'independent')

        form.addParam('refine_FSC_inflate_factor', IntParam, default=1,
                      validators=[Positive],
                      expertLevel=LEVEL_ADVANCED,
                      label="FSC Inflate Factor")

        form.addParam('refine_clip', BooleanParam, default=True,
                      label="Enforce non-negativity",
                      help='Clip negative density. Probably should be false')

        form.addParam('refine_window', BooleanParam, default=True,
                      label="Skip interpolant premult",
                      help='Softly window the structure in real space with a '
                           'spherical window. Should be true')

        form.addParam('refine_skip_premult', BooleanParam, default=True,
                      label="Window structure in real space",
                      help='Leave this as true')

        form.addParam('refine_ignore_dc', BooleanParam, default=True,
                      label="Ignore DC component",
                      help='Ignore the DC component of images. Should be true')

        form.addParam('refine_batchsize_init', IntParam, default=0,
                      expertLevel=LEVEL_ADVANCED,
                      label="Initial batchsize",
                      help='Number of images used in the initial iteration. '
                           'Set to zero to autotune')

        form.addParam('refine_batchsize_epsilon', FloatParam, default=0.001,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label="Batchsize epsilon",
                      help='Controls batch size when autotuning batchsizes. '
                           'Set closer to zero for larger batches')

        form.addParam('refine_batchsize_snrfactor', FloatParam, default=40.0,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label="Batchsize snrfactor",
                      help='Specifies the desired improvement in SNR from the '
                           'images when autotuning batchsizes. Directly '
                           'multiplies the number of images in the batch')

        form.addParam('refine_scale_min', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Minimize over per-particle scale")

        form.addParam('refine_scale_align_use_prev', BooleanParam,
                      default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Use scales from previous iteration during "
                            "alignment")

        form.addParam('refine_scale_ctf_use_current', BooleanParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=True,
                      label="Use scales from current alignment in reconstruction",
                      help='Use scales from current alignment in reconstruction')

        form.addParam('refine_scale_start_iter', IntParam, default=0,
                      label="Scale min/use start iter",
                      help='Iteration to start minimizing over per-particle scale')

        form.addParam('refine_noise_model', EnumParam,
                      choices=['symmetric', 'white', 'coloured'],
                      default=0,
                      label="Noise model:",
                      help='Noise model to be used. Valid options are white, '
                           'coloured or symmetric. Symmetric is the default, '
                           'meaning coloured with radial symmetry')

        form.addParam('refine_noise_priorw', IntParam, default=50,
                      validators=[Positive],
                      expertLevel=LEVEL_ADVANCED,
                      label="Noise priorw",
                      help='Weight of the prior for estimating noise (units of '
                           '# of images)')

        form.addParam('refine_noise_initw', IntParam, default=200,
                      validators=[Positive],
                      expertLevel=LEVEL_ADVANCED,
                      label="Noise initw",
                      help='Weight of the initial noise estimate (units of # '
                           'of images)')

        form.addParam('refine_noise_init_sigmascale', IntParam, default=3,
                      validators=[Positive],
                      expertLevel=LEVEL_ADVANCED,
                      label="Noise initial sigma-scale",
                      help='Scale factor initially applied to the base noise '
                           'estimate')

        form.addParam('refine_mask', EnumParam,
                      choices=['dynamic', 'static', 'null'],
                      default=0,
                      label="Mask:",
                      help='Type of masking to use. Either "dynamic", '
                           '"static", or "null"')

        form.addParam('refine_dynamic_mask_thresh_factor', FloatParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=0.2,
                      validators=[Positive],
                      label="Dynamic mask threshold (0-1)",
                      help='Level set threshold for selecting regions that are '
                           'included in the dynamic mask. Probably don\'t need '
                           'to change this')

        form.addParam('refine_dynamic_mask_near_ang', FloatParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=3.0,
                      validators=[Positive],
                      label="Dynamic mask near (A)",
                      help='Controls extent to which mask is expanded. At the '
                           'near distance, the mask value is 1.0 (in A)')

        form.addParam('refine_dynamic_mask_far_ang', FloatParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=6,
                      validators=[Positive],
                      label="Dynamic mask far (A)",
                      help='Controls extent to which mask is expanded. At the '
                           'far distance the mask value becomes 0.0 (in A)')

        # --------------[Compute settings]---------------------------
        form.addSection(label="Compute settings")
        addComputeSectionParams(form, allowMultipleGPUs=False)

    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._defineFileNames()
        self._defineParamsName()
        self._initializeCryosparcProject()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.processStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions ------------------------------
    def processStep(self):
        print(pwutils.yellowStr("Local Refinement started..."), flush=True)
        self.doLocalRefine()

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to Relion file
        """
        self._initializeUtilsVariables()
        idd, itera = self.findLastIteration(self.runLocalRefinement.get())
        csOutputFolder = os.path.join(self.projectPath, self.projectName.get(),
                                      self.runLocalRefinement.get())
        csOutputPattern = "cryosparc_%s_%s_%s" % (self.projectName.get(),
                                                  self.runLocalRefinement.get(),
                                                  itera)
        csParticlesName = csOutputPattern + "_particles.cs"

        fnVolName = csOutputPattern + "_volume_map.mrc"
        half1Name = csOutputPattern + "_volume_map_half_A.mrc"
        half2Name = csOutputPattern + "_volume_map_half_B.mrc"

        # Copy the CS output volume and half to extra folder
        copyFiles(csOutputFolder, self._getExtraPath(), files=[csParticlesName,
                                                               fnVolName,
                                                               half1Name,
                                                               half2Name])

        csFile = os.path.join(self._getExtraPath(), csParticlesName)

        outputStarFn = self._getFileName('out_particles')
        argsList = [csFile, outputStarFn]

        parser = defineArgs()
        args = parser.parse_args(argsList)
        convertCs2Star(args)

        fnVol = os.path.join(self._getExtraPath(), fnVolName)
        half1 = os.path.join(self._getExtraPath(), half1Name)
        half2 = os.path.join(self._getExtraPath(), half2Name)

        imgSet = self._getInputParticles()
        vol = Volume()
        fixVolume([fnVol, half1, half2])
        vol.setFileName(fnVol)
        vol.setSamplingRate(calculateNewSamplingRate(vol.getDim(),
                                                     imgSet.getSamplingRate(),
                                                     imgSet.getDim()))
        vol.setHalfMaps([half1, half2])

        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(imgSet)
        self._fillDataFromIter(outImgSet)

        self._defineOutputs(outputVolume=vol)
        self._defineSourceRelation(self.inputParticles.get(), vol)
        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(self.inputParticles.get(), outImgSet)
        self.createFSC(idd, imgSet, vol)

    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        """ Should be overwritten in subclasses to
               return summary message for NORMAL EXECUTION.
               """
        validateMsgs = cryosparcValidate()
        if not validateMsgs:
            validateMsgs = gpusValidate(self.getGpuList(), checkSingleGPU=True)
            if not validateMsgs:
                self._validateDim(self._getInputParticles(), self.refVolume.get(),
                                  validateMsgs, 'Input particles', 'Input volume')
        return validateMsgs

    def _summary(self):
        summary = []
        if (not hasattr(self, 'outputVolume') or
                not hasattr(self, 'outputParticles') or
                not hasattr(self, 'outputFSC')):
            summary.append("Output objects not ready yet.")
        else:
            summary.append("Input Particles: %s" %
                           self.getObjectTag('inputParticles'))
            summary.append("Input Volume: %s" %
                           self.getObjectTag('refVolume'))
            summary.append("Input Mask: %s" %
                           self.getObjectTag('refMask'))
            summary.append("------------------------------------------")
            summary.append("Output particles %s" %
                           self.getObjectTag('outputParticles'))
            summary.append("Output volume %s" %
                           self.getObjectTag('outputVolume'))
            if self.hasAttribute('mapResolution'):
                summary.append("\nMap Resolution: %s" % self.mapResolution.get())
            if self.hasAttribute('estBFactor'):
                summary.append('\nEstimated Bfactor: %s' % self.estBFactor.get())
        return summary

    # ---------------Utils Functions-------------------------------------------

    def _fillDataFromIter(self, imgSet):
        outImgsFn = 'particles@' + self._getFileName('out_particles')
        imgSet.setAlignmentProj()
        imgSet.copyItems(self._getInputParticles(),
                         updateItemCallback=self._createItemMatrix,
                         itemDataIterator=md.iterRows(outImgsFn,
                                                      sortByLabel=md.RLN_IMAGE_ID))

    def _createItemMatrix(self, particle, row):
        createItemMatrix(particle, row, align=ALIGN_PROJ)
        setCryosparcAttributes(particle, row, md.RLN_PARTICLE_RANDOM_SUBSET)

    def _defineParamsName(self):
        """ Define a list with all protocol parameters names"""
        self._paramsName = ['local_align_extent_pix', 'local_align_extent_deg',
                            'local_align_max_align', 'local_align_grid_r',
                            'local_align_grid_t', 'override_final_radwn',
                            'n_iterations',
                            'refine_num_final_iterations',
                            'refine_res_init',
                            'refine_res_gsfsc_split',
                            'refine_clip',
                            'refine_window', 'refine_skip_premult',
                            'refine_ignore_dc',
                            'refine_batchsize_init',
                            'refine_batchsize_snrfactor',
                            'refine_batchsize_epsilon',
                            'refine_scale_min', 'refine_scale_align_use_prev',
                            'refine_scale_ctf_use_current',
                            'refine_scale_start_iter',
                            'refine_noise_model', 'refine_noise_priorw',
                            'refine_noise_initw',
                            'refine_mask',
                            'refine_dynamic_mask_thresh_factor',
                            'refine_dynamic_mask_near_ang',
                            'refine_dynamic_mask_far_ang',
                            'compute_use_ssd']
        self.lane = str(self.getAttributeValue('compute_lane'))

    def doLocalRefine(self):
        """
        :return:
        """
        if self.mask.get() is not None:
            input_group_connect = {"particles": self.particles.get(),
                                   "volume": self.volume.get(),
                                   "mask": self.mask.get()}
        else:
            input_group_connect = {"particles": self.particles.get(),
                                   "volume": self.volume.get()}

        input_result_connect = None
        if self._getInputVolume().hasHalfMaps():
            input_result_connect = {"volume.0.map_half_A": self.importVolumeHalfA.get(),
                                    "volume.0.map_half_B": self.importVolumeHalfB.get()}

        params = {}

        for paramName in self._paramsName:
            if (paramName != 'refine_noise_model' and
                    paramName != 'refine_mask'):
                params[str(paramName)] = str(self.getAttributeValue(paramName))

            elif paramName == 'refine_noise_model':
                params[str(paramName)] = str(
                    NOISE_MODEL_CHOICES[self.refine_noise_model.get()])
            elif paramName == 'refine_mask':
                params[str(paramName)] = str(
                    REFINE_MASK_CHOICES[self.refine_mask.get()])

        # Determinate the GPUs to use (in dependence of
        # the cryosparc version)
        try:
            gpusToUse = self.getGpuList()
        except Exception:
            gpusToUse = False

        runLocalRefinementJob = enqueueJob(self._className, self.projectName.get(),
                                             self.workSpaceName.get(),
                                             str(params).replace('\'', '"'),
                                             str(input_group_connect).replace('\'', '"'),
                                             self.lane, gpusToUse,
                                             result_connect=input_result_connect)

        self.runLocalRefinement = String(runLocalRefinementJob.get())
        self.currenJob.set(runLocalRefinementJob.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(), self.runLocalRefinement.get(),
                         "An error occurred in the local refinement process. "
                         "Please, go to cryosPARC software for more "
                         "details.")
        clearIntermediateResults(self.projectName.get(), self.runLocalRefinement.get())
