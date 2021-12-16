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

from pkg_resources import parse_version

from pwem import ALIGN_PROJ
from pwem.protocols import ProtOperateParticles

import pyworkflow.utils as pwutils
from pyworkflow import BETA
from pyworkflow.object import String
from pyworkflow.protocol.params import (PointerParam, FloatParam, IntParam,
                                        Positive, BooleanParam, EnumParam)
from pwem.objects import Volume

from .protocol_base import ProtCryosparcBase
from ..convert import (defineArgs, convertCs2Star, createItemMatrix,
                       setCryosparcAttributes)
from ..utils import (addComputeSectionParams, calculateNewSamplingRate,
                     cryosparcValidate, gpusValidate, enqueueJob,
                     waitForCryosparc, clearIntermediateResults,
                     addSymmetryParam, getCryosparcVersion, getSymmetry,
                     fixVolume, copyFiles)
from ..constants import *


class ProtCryoSparcLocalRefine(ProtCryosparcBase, ProtOperateParticles):
    """ Signal subtraction protocol of cryoSPARC.
        Subtract projections of a masked volume from particles.
        """
    _label = 'local refinement'
    _devStatus = BETA
    _protCompatibility = [V3_0_0, V3_1_0, V3_2_0, V3_3_0, V3_3_1]
    _className = "new_local_refine"
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
                      help="Provide a soft mask where the protein density "
                           "you wish to subtract from the experimental "
                           "particles is white (1) and the rest of the "
                           "protein and the solvent is black (0). "
                           "That is: *the mask should INCLUDE the part of the "
                           "volume that you wish to SUBTRACT.*")

        # -----------[Alignment Parameters]------------------------
        form.addSection(label="Alignment Parameters")

        form.addParam('use_alignment_prior', BooleanParam, default=False,
                      label='Use pose/shift gaussian prior during alignment',
                      help='This can help softly penalise rotations/shifts far '
                           'away from the known initial pose, hence increasing '
                           'stability.')

        form.addParam('sigma_prior_r', IntParam, default=15,
                  validators=[Positive],
                  condition="use_alignment_prior == True",
                  label="Standard deviation (deg) of prior over rotation",
                  help='Standard deviation of gaussian prior over rotation magnitude in degrees.')

        form.addParam('sigma_prior_s', IntParam, default=7,
                      validators=[Positive],
                      condition="use_alignment_prior == True",
                      label="Standard deviation (A) of prior over shifts",
                      help='Standard deviation of gaussian prior over shift magnitude in Angstroms.')

        form.addParam('init_r_extent', IntParam, default=20,
                      validators=[Positive],
                      label="Rotation search extent (deg)",
                      help='Rotation search extent in degrees.')

        form.addParam('init_s_extent', IntParam, default=10,
                      validators=[Positive],
                      label="Shift search extent (A)",
                      help='Shift search extent in Angstroms.')

        form.addParam('fulcrum', EnumParam,
                      choices=['mask_center', 'box_center'],
                      default=0,
                      label="Default fulcrum location",
                      help="Where to place the fulcrum by default. Can be set "
                           "to the center of mass of the mask, or the "
                           "box center.")

        # -----------[Homogeneous Refinement]------------------------
        form.addSection(label="Homogeneous Refinement")

        addSymmetryParam(form, help="Symmetry String (C, D, I, O, T). E.g. C1, "
                                    "D7, C4, etc")

        form.addParam('refine_res_align_max', FloatParam, default=None,
                      allowsNull=True,
                      label="Maximum align resolution (A)",
                      help='Manual override for maximum resolution that is '
                           'used for alignment. This value is normally '
                           'set by the GS-FSC')

        form.addParam('refine_res_init', IntParam, default=12,
                      validators=[Positive],
                      label="Initial lowpass resolution (A)",
                      help='Applied to input structure')

        form.addParam('refine_gs_resplit', BooleanParam, default=False,
                      label='Force re-do GS split',
                      help='Force re-splitting the particles into two random '
                           'gold-standard halves. If this is not set, split '
                           'is preserved from input alignments (if connected). '
                           'Note: if particles are coming directly from an '
                           'ab-initio job, this must be True.')

        form.addParam('refine_do_marg', BooleanParam, default=True,
                      label='Marginalization',
                      help='Efficiently marginalize over poses and shifts. '
                           'Can improve results on small molecules..')

        form.addParam('refine_nu_enable', BooleanParam, default=True,
                      label='Non-uniform refine enable',
                      help='Enable cross-validation-optimal non-uniform '
                           'regularization during refinement.')

        form.addParam('refine_clip', BooleanParam, default=False,
                      label='Enforce non-negativity',
                      help='Bring negative density up to 0, prior to alignment. '
                           'May help in some cases, but recommended to leave '
                           'off in most cases.')

        form.addParam('refine_mask', EnumParam,
                      choices=['dynamic', 'static', 'null'],
                      default=0,
                      label="Mask:",
                      help='Type of masking to use. Either "dynamic", '
                           '"static", or "null"')

        form.addParam('refine_dynamic_mask_near_ang', FloatParam, default=3.0,
                      validators=[Positive],
                      label="Dynamic mask near (A)",
                      help='Controls extent to which mask is expanded. At the '
                           'near distance, the mask value is 1.0 (in A)')

        form.addParam('refine_dynamic_mask_far_ang', FloatParam, default=12.0,
                      validators=[Positive],
                      label="Dynamic mask far  (A)",
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

        # Copy the CS output to extra folder
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
            csVersion = getCryosparcVersion()
            if [version for version in self._protCompatibility
                if parse_version(version) >= parse_version(csVersion)]:
                validateMsgs = gpusValidate(self.getGpuList(),
                                            checkSingleGPU=True)
                if not validateMsgs:
                    particles = self._getInputParticles()
                    self._validateDim(particles,
                                      self.refVolume.get(),
                                      validateMsgs, 'Input particles',
                                      'Input volume')
                    if not particles.hasCTF():
                        validateMsgs.append("The Particles has not associated a "
                                            "CTF model")
                        if not validateMsgs and not particles.hasAlignment3D():
                            validateMsgs.append("The Particles has not a 3D "
                                                "alignment")
            else:
                validateMsgs.append("The protocol is not compatible with the "
                                    "cryoSPARC version %s" % csVersion)

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
                summary.append(
                    "\nMap Resolution: %s" % self.mapResolution.get())
            if self.hasAttribute('estBFactor'):
                summary.append(
                    '\nEstimated Bfactor: %s' % self.estBFactor.get())
        return summary

    # ---------------Utils Functions------------------------------------

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
        self._paramsName = ['use_alignment_prior',
                            'init_r_extent',
                            'init_s_extent',
                            'fulcrum',
                            'refine_res_align_max',
                            'refine_res_init',
                            'refine_gs_resplit',
                            'refine_do_marg',
                            'refine_nu_enable',
                            'refine_clip',
                            'refine_mask',
                            'refine_dynamic_mask_near_ang',
                            'refine_dynamic_mask_far_ang',
                            'intermediate_plots',
                            'sigma_prior_r',
                            'sigma_prior_s',
                            'compute_use_ssd',
                            'refine_symmetry']
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

        params = {}

        for paramName in self._paramsName:
            if (paramName != 'refine_mask' and
                    paramName != 'refine_res_align_max' and
                    paramName != 'fulcrum' and
                    paramName != 'intermediate_plots' and
                    paramName != 'sigma_prior_r' and
                    paramName != 'sigma_prior_s' and
                    paramName != 'refine_symmetry'):
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif paramName == 'refine_mask':
                params[str(paramName)] = str(
                    REFINE_MASK_CHOICES[self.refine_mask.get()])
            elif paramName == 'fulcrum':
                params[str(paramName)] = str(
                    REFINE_FULCRUM_LOCATION[self.fulcrum.get()])
            elif (paramName == 'refine_res_align_max' and
                  self.getAttributeValue(paramName) is not None and
                  float(self.getAttributeValue(paramName)) > 0):
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif paramName == 'intermediate_plots':
                params[str(paramName)] = 'False'
            elif (paramName == 'sigma_prior_r' or
                  paramName == 'sigma_prior_s') and self.getAttributeValue('use_alignment_prior'):
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif paramName == 'refine_symmetry':
                symetryValue = getSymmetry(self.symmetryGroup.get(),
                                           self.symmetryOrder.get())
                params[str(paramName)] = symetryValue
            params['refine_dynamic_mask_start_res'] = str(1000)

        # Determinate the GPUs to use (in dependence of
        # the cryosparc version)
        try:
            gpusToUse = self.getGpuList()
        except Exception:
            gpusToUse = False

        runLocalRefinementJob = enqueueJob(self._className, self.projectName.get(),
                                             self.workSpaceName.get(),
                                             str(params).replace('\'', '"'),
                                             str(input_group_connect).replace('\'',
                                                                             '"'),
                                             self.lane, gpusToUse)

        self.runLocalRefinement = String(runLocalRefinementJob.get())
        self.currenJob.set(runLocalRefinementJob.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(), self.runLocalRefinement.get(),
                         "An error occurred in the local refinement process. "
                         "Please, go to cryosPARC software for more "
                         "details.")
        clearIntermediateResults(self.projectName.get(), self.runLocalRefinement.get())
