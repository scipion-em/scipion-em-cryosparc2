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

from pwem import ALIGN_PROJ
from pwem.protocols import ProtOperateParticles

import pyworkflow.utils as pwutils
from pyworkflow.object import String
from pyworkflow.protocol.params import (PointerParam, FloatParam, IntParam,
                                        Positive, BooleanParam, EnumParam)
from pwem.objects import Volume

from .protocol_base import ProtCryosparcBase
from ..convert import (convertCs2Star, createItemMatrix,
                       setCryosparcAttributes)
from ..utils import (addComputeSectionParams, calculateNewSamplingRate,
                     cryosparcValidate, gpusValidate, enqueueJob,
                     waitForCryosparc, clearIntermediateResults,
                     addSymmetryParam, getSymmetry,
                     fixVolume, copyFiles, getOutputPreffix)
from ..constants import *


class ProtCryoSparcLocalRefine(ProtCryosparcBase, ProtOperateParticles):
    """ Signal subtraction protocol of cryoSPARC.
        Subtract projections of a masked volume from particles.
        """
    """
        Performs local refinement of cryo-EM particles in cryoSPARC using a
        reference volume and an optional mask to improve local structural
        features while preserving previously estimated particle orientations.

        AI Generated:

        Local Refinement (ProtCryoSparcLocalRefine) — User Manual
            Overview

            The Local Refinement protocol performs focused refinement of
            cryo-EM particle datasets using cryoSPARC local refinement
            algorithms. Its main purpose is to improve the quality and local
            resolution of a specific structural region while maintaining the
            overall orientation information obtained from previous refinement
            steps. In practical cryo-EM workflows, this protocol is commonly
            used after consensus refinement when flexible domains, peripheral
            regions, or locally heterogeneous areas require additional
            optimization.

            From a biological perspective, local refinement is especially
            useful for improving densities corresponding to mobile domains,
            ligand-binding regions, membrane-associated components, or
            conformationally variable substructures. By restricting alignment
            and refinement to a masked region, the protocol can recover local
            high-resolution information that may otherwise remain blurred in
            global refinements.

            Inputs and General Workflow

            The protocol requires a set of input particles with existing 3D
            alignment information, a reference volume, and a mask defining
            the region of interest. The input particles are refined locally
            against the reference volume while the mask constrains the
            alignment focus to biologically relevant regions.

            A good practical strategy is to use as reference a previously
            refined consensus map generated from the same particle dataset.
            The mask should isolate the structural region intended for local
            optimization while excluding unrelated density and solvent areas.
            Using particles and references originating from the same dataset
            is essential to ensure consistency in scale, alignment, and noise
            statistics.

            Alignment Parameters and Local Search

            The protocol provides several alignment controls that determine
            how particle orientations and shifts are refined during local
            optimization. Rotation and shift search extents define the
            exploration range around the initial particle poses. Smaller
            ranges are generally preferred because local refinement assumes
            that particles are already approximately aligned from previous
            refinement stages.

            An optional Gaussian prior over rotations and shifts can be
            enabled to stabilize refinement. This prior softly penalizes
            orientations that deviate excessively from the initial alignment.
            Biologically, this is particularly useful when refining flexible
            regions with weak signal, where unconstrained searches may lead
            to unstable alignments or overfitting.

            The protocol also allows iterative re-centering of rotations and
            shifts. These options can improve convergence in difficult cases,
            although they are typically recommended only when alignment
            priors are enabled.

            Masking and Focused Refinement

            Masking is one of the central elements of local refinement
            because it determines which structural regions contribute to the
            alignment and reconstruction process. The input mask should
            include the density intended for refinement while excluding
            unrelated regions that may introduce alignment ambiguity.

            In biological applications, focused masks are commonly used to
            refine flexible domains, membrane proteins, ligand interaction
            sites, or compositional variants within large assemblies.
            Accurate masking often produces substantial improvements in local
            map quality and interpretability.

            The protocol supports different masking strategies, including
            dynamic masking, static masking, or no masking. Dynamic masking
            adapts during refinement and is generally useful for flexible
            structures or regions with varying density boundaries. Parameters
            controlling near and far mask expansion determine how smoothly
            the mask transitions into surrounding solvent regions.

            Homogeneous and Non-Uniform Refinement

            The refinement stage includes several cryoSPARC optimization
            strategies designed to improve reconstruction quality. Non-uniform
            refinement can be enabled to apply adaptive regularization during
            optimization, improving reconstruction quality in heterogeneous
            or flexible regions.

            Marginalization over poses and shifts can also be activated to
            improve stability for smaller particles or noisy datasets. This
            approach efficiently integrates uncertainty in alignment
            parameters and may improve convergence under difficult conditions.

            Additional controls include enforcing non-negative density values,
            manually limiting alignment resolution, and forcing a new
            gold-standard split of particles. Re-splitting is particularly
            important when particles originate from ab-initio reconstructions
            or workflows where previous half-set consistency is uncertain.

            Symmetry and Resolution Control

            The protocol supports standard symmetry definitions such as
            cyclic, dihedral, tetrahedral, octahedral, and icosahedral
            symmetries. Proper symmetry assignment is biologically important
            because incorrect symmetry can introduce structural artifacts or
            mask meaningful asymmetry.

            Users may optionally define a maximum alignment resolution to
            restrict high-frequency information during alignment. This can
            improve stability in cases where local high-resolution features
            are unreliable or dominated by noise.

            Outputs and Their Interpretation

            After execution, the protocol generates a refined 3D volume,
            updated particle alignments, and associated half maps used for
            FSC-based resolution estimation. The refined particles preserve
            their original identities but contain updated local alignment
            parameters optimized for the selected region.

            The resulting refined map should be interpreted in the context of
            the applied mask and refinement strategy. Improvements in local
            resolution often indicate successful focusing of structurally
            stable regions, while limited improvement may reflect intrinsic
            flexibility or insufficient particle signal.

            FSC curves and estimated map resolution provide quantitative
            assessment of reconstruction quality. However, biological
            interpretation should also include careful visual inspection of
            density continuity, side-chain definition, and consistency with
            known structural features.

            Practical Recommendations

            In most cryo-EM workflows, local refinement is most effective
            when starting from a high-quality consensus refinement and using
            a carefully designed soft mask. Excessively large masks may
            reduce refinement specificity, while overly restrictive masks may
            introduce edge artifacts or unstable alignments.

            For flexible assemblies, enabling alignment priors and
            non-uniform refinement often improves robustness. Dynamic masking
            is typically beneficial for regions with poorly defined density
            boundaries, whereas static masks may provide more stable behavior
            for rigid domains.

            Users should begin with moderate rotational and translational
            search ranges and only expand them if convergence problems are
            observed. Excessively broad searches may increase runtime and
            reduce alignment stability without improving biological results.

            Final Perspective

            Local refinement is one of the most important strategies for
            improving biologically relevant structural details in cryo-EM
            reconstructions. Rather than refining the entire particle
            globally, this protocol focuses computational effort on specific
            regions of interest, enabling improved interpretation of flexible
            domains, interaction interfaces, and localized conformational
            variability.

            Successful application depends strongly on the quality of the
            initial consensus refinement, the biological relevance of the
            selected mask, and the careful tuning of alignment constraints.
            When applied appropriately, local refinement can substantially
            improve structural interpretability and downstream biological
            analysis.
        """
    _label = 'local refinement'
    _protCompatibility = [V3_3_1, V3_3_2, V4_0_0,  V4_0_1, V4_0_2, V4_0_3,
                          V4_1_0, V4_1_1, V4_1_2, V4_2_0, V4_2_1, V4_3_1, V4_4_0, V4_4_1, V4_5_1,
                          V4_5_3, V4_6_0, V4_6_1, V4_6_2, V4_7_0, V4_7_1]
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

        form.addParam('reinitialize_rs', BooleanParam, default=False,
                      label='Re-center rotations each iteration?',
                      help='If true, strongly recommended to use prior.')

        form.addParam('reinitialize_ss', BooleanParam, default=False,
                      label='Re-center shifts each iteration?',
                      help='If true, strongly recommended to use prior.')

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
                      help='Controls extent to which mask is expanded. At the ')

        form.addParam('refine_dynamic_mask_start_res', FloatParam, default=12.0,
                      validators=[Positive],
                      label="Dynamic mask start resolution (A)",
                      help='Map resolution at which to start dynamic masking (in A)')

        form.addParam('refine_dynamic_mask_use_abs', BooleanParam, default=False,
                      label='Dynamic mask use absolute value',
                      help='Include negative regions if they are more negative than the threshold')


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
        self.info(pwutils.yellowStr("Local Refinement started..."))
        self.doLocalRefine()

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to Relion file
        """
        self._initializeUtilsVariables()
        idd, itera = self.findLastIteration(self.runLocalRefinement.get())

        csOutputFolder = os.path.join(self.projectDir.get(),
                                      self.runLocalRefinement.get())

        csOutputPattern = "%s%s_%s" % (getOutputPreffix(self.projectName.get()),
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
        convertCs2Star(argsList)

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

        return validateMsgs

    def _summary(self):
        summary = []
        if (not hasattr(self, 'outputVolume') or
                not hasattr(self, 'outputParticles')):
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
                         itemDataIterator=emtable.Table.iterRows(outImgsFn))

    def _createItemMatrix(self, particle, row):
        createItemMatrix(particle, row, align=ALIGN_PROJ)
        setCryosparcAttributes(particle, row,
                               RELIONCOLUMNS.rlnRandomSubset.value)

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
                            'refine_symmetry',
                            'reinitialize_ss',
                            'reinitialize_rs',
                            'refine_dynamic_mask_start_res',
                            'refine_dynamic_mask_use_abs']
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

        # Determinate the GPUs to use (in dependence of
        # the cryosparc version)
        try:
            if not self.useQueueForSteps() and not self.useQueue():  # not using queue system
                gpusToUse = self.getGpuList()
            else:  # using queue system
                gpusToUse = False
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
                         "Please, go to cryoSPARC software for more "
                         "details.", self)
        clearIntermediateResults(self.projectName.get(), self.runLocalRefinement.get())
