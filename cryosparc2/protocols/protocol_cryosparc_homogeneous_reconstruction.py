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
import ast
import json
import os

import emtable

from pwem import ALIGN_PROJ

import pyworkflow.utils as pwutils
from pyworkflow import NEW
from pyworkflow.object import String
from pyworkflow.protocol.params import (PointerParam, LEVEL_ADVANCED, IntParam,
                                        BooleanParam, EnumParam, FloatParam)
from pwem.objects import Volume

from .protocol_base import ProtCryosparcBase
from ..convert import (convertCs2Star, createItemMatrix,
                       setCryosparcAttributes)
from ..utils import (addComputeSectionParams, calculateNewSamplingRate,
                     cryosparcValidate, gpusValidate, enqueueJob,
                     waitForCryosparc, clearIntermediateResults, fixVolume,
                     copyFiles, addSymmetryParam, getSymmetry,
                     getCryosparcVersion, get_job_streamlog, getOutputPreffix, parse_version, getVersionedEnumValue)
from ..constants import *


class ProtCryoSparcHomogeneousReconstruct(ProtCryosparcBase):
    """
    Creates a 3D reconstruction from particles that already contain
    valid 3D alignment information. The protocol generates a
    homogeneous density map together with gold-standard half maps
    suitable for resolution estimation and downstream structural
    interpretation.

    AI Generated:

    Homogeneous Reconstruction (ProtCryoSparcHomogeneousReconstruct) — User Manual
        Overview

        The Homogeneous Reconstruction protocol generates a 3D cryo-EM
        density map from particles that already possess assigned viewing
        directions and alignment parameters. Its primary goal is to
        reconstruct a consensus structure representing a single dominant
        conformational state using cryoSPARC reconstruction methods.

        In practical biological workflows, this protocol is commonly used
        after particle alignment has already been completed through
        refinement or classification procedures. Because particle poses
        are already known, the protocol focuses on reconstructing the
        highest-quality density map possible from the aligned particle
        population. This makes it especially useful for producing final
        consensus maps, validating alignment quality, or generating maps
        for model building and interpretation.

        Inputs and Reconstruction Strategy

        The protocol requires a set of particles with valid 3D alignment
        parameters and associated CTF information. These alignments define
        how each particle contributes to the final reconstruction. Without
        reliable alignment information, the resulting map may become
        blurred, distorted, or biologically misleading.

        An optional mask may also be provided for FSC calculations. From a
        biological perspective, masking is important because it restricts
        resolution estimation to structurally meaningful regions of the
        reconstruction while excluding solvent noise. Proper masking
        generally produces more reliable FSC curves and more realistic
        resolution estimates.

        The protocol reconstructs both the final consensus map and the
        independent half maps used for gold-standard FSC validation.
        These outputs are critical for assessing map quality and for
        determining whether the reconstruction is suitable for atomic
        interpretation.

        Symmetry and Biological Interpretation

        The protocol supports cyclic, dihedral, tetrahedral, octahedral,
        and icosahedral symmetry as well as helical symmetry options.
        Applying the correct symmetry can dramatically improve map quality
        by increasing the effective signal-to-noise ratio. However,
        biological caution is essential because imposing incorrect
        symmetry may artificially distort structural features or hide
        biologically important asymmetry.

        For globular complexes with known symmetry, standard point-group
        symmetry is often appropriate. In contrast, asymmetric assemblies,
        flexible complexes, or particles containing partially occupied
        regions may require reconstruction without imposed symmetry in
        order to preserve biologically meaningful heterogeneity.

        The protocol also supports helical reconstruction parameters,
        including twist, rise, and helical symmetry order. These options
        are particularly important for filamentous assemblies such as
        cytoskeletal filaments, amyloid fibrils, or helical viral
        structures. Accurate helical parameters are biologically critical
        because even small errors may produce severe artifacts or incorrect
        filament geometry.

        Gold-Standard Reconstruction and Validation

        The reconstruction process preserves or regenerates independent
        gold-standard half sets used for FSC validation. This separation
        is essential in cryo-EM because it minimizes overfitting and
        provides a more realistic estimate of structural resolution.

        In some workflows, particles may already contain an existing
        gold-standard split inherited from previous refinements. In other
        situations, users may choose to regenerate the split in order to
        ensure independence or adapt the reconstruction to a new analysis
        strategy.

        The resulting FSC curves and half maps are especially important
        when maps will be interpreted biologically, deposited in public
        databases, or used for atomic modeling. Reliable validation is
        essential for distinguishing genuine structural detail from noise
        amplification.

        Handedness and Structural Orientation

        The protocol allows inversion of reconstruction handedness. This
        option becomes relevant when the map appears mirrored relative to
        known biological structures or when correcting handedness after
        reconstruction. Biological interpretation strongly depends on
        correct handedness because incorrect chirality can invalidate
        structural conclusions and atomic models.

        When helical symmetry is used together with handedness inversion,
        the helical twist convention must remain consistent with the new
        orientation. Careful verification against known structural motifs
        or external biological information is strongly recommended.

        Aberration Handling and Optical Corrections

        Advanced options allow selective handling of higher-order optical
        aberrations such as beam tilt, trefoil, tetrafoil, and anisotropic
        magnification effects. These parameters are particularly relevant
        for high-resolution datasets where subtle optical distortions may
        limit achievable detail.

        For most routine biological reconstructions, default settings are
        generally sufficient. However, in near-atomic or atomic resolution
        projects, proper correction of these effects may substantially
        improve map sharpness and interpretability.

        The protocol also supports Ewald sphere correction, which becomes
        increasingly important for large particles, very high resolutions,
        or thick specimens. From a biological standpoint, these corrections
        can improve the accuracy of fine structural features that are
        critical for mechanistic interpretation.

        Filtering and Sharpening

        The reconstruction workflow includes optional control over FSC-based
        filtering and sharpening behavior. In standard workflows, automatic
        FSC filtering provides a balanced approach that suppresses noise
        while preserving meaningful signal.

        Advanced users may override filtering behavior to apply manual
        filtering resolutions or sharpening factors. This flexibility is
        useful for specialized visualization goals or experimental map
        interpretation, although aggressive sharpening can easily amplify
        noise and create misleading features.

        Biologically meaningful sharpening should improve visibility of
        secondary structure and side-chain detail without introducing
        artificial fragmentation or discontinuities in the density.

        GPU and Computational Considerations

        The protocol is optimized for GPU execution and supports control
        over reconstruction batch sizes. These settings primarily affect
        computational efficiency and memory usage rather than biological
        interpretation. For very large datasets or limited GPU memory,
        adjusting batch size parameters may improve execution stability.

        Intermediate plotting and diagnostic options are also available
        for users who wish to inspect reconstruction progress in more
        detail. While these outputs are not always necessary for routine
        processing, they can be useful during optimization or troubleshooting
        of challenging datasets.

        Outputs and Their Interpretation

        The protocol produces a reconstructed 3D volume, two independent
        half maps, FSC validation information, and an updated particle set
        containing the reconstruction-associated metadata. These outputs
        form the basis for downstream cryo-EM analysis including map
        visualization, local resolution analysis, atomic modeling, and
        biological interpretation.

        The reconstructed map should always be interpreted together with
        its FSC validation and visual inspection. Strong FSC values alone
        do not guarantee biological correctness if masking, symmetry, or
        alignment assumptions are inappropriate.

        Practical Recommendations

        For most biological applications, it is advisable to begin with
        conservative settings and only introduce advanced corrections when
        justified by data quality or resolution goals. Correct symmetry
        assignment and reliable input alignments are usually the most
        important factors controlling reconstruction quality.

        When processing filamentous assemblies, special attention should
        be given to helical twist and rise parameters because small errors
        may strongly affect map continuity and interpretability. For highly
        symmetric particles, applying the correct point-group symmetry can
        significantly improve reconstruction quality.

        Users should visually inspect half maps, FSC behavior, and map
        features before proceeding to model building. Over-sharpened maps,
        unstable FSC curves, or unrealistic density continuity often
        indicate issues requiring further refinement or parameter adjustment.

        Final Perspective

        Homogeneous reconstruction represents one of the central stages of
        cryo-EM structural analysis because it transforms aligned particle
        images into an interpretable three-dimensional biological map.
        Reliable alignments, appropriate symmetry handling, careful
        validation, and biologically sensible interpretation are all
        essential for producing structurally meaningful reconstructions
        suitable for downstream scientific conclusions.
    """
    _label = 'homogeneous reconstruction'
    _className = "homo_reconstruct"
    _devStatus = NEW
    _fscColumns = 6
    ewsParamsName = []

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
                      label="Input particles", important=True,
                      help='Select the experimental particles.')
        form.addParam('refMask', PointerParam, pointerClass='VolumeMask',
                      label='Mask used for computing FSC',
                      allowsNull=True,
                      help="The mask is used for computing"
                           "the FSC between the reconstructed half-maps")

        form.addSection(label='Homogeneous Reconstruction')

        form.addParam('refine_N_cmp', IntParam, default=None,
                      allowsNull=True,
                      label="Reconstruction box size (voxels)",
                      help='The volume size to use for refinement. Default: raw particle size')

        addSymmetryParam(form, help="Symmetry String (C, D, I, O, T). E.g. C1, "
                                    "D7, C4, etc")

        form.addParam('refine_helical_twist', IntParam, default=None,
                      allowsNull=True,
                      label="Helical twist (degrees)",
                      help='Helical twist to apply in degrees')

        form.addParam('refine_helical_shift', IntParam, default=None,
                      allowsNull=True,
                      label="Helical rise (Angstroms)",
                      help='Helical rise to apply in Angstroms')

        form.addParam('refine_hsym_order', IntParam, default=None,
                      allowsNull=True,
                      label="Helical symmetry order to apply",
                      help='The maximum amount of helical symmetry to impose '
                           'during reconstruction; i.e., the maximum number '
                           'of (twist, rise) pairs to backproject each '
                           'particle image with. For particles picked outside '
                           'of the filament tracer or template picker, this '
                           'should be set to the distance between extracted '
                           'boxes, divided by the helical rise. If left as '
                           'None, will be calculated based on the inter-box '
                           'distance (if available). Set to 1 for no helical '
                           'symmetry during reconstruction.')

        form.addParam('refine_gs_resplit', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Force re-do GS split",
                      help='Force re-splitting the particles into two random '
                           'gold-standard halves. If this is not set, split is '
                           'preserved from input alignments')

        form.addParam('refine_do_flip', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Flip the reconstruction hand",
                      help='Flip the hand of the reconstruction, and correct '
                           'the alignments. If flipping the hand and applying '
                           'helical symmetry, ensure that the sign of the '
                           'helical twist is inverted to match the hand of '
                           'the new symmetry.')

        form.addParam('refine_ignore_tilt', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Ignore the tilt",
                      help='Ignore the tilt')


        form.addParam('refine_ignore_trefoil', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Ignore trefoil",
                      help='Ignore trefoil')

        form.addParam('refine_ignore_tetra', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Ignore tetra",
                      help='Ignore tetra')

        form.addParam('refine_ignore_anisomag', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Ignore anisomag",
                      help='Ignore the anisomag')

        csVersion = getCryosparcVersion()
        if parse_version(csVersion) >= parse_version(V3_3_1):

            form.addParam('refine_do_ews_correct', BooleanParam, default=False,
                          expertLevel=LEVEL_ADVANCED,
                          label="Do EWS correction",
                          help='Whether or not to correct for the curvature of the '
                               'Ewald Sphere.')

            form.addParam('refine_ews_zsign', EnumParam,
                          choices=['positive', 'negative'],
                          default=0,
                          expertLevel=LEVEL_ADVANCED,
                          label="EWS curvature sign",
                          help='Whether to use positive or negative curvature in '
                               'Ewald Sphere correction.')

            form.addParam('refine_ews_simple', EnumParam,
                          choices=['simple', 'iterative'],
                          default=0,
                          expertLevel=LEVEL_ADVANCED,
                          label="EWS correction method",
                          help='Whether to use the simple insertion method, or to '
                               'use an iterative optimization method, for Ewald '
                               'Sphere correction.')

            form.addParam('refine_fsc_mask_opt', BooleanParam, default=False,
                          expertLevel=LEVEL_ADVANCED,
                          label="Optimize FSC mask",
                          help='Whether or not to optimize the mask used for '
                               'calculating FSCs')

            form.addParam('refine_override_filter', BooleanParam, default=False,
                          expertLevel=LEVEL_ADVANCED,
                          label="Override FSC Filtering",
                          help='Whether to override the FSC-filtering and use '
                               'manual filtering and sharpening (set below) '
                               'instead. FSC is still computed, just not used '
                               'for filtering.')

            form.addParam('refine_override_filter_res', FloatParam, default=None,
                          allowsNull=True,
                          expertLevel=LEVEL_ADVANCED,
                          label="Override Filtering Resolution",
                          help='Override filter corner resolution (A)')

            form.addParam('refine_override_filter_order', IntParam, default=8,
                          expertLevel=LEVEL_ADVANCED,
                          label="Override Filtering Order",
                          help='Override filter order (Butterworth)')

            form.addParam('refine_override_filter_bfactor', FloatParam, default=0,
                          expertLevel=LEVEL_ADVANCED,
                          label="Override Filtering Bfactor",
                          help='Override sharpening B-factor. Negative to sharpen')

            self.ewsParamsName = ['refine_do_ews_correct',
                                  'refine_ews_zsign', 'refine_fsc_mask_opt',
                                  'refine_override_filter',
                                  'refine_override_filter_res',
                                  'refine_ews_simple',
                                  'refine_override_filter_order',
                                   'refine_override_filter_bfactor']

        form.addParam('recon_do_expand_mask', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Wide mask final output",
                      help='Apply a wide mask to the final map, to reduce file '
                           'size after compression.')

        form.addParam('intermediate_plots', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Show plots from intermediate steps",
                      help='Hide plots from intermediate steps to speed up '
                           'processing.')

        form.addParam('refine_compute_batch_size', IntParam, default=None,
                      expertLevel=LEVEL_ADVANCED,
                      allowsNull=True,
                      label="GPU batch size of images",
                      help='Batch size of images to process at a time on the '
                           'GPU. If you run out of GPU memory, try setting '
                           'this to a small number to override the '
                           'auto-detect procedure.')

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

    def processStep(self):
        self.info(pwutils.yellowStr("Homogeneous Reconstruction started..."))
        self.doHomogeneousReconstruction()

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to Relion file
        """
        self._initializeUtilsVariables()
        idd = self.findLastIteration(self.runHomogeneousReconstruction.get())

        csOutputFolder = os.path.join(self.projectDir.get(),
                                       self.runHomogeneousReconstruction.get())

        csOutputPattern = "%s%s" % (getOutputPreffix(self.projectName.get()),
                                    self.runHomogeneousReconstruction.get())
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
        self._defineSourceRelation(self.inputParticles, vol)
        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(self.inputParticles, outImgSet)
        self.createFSC(idd, imgSet, vol)

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

    def findLastIteration(self, jobName):
        get_job_streamlog(self.projectName.get(),
                          jobName,
                          self._getFileName('stream_log'))

        # Get the metadata information from stream.log
        with open(self._getFileName('stream_log')) as f:
            data = f.readlines()

        x = ast.literal_eval(data[0])
        if isinstance(x, str):
            x = json.loads(x)

        # Find the ID of last iteration and the map resolution
        for y in x:
            if 'text' in y:
                z = str(y['text'])

                if z.startswith('FSC, after mask auto-tightening'):
                    idd = y['imgfiles'][2]['fileid']

        return idd

    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        """ Should be overwritten in subclasses to
               return summary message for NORMAL EXECUTION.
               """
        validateMsgs = cryosparcValidate()
        if not validateMsgs:
            validateMsgs = gpusValidate(self.getGpuList(), checkSingleGPU=True)
            if not validateMsgs:
                particles = self._getInputParticles()
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
            summary.append("Input Mask: %s" %
                           self.getObjectTag('refMask'))
            summary.append("Symmetry: %s" %
                           getSymmetry(self.symmetryGroup.get(),
                                       self.symmetryOrder.get())
                           )
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

    def _defineParamsName(self):
        """ Define a list with all protocol parameters names"""
        self._paramsName = ['refine_N_cmp', 'refine_gs_resplit',
                            'refine_do_flip', 'refine_ignore_tilt',
                            'refine_ignore_trefoil', 'refine_ignore_tetra',
                            'refine_ignore_anisomag', 'refine_fsc_mask_opt',
                            'recon_do_expand_mask',
                            'intermediate_plots', 'refine_symmetry',
                            'refine_compute_batch_size', 'refine_helical_twist',
                            'refine_helical_shift', 'refine_hsym_order',
                            'compute_use_ssd'] + self.ewsParamsName

        self.lane = str(self.getAttributeValue('compute_lane'))

    def doHomogeneousReconstruction(self):
        input_group_connect = {"particles": self.particles.get()}
        if self.mask.get() is not None:
            input_group_connect["mask"] = self.mask.get()

        params = {}

        for paramName in self._paramsName:
            if (paramName != 'refine_N_cmp' and
                    paramName != 'refine_compute_batch_size' and
                    paramName != 'refine_symmetry' and
                    paramName != 'refine_helical_twist' and
                    paramName != 'refine_helical_shift' and
                    paramName != 'refine_hsym_order' and
                    paramName != 'refine_fsc_mask_opt' and
                    paramName != 'refine_ews_zsign' and
                    paramName != 'refine_ews_simple' and
                    paramName != 'refine_override_filter_res'):
                params[str(paramName)] = str(self.getAttributeValue(paramName))

            elif paramName == 'refine_fsc_mask_opt':
                params[str(paramName)] = str("True")
            elif paramName == 'refine_N_cmp' and self.refine_N_cmp.get() is not None:
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif paramName == 'refine_compute_batch_size' and self.refine_compute_batch_size.get() is not None:
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif paramName == 'refine_symmetry':
                symetryValue = getSymmetry(self.symmetryGroup.get(),
                                           self.symmetryOrder.get())
                params[str(paramName)] = symetryValue
            elif paramName == 'refine_helical_twist' and self.refine_helical_twist.get() is not None:
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif paramName == 'refine_helical_shift' and self.refine_helical_shift.get() is not None:
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif paramName == 'refine_hsym_order' and self.refine_hsym_order.get() is not None:
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif paramName == 'refine_ews_zsign':
                params[str(paramName)] = str(getVersionedEnumValue(self.refine_ews_zsign.get(), EWS_CURVATURE_SIGN, EWS_CURVATURE_SIGN_V5))
            elif paramName == 'refine_ews_simple':
                params[str(paramName)] = (getVersionedEnumValue(self.refine_ews_simple.get(), EWS_CORRECTION_METHOD, EWS_CORRECTION_METHOD_V5))
            elif paramName == 'refine_override_filter_res' and self.refine_override_filter_res.get() is not None:
                params[str(paramName)] = str(self.getAttributeValue(paramName))

        # Determinate the GPUs to use (in dependence of
        # the cryosparc version)
        try:
            if not self.useQueueForSteps() and not self.useQueue():  # not using queue system
                gpusToUse = self.getGpuList()
            else:  # using queue system
                gpusToUse = False
        except Exception:
            gpusToUse = False

        runHomogeneousReconstructionJob = enqueueJob(self._className, self.projectName.get(),
                                                       self.workSpaceName.get(),
                                                       str(params).replace('\'', '"'),
                                                       str(input_group_connect).replace('\'', '"'),
                                                       self.lane, gpusToUse)
        self.runHomogeneousReconstruction = String(runHomogeneousReconstructionJob.get())
        self.currenJob.set(runHomogeneousReconstructionJob.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(), self.runHomogeneousReconstruction.get(),
                         "An error occurred in the homogeneous reconstruction process. "
                         "Please, go to cryoSPARC software for more "
                         "details.", self)
        clearIntermediateResults(self.projectName.get(), self.runHomogeneousReconstruction.get())








