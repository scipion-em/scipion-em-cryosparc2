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
import pwem.protocols as pwprot

import pyworkflow.utils as pwutils
from pyworkflow.object import String
from pyworkflow.protocol.params import (PointerParam, FloatParam, IntParam,
                                        LEVEL_ADVANCED, Positive, BooleanParam,
                                        EnumParam)

from .protocol_base import ProtCryosparcBase
from ..convert import (convertCs2Star, createItemMatrix,
                       setCryosparcAttributes)
from ..utils import (addComputeSectionParams, cryosparcValidate, gpusValidate,
                     enqueueJob, waitForCryosparc, copyFiles,
                     getCryosparcVersion, parse_version)

from ..constants import *


class ProtCryoSparcGlobalCtfRefinement(ProtCryosparcBase, pwprot.ProtParticles):
    """
    Refines global Contrast Transfer Function parameters for cryo-EM particle
    datasets using a high-quality 3D reference and exposure-based aberration
    estimation. The protocol improves the optical consistency of particle data
    by correcting higher-order microscope aberrations and refining beam-related
    effects that can limit the final achievable resolution.

    AI Generated:

    Global CTF Refinement (ProtCryoSparcGlobalCtfRefinement) — User Manual
        Overview

        The Global CTF Refinement protocol is designed to improve the optical
        accuracy of cryo-EM datasets by refining microscope aberration
        parameters against a reconstructed three-dimensional reference map.
        In practical cryo-EM workflows, this stage is commonly performed after
        an initial refinement has already produced a reliable volume and stable
        particle alignments. The protocol uses the information present across
        groups of exposures to estimate and correct systematic imaging errors
        that are not fully captured during standard CTF estimation.

        From a biological perspective, this refinement step is particularly
        important when pursuing high-resolution reconstructions. Small optical
        distortions such as beam tilt, anisotropic magnification, trefoil, or
        tetrafoil aberrations may only produce subtle visual artifacts in early
        maps, yet they can strongly reduce the interpretability of atomic
        features in the final reconstruction. Correcting these effects often
        leads to sharper density, improved side-chain visibility, and more
        reliable structural interpretation.

        Inputs and General Workflow

        The protocol requires three principal inputs: a set of aligned
        particles, a reference volume, and an appropriate soft mask. The
        particle set should already contain stable projection alignments,
        because the refinement assumes that the orientations are sufficiently
        accurate for optical parameter estimation. The reference volume serves
        as the structural target against which aberrations are measured and
        corrected.

        The soft mask is biologically important because it defines which
        regions of the reconstruction contribute most strongly to the
        refinement process. In most cases, the mask should include the ordered
        molecular density while excluding solvent noise and flexible peripheral
        regions. Poor masking may reduce the stability of aberration
        estimation, especially when the reconstruction contains heterogeneous
        or highly mobile domains.

        In workflows where half maps are available, the protocol can use this
        information to preserve proper gold-standard refinement behavior. This
        is particularly important for preventing overfitting in
        high-resolution reconstructions.

        Refinement of Optical Aberrations

        The protocol supports refinement of several higher-order aberration
        terms that commonly affect cryo-EM data quality. Beam tilt correction
        compensates for systematic angular deviations in the electron beam,
        often improving high-frequency information consistency across the
        dataset. Trefoil and tetrafoil refinement target more complex optical
        distortions that may arise from microscope alignment imperfections or
        beam-shaping effects.

        Spherical aberration refinement allows the protocol to compensate for
        deviations in the microscope optical system that can influence phase
        accuracy at high resolution. In newer workflows, anisotropic
        magnification correction is also available, helping account for slight
        directional differences in magnification that can distort structural
        dimensions.

        For many biological datasets, enabling all major aberration refinement
        terms provides the best starting point. However, very noisy datasets,
        low-resolution reconstructions, or small particle populations may not
        support reliable estimation of every parameter. In such cases, more
        conservative refinement strategies are often preferable.

        Iterative Refinement Strategy

        The protocol allows multiple refinement iterations. Iterative
        refinement can be biologically useful because improvements in one
        optical parameter may reveal or stabilize corrections in others. For
        example, correcting anisotropic magnification can improve the accuracy
        of subsequent trefoil or beam tilt estimation.

        In many routine datasets, a single iteration is sufficient to obtain
        measurable improvements. For demanding near-atomic-resolution studies,
        multiple iterations may provide additional gains, especially when
        combined with later rounds of particle polishing or high-resolution
        refinement.

        Ewald Sphere Correction

        Advanced workflows may optionally account for Ewald Sphere curvature.
        This correction becomes increasingly relevant for large particles,
        thick specimens, or very high-resolution reconstructions where the
        flat projection approximation begins to break down.

        From a biological interpretation standpoint, Ewald Sphere correction
        may improve the consistency of fine structural features and enhance map
        interpretability in challenging datasets. However, it is computationally
        more demanding and is usually most beneficial only when the underlying
        data quality already supports near-atomic detail.

        Resetting Previous Aberration Estimates

        The protocol provides options to reset previously estimated aberration
        terms before refinement. This functionality is useful when earlier
        processing stages introduced unstable or biologically implausible
        optical corrections. Resetting allows the refinement to begin again
        from neutral optical conditions.

        In practice, resetting parameters may help when combining datasets from
        different acquisition sessions, changing processing strategies, or
        troubleshooting refinement instability. However, unnecessary resets may
        discard valid corrections that were already improving the data quality.

        Outputs and Their Interpretation

        The main output of the protocol is a refined particle set containing
        updated optical and alignment-related metadata. These corrected
        particles are intended for subsequent refinement stages, where the
        benefits of improved aberration estimation become visible in the final
        map quality.

        Biologically, successful refinement often manifests as clearer
        secondary structure boundaries, improved side-chain density, sharper
        ligand features, and more stable local resolution estimates. The impact
        may appear subtle in intermediate maps but can become highly important
        during atomic model building and validation.

        Practical Recommendations

        In most biological workflows, it is advisable to perform global CTF
        refinement only after obtaining a reliable consensus reconstruction.
        Attempting aberration refinement too early may lead to unstable
        estimates because the structural reference does not yet contain enough
        high-resolution information.

        A good practical strategy is to begin with beam tilt and spherical
        aberration refinement enabled, then evaluate whether additional terms
        such as trefoil, tetrafoil, or anisotropic magnification provide
        measurable improvements. Careful inspection of FSC behavior and map
        interpretability remains essential.

        Datasets collected on well-aligned microscopes may show only modest
        gains, whereas older datasets or difficult acquisition conditions may
        benefit substantially from these corrections. Users pursuing
        publication-level reconstructions should generally consider at least
        one round of global CTF refinement as part of a modern high-resolution
        cryo-EM workflow.

        Final Perspective

        Global optical refinement is a critical component of modern cryo-EM
        image processing because it addresses subtle microscope imperfections
        that directly affect structural accuracy. Although these corrections
        may seem highly technical, they ultimately improve the biological
        interpretability of reconstructed maps. Careful selection of the
        reference volume, proper masking, and realistic refinement settings are
        essential for obtaining reliable and meaningful improvements in final
        reconstruction quality.
    """
    _label = 'global ctf refinement'
    _className = "ctf_refine_global"
    _protCompatibility = [V3_3_1, V3_3_2, V4_0_0, V4_0_1, V4_0_2, V4_0_3, V4_1_0,
                          V4_1_1, V4_1_2, V4_2_0, V4_2_1, V4_3_1, V4_4_0, V4_4_1, V4_5_1,
                          V4_5_3, V4_6_0, V4_6_1, V4_6_2, V4_7_0, V4_7_1]
    newParamsName = []

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
        form.addParam('refVolume', PointerParam, pointerClass='Volume',
                      important=True,
                      label="Input volume",
                      help='Provide a reference volume for global '
                           'CTF refinement.')
        form.addParam('refMask', PointerParam, pointerClass='VolumeMask',
                      label='Mask to be applied to this map',
                      important=True,
                      help="Provide a soft mask. if mask is present, use that, "
                           "otherwise use mask_refine if present, otherwise "
                           "fail")

        # -----------[Global CTF Refinement]------------------------
        form.addSection(label="Global CTF Refinement")
        form.addParam('crg_num_iters', IntParam, default=1,
                      validators=[Positive],
                      label='Number of iterations',
                      help='Number of times refinement of the various CTF '
                           'parameters is done. Using 2 or more iterations '
                           'allows changes in one parameter (eg. anisomag) to '
                           'affect the estimation of other parameters '
                           '(eg. tetrafoil).')
        form.addParam('crg_num_plots', IntParam, default=50,
                      validators=[Positive],
                      label='Num. groups to plot',
                      help='Number of exposure groups to make plots for. '
                           'After this many, stop plotting to save time.')

        form.addParam('crg_plot_binfactor', IntParam, default=1,
                      validators=[Positive],
                      expertLevel=LEVEL_ADVANCED,
                      label='Binning to apply to plots',
                      help='Binning makes it easier to see tilt/trefoil/tetrafoil '
                           'etc in the data plots, but does not change the '
                           'results')

        form.addParam('crg_min_res_A', FloatParam, default=10,
                      validators=[Positive],
                      label='Minimum Fit Res (A)',
                      help='The minimum resolution to use during refinement of '
                           'image aberrations.')

        form.addParam('crg_do_tilt', BooleanParam, default=True,
                      label="Fit Tilt",
                      help="Whether to fit beam tilt.")

        form.addParam('crg_do_trefoil', BooleanParam, default=True,
                      label="Fit Trefoil",
                      help="Whether to fit beam trefoil.")

        form.addParam('crg_do_spherical', BooleanParam, default=True,
                      label="Fit Spherical Aberration",
                      help="Whether to fit spherical aberration.")

        form.addParam('crg_do_tetrafoil', BooleanParam, default=True,
                      label="Fit Tetrafoil",
                      help="Whether to fit beam tetrafoil.")

        # new parameter to V3.3.1
        csVersion = getCryosparcVersion()
        if parse_version(csVersion) >= parse_version(V3_3_1):

            form.addParam('crg_do_anisomag', BooleanParam, default=False,
                          label="Fit Anisotropic Mag.",
                          help="Whether to fit beam anisotropic magnification.")

            form.addParam('crg_do_ews_correct', BooleanParam, default=False,
                          label="Account for EWS curvature",
                          expertLevel=LEVEL_ADVANCED,
                          help="Whether or not to correct for the curvature of "
                               "the Ewald Sphere")
            form.addParam('crg_ews_zsign', EnumParam,
                          choices=['positive', 'negative'],
                          expertLevel=LEVEL_ADVANCED,
                          default=0,
                          label="EWS curvature sign",
                          help='Whether to use positive or negative curvature in '
                               'Ewald Sphere correction.')

            form.addParam('ctf_reset_tilt', BooleanParam, default=False,
                          label="Reset Tilt to default",
                          expertLevel=LEVEL_ADVANCED,
                          help="Reset tilt and shift CTF parameters to 0 "
                               "before refining.")

            form.addParam('ctf_reset_trefoil', BooleanParam, default=False,
                          label="Reset Trefoil to default",
                          expertLevel=LEVEL_ADVANCED,
                          help="Reset trefoil CTF parameters to 0 before "
                               "refining.")

            form.addParam('ctf_reset_tetra', BooleanParam, default=False,
                          label="Reset Tetrafoil to default",
                          expertLevel=LEVEL_ADVANCED,
                          help="Reset tetrafoil CTF parameters to 0 "
                               "before refining.")

            form.addParam('ctf_reset_anisomag', BooleanParam, default=False,
                          label="Reset Anisotropic Magnification to default",
                          expertLevel=LEVEL_ADVANCED,
                          help="Reset anisotropic magnification parameters to "
                               "0 before refining.")

            self.newParamsName = ['ctf_reset_anisomag', 'ctf_reset_tetra',
                                  'ctf_reset_trefoil', 'crg_do_ews_correct',
                                  'crg_ews_zsign', 'crg_do_anisomag',
                                  'ctf_reset_tilt']


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
        self._paramsName = ['crg_num_iters',
                            'crg_num_plots',
                            'crg_plot_binfactor',
                            'crg_min_res_A',
                            'crg_do_tilt',
                            'crg_do_trefoil',
                            'crg_do_spherical',
                            'crg_do_tetrafoil',
                            'compute_use_ssd'] + self.newParamsName
        self.lane = str(self.getAttributeValue('compute_lane'))

    # --------------------------- STEPS functions ------------------------------
    def processStep(self):
        self.info(pwutils.yellowStr("Ctf Refinement started..."))
        self.doGlobalCtfRefinement()

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to Relion file
        """
        self._initializeUtilsVariables()
        outputStarFn = self._getFileName('out_particles')
        csOutputFolder = os.path.join(self.projectDir.get(),
                                      self.runGlobalCtfRefinement.get())
        csFileName = "particles.cs"

        # Copy the CS output particles to extra folder
        copyFiles(csOutputFolder, self._getExtraPath(), files=[csFileName])

        csFile = os.path.join(self._getExtraPath(), csFileName)

        argsList = [csFile, outputStarFn]

        convertCs2Star(argsList)

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
                         itemDataIterator=emtable.Table.iterRows(outImgsFn))

    def _createItemMatrix(self, particle, row):
        createItemMatrix(particle, row, align=ALIGN_PROJ)
        setCryosparcAttributes(particle, row,
                               RELIONCOLUMNS.rlnRandomSubset.value)

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
            summary.append("Number of Iterations: %s" %
                           self.crg_num_iters.get())
            summary.append("--------------------------------------------------")
            summary.append("Output particles %s" %
                           self.getObjectTag('outputParticles'))
        return summary

    def doGlobalCtfRefinement(self):
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
            if paramName != 'crg_ews_zsign':
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            else:
                params[str(paramName)] = str(EWS_CURVATURE_SIGN[self.crg_ews_zsign.get()])

        # Determinate the GPUs to use (in dependence of
        # the cryosparc version)
        try:
            if not self.useQueueForSteps() and not self.useQueue():  # not using queue system
                gpusToUse = self.getGpuList()
            else:  # using queue system
                gpusToUse = False
        except Exception:
            gpusToUse = False

        runGlobalCtfRefinementJob = enqueueJob(self._className, self.projectName.get(),
                                                 self.workSpaceName.get(),
                                                 str(params).replace('\'', '"'),
                                                 str(input_group_connect).replace('\'', '"'),
                                                 self.lane, gpusToUse,
                                                 result_connect=input_result_connect)

        self.runGlobalCtfRefinement = String(runGlobalCtfRefinementJob)
        self.currenJob.set(runGlobalCtfRefinementJob.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(), self.runGlobalCtfRefinement.get(),
                         "An error occurred in the particles subtraction process. "
                         "Please, go to cryoSPARC software for more "
                         "details.", self)