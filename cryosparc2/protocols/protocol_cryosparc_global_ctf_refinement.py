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
from pkg_resources import parse_version

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
from ..utils import (addComputeSectionParams, addPreprocessLaneParam, cryosparcValidate, gpusValidate,
                     enqueueJob, waitForCryosparc, copyFiles,
                     getCryosparcVersion)

from ..constants import *


class ProtCryoSparcGlobalCtfRefinement(ProtCryosparcBase, pwprot.ProtParticles):
    """
    Wrapper protocol for the Cryosparc's per-particle Global CTF refinement.
    Performs per-exposure-group CTF parameter refinement of higher-order
    aberrations, against a given 3D reference
    """

    """
        Global CTF Refinement (ProtCryoSparcGlobalCtfRefinement) — User Manual

        Overview

        The Global CTF Refinement protocol performs per-exposure-group refinement
        of higher-order optical aberrations using CryoSPARC refinement strategies.
        The protocol operates on an existing set of aligned particles together with
        a reference 3D volume and improves the accuracy of Contrast Transfer Function
        (CTF) parameters across the dataset. The main objective is to correct optical
        imperfections that may limit high-resolution reconstruction quality.

        In practical cryo-EM workflows, this refinement stage is commonly applied
        after obtaining a reliable consensus reconstruction and accurate particle
        alignments. By refining aberration parameters such as beam tilt, trefoil,
        tetrafoil, spherical aberration, or anisotropic magnification, the protocol
        improves consistency between experimental images and the reference projection
        model, which can significantly enhance downstream map resolution and structural
        interpretability.

        Inputs and General Workflow

        The protocol requires three main inputs: a set of particles with projection
        alignment information, a reference volume, and a soft mask associated with
        the reference map. The particles must already contain valid alignment
        parameters because the refinement process depends on comparing experimental
        particle images against projections derived from the reference structure.

        The reference volume defines the structural model used during refinement.
        Biologically, the quality of this reference strongly influences the quality
        of the refined optical parameters. A high-resolution and well-converged map
        generally produces more stable refinement results.

        The soft mask is also an important component of the workflow because it
        restricts the refinement to biologically relevant regions of the structure.
        Applying an appropriate mask reduces solvent influence and minimizes noise
        contributions from flexible or poorly resolved regions.

        During execution, the protocol initializes a CryoSPARC project environment,
        converts the input metadata into CryoSPARC-compatible formats, launches the
        refinement job, waits for completion, and finally converts the resulting
        particle metadata back into STAR format for Scipion integration.

        Higher-Order Aberration Refinement

        The protocol provides refinement of several higher-order aberration terms.
        Beam tilt correction compensates for systematic angular deviations in the
        electron beam that can introduce phase errors into the reconstruction.
        Trefoil and tetrafoil refinement model more complex optical distortions
        originating from microscope imperfections. Spherical aberration refinement
        further improves phase consistency at high spatial frequencies.

        For CryoSPARC versions supporting anisotropic magnification refinement,
        the protocol can additionally estimate directional magnification distortions.
        This correction becomes particularly relevant in high-resolution datasets
        where small calibration inaccuracies may otherwise reduce reconstruction quality.

        The protocol also supports Ewald Sphere curvature correction. This option
        becomes increasingly important for large particles or near-atomic-resolution
        reconstructions where curvature effects can no longer be neglected. The user
        may define whether positive or negative curvature should be applied during
        correction.

        Iterative Refinement Strategy

        The refinement process can be repeated for multiple iterations. Performing
        several iterations allows corrections estimated in one cycle to influence
        the estimation of other aberration parameters in subsequent cycles. For
        example, anisotropic magnification refinement may improve the stability of
        trefoil or tetrafoil estimation during later iterations.

        The protocol also provides reset options for tilt, trefoil, tetrafoil,
        and anisotropic magnification parameters. These settings restore selected
        aberration terms to their default values before refinement begins. Such
        resets can be useful when previous refinement attempts introduced unstable
        or biologically implausible parameter estimates.

        Resolution and Plotting Parameters

        The minimum fitting resolution parameter determines the lowest spatial
        frequency considered during aberration estimation. Restricting refinement
        to higher-resolution information can improve sensitivity to subtle optical
        effects, although excessively aggressive settings may reduce robustness
        in noisy datasets.

        The protocol additionally allows generation of diagnostic plots for a
        selected number of exposure groups. These plots help evaluate refinement
        quality and visualize optical distortions across the dataset. Optional
        plot binning improves visual interpretation of trefoil, tilt, and related
        aberration patterns without altering the refinement results themselves.

        GPU and CryoSPARC Integration

        The implementation automatically integrates with CryoSPARC job management.
        Depending on the execution environment, GPU resources are assigned either
        directly or through a queue system. The protocol constructs parameter
        dictionaries dynamically, launches the CryoSPARC refinement job, monitors
        execution, and waits until the refinement process completes successfully.

        If the input reference contains half maps, the protocol automatically links
        them to the CryoSPARC refinement job. This behavior ensures compatibility
        with advanced refinement workflows that rely on independent half-map
        information.

        Outputs and Interpretation

        After completion, the protocol generates a refined set of particles with
        updated projection alignment and refined optical parameters. The output
        particles preserve the original dataset structure while incorporating the
        corrected CTF information estimated during refinement.

        Internally, CryoSPARC output metadata is converted into STAR format and
        copied into the Scipion project structure. Particle transformation matrices
        and CryoSPARC-specific attributes are reconstructed and assigned to each
        particle object to maintain compatibility with downstream processing steps.

        From a biological perspective, successful global CTF refinement often leads
        to improved map sharpness, enhanced high-resolution features, and better
        interpretability of structural details such as side chains, ligand density,
        or secondary-structure elements.

        Validation and Compatibility

        The protocol validates CryoSPARC installation compatibility, GPU availability,
        and consistency between input particles and the reference volume dimensions.
        Only particles with projection alignment information are accepted because
        accurate alignment parameters are essential for reliable optical refinement.

        The implementation is compatible with multiple CryoSPARC versions ranging
        from v3.3.1 up to v4.7.1, with additional refinement features automatically
        enabled depending on the detected software version.

        Practical Recommendations

        In routine cryo-EM processing, global CTF refinement is typically most useful
        after achieving stable particle alignments and an accurate consensus map.
        Applying the protocol too early in the workflow may produce unstable results
        because aberration estimation depends strongly on the quality of the reference.

        Beam tilt and spherical aberration refinement are usually safe starting
        points for most datasets. Trefoil, tetrafoil, and anisotropic magnification
        refinement become increasingly important for high-resolution studies where
        subtle optical distortions limit reconstruction quality.

        Careful mask design remains essential. Masks should include the stable and
        biologically meaningful regions of the structure while excluding solvent and
        highly flexible domains. Poor masking may bias aberration estimation and
        reduce refinement stability.

        Final Perspective

        Global CTF refinement is a critical high-resolution optimization step in
        modern cryo-EM workflows. Beyond simple parameter correction, it improves
        the physical consistency between experimental images and the reconstruction
        model, enabling more accurate structural interpretation. When combined with
        reliable alignments, appropriate masking, and a high-quality reference map,
        this protocol can substantially improve the final reconstruction quality and
        enhance confidence in downstream biological conclusions.
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
        addPreprocessLaneParam(form)

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