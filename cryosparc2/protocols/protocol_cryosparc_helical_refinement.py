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

from pwem import SCIPION_SYM_NAME

import pyworkflow.utils as pwutils
from pyworkflow.object import String
from pyworkflow.protocol.params import (FloatParam, Positive, IntParam,
                                        BooleanParam, EnumParam, PointerParam)

from .protocol_cryosparc_homogeneous_refine import ProtCryoSparc3DHomogeneousRefine
from ..utils import (getSymmetry, enqueueJob, waitForCryosparc,
                     clearIntermediateResults, addComputeSectionParams, addPreprocessLaneParam,
                     getVersionedEnumValue, cryosparcValidate, gpusValidate)
from ..constants import *


class ProtCryoSparcHelicalRefine3D(ProtCryoSparc3DHomogeneousRefine):
    """
    Refines helical cryo-EM structures by reconstructing filamentous or
    helical assemblies while optionally optimizing helical symmetry
    parameters and improving map quality through iterative refinement.

    AI Generated:

    Helical Refinement (ProtCryoSparcHelicalRefine3D) — User Manual
        Overview

        The Helical Refinement protocol is designed for the reconstruction
        and refinement of filamentous or helical biological assemblies in
        cryo-electron microscopy workflows. Typical applications include
        actin filaments, microtubules, amyloid fibrils, viral helices,
        nucleoprotein assemblies, and other repeating biological polymers
        that exhibit helical organization.

        The protocol combines particle alignment, three-dimensional
        reconstruction, and symmetry-guided refinement into a unified
        workflow capable of producing high-resolution maps from segmented
        filament particles. It is especially useful when the biological
        system contains repeating subunits arranged along a helical axis,
        where exploiting the inherent symmetry can significantly improve
        signal-to-noise ratio and structural resolution.

        In practical cryo-EM analysis, this protocol is often applied after
        particle extraction and preliminary processing steps have already
        generated aligned filament segments or an approximate initial model.
        The refinement process progressively improves both the particle
        alignments and the reconstructed density while enforcing helical
        constraints consistent with the underlying biological assembly.

        Inputs and Initial Requirements

        The protocol requires a set of input particles corresponding to
        extracted filament segments. These particles should ideally contain
        reliable contrast transfer function information and sufficient
        angular diversity for stable refinement.

        An initial three-dimensional reference volume can optionally be
        provided. In most biological workflows, supplying a reasonable
        starting map greatly improves convergence and reduces the risk of
        refinement bias or instability. Suitable initial references often
        originate from ab initio reconstruction, previous refinements, or
        related structures.

        When no initial reference is available, the protocol can generate a
        simplified cylindrical starting model. This approach is particularly
        useful for long helical assemblies with approximately uniform
        diameters. However, cylindrical models provide only coarse initial
        information and are generally less reliable than experimentally
        derived references.

        Helical Symmetry Parameters

        A central aspect of helical reconstruction is the definition of the
        helical rise and twist. The rise corresponds to the translational
        displacement between adjacent subunits along the filament axis,
        whereas the twist defines the angular rotation between neighboring
        asymmetric units.

        Accurate estimates of these parameters are biologically important
        because they determine how the repeating units are organized within
        the filament. Incorrect values may produce blurred densities,
        distorted symmetry relationships, or convergence toward incorrect
        structural solutions.

        Positive and negative twist values correspond to right-handed and
        left-handed helices, respectively. Users should therefore verify the
        handedness of the biological assembly whenever possible using
        complementary experimental evidence or previously established
        structural information.

        The protocol also supports limiting the maximum helical symmetry
        order during reconstruction. This parameter determines how many
        neighboring asymmetric units contribute during symmetry averaging.
        For rigid and highly ordered helices, larger symmetry expansion may
        improve map quality substantially. In flexible or heterogeneous
        assemblies, however, excessive symmetry averaging may introduce
        artifacts or smear biologically relevant variability.

        Point Group Symmetry

        In addition to helical symmetry, the protocol supports cyclic and
        dihedral point group symmetries. These are commonly encountered in
        decorated filaments, tubular assemblies, and other higher-order
        helical systems.

        Applying point group symmetry can improve reconstruction quality
        when the biological specimen truly obeys the assumed symmetry.
        However, incorrect symmetry assignment can artificially distort the
        reconstructed density and hide meaningful asymmetry. Biological
        validation of the imposed symmetry is therefore essential before
        interpreting fine structural details.

        Real-Space Symmetry Enforcement

        The protocol allows symmetry enforcement in real space during the
        refinement process. This strategy can improve alignment stability
        and convergence, especially in noisy datasets or during early
        refinement stages.

        From a biological perspective, enforcing symmetry at relatively low
        resolution may help stabilize the reconstruction while preserving
        flexibility at higher resolution. Excessively aggressive symmetry
        enforcement, however, may suppress genuine structural heterogeneity
        or conformational variability that could be functionally important.

        Non-Uniform Refinement

        An optional non-uniform refinement strategy can be enabled to
        improve local map quality and enhance high-resolution structural
        features. This approach is especially beneficial for biological
        assemblies containing flexible domains, variable occupancy, or
        regions with uneven local resolution.

        In practice, non-uniform refinement often improves interpretability
        of side chains, secondary structure elements, and flexible filament
        interfaces. Nevertheless, users should carefully validate the final
        maps to distinguish genuine structural features from refinement
        artifacts.

        Initial Model Generation

        When generating a cylindrical initial model, the filament geometry
        becomes particularly important. Parameters such as outer diameter,
        inner diameter, and padding distances define the approximate
        physical dimensions of the starting density.

        Biologically realistic diameter estimates are critical because they
        influence the early alignment process and determine how rapidly the
        refinement converges toward a meaningful structure. Overly narrow or
        excessively large cylinders may delay convergence or bias the
        resulting reconstruction.

        The initial lowpass filtering stage is also biologically important.
        Starting with low-resolution information reduces the influence of
        noise and helps prevent overfitting during the earliest refinement
        iterations. This is particularly valuable for challenging filament
        datasets with preferred orientation or limited contrast.

        Alignment and Masking Strategies

        The protocol provides several controls for alignment resolution and
        masking behavior. Restricting the maximum alignment resolution may
        improve stability during difficult refinements, especially when the
        current reconstruction quality is still limited.

        Dynamic masking automatically adapts to the evolving density during
        refinement and is often the preferred option for flexible or
        partially disordered helices. Static masking, by contrast, provides
        a fixed region of interest and may be preferable for highly stable
        filament cores or well-characterized assemblies.

        Careful masking is biologically important because it determines
        which regions dominate the alignment process. Masks that exclude
        disordered solvent regions while preserving the structurally stable
        filament core usually produce more reliable results.

        Outputs and Interpretation

        The protocol produces refined three-dimensional maps together with
        updated particle alignments and associated refinement statistics.
        The final reconstruction represents the consensus structure of the
        helical assembly under the imposed symmetry assumptions.

        The resulting maps can be used for atomic modeling, structural
        interpretation, variability analysis, or downstream comparative
        studies. As with all symmetry-based reconstructions, biological
        interpretation should consider the possibility that local
        heterogeneity or flexibility may be partially averaged during the
        reconstruction process.

        Practical Recommendations

        In routine cryo-EM workflows, it is generally advisable to begin
        with conservative symmetry assumptions and gradually introduce more
        aggressive symmetry enforcement only after stable refinement has
        been achieved. Reliable initial estimates of twist and rise often
        determine whether refinement converges successfully.

        For flexible filaments, limiting symmetry averaging and enabling
        non-uniform refinement usually improves interpretability. For rigid
        and highly ordered helices, stronger symmetry enforcement may
        significantly enhance achievable resolution.

        When using cylindrical initial models, users should carefully
        inspect intermediate reconstructions to ensure that biologically
        meaningful features are emerging and that the refinement has not
        converged toward an incorrect helical solution.

        Final Perspective

        Helical refinement is one of the most biologically specialized
        stages in cryo-EM image processing because it combines structural
        averaging with strict geometric constraints imposed by filament
        organization. Successful reconstruction depends not only on
        computational refinement but also on accurate biological knowledge
        of filament architecture, symmetry, flexibility, and handedness.

        Careful selection of symmetry parameters, realistic initial models,
        and appropriate masking strategies are essential for obtaining
        reliable and biologically interpretable helical reconstructions.
    """
    _label = '3D helical refinement'
    _fscColumns = 4
    _className = "helix_refine"

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      label="Input particles", important=True,
                      validators=[Positive],
                      help='Particle stacks to use. Multiple stacks will '
                           'be concatenated.')
        form.addParam('referenceVolume', PointerParam, pointerClass='Volume',
                      default=None,
                      allowsNull=True,
                      label="Initial volume",
                      help='Initial volume to use for helical refinement.')
        form.addParam('refMask', PointerParam, pointerClass='VolumeMask',
                      default=None,
                      label='Mask to be applied to this map(Optional)',
                      allowsNull=True,
                      help='Initial volume mask raw data.')

        form.addSection(label='Helical Refinement')

        form.addParam('refine_init_twist', FloatParam,
                      default=None,
                      allowsNull=True,
                      label="Helical twist estimate (degrees)",
                      help='Angular distance between adjacent subunits, in the '
                           'range (-180,180). Positive and negative values '
                           'correspond to right and left-handed helices, '
                           'respectively.')

        form.addParam('refine_init_shift', FloatParam,
                      default=None,
                      allowsNull=True,
                      label="Helical rise estimate (A)",
                      help='Positive non-zero translation distance (along '
                           'helical axis) between adjacent subunits.')

        form.addParam('refine_limit_shifts', BooleanParam,
                      default=True,
                      label="Limit shifts along the helical axis",
                      help='Limit alignment shifts along the meridian to +/- '
                           '0.5 helical rises, once all particles are seen. '
                           'For fairly rigid helices, can improve resolution'
                           ' by ensuring symmetry averaging occurs over only '
                           'central asymmetric units. Not recommended for '
                           'helices with significant flexibility.')

        form.addParam('refine_hsym_order', FloatParam,
                      default=None,
                      allowsNull=True,
                      label="Maximum symmetry order to apply during reconstruction",
                      help='The maximum amount of helical symmetry to impose '
                           'during reconstruction; i.e., the maximum number of '
                           '(twist, rise) pairs to backproject each particle '
                           'image with. For particles picked outside of the '
                           'filament tracer or template picker, this should be '
                           'set to the distance between extracted boxes, '
                           'divided by the helical rise. If left as None, '
                           'will be calculated based on the inter-box distance '
                           '(if available). Set to 1 for no helical Fourier'
                           ' space symmetrization.')

        form.addParam('refine_sym_enforce_r', FloatParam,
                      default=8,
                      validators=[Positive],
                      label="Resolution to begin real-space symmetrization",
                      help='At what resolution (A) to begin enforce symmetry '
                           'in real-space prior to alignment; may improve'
                           ' particle alignments. Set to 0 to disable '
                           'real-space symmetry enforcement. If symmetry '
                           'parameters are being searched, it\'s recommended '
                           'to set this to a fairly high resolution '
                           '(e.g. ~5 A).')

        form.addParam('symmetryGroup', EnumParam,
                      choices=[CS_SYM_NAME[SYM_CYCLIC] +
                               " (" + SCIPION_SYM_NAME[SYM_CYCLIC] + ")",
                               CS_SYM_NAME[SYM_DIHEDRAL_Y] +
                               " (" + SCIPION_SYM_NAME[SYM_DIHEDRAL_Y] + ")"],
                      default=SYM_CYCLIC,
                      label="Symmetry",
                      help="Symmetry String (C, D). E.g. C1, D7, C4, etc. "
                           "Only cyclic and dihedral symmetries are supported."
                      )

        form.addParam('symmetryOrder', IntParam, default=1,
                      condition='symmetryGroup==%d or symmetryGroup==%d' %
                                (SYM_DIHEDRAL_Y - 1, SYM_CYCLIC),
                      label='Point group symmetry',
                      validators=[Positive],
                      help='Order of symmetry.')

        form.addSection(label='Non-Uniform Refinement')

        form.addParam('nu_refine', BooleanParam,
                      default=False,
                      label="Use Non-Uniform Refinement?",
                      help='Use Non-Uniform regularization during refinement '
                           'to achieve higher resolution and map quality.')

        form.addSection(label='Initial Model')

        form.addParam('refine_res_init', FloatParam,
                      default=20,
                      validators=[Positive],
                      label="Initial lowpass resolution (A)",
                      help='Lowpass filter resolution applied to input '
                           'structure. Values between 15 and 35 Angstroms '
                           'may produce best results.')

        form.addParam('refine_initmodel_numimages', FloatParam,
                      default=5000,
                      validators=[Positive],
                      label="Number of images for initial density generation",
                      help='Number of images used in generating the initial '
                           'density.')

        form.addParam('use_cylindrical_model', BooleanParam,
                      default=False,
                      label="Generate a cylindrical initial model?",
                      help='Whether or not to generate a cylindrical initial model')

        form.addParam('filament_outer_diameter', FloatParam,
                      default=None,
                      allowsNull=True,
                      label="Filament Outer Diameter (Angstrom)",
                      help='Approximate outer diameter of the filament in Angstroms')

        form.addParam('filament_inner_diameter', FloatParam,
                      default=0,
                      label="Filament Inner Diameter (Angstrom)",
                      help='Approximate inner diameter of the filament in Angstroms')

        form.addParam('filament_far_dist_A', FloatParam,
                      default=6,
                      label="Far distance (Angstrom)",
                      help='Distance over which the model is padded, with voxel values fading to 0.')

        form.addSection(label='Refinement')

        form.addParam('refine_res_align_max', FloatParam,
                      default=None,
                      allowsNull=True,
                      label="Maximum align resolution (A)",
                      help='Manual override for maximum resolution that is '
                           'used for alignment. This value is normally '
                           'set by the GS-FSC')

        form.addParam('refine_res_gsfsc_split', FloatParam,
                      default=20,
                      validators=[Positive],
                      label="GSFSC split resolution (A)",
                      help='Resolution beyond which two GS-FSC halves are '
                           'independent')

        form.addParam('refine_mask', EnumParam,
                      choices=["static", "dynamic"],
                      default=0,
                      label="Mask",
                      help='Type of masking to use. Either "dynamic", or '
                           '"static".')

        form.addParam('refine_dynamic_mask_thresh_factor', FloatParam,
                      default=0.3,
                      validators=[Positive],
                      label="Dynamic mask threshold (0-1)",
                      help='Level set threshold for selecting regions that are '
                           'included in the dynamic mask.')

        # --------------[Compute settings]---------------------------
        form.addSection(label="Compute settings")
        addComputeSectionParams(form, allowMultipleGPUs=False)
        addPreprocessLaneParam(form)


    def _insertAllSteps(self):
        ProtCryoSparc3DHomogeneousRefine._insertAllSteps(self)

    def processStep(self):
        self.info(pwutils.yellowStr("Refinement started..."))
        self.doRunRefine()

    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        validateMsgs = cryosparcValidate()
        if not validateMsgs:
            validateMsgs = gpusValidate(self.getGpuList(), checkSingleGPU=True)
            if not validateMsgs:
                if self.referenceVolume.get() is None and not self.use_cylindrical_model.get():
                    validateMsgs.append("Cannot generate initial model "
                                        "without in-plane rotation "
                                        "information. Please input an "
                                        "initial model from a previous ab-initio or "
                                        "refinement protocol, or activate the "
                                        "'Generate a cylindrical initial "
                                        "model?' parameter")

                if self.use_cylindrical_model.get() and self.filament_outer_diameter.get() is None:
                    validateMsgs.append("Must set the filament outer diameter to use a cylindrical model")

        return validateMsgs

    def _defineParamsName(self):
        """ Define a list with all protocol parameters names"""
        self._paramsName = ['refine_init_twist',
                            'refine_init_shift',
                            'refine_hsym_order',
                            'refine_limit_shifts',
                            'refine_sym_enforce_r',
                            'refine_pg_symmetry',
                            'nu_refine',
                            'refine_res_init',
                            'refine_initmodel_numimages',
                            'use_cylindrical_model',
                            'refine_res_align_max',
                            'refine_res_gsfsc_split',
                            'refine_mask',
                            'refine_dynamic_mask_thresh_factor',
                            'filament_outer_diameter',
                            'filament_inner_diameter',
                            'filament_far_dist_A',
                            'compute_use_ssd']
        self.lane = str(self.getAttributeValue('compute_lane'))

    def doRunRefine(self):
        input_group_connect = {"particles": self.particles.get()}
        if self.volume.get() is not None:
            input_group_connect["volume"] = self.volume.get()
        if self.mask.get() is not None:
            input_group_connect["mask"] = self.mask.get()
        params = {}

        for paramName in self._paramsName:
            if (paramName != 'refine_hsym_order' and
                    paramName != 'refine_pg_symmetry' and
                    paramName != "refine_init_twist" and
                    paramName != "refine_init_shift" and
                    paramName != 'filament_outer_diameter' and
                    paramName != 'refine_mask' and
                    paramName != 'refine_res_align_max'):
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif paramName == 'refine_pg_symmetry':
                symetryValue = getSymmetry(self.symmetryGroup.get(),
                                           self.symmetryOrder.get())
                params[str(paramName)] = symetryValue
            elif self.getAttributeValue(paramName) is not None and float(self.getAttributeValue(paramName)) > 0:
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif paramName == 'refine_mask':
                params[str(paramName)] = str(
                    getVersionedEnumValue(
                        self.refine_mask.get(),
                        HELIXREFINE_MASK_CHOICES,
                        HELIXREFINE_MASK_CHOICES_V5
                    )
                )


        # Determinate the GPUs to use (in dependence of
        # the cryosparc version)
        try:
            if not self.useQueueForSteps() and not self.useQueue():  # not using queue system
                gpusToUse = self.getGpuList()
            else:  # using queue system
                gpusToUse = False
        except Exception:
            gpusToUse = False

        runRefineJob = enqueueJob(self._className, self.projectName.get(),
                                    self.workSpaceName.get(),
                                    str(params).replace('\'', '"'),
                                    str(input_group_connect).replace('\'', '"'),
                                    self.lane, gpusToUse)

        self.runRefine = String(runRefineJob.get())
        self.currenJob.set(runRefineJob.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(), self.runRefine.get(),
                         "An error occurred in the Refinement process. "
                         "Please, go to cryoSPARC software for more "
                         "details.", self)
        clearIntermediateResults(self.projectName.get(), self.runRefine.get())