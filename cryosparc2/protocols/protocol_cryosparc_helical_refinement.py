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
                     cryosparcValidate, gpusValidate)
from ..constants import *


class ProtCryoSparcHelicalRefine3D(ProtCryoSparc3DHomogeneousRefine):
    """ Reconstruct and refine a homogeneous helical assembly, with or without
    imposition and refinement of symmetry parameters. Helical Refinement (BETA)
    uses an algorithm that is conceptually similar to Egelman's Iterative
    Helical Real Space Reconstruction (IHRSR) algorithm, while incorporating
    the same maximum likelihood framework, accelerated branch-and-bound
    alignment algorithm, and optional Non-Uniform regularization as used in
    other cryoSPARC refinement jobs.
    """
    """
        Performs homogeneous 3D refinement of helical assemblies using cryoSPARC’s
        helical refinement framework. The protocol is designed for the reconstruction
        and optimization of filamentous or helical structures while optionally applying
        symmetry constraints and Non-Uniform refinement strategies. The implementation
        follows concepts similar to the IHRSR methodology while integrating cryoSPARC’s
        accelerated alignment and maximum-likelihood refinement algorithms.

        AI Generated:

        Helical Refinement 3D (ProtCryoSparcHelicalRefine3D) — User Manual
            Overview

            The Helical Refinement protocol reconstructs and refines helical particles
            from cryo-EM datasets using cryoSPARC’s dedicated helical refinement engine.
            Its primary objective is to optimize the alignment and reconstruction of
            filamentous assemblies while preserving biologically meaningful helical
            symmetry. The protocol is especially useful for structures such as amyloid
            fibrils, cytoskeletal filaments, membrane-associated helices, and other
            repetitive macromolecular assemblies.

            In practical cryo-EM workflows, this protocol is commonly used after particle
            extraction and initial model generation. It refines the angular alignment,
            translational positioning, and helical symmetry parameters in order to improve
            map resolution and structural interpretability. The refinement process can
            optionally incorporate Non-Uniform regularization to improve local map quality
            in flexible or heterogeneous regions.

            Inputs and General Workflow

            The protocol requires a set of input particles representing filament segments.
            An initial reference volume can optionally be provided to guide the refinement
            process. If no reference volume is available, the protocol supports generation
            of a cylindrical initial model, which is particularly useful during early-stage
            refinement of filamentous datasets.

            An optional soft mask may also be supplied. The mask defines which regions of
            the reconstruction contribute most strongly during refinement and can improve
            convergence by reducing solvent noise and excluding unstable peripheral regions.
            For biological assemblies with flexible domains or variable density regions,
            careful masking frequently improves the robustness of refinement.

            Helical Parameters and Symmetry

            One of the central components of this protocol is the definition and refinement
            of helical geometry. The user may define the initial helical twist and rise,
            corresponding respectively to the angular rotation and translational displacement
            between adjacent subunits along the filament axis. Positive and negative twist
            values represent right-handed and left-handed helices.

            The protocol also supports limiting particle shifts along the helical axis.
            This option is useful for relatively rigid filaments because it constrains the
            alignment search space and improves consistency between symmetry-related regions.
            However, for flexible helices, excessive constraints may reduce refinement quality.

            Helical symmetry can additionally be enforced during reconstruction through the
            helical symmetry order parameter. This controls how many symmetry-related copies
            are incorporated during backprojection and reconstruction. Lower symmetry orders
            reduce averaging effects, while larger values increase symmetry enforcement and
            may improve signal-to-noise ratio for highly regular filaments.

            Point Group Symmetry

            Besides helical symmetry, the protocol supports cyclic and dihedral point-group
            symmetry. These symmetry definitions are integrated into the refinement process
            and are particularly important for assemblies with additional rotational symmetry
            components beyond the helical arrangement itself.

            Correct symmetry assignment is biologically important because inappropriate
            symmetry enforcement may introduce structural artifacts or artificially distort
            asymmetric regions. In most practical workflows, symmetry parameters are selected
            based on prior structural knowledge or previous reconstruction experiments.

            Initial Model Generation

            The protocol includes tools for generating an initial cylindrical density model.
            This option is especially useful when no reliable starting reconstruction exists.
            The cylindrical model can be customized through parameters such as filament outer
            diameter, inner diameter, and padding distance. These values determine the basic
            geometric properties of the initial density distribution.

            Initial lowpass filtering is also available and is commonly used to suppress
            high-frequency noise during the early refinement stages. In practice, moderate
            lowpass filtering improves convergence stability and reduces alignment bias
            toward noisy features.

            Refinement Strategy

            During refinement, the protocol iteratively optimizes particle orientations,
            translational shifts, and symmetry-related parameters against the evolving 3D
            reconstruction. The refinement process can optionally use Non-Uniform refinement,
            which applies adaptive regularization to improve local density quality and
            reduce overfitting in heterogeneous regions.

            The protocol also supports control over alignment resolution limits and GS-FSC
            split resolution thresholds. These parameters influence how resolution-dependent
            information contributes during refinement and can significantly affect the
            balance between convergence speed and reconstruction quality.

            Dynamic and static masking strategies are available during refinement. Dynamic
            masks adapt automatically to the evolving density map, while static masks apply
            fixed masking regions. Dynamic masking is generally preferred for flexible or
            partially disordered assemblies because it better follows structural variability
            during iterative refinement.

            Validation and Parameter Control

            The implementation performs several validation checks before execution. The
            protocol verifies GPU compatibility, validates the presence of required input
            data, and ensures that cylindrical model parameters are correctly defined when
            no initial reference map is provided.

            During execution, refinement parameters are collected and converted into the
            format expected by cryoSPARC jobs. The protocol then launches the corresponding
            cryoSPARC refinement task, manages GPU allocation, monitors job execution, and
            waits for successful completion before continuing.

            Outputs and Interpretation

            After completion, the protocol generates refined helical reconstructions together
            with updated alignment parameters for all particles. These outputs represent the
            optimized spatial arrangement of filament segments within the reconstructed
            helical assembly.

            The resulting maps can be used for downstream structural interpretation,
            atomic model building, flexibility analysis, or comparative structural studies.
            Biologically, the final reconstruction quality strongly depends on the accuracy
            of the initial symmetry estimates, masking strategy, and refinement constraints.

            Practical Recommendations

            In routine cryo-EM workflows, refinement is often started using moderate
            lowpass filtering and conservative symmetry assumptions. If convergence is
            unstable, introducing a carefully designed soft mask or enabling Non-Uniform
            refinement frequently improves results.

            For highly regular filaments, enforcing helical symmetry and limiting axial
            shifts can substantially improve resolution. Conversely, flexible or polymorphic
            assemblies may require weaker symmetry enforcement and broader alignment freedom
            to preserve biologically relevant variability.

            When no reliable reference map exists, the cylindrical initial model provides
            a practical starting point, although refinement quality should always be
            evaluated visually and through independent resolution metrics.

            Final Perspective

            Helical refinement is a biologically significant stage in cryo-EM analysis
            because the imposed symmetry and alignment strategy directly influence the
            interpretability of filamentous assemblies. Reliable results depend on accurate
            estimation of helical parameters, appropriate masking, and careful refinement
            configuration. When applied correctly, the protocol enables high-resolution
            structural characterization of complex helical systems while preserving the
            essential biological organization of the assembly.
        """
    _label = '3D helical refinement'
    _fscColumns = 4
    _protCompatibility = [V3_3_1, V3_3_2, V4_0_0, V4_0_1, V4_0_2, V4_0_3,
                          V4_1_0, V4_1_1, V4_1_2, V4_2_0, V4_2_1, V4_3_1, V4_4_0, V4_4_1, V4_5_1,
                          V4_5_3, V4_6_0, V4_6_1, V4_6_2, V4_7_0, V4_7_1]
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
                      choices=["dynamic", "static"],
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
                    paramName != 'refine_res_align_max'):
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif paramName == 'refine_pg_symmetry':
                symetryValue = getSymmetry(self.symmetryGroup.get(),
                                           self.symmetryOrder.get())
                params[str(paramName)] = symetryValue
            elif self.getAttributeValue(paramName) is not None and float(self.getAttributeValue(paramName)) > 0:
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