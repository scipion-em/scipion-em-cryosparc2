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
                     clearIntermediateResults, addComputeSectionParams,
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
    ProtCryoSparcHelicalRefine3D — Helical 3D Refinement Protocol

    Overview
    --------
    Performs high-resolution refinement of helical assemblies using cryoSPARC’s
    helical refinement framework. The protocol is designed for filamentous or
    helical biological structures and combines iterative helical reconstruction
    strategies with cryoSPARC’s maximum-likelihood refinement algorithms.

    The method incorporates:
        - Helical symmetry refinement.
        - Branch-and-bound alignment optimization.
        - Optional Non-Uniform (NU) regularization.
        - Real-space and Fourier-space symmetry enforcement.
        - Flexible refinement of twist and rise parameters.

    Typical applications include:
        - Refinement of amyloid fibrils.
        - Helical virus reconstruction.
        - Cytoskeletal filament analysis.
        - High-resolution reconstruction of filamentous protein assemblies.

    Inputs and Workflow
    -------------------
    - Input Particles:
        Particle stacks corresponding to extracted helical segments.
        Multiple particle sets may be concatenated.

    - Initial Volume:
        A starting 3D reconstruction used to initialize refinement.
        Recommended for stable convergence and orientation estimation.

    - Optional Reference Mask:
        Soft mask applied during refinement to focus optimization on
        biologically relevant regions and reduce solvent noise.

    Workflow summary:
        1. Load particle stacks and optional reference volume.
        2. Configure helical symmetry and refinement parameters.
        3. Launch cryoSPARC helical refinement job.
        4. Monitor execution and GPU allocation.
        5. Generate refined helical reconstruction outputs.

    Helical Symmetry Parameters
    ---------------------------
    The protocol supports explicit refinement of helical geometry:

    - Helical Twist:
        Angular rotation between adjacent subunits.
        Positive and negative values represent right-handed and
        left-handed helices respectively.

    - Helical Rise:
        Axial translation between neighboring asymmetric units.

    - Symmetry Order:
        Defines the amount of symmetry applied during reconstruction.
        Higher values increase averaging but may introduce artifacts if
        symmetry assumptions are incorrect.

    - Point Group Symmetry:
        Supports cyclic (Cn) and dihedral (Dn) symmetry groups.

    Best practices:
        * Use experimentally estimated twist/rise values whenever possible.
        * Incorrect symmetry values may reduce resolution or distort features.
        * Flexible helices may require reduced symmetry enforcement.

    Real-Space Symmetrization
    -------------------------
    The protocol can enforce symmetry directly in real space before alignment.

    Features:
        - Improves particle alignment consistency.
        - Enhances signal for rigid helical assemblies.
        - Can stabilize refinement at intermediate resolutions.

    Recommendations:
        * Use moderate enforcement resolution (~5–8 Å) for stable helices.
        * Disable or reduce enforcement for highly flexible filaments.

    Initial Model Generation
    ------------------------
    Two initialization strategies are supported:

    - External Initial Volume:
        Uses a previously reconstructed map from ab-initio or refinement.

    - Cylindrical Initial Model:
        Generates a simple cylindrical density approximation when no
        starting volume is available.

    Cylindrical model parameters include:
        * Filament outer diameter.
        * Filament inner diameter.
        * Padding/falloff distance.

    Practical considerations:
        * Cylindrical initialization is useful for unknown structures.
        * Accurate filament diameter estimates improve convergence.

    Non-Uniform Refinement
    ----------------------
    Optional Non-Uniform (NU) refinement can be enabled to improve:

        - Local map quality.
        - High-resolution recovery.
        - Density consistency in flexible regions.

    NU refinement is especially useful for:
        * Structurally heterogeneous helices.
        * Flexible filament assemblies.
        * Maps with uneven local resolution.

    Refinement Controls
    -------------------
    The protocol exposes several refinement tuning parameters:

    - Initial lowpass filtering of the reference volume.
    - Maximum alignment resolution.
    - GS-FSC split resolution.
    - Dynamic or static masking modes.
    - Dynamic mask threshold control.

    Dynamic masking:
        * Adapts mask boundaries during refinement.
        * Helps preserve flexible or variable-density regions.

    Static masking:
        * Uses a fixed mask throughout refinement.
        * Often preferred for rigid, well-defined assemblies.

    Shift Limitation Along Helical Axis
    ----------------------------------
    The protocol optionally restricts translational shifts along the
    filament axis.

    Advantages:
        - Prevents alignment drift.
        - Improves averaging consistency.
        - Enhances resolution for rigid helices.

    Limitations:
        - May negatively affect flexible filaments.
        - Should be disabled when strong conformational variability exists.

    Outputs
    -------
    - Refined helical 3D reconstruction.
    - Updated particle alignment parameters.
    - Symmetry-refined helical geometry.
    - FSC and refinement statistics.
    - Optional Non-Uniform refined map.

    Validation and Safety Checks
    ----------------------------
    The protocol validates:
        - cryoSPARC environment compatibility.
        - GPU availability.
        - Presence of a valid initial model or cylindrical model option.
        - Filament diameter requirements for cylindrical initialization.

    These checks help avoid unstable refinements or incomplete job execution.

    Computational Workflow
    ----------------------
    Internally, the protocol:
        - Builds parameter dictionaries dynamically.
        - Converts symmetry definitions into cryoSPARC-compatible format.
        - Manages GPU execution and queue handling.
        - Launches refinement jobs using cryoSPARC job scheduling.
        - Waits for completion and cleans intermediate files.

    Practical Recommendations
    -------------------------
    - Start with reliable twist and rise estimates whenever possible.
    - Use soft masks to suppress solvent noise.
    - Enable Non-Uniform refinement for flexible or heterogeneous filaments.
    - Use cylindrical initialization only when no prior map exists.
    - Carefully inspect symmetry assumptions before high-resolution refinement.
    - Avoid excessive symmetry enforcement for polymorphic helices.

    Biological Perspective
    ----------------------
    Helical refinement is critical for resolving filamentous biological
    assemblies at near-atomic resolution.

    Reliable reconstructions depend on:
        * Accurate helical symmetry estimation.
        * Appropriate masking strategies.
        * Careful handling of filament flexibility.
        * Stable initialization and refinement settings.

    Proper refinement enables structural interpretation of complex
    biological filaments, including molecular packing, symmetry organization,
    and conformational variability.

    """