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
    ProtCryoSparcLocalRefine — Local Refinement Protocol

    Overview
    --------
    Performs localized 3D refinement of cryo-EM particles using cryoSPARC’s 
    Local Refinement workflow. This protocol refines particle orientations and 
    shifts around a previously known alignment while focusing on a specific 
    structural region defined by a mask.

    Local refinement is especially useful for:
        - Improving flexible or mobile regions of a complex.
        - Refining local domains without disturbing the global structure.
        - Enhancing high-resolution features in targeted areas.
        - Performing focused refinement after consensus reconstruction.

    Inputs and Workflow
    -------------------
    - Input Particles:
        Particle stack with existing projection alignment information.
        The protocol requires particles with associated CTF models and 3D alignment.

    - Reference Volume:
        Initial 3D map used as refinement reference.
        Best practices:
            * Use a high-quality consensus refinement.
            * Ensure particles and volume originate from the same dataset.
            * Prefer maps already aligned to the input particles.

    - Reference Mask:
        Defines the region to refine locally.
        Important considerations:
            * The mask should include the region of interest.
            * Soft masks are recommended to avoid edge artefacts.
            * Flexible regions benefit most from focused masking.

    Alignment Strategy
    ------------------
    The protocol performs restrained local searches around existing particle poses.

    Main alignment controls:
        - Rotation search extent:
            Defines angular exploration around current orientations.

        - Shift search extent:
            Defines translational search range in Angstroms.

        - Gaussian pose priors:
            Optional soft restraints on rotations and shifts.
            Useful for:
                * Small flexible motions.
                * Preventing unstable refinements.
                * Maintaining consistency with prior alignments.

    Fulcrum and Re-centering
    ------------------------
    The refinement can define the rotational fulcrum using:
        - Mask center of mass.
        - Box center.

    Optional re-centering:
        - Re-center rotations every iteration.
        - Re-center shifts every iteration.

    These options may improve convergence in highly flexible regions, especially
    when pose priors are enabled.

    Refinement and Regularization
    -----------------------------
    The protocol includes cryoSPARC homogeneous refinement features adapted for 
    local refinement.

    Available refinement options:
        - Non-uniform refinement:
            Improves reconstruction quality for heterogeneous regions.

        - Marginalization:
            Efficiently integrates pose uncertainty and can improve results for
            small or noisy particles.

        - Dynamic masking:
            Automatically adapts masking during refinement.

        - Non-negativity enforcement:
            Removes negative density before alignment.

        - Gold-standard splitting:
            Optionally re-generates independent half-sets for FSC validation.

    Dynamic Masking
    ----------------
    Dynamic masking expands the refinement region adaptively.

    Main parameters:
        - Near distance:
            Region fully included in the mask.

        - Far distance:
            Region gradually tapered to zero.

        - Start resolution:
            Resolution at which dynamic masking begins.

        - Absolute value masking:
            Allows inclusion of strong negative density regions.

    Symmetry Handling
    -----------------
    Supports standard cryo-EM symmetry groups:
        - Cyclic (Cn)
        - Dihedral (Dn)
        - Icosahedral (I)
        - Octahedral (O)
        - Tetrahedral (T)

    Proper symmetry selection is critical for:
        - Maximizing signal.
        - Avoiding reconstruction artefacts.
        - Achieving correct local refinement behaviour.

    Outputs
    -------
    The protocol generates:
        - Refined particle set with updated alignments.
        - Refined 3D volume.
        - Gold-standard half maps.
        - FSC information for resolution estimation.

    Output particle metadata includes:
        * Updated projection matrices.
        * cryoSPARC refinement attributes.
        * Random subset assignments.

    Validation and Safety Checks
    ----------------------------
    Before execution, the protocol validates:
        - GPU configuration.
        - Input particle dimensional consistency.
        - Presence of CTF information.
        - Availability of 3D alignment parameters.
        - Compatibility between particles and reference volume.

    Execution Workflow
    ------------------
    The refinement pipeline follows these stages:
        1. Initialize cryoSPARC project.
        2. Convert and prepare input data.
        3. Configure refinement parameters.
        4. Launch cryoSPARC local refinement job.
        5. Monitor execution and wait for completion.
        6. Import refined particles and volumes.
        7. Generate FSC and final outputs.

    GPU and Compute Management
    --------------------------
    The protocol supports:
        - Single GPU execution.
        - Queue-system integration.
        - SSD particle caching for improved performance.

    GPU allocation is automatically adapted depending on whether the protocol
    runs interactively or through a scheduler queue.

    Practical Recommendations
    -------------------------
    - Use consensus refinements as starting references.
    - Apply soft masks tightly around the flexible region.
    - Enable pose priors for highly mobile domains.
    - Use non-uniform refinement for heterogeneous particles.
    - Avoid excessively large angular searches.
    - Inspect FSC curves and refined maps carefully.

    Biological Perspective
    ----------------------
    Local refinement is essential for studying structural heterogeneity in
    macromolecular complexes.

    Typical biological applications include:
        - Flexible domain analysis.
        - Ligand-binding region refinement.
        - Membrane protein conformational variability.
        - Ribosome and spliceosome substructure refinement.
        - Multi-body and continuous flexibility studies.

    Accurate local refinement can significantly improve interpretability of
    biologically relevant regions while preserving global structural consistency.
    """