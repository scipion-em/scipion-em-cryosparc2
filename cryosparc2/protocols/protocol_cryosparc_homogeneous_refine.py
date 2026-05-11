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


import pwem.objects as pwobj
import pyworkflow.utils as pwutils
from pwem.convert import getSymmetryMatrices, getUnitCell
from pwem.convert.symmetry import moveParticleInsideUnitCell
from pyworkflow.protocol.params import *

from .protocol_base import ProtCryosparcBase
from ..convert import (convertCs2Star, createItemMatrix,
                       setCryosparcAttributes)
from ..utils import (addSymmetryParam, addComputeSectionParams,
                     calculateNewSamplingRate,
                     cryosparcValidate, gpusValidate, getSymmetry,
                     waitForCryosparc, clearIntermediateResults, enqueueJob,
                     getCryosparcVersion, fixVolume, copyFiles,
                     getOutputPreffix, parse_version)
from ..constants import *


class ProtCryoSparc3DHomogeneousRefine(ProtCryosparcBase):
    """ Protocol to refine a 3D map using cryosparc.
        Rapidly refine a single homogeneous structure to high-resolution and
        validate using the gold-standard FSC. Using new faster GPU code, and
        support for higher-order aberration (beam tilt, spherical aberration,
        trefoil, tetrafoil) correction and per-particle defocus refinement on
        the fly.
    """


    """
    ProtCryoSparc3DHomogeneousRefine — 3D Homogeneous Refinement Protocol

    Overview
    --------
    Performs high-resolution homogeneous 3D refinement of cryo-EM particles 
    using cryoSPARC refinement algorithms. The protocol refines a single 
    consensus structure while optimizing particle alignments, CTF parameters, 
    masking strategies, and reconstruction quality through iterative 
    gold-standard FSC validation.

    The protocol supports:
        - High-resolution homogeneous refinement.
        - Symmetry-aware reconstruction and refinement.
        - Per-particle defocus optimization.
        - Global CTF aberration refinement.
        - Dynamic masking strategies.
        - Ewald Sphere correction.
        - GPU-accelerated processing.

    Inputs and Workflow
    -------------------
    - Input Particles:
        Experimental particles with associated CTF information.
    - Reference Volume:
        Initial 3D reference map used as starting model.
    - Optional Mask:
        Soft mask defining protein and solvent regions.

    Typical workflow:
        1. Import particles and reference volume.
        2. Define symmetry and refinement parameters.
        3. Configure masking and noise models.
        4. Run iterative refinement in cryoSPARC.
        5. Generate refined maps, half-maps, and updated particle alignments.

    Practical recommendations:
        * Use a clean and properly centered initial volume.
        * Ensure particles and reference share identical pixel size and box size.
        * Use masks to focus refinement on stable structural regions.
        * Validate refinement quality with FSC and visual inspection.

    Symmetry Handling
    -----------------
    - Supports standard point-group symmetries:
        * Cyclic (Cn)
        * Dihedral (Dn)
        * Icosahedral (I)
        * Octahedral (O)
        * Tetrahedral (T)

    Additional symmetry features:
        * Automatic symmetry alignment.
        * Symmetry relaxation modes:
            - None
            - Maximization
            - Marginalization

    Notes:
        * Symmetry relaxation cannot be used with C1 symmetry.
        * Correct symmetry selection is critical for accurate refinement.

    Refinement Configuration
    ------------------------
    The protocol exposes advanced refinement controls including:

    - Initial lowpass filtering.
    - GS-FSC split resolution.
    - Highpass filtering.
    - Non-negativity constraints.
    - Batch size autotuning.
    - Per-particle scale minimization.
    - Noise model estimation.

    Supported noise models:
        * Symmetric
        * White
        * Coloured

    These parameters allow balancing:
        - Computational efficiency.
        - Refinement stability.
        - High-resolution recovery.

    Masking Strategies
    ------------------
    Supports multiple masking approaches:

    - Dynamic Mask:
        Automatically adapts during refinement.
    - Static Mask:
        Uses predefined masking regions.
    - Null Mask:
        No masking applied.

    Dynamic masking parameters include:
        * Threshold factor.
        * Near/Far expansion distances.
        * Start resolution.
        * Absolute-value masking.

    Best practices:
        * Dynamic masking is recommended for most datasets.
        * Static masks are useful for highly flexible systems.
        * Overly tight masks may introduce artifacts.

    Defocus Refinement
    ------------------
    Optimizes per-particle defocus values during refinement.

    Features:
        * Iterative defocus estimation.
        * Adjustable search ranges.
        * GPU batch optimization.
        * Diagnostic plotting support.

    Benefits:
        * Improved high-resolution signal recovery.
        * Better particle consistency.

    Caution:
        * Small or noisy datasets may become unstable with aggressive refinement.

    Global CTF Refinement
    ---------------------
    Refines higher-order optical aberrations per exposure group.

    Supported aberration corrections:
        * Beam tilt
        * Trefoil
        * Tetrafoil
        * Spherical aberration

    This refinement improves:
        - Map sharpness.
        - High-frequency accuracy.
        - Optical consistency across micrographs.

    Ewald Sphere Correction
    -----------------------
    Supports correction for Ewald Sphere curvature effects.

    Available options:
        * Simple correction.
        * Iterative correction.
        * Positive or negative curvature modes.
        * Alignment-aware correction.

    Particularly useful for:
        - Large particles.
        - Near-atomic resolution reconstructions.

    Output Generation
    -----------------
    After refinement, the protocol:

    - Retrieves cryoSPARC output maps.
    - Converts particle metadata into STAR format.
    - Generates:
        * Refined volume map.
        * Half-maps.
        * Updated particle alignments.
        * FSC estimation.

    Output particles are updated with:
        - Projection alignment matrices.
        - CryoSPARC refinement metadata.
        - Unit-cell corrected orientations.

    GPU and Execution Management
    ----------------------------
    - Supports multi-GPU execution.
    - Compatible with queue systems.
    - Automatically handles GPU allocation when possible.
    - Supports SSD acceleration.

    Internal workflow:
        1. Convert input data.
        2. Submit cryoSPARC refinement job.
        3. Monitor execution.
        4. Collect and convert outputs.
        5. Clean intermediate files.

    Validation and Safety Checks
    ----------------------------
    Before execution, the protocol validates:

    - cryoSPARC installation compatibility.
    - GPU availability.
    - Presence of CTF information.
    - Symmetry relaxation consistency.

    Invalid configurations are rejected before refinement starts.

    Practical Recommendations
    -------------------------
    - Begin with conservative refinement parameters.
    - Use dynamic masking for most datasets.
    - Enable CTF refinement only when data quality supports it.
    - Use Ewald Sphere correction for high-resolution datasets.
    - Monitor FSC evolution and map interpretability carefully.
    - Validate symmetry choices experimentally when possible.

    Biological Perspective
    ----------------------
    Homogeneous refinement is a critical step in single-particle cryo-EM
    reconstruction pipelines. Accurate refinement enables:

        * High-resolution structural interpretation.
        * Improved atomic model building.
        * Better visualization of functional states.
        * Reliable biological conclusions.

    The combination of alignment optimization, masking strategies,
    aberration correction, and symmetry handling directly impacts the
    interpretability and biological relevance of the final structure.

    """