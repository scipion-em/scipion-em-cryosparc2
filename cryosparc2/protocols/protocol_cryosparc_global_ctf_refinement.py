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
    Wrapper protocol for the Cryosparc's per-particle Global CTF refinement.
    Performs per-exposure-group CTF parameter refinement of higher-order
    aberrations, against a given 3D reference
    """

    """
    ProtCryoSparcGlobalCtfRefinement — Global CTF Refinement Protocol

    Overview
    --------
    Performs global per-particle CTF refinement using cryoSPARC refinement
    algorithms against a reference 3D volume. The protocol estimates and refines
    higher-order optical aberrations and imaging parameters across exposure groups
    to improve particle alignment accuracy and final reconstruction quality.

    This refinement stage is typically applied after consensus refinement,
    once reliable particle orientations and a high-quality reference map
    are available.

    Inputs and Workflow
    -------------------
    - Input Particles:
        * Requires particles with projection alignment information.
        * Particles should originate from a refined reconstruction workflow.

    - Reference Volume:
        * Provides the structural model used during CTF refinement.
        * High-resolution and well-refined maps improve parameter estimation.

    - Reference Mask:
        * Soft mask defining the region used for refinement.
        * Helps focus refinement on meaningful structural signal while
          excluding solvent noise and flexible regions.

    Typical workflow:
        1. Load aligned particles and reference volume.
        2. Apply optional masking.
        3. Refine optical aberration parameters iteratively.
        4. Update particle metadata and alignment information.
        5. Export refined particles for downstream processing.

    Higher-Order Aberration Refinement
    ----------------------------------
    The protocol supports refinement of several advanced optical parameters:

    - Beam Tilt:
        * Corrects systematic angular deviations in the electron beam.
        * Important for high-resolution reconstructions.

    - Trefoil Aberration:
        * Refines third-order asymmetric optical distortions.

    - Tetrafoil Aberration:
        * Corrects fourth-order beam aberrations.

    - Spherical Aberration:
        * Refines deviations caused by lens imperfections.

    - Anisotropic Magnification:
        * Corrects directional magnification distortions.
        * Particularly useful for datasets collected on imperfectly calibrated systems.

    Each aberration type can be independently enabled or disabled depending
    on dataset quality and refinement goals.

    Ewald Sphere Correction
    -----------------------
    The protocol optionally supports Ewald Sphere curvature correction.

    Features:
        * Accounts for curvature effects in high-resolution datasets.
        * Allows positive or negative curvature sign selection.
        * Particularly relevant for large particles or near-atomic resolution studies.

    This option is considered advanced and is generally recommended only
    for high-resolution refinement workflows.

    Iterative Refinement Strategy
    -----------------------------
    - Multiple refinement iterations can be performed.
    - Iterative refinement allows improvements in one parameter to influence
      estimation of others.
    - Additional plotting controls are available for diagnostic visualization
      of exposure group aberrations.

    Practical considerations:
        * One iteration is usually sufficient for standard refinement.
        * Two or more iterations may improve convergence for challenging datasets.
        * Plot binning can improve visualization of aberration patterns.

    Output Generation
    -----------------
    The protocol converts cryoSPARC outputs into Scipion-compatible particle sets.

    Generated outputs include:
        - Refined particle set with updated alignment information.
        - Updated CTF-related metadata.
        - Refined optical parameters associated with particles.

    cryoSPARC .cs files are automatically converted into STAR-compatible formats
    for downstream interoperability.

    GPU and Job Management
    ----------------------
    - Supports GPU-accelerated execution through cryoSPARC scheduling.
    - Automatically handles queue-based or direct GPU execution modes.
    - Supports compatibility across multiple cryoSPARC versions.
    - Intermediate files and refinement jobs are managed automatically.

    Validation and Compatibility
    ----------------------------
    The protocol validates:
        * cryoSPARC environment configuration.
        * GPU availability and compatibility.
        * Dimensional consistency between particles and reference volume.

    Additional support exists for:
        * Half-map connections.
        * Advanced refinement parameters introduced in newer cryoSPARC versions.

    Practical Recommendations
    -------------------------
    - Use high-quality consensus refinements before running global CTF refinement.
    - Apply soft masks focused on stable structural regions.
    - Enable anisotropic magnification refinement for high-resolution datasets.
    - Use Ewald Sphere correction only when resolution and particle size justify it.
    - Start with a single iteration before attempting multi-iteration refinement.
    - Inspect diagnostic plots to verify aberration convergence.

    Biological Perspective
    ----------------------
    Accurate CTF refinement is essential for extracting high-resolution structural
    information from cryo-EM datasets.

    Key elements for reliable refinement:
        * High-quality reference maps.
        * Accurate particle alignments.
        * Appropriate masking strategies.
        * Careful use of advanced aberration corrections.

    Proper global CTF refinement can significantly improve map interpretability,
    reduce optical distortions, and enhance the reliability of downstream
    structural and biological analyses.
    """