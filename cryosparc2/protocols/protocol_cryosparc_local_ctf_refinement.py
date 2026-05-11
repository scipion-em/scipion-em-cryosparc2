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
from pwem.protocols import ProtParticles
import pyworkflow.utils as pwutils
from pyworkflow.object import String
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        LEVEL_ADVANCED, IntParam, Positive)

from .protocol_base import ProtCryosparcBase
from .. import RELIONCOLUMNS
from ..convert import (convertCs2Star, createItemMatrix,
                       setCryosparcAttributes)
from ..utils import (addComputeSectionParams, cryosparcValidate, gpusValidate,
                     enqueueJob, waitForCryosparc, copyFiles)


class ProtCryoSparcLocalCtfRefinement(ProtCryosparcBase, ProtParticles):
    """
    Wrapper protocol for the Cryosparc's per-particle Local CTF refinement.
    Performs per-particle defocus estimation for each particle in a dataset,
    against a given 3D reference structure.
    """
    """
    ProtCryoSparcLocalCtfRefinement — Local CTF Refinement Protocol

    Overview
    --------
    Performs per-particle local CTF refinement using cryoSPARC algorithms.
    The protocol estimates and refines local defocus parameters for each
    particle independently against a reference 3D reconstruction.

    This refinement improves the consistency between experimental particle
    images and the reconstructed volume, helping recover higher-resolution
    structural information.

    Main objectives include:
        - Refining per-particle defocus values.
        - Improving local CTF accuracy.
        - Enhancing high-resolution reconstruction quality.
        - Correcting local optical variations across particles.

    Inputs and Workflow
    -------------------
    - Input Particles:
        Particles with projection alignment information already assigned.
    - Reference Volume:
        3D map used as the structural reference during refinement.
    - Reference Mask:
        Soft mask defining relevant structural regions for refinement.

    Typical workflow:
        1. Import aligned particles.
        2. Provide a refined reference volume.
        3. Define masking and refinement parameters.
        4. Run local CTF optimization in cryoSPARC.
        5. Export updated particles with refined CTF metadata.

    Practical recommendations:
        * Use particles from a stable and well-refined reconstruction.
        * Ensure accurate alignment information before refinement.
        * Use a biologically meaningful soft mask.
        * High-quality half-maps improve FSC-guided refinement.

    Local CTF Refinement Strategy
    -----------------------------
    The protocol estimates local defocus deviations for each particle by
    comparing experimental images against reprojections of the reference map.

    Refinement parameters include:
        * Refinement box size.
        * Minimum fitting resolution.
        * Maximum fitting resolution.
        * Defocus search range.
        * Number of diagnostic plots.

    The search is performed around the current defocus estimate using
    user-defined defocus ranges.

    This process helps compensate for:
        - Ice thickness variations.
        - Local beam-induced effects.
        - Defocus inaccuracies across the dataset.

    Resolution Control
    ------------------
    The refinement can be constrained within specific resolution ranges.

    Parameters:
        * Minimum Fit Resolution:
            Defines the lowest resolution used during fitting.
        * Maximum Fit Resolution:
            Defines the highest resolution used during fitting.
            If omitted, cryoSPARC automatically estimates the limit
            using FSC information from half-maps.

    Best practices:
        * Conservative resolution limits improve stability.
        * Very high-resolution fitting may become unstable on noisy datasets.
        * FSC-driven automatic estimation is recommended for most workflows.

    Masking Strategy
    ----------------
    The protocol requires a soft mask to focus refinement on meaningful
    structural regions.

    Mask priority:
        1. User-provided mask.
        2. Refinement mask from previous processing steps.
        3. Validation failure if no valid mask exists.

    Biological rationale:
        * Excluding solvent regions improves fitting stability.
        * Flexible or noisy regions can be minimized during optimization.
        * Proper masking improves local defocus estimation accuracy.

    Half-Map Support
    ----------------
    If the input reference volume contains half-maps, the protocol
    automatically connects them to cryoSPARC refinement jobs.

    Benefits:
        * FSC-guided resolution estimation.
        * Improved local refinement reliability.
        * Better high-resolution validation.

    Output Generation
    -----------------
    After refinement, the protocol:

        - Retrieves cryoSPARC particle outputs.
        - Converts cryoSPARC metadata into STAR format.
        - Updates particle alignment and CTF information.
        - Generates a refined particle dataset.

    Output particles preserve:
        * Projection alignments.
        * CryoSPARC metadata.
        * Random subset assignments.

    The resulting particle set can be used in:
        - Further refinement.
        - High-resolution reconstruction.
        - Classification workflows.
        - Validation procedures.

    GPU and Execution Management
    ----------------------------
    - Supports GPU-accelerated execution.
    - Compatible with cryoSPARC queue systems.
    - Automatically determines GPU allocation when possible.
    - Supports SSD acceleration options.

    Internal execution workflow:
        1. Convert input metadata.
        2. Prepare cryoSPARC job parameters.
        3. Launch local CTF refinement job.
        4. Monitor execution status.
        5. Convert and store refined outputs.

    Validation and Safety Checks
    ----------------------------
    Before execution, the protocol validates:

        - cryoSPARC installation compatibility.
        - GPU availability.
        - Dimensional consistency between particles and volume.
        - Presence of valid alignment information.
        - Presence of valid masking information.

    Incorrect configurations are rejected before job submission.

    Practical Recommendations
    -------------------------
    - Perform homogeneous refinement before local CTF refinement.
    - Use accurate masks focused on stable protein density.
    - Avoid overly aggressive resolution limits.
    - Use half-maps whenever available.
    - Verify refinement improvements through FSC and map inspection.
    - Re-run refinement if substantial alignment updates occur later.

    Biological Perspective
    ----------------------
    Local CTF refinement is an important step in modern cryo-EM processing
    pipelines, especially for near-atomic resolution studies.

    Accurate local defocus estimation improves:

        * Structural sharpness.
        * High-frequency signal recovery.
        * Atomic model interpretability.
        * Reliability of biological conclusions.

    By refining optical parameters at the particle level, the protocol helps
    compensate for experimental variability and increases the overall quality
    of the final reconstruction.

    """