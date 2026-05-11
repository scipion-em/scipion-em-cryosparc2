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
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        LEVEL_ADVANCED, Positive, BooleanParam)

from .protocol_base import ProtCryosparcBase
from ..convert import (convertCs2Star, cryosparcToLocation, rowToCtfModel, rowToAlignment)
from ..utils import (addComputeSectionParams, calculateNewSamplingRate,
                     cryosparcValidate, gpusValidate, enqueueJob,
                     waitForCryosparc, clearIntermediateResults, copyFiles)
from ..constants import *


class ProtCryoSparcSubtract(ProtCryosparcBase, ProtOperateParticles):
    """ Signal subtraction protocol of cryoSPARC.
        Subtract projections of a masked volume from particles.
        """

    """
    ProtCryoSparcSubtract — Particle Signal Subtraction Protocol

    Overview
    --------
    Performs particle signal subtraction using cryoSPARC by projecting a masked
    reference volume and subtracting it from experimental particles.
    This protocol is commonly used to isolate flexible regions, remove dominant
    structural components, or focus downstream refinements on specific domains.

    Signal subtraction is especially useful for:
        - Focused classification of flexible regions.
        - Removing stable scaffold densities.
        - Improving analysis of conformational heterogeneity.
        - Enhancing weak or partially occupied structural features.

    Inputs and Workflow
    -------------------
    - Input Particles:
        Experimental particles with valid projection alignments.
        The protocol assumes particles already contain reliable orientation
        parameters.

    - Reference Volume:
        3D map used to generate projections that will be subtracted from the
        particles.
        Practical recommendations:
            * Use a volume reconstructed from the same particle dataset.
            * Ensure matching pixel size and box dimensions.
            * Maintain consistent grayscale normalization.

    - Reference Mask:
        Defines the region to subtract from particles.
        Important:
            * White regions (1) represent densities to REMOVE.
            * Black regions (0) preserve densities in the particles.
            * Soft masks are recommended to avoid subtraction artifacts.

    Particle Subtraction Strategy
    -----------------------------
    The protocol projects the masked reference map according to each particle
    orientation and subtracts the resulting signal from the experimental images.

    This enables:
        - Isolation of localized conformational changes.
        - Improved focused refinements.
        - Reduction of dominant structural contributions.
        - Better characterization of dynamic assemblies.

    Windowing and Scaling Parameters
    --------------------------------
    - Inner Radius:
        Defines the fully preserved central region during windowing.

    - Outer Radius:
        Controls the transition region for soft masking/windowing.

    - Premultiplier Scaling:
        Improves subtraction stability by adjusting real-space scaling before
        subtraction.

    Practical guidance:
        * Smaller radii may remove useful signal.
        * Larger radii may preserve unwanted densities.
        * Default values are generally appropriate for most datasets.

    Gold-Standard Subtraction
    -------------------------
    The protocol supports half-map aware subtraction for gold-standard workflows.

    - Half-map subtraction:
        Each particle subset is subtracted using its corresponding half-map.

    Benefits:
        * Preserves FSC independence.
        * Reduces overfitting risk.
        * Maintains proper resolution validation.

    Recommendation:
        Keep half-map subtraction enabled whenever half-maps are available.

    Volume and Mask Processing
    --------------------------
    Optional preprocessing steps can improve subtraction quality:

    - Low-pass Filtering:
        Applies resolution filtering to the input structure before subtraction.
        Useful when high-resolution noise affects subtraction stability.

    - Mask Thresholding:
        Converts masks into binarized regions for cleaner subtraction boundaries.

    - Hole Filling:
        Removes discontinuities inside binary masks.

    Best practices:
        * Use smooth masks whenever possible.
        * Avoid aggressive thresholding.
        * Apply low-pass filtering cautiously.

    Outputs
    -------
    - Subtracted Particles:
        New particle stack containing signal-subtracted images.

    - Updated Metadata:
        Includes:
            * Particle locations.
            * CTF information.
            * Projection alignments.
            * Updated sampling rates.

    Output particles remain compatible with downstream cryoSPARC and Scipion
    refinement or classification workflows.

    Validation and Compatibility
    ----------------------------
    The protocol validates:
        - GPU availability.
        - Consistency between particles and reference volume dimensions.
        - Presence of projection alignment metadata.

    Additional support:
        - Compatible with multiple cryoSPARC versions.
        - Supports SSD caching and GPU acceleration.
        - Preserves alignment information during conversion.

    Practical Recommendations
    -------------------------
    - Use accurate masks focused on the region to subtract.
    - Ensure reference maps originate from the same dataset.
    - Preserve half-map workflows for reliable FSC estimation.
    - Apply subtraction before focused classification or local refinement.
    - Visually inspect subtracted particles for residual artifacts.

    Biological Perspective
    ----------------------
    Particle subtraction is a powerful strategy for studying structural
    heterogeneity in cryo-EM datasets.

    Typical biological applications include:
        - Flexible domain analysis.
        - Ligand occupancy studies.
        - Membrane protein conformational variability.
        - Multi-body and focused refinement workflows.

    Reliable subtraction depends on:
        * Accurate alignments.
        * High-quality masks.
        * Properly normalized reference maps.
        * Conservative preprocessing choices.

    """