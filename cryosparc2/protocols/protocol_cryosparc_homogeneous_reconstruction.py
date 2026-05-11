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
import ast
import os

import emtable

from pwem import ALIGN_PROJ

import pyworkflow.utils as pwutils
from pyworkflow import NEW
from pyworkflow.object import String
from pyworkflow.protocol.params import (PointerParam, LEVEL_ADVANCED, IntParam,
                                        BooleanParam, EnumParam, FloatParam)
from pwem.objects import Volume

from .protocol_base import ProtCryosparcBase
from ..convert import (convertCs2Star, createItemMatrix,
                       setCryosparcAttributes)
from ..utils import (addComputeSectionParams, calculateNewSamplingRate,
                     cryosparcValidate, gpusValidate, enqueueJob,
                     waitForCryosparc, clearIntermediateResults, fixVolume,
                     copyFiles, addSymmetryParam, getSymmetry,
                     getCryosparcVersion, get_job_streamlog, getOutputPreffix, parse_version)
from ..constants import *


class ProtCryoSparcHomogeneousReconstruct(ProtCryosparcBase):
    """ Create a 3D reconstruction from input particles that already have alignments in 3D.
    """

    """
    ProtCryoSparcHomogeneousReconstruct — Homogeneous 3D Reconstruction Protocol

    Overview
    --------
    Generates a homogeneous 3D reconstruction from particles that already 
    contain valid 3D alignment parameters. The protocol wraps cryoSPARC’s 
    homogeneous reconstruction job and produces a refined consensus map 
    together with updated particle alignment information.

    The protocol is designed for:
        - Generating consensus reconstructions from aligned particles.
        - Producing sharpened half-maps and FSC statistics.
        - Applying symmetry and optional helical constraints.
        - Refining maps using Ewald Sphere correction and advanced filtering.

    Inputs and Workflow
    -------------------
    - Input Particles:
        * Requires particles with associated CTF information.
        * Requires valid 3D alignment parameters.
        * Particle orientations are reused during reconstruction.

    - Optional FSC Mask:
        * Used for FSC computation between half-maps.
        * Helps obtain more realistic resolution estimation.

    Reconstruction Parameters
    -------------------------
    - Reconstruction Box Size:
        * Allows overriding the reconstruction volume dimensions.
        * Default behavior uses the original particle box size.

    - Symmetry Support:
        * Supports cyclic, dihedral, tetrahedral, octahedral, and icosahedral symmetry.
        * Symmetry can significantly improve SNR and map quality.
        * Incorrect symmetry assignment may introduce reconstruction artifacts.

    - Helical Reconstruction:
        * Optional helical twist and rise parameters.
        * Helical symmetry order can be specified.
        * Useful for filamentous or helical assemblies.

    Gold-Standard and Handedness Options
    ------------------------------------
    - Gold-Standard Re-Splitting:
        * Forces generation of new random half-sets.
        * Useful for validation or reconstruction reproducibility.

    - Hand Flipping:
        * Allows inversion of map handedness.
        * Automatically adjusts particle alignments accordingly.
        * Important when correcting incorrect chirality assignments.

    Aberration and Optical Corrections
    ----------------------------------
    - Tilt, Trefoil, Tetrafoil, and Anisotropic Magnification:
        * Individual optical aberrations can be ignored if necessary.
        * Useful for testing reconstruction stability or debugging refinements.

    - Ewald Sphere Correction (EWS):
        * Optional correction for Ewald Sphere curvature effects.
        * Supports:
            - Positive or negative curvature signs.
            - Simple or iterative correction modes.
        * Particularly relevant for high-resolution reconstructions.

    FSC Filtering and Sharpening
    ----------------------------
    - FSC Mask Optimization:
        * Automatically optimizes FSC masking during refinement.

    - Manual Filtering Override:
        * Allows disabling automatic FSC filtering.
        * Users may define:
            - Manual filtering resolution.
            - Butterworth filter order.
            - Sharpening B-factor.

    - Wide Mask Output:
        * Produces compressed final maps with expanded masks.
        * Helps reduce output file size.

    GPU and Performance Settings
    ----------------------------
    - GPU Batch Size:
        * Controls the number of images processed simultaneously.
        * Lower values help reduce GPU memory usage.

    - Intermediate Plot Generation:
        * Optional visualization of intermediate processing steps.
        * Can be disabled for faster execution.

    Outputs
    -------
    - Reconstructed 3D Volume:
        * Final consensus reconstruction map.
        * Includes associated half-maps for validation.

    - Updated Particle Set:
        * Particles with updated alignment metadata.
        * Preserves projection alignment information.

    - FSC Information:
        * Gold-standard FSC curves.
        * Resolution estimation derived from half-map comparison.

    Validation and Safety Checks
    ----------------------------
    - Ensures:
        * Input particles contain valid CTF models.
        * Input particles include 3D alignment information.
        * GPU configuration is compatible with execution requirements.

    Practical Recommendations
    -------------------------
    - Use accurate symmetry whenever biologically justified.
    - Enable Ewald Sphere correction for near-atomic resolution goals.
    - Use FSC mask optimization for more reliable resolution estimation.
    - Adjust GPU batch size if memory limitations occur.
    - Verify handedness carefully before downstream interpretation.

    Biological Perspective
    ----------------------
    - Homogeneous reconstruction aims to recover a single consensus structure
      from aligned particle images.
    - Proper symmetry handling and optical corrections are essential for
      achieving high-resolution structural information.
    - Half-map validation and FSC analysis are critical for assessing map quality
      and avoiding overfitting.

    """