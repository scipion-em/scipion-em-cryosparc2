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
import numpy

from pwem.objects import CTFModel
import pyworkflow.utils as pwutils
from pyworkflow import NEW
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        BooleanParam, IntParam,
                                        String)

from .protocol_base import ProtCryosparcBase
from .. import RELIONCOLUMNS
from ..convert import convertCs2Star
from ..utils import (addComputeSectionParams, cryosparcValidate,  enqueueJob, waitForCryosparc,
                     copyFiles)


class ProtCryoSparcPatchCTFEstimate(ProtCryosparcBase):
    """
    Patch-based CTF estimation automatically estimates defocus variation for tilted, bent,
    deformed samples and is accurate for all particle sizes and types including flexible and membrane proteins.
    """

    """
    ProtCryoSparcPatchCTFEstimate — Patch-Based CTF Estimation Protocol

    Overview
    --------
    Estimates Contrast Transfer Function (CTF) parameters from cryo-EM
    micrographs using cryoSPARC’s patch-based estimation strategy.
    The protocol is designed to accurately model local defocus variations
    across the micrograph, making it especially effective for:
        - Tilted or bent samples.
        - Deformed ice regions.
        - Flexible macromolecular complexes.
        - Membrane proteins and heterogeneous specimens.

    Compared to traditional global CTF estimation methods, the patch-based
    approach improves robustness by fitting local regions independently,
    leading to more reliable defocus estimation under challenging imaging
    conditions.

    Inputs and Workflow
    -------------------
    - Input Micrographs:
        * Motion-corrected cryo-EM micrographs.
        * Micrographs should contain valid acquisition metadata.
        * Consistent pixel size and microscope parameters are recommended.

    Workflow:
        1. Import micrographs into cryoSPARC.
        2. Perform patch-based CTF estimation.
        3. Search for optimal defocus and phase-shift parameters.
        4. Convert cryoSPARC outputs into Scipion-compatible CTF objects.
        5. Export estimated CTF information for downstream processing.

    CTF Estimation Parameters
    -------------------------
    - Amplitude Contrast:
        * Defines the amplitude contribution in the CTF model.
        * Typical values:
            - 0.07 for standard cryo-EM datasets.
            - 0.10 for membrane proteins or thick ice conditions.

    - Resolution Limits:
        * Minimum Resolution:
            Defines the low-frequency boundary used during fitting.
        * Maximum Resolution:
            Defines the high-frequency limit for CTF estimation.
        * Proper selection improves fitting stability and prevents
          overfitting noise at high frequencies.

    Defocus and Phase Shift Search
    ------------------------------
    - Defocus Search Range:
        * Defines the minimum and maximum defocus values explored
          during grid search.
        * Wide ranges are useful for tilted specimens or variable ice thickness.

    - Phase Shift Search:
        * Supports phase plate datasets.
        * Searches phase-shift values in radians.

    - Phase Shift Refinement Only:
        * Restricts optimization to phase-shift estimation while
          preserving existing defocus values.
        * Useful for datasets with already validated defocus estimates.

    Patch-Based Local Estimation
    ----------------------------
    - The micrograph is divided into smaller patches.
    - Local CTF parameters are estimated independently for each region.
    - Particularly effective for:
        * Spatial defocus gradients.
        * Uneven ice thickness.
        * Beam-induced sample deformation.
        * Large micrographs with local optical variations.

    Outputs
    -------
    - Set of estimated CTF models associated with input micrographs.
    - Defocus U and Defocus V values.
    - Defocus angle estimation.
    - Phase-shift estimation.
    - Estimated maximum CTF resolution.

    Output CTF models are fully compatible with downstream Scipion workflows,
    including:
        - Particle picking.
        - Particle extraction.
        - 2D classification.
        - 3D refinement pipelines.

    Practical Recommendations
    -------------------------
    - Use motion-corrected micrographs before CTF estimation.
    - Start with default resolution limits for most datasets.
    - Increase the defocus search range for tilted or heterogeneous samples.
    - Enable phase-shift estimation for Volta phase plate acquisitions.
    - Visually inspect CTF fits before continuing downstream processing.

    Performance Considerations
    --------------------------
    - GPU acceleration is supported.
    - This implementation restricts execution to a single GPU.
    - Patch-based estimation is computationally more demanding than
      global CTF fitting, but generally provides higher accuracy
      for difficult datasets.

    Biological Perspective
    ----------------------
    Accurate CTF estimation is a critical step in cryo-EM processing because it
    directly affects:
        * Particle alignment precision.
        * High-resolution signal recovery.
        * Final reconstruction quality.
        * Structural interpretability.

    Patch-based CTF estimation is particularly valuable for modern cryo-EM
    datasets containing flexible proteins, membrane complexes, or
    spatially heterogeneous ice conditions.

    """
