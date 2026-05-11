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

from pwem.objects import Coordinate, CTFModel
import pyworkflow.utils as pwutils
from pyworkflow import NEW
from pyworkflow.object import String
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        BooleanParam, IntParam)

from .protocol_base import ProtCryosparcBase
from .. import RELIONCOLUMNS
from ..convert import convertCs2Star
from ..utils import (addComputeSectionParams, cryosparcValidate,  enqueueJob, waitForCryosparc, clearIntermediateResults,
                     copyFiles)


class ProtCryoSparcBlobPicker(ProtCryosparcBase):
    """
    Automatically picks particles by searching for Gaussian signals.
    """

    """
    ProtCryoSparcBlobPicker — Automatic Particle Picking with cryoSPARC

    Overview
    --------
    Automatically detects and picks particles from cryo-EM micrographs using
    cryoSPARC blob-based particle picking algorithms. The protocol identifies
    particle-like Gaussian signals without requiring templates, making it useful
    for early-stage processing or datasets without prior structural information.

    The protocol supports multiple blob geometries and optional CTF estimation,
    providing a flexible workflow for rapid particle detection.

    Inputs and Workflow
    -------------------
    - Input Micrographs: Raw or preprocessed cryo-EM micrographs used for
      particle detection.
    - Particle Diameter Range:
        * Minimum diameter defines the smallest expected particle size.
        * Maximum diameter defines the largest expected particle size.
        * These values guide blob generation and coordinate extraction.
    - Optional CTF estimation using cryoSPARC Patch CTF before particle picking.

    Typical workflow:
        1. Load input micrographs.
        2. Optionally estimate CTF parameters.
        3. Run blob-based particle picking.
        4. Convert cryoSPARC outputs into STAR format.
        5. Generate coordinate and optional CTF output sets.

    Blob Detection Modes
    --------------------
    The protocol supports several blob geometries for particle detection:

    - Circular Blob:
        * Default and recommended mode.
        * Uses circular Gaussian blobs across the specified diameter range.
        * Suitable for approximately spherical or compact particles.

    - Elliptical Blob:
        * Uses elongated blobs defined by minimum and maximum diameters.
        * Useful for anisotropic or elongated particles.

    - Ring Blob:
        * Detects ring-like intensity patterns.
        * Can improve picking for hollow or membrane-associated particles.

    Multiple blob modes may be combined, although enabling too many modes
    simultaneously can increase false positives.

    Particle Separation and Peak Detection
    --------------------------------------
    - Minimum Separation Distance:
        * Prevents overlapping particle picks.
        * Expressed in units of particle diameter.
        * Lower values increase particle density but may introduce duplicates.

    - Maximum Number of Peaks:
        * Controls the number of local maxima evaluated during picking.
        * Higher values improve sensitivity but increase computational cost.

    - Number of Micrographs to Process:
        * Allows quick testing on subsets of the dataset.
        * Useful for parameter optimization before full processing.

    CTF Estimation
    --------------
    Optional Patch CTF estimation can be performed before picking.

    Features:
        * Estimates defocus and phase shift parameters.
        * Produces CTF models associated with each micrograph.
        * Integrates directly into the cryoSPARC workflow.

    This step is particularly useful when downstream processing requires
    validated CTF metadata.

    Output Generation
    -----------------
    The protocol converts cryoSPARC outputs into Scipion-compatible objects:

    - Coordinate Set:
        * Particle coordinates extracted from picked particles.
        * Coordinates are mapped back to their corresponding micrographs.
        * Box size is estimated from the average particle diameter.

    - Optional CTF Set:
        * Contains estimated CTF parameters per micrograph.
        * Includes defocus, phase shift, and resolution estimates.

    Output STAR files are automatically generated from cryoSPARC .cs files.

    GPU and Processing Management
    -----------------------------
    - Supports GPU execution through cryoSPARC job scheduling.
    - Automatically detects GPU usage depending on queue configuration.
    - Intermediate cryoSPARC files are copied and managed internally.
    - Temporary results can be cleaned automatically after execution.

    Practical Recommendations
    -------------------------
    - Start with circular blobs and conservative diameter ranges.
    - Use a subset of micrographs for parameter tuning.
    - Increase minimum separation distance to reduce duplicate picks.
    - Use elliptical blobs for filamentous or elongated particles.
    - Enable CTF estimation when downstream refinement requires accurate optics.
    - Visually inspect coordinates to validate picking quality.

    Biological Perspective
    ----------------------
    Blob-based particle picking is a rapid and unbiased strategy for detecting
    particles in cryo-EM datasets.

    Key considerations for reliable results:
        * Accurate particle diameter estimation.
        * Appropriate blob geometry selection.
        * Careful validation of false positives and contaminants.
        * Consistency between picking parameters and particle morphology.

    This protocol is particularly valuable during initial dataset assessment,
    rapid screening, and early ab-initio reconstruction workflows.
    """