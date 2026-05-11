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

from pwem.objects import Volume
from pwem.protocols import ProtAnalysis3D
import pyworkflow.utils as pwutils
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        LEVEL_ADVANCED, BooleanParam, IntParam,
                                        String)

from .protocol_base import ProtCryosparcBase
from ..utils import (addComputeSectionParams, cryosparcValidate, gpusValidate,
                     enqueueJob, waitForCryosparc, clearIntermediateResults,
                     fixVolume, copyFiles, getOutputPreffix)


class ProtCryoSparcSharppening(ProtCryosparcBase, ProtAnalysis3D):
    """
    Wrapper protocol for the Cryosparc's to calculate the sharpened map.
    """

    """
    ProtCryoSparcSharppening — CryoSPARC Map Sharpening Protocol

    Overview
    --------
    Applies post-processing sharpening to cryo-EM reconstructed volumes
    using cryoSPARC sharpening tools. The protocol enhances high-frequency
    structural details by applying a negative B-factor correction and
    optional masking strategies.

    Sharpening improves the visual interpretability of density maps and can
    facilitate downstream analysis such as:
        - Atomic model building.
        - Secondary structure interpretation.
        - Visualization of side-chain densities.
        - Structural validation and presentation.

    Inputs and Workflow
    -------------------
    - Input Volume:
        * A reconstructed cryo-EM map obtained from refinement.
        * Half-maps can optionally be provided for FSC-based sharpening.
        * Input maps should already be reasonably refined and masked.

    Workflow:
        1. Import reconstructed volume.
        2. Optionally import associated half-maps.
        3. Apply B-factor sharpening.
        4. Optionally generate new masks for FSC estimation.
        5. Produce sharpened output volume.

    Sharpening Parameters
    ---------------------
    - B-Factor:
        * Controls enhancement of high-frequency signal.
        * Negative values increase sharpness.
        * Strong sharpening may amplify noise and artifacts.
        * Typical values depend on map quality and resolution.

    - Full FSC Usage:
        * Uses the full-map FSC instead of half-map FSC.
        * Half-map FSC is generally more conservative and preferred
          for avoiding overestimation of resolution.

    Low-Pass Filtering
    ------------------
    - Falloff Order:
        * Controls the steepness of the low-pass filter transition.
        * Higher values produce sharper cutoff behavior.

    - Falloff Offset:
        * Adjusts the filter cutoff relative to the FSC-derived resolution shell.

    Filtering helps suppress amplified high-frequency noise after sharpening.

    Mask Generation and FSC Control
    -------------------------------
    - Generate New FSC Mask:
        * Creates a new mask specifically optimized for sharpening.
        * Prevents dependence on masks generated during refinement.

    - Mask Threshold:
        * Defines density threshold for mask creation.

    - Near and Far Mask Expansion:
        * Controls dilation of the generated mask.
        * Helps include relevant structural regions while excluding solvent.

    Final Map Masking
    -----------------
    - Spherical Mask:
        * Applies a spherical mask to suppress noisy corners.
        * Reduces file size after compression.
        * Commonly recommended for visualization purposes.

    - Wide Mask:
        * Applies a broader masking strategy to preserve additional density.
        * Useful for large or elongated complexes.

    Outputs
    -------
    - Sharpened cryo-EM volume in MRC format.
    - Optional FSC-optimized masking effects applied to the final map.
    - Output volume compatible with downstream visualization and modeling tools.

    The sharpened map can be used in:
        - Model building software.
        - Validation workflows.
        - Structural interpretation pipelines.
        - Publication-quality visualization.

    Practical Recommendations
    -------------------------
    - Start with moderate negative B-factor values.
    - Avoid excessive sharpening, which may introduce artifacts.
    - Use half-map FSC when possible for more reliable sharpening.
    - Generate a new FSC mask for heterogeneous or flexible structures.
    - Visually inspect the sharpened map before downstream interpretation.

    Performance Considerations
    --------------------------
    - GPU acceleration is supported.
    - Execution is restricted to a single GPU in this implementation.
    - Sharpening is computationally lightweight compared to refinement steps.

    Biological Perspective
    ----------------------
    Map sharpening is an essential post-processing step in cryo-EM because it:
        * Enhances structural interpretability.
        * Improves visibility of fine structural features.
        * Facilitates atomic model fitting.
        * Helps distinguish biologically relevant densities from noise.

    Proper sharpening requires balancing signal enhancement and noise
    amplification to preserve biologically meaningful structural information.

    """