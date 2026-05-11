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
from pwem.objects import Volume
from pyworkflow import BETA
from pyworkflow.protocol.params import (PointerParam,  BooleanParam)
from . import ProtCryosparcBase
from ..convert import *
from ..utils import *
from ..constants import *


class ProtCryoSparc3DFlexReconstruction(ProtCryosparcBase):
    """
    Takes in a checkpoint from training as well as prepared high-resolution
    particles and performs high-resolution refinement using L-BFGS under
    the 3DFlex model. This is the stage at which improvements to density
    in high-res regions are computed. Outputs two half-maps that can be used
    for FSC validation, sharpening, and other downstream tasks.
    """


    """
    ProtCryoSparc3DFlexReconstruction — 3D Flex High-Resolution Reconstruction Protocol

    Overview
    --------
    Performs high-resolution cryo-EM reconstruction using a trained cryoSPARC
    3D Flex deformation model and prepared particle datasets.

    The protocol applies flexible refinement under the learned deformation
    framework in order to improve density quality in structurally variable
    regions while preserving high-resolution information.

    Reconstruction is performed using L-BFGS optimization and produces
    independently refined half-maps suitable for gold-standard FSC validation,
    map sharpening, and downstream structural analysis.

    The protocol can also generate a rigid baseline reconstruction for direct
    comparison against the flexible refinement results.

    Main Objectives
    ---------------
    - Perform high-resolution flexible reconstruction.
    - Refine structural variability using trained deformation models.
    - Improve density quality in flexible regions.
    - Generate gold-standard half-maps for validation.
    - Produce rigid and flexible reconstruction comparisons.

    Inputs and Workflow
    -------------------
    - 3D Flex Training Protocol:
        * Provides the trained Flex deformation model.
        * Supplies prepared particles and latent-space information.
        * Defines the learned continuous conformational landscape.

    The workflow follows these stages:
        1. Initialize cryoSPARC project environment.
        2. Connect particles and trained Flex model.
        3. Configure reconstruction parameters.
        4. Launch cryoSPARC high-resolution refinement job.
        5. Wait for reconstruction completion.
        6. Copy generated maps and half-maps.
        7. Build Scipion-compatible output volumes.

    Flexible Reconstruction Strategy
    --------------------------------
    The protocol performs refinement under the 3D Flex deformation model.

    Key features include:
        - Flexible deformation-aware refinement.
        - Continuous heterogeneity modeling.
        - High-resolution density optimization.
        - Latent-space-guided reconstruction.

    This allows structural regions undergoing continuous motion to be
    reconstructed more accurately than with traditional rigid refinement.

    L-BFGS Optimization
    -------------------
    Reconstruction is optimized using the L-BFGS algorithm.

    - Maximum BFGS Iterations:
        * Controls the refinement convergence process.
        * Higher values may improve very high-resolution reconstructions.
        * Larger volumes may require additional iterations.

    The default iteration count is generally sufficient for most datasets,
    balancing reconstruction quality and computational cost.

    Gold-Standard Refinement
    ------------------------
    The protocol supports gold-standard half-set reconstruction workflows.

    - Force Gold-Standard Resplitting:
        * Reassigns particles into new random half-sets.
        * Useful when input alignments lack balanced half-set distributions.
        * Helps maintain robust FSC validation procedures.

    If resplitting is disabled, original half-set assignments are preserved.

    Rigid Reconstruction Baseline
    -----------------------------
    An optional rigid reconstruction can be generated alongside the
    flexible refinement.

    This baseline reconstruction:
        - Uses the same L-BFGS reconstruction framework.
        - Excludes deformation modeling.
        - Enables direct comparison between rigid and flexible refinement.

    Comparing both reconstructions helps evaluate the benefits of
    flexibility-aware modeling for a given dataset.

    cryoSPARC Integration
    ---------------------
    The protocol interfaces directly with cryoSPARC through:
        - Job enqueueing.
        - Input connection mapping.
        - GPU allocation.
        - Execution monitoring.
        - Automatic synchronization with project workflows.

    The protocol supports:
        - Queue-based execution systems.
        - GPU-enabled refinement.
        - Integration with ongoing training workflows.

    Output Generation
    -----------------
    The protocol generates two reconstruction families:

    - Flexible Reconstruction Outputs:
        * Final flexible refined map.
        * Gold-standard half-map A.
        * Gold-standard half-map B.

    - Rigid Reconstruction Outputs:
        * Final rigid refined map.
        * Rigid half-map A.
        * Rigid half-map B.

    Output maps are:
        - Converted into Scipion-compatible Volume objects.
        - Corrected for CCP4 header consistency.
        - Assigned proper sampling rates.
        - Linked with corresponding half-maps.

    Validation and Downstream Analysis
    ----------------------------------
    Generated half-maps can be used for:
        - FSC resolution estimation.
        - Post-processing and sharpening.
        - Local resolution analysis.
        - Model validation.
        - Structural interpretation.

    The flexible and rigid reconstructions can also be compared to assess:
        - Improvement in flexible regions.
        - Recovery of dynamic features.
        - Reduction of conformational blurring.

    Compute Configuration
    ---------------------
    - GPU acceleration is supported and recommended.
    - Compatible with cryoSPARC compute lanes.
    - Supports queue-managed execution environments.
    - Optional SSD caching can improve I/O performance.

    Practical Recommendations
    -------------------------
    - Use well-trained 3D Flex models before reconstruction.
    - Keep default BFGS iterations for initial tests.
    - Increase iterations only for very high-resolution targets.
    - Enable rigid reconstruction for benchmarking purposes.
    - Verify FSC curves using generated half-maps.
    - Inspect flexible regions carefully for biologically meaningful improvements.

    Biological Perspective
    ----------------------
    Flexible reconstruction enables recovery of structural information that
    may be blurred or lost in traditional rigid refinement approaches.

    By incorporating continuous conformational variability directly into the
    refinement process, the protocol improves:
        - Representation of molecular motions.
        - Density quality in flexible domains.
        - Interpretation of dynamic assemblies.
        - Resolution of heterogeneous conformations.

    This approach is particularly valuable for studying molecular machines,
    flexible complexes, and proteins exhibiting continuous structural dynamics.

    """