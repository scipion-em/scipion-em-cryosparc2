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
from pyworkflow.protocol.params import PointerParam
from . import ProtCryosparcBase
from ..convert import *
from ..utils import *
from ..constants import *


class ProtCryoSparc3DFlexDataPrepare(ProtCryosparcBase):
    """
    Prepares particles for use in 3DFlex training and reconstruction. At the same
    way, takes in a consensus (rigid) refinement density map, plus optionally
    a segmentation and generates a tetrahedral mesh for 3DFlex.
    """

    """
    ProtCryoSparc3DFlexDataPrepare — 3D Flex Data Preparation Protocol

    Overview
    --------
    Prepares particle datasets and consensus volumes for use in cryoSPARC 3D Flex 
    training and reconstruction workflows. The protocol performs preprocessing of 
    particle stacks, adapts volumes to compatible dimensions, and generates the 
    required metadata and intermediate files needed by the 3D Flex framework.

    The protocol is designed to standardize datasets before flexible reconstruction,
    ensuring that particles and reference maps are correctly formatted and optimized
    for downstream neural-network-based heterogeneity analysis.

    Main Objectives
    ---------------
    - Prepare particles for 3D Flex training.
    - Crop and downsample particles to appropriate box sizes.
    - Filter particles based on scale quality metrics.
    - Ensure compatibility with cryoSPARC 3D Flex requirements.
    - Generate consensus maps and metadata outputs.
    - Convert cryoSPARC outputs into RELION-compatible STAR files.

    Inputs and Workflow
    -------------------
    - Input Particles:
        * A SetOfParticles object containing aligned particle stacks.
        * Particles must include:
            - CTF information.
            - 3D alignment parameters.

    - Reference Volume:
        * Consensus rigid refinement volume used as structural reference.
        * Serves as the initial density map for 3D Flex preparation.

    The workflow follows these stages:
        1. Initialize protocol paths and cryoSPARC project environment.
        2. Convert input particle metadata.
        3. Launch cryoSPARC 3D Flex preparation job.
        4. Wait for cryoSPARC execution to complete.
        5. Convert generated outputs into STAR-compatible metadata.
        6. Build final particle and volume outputs.

    Particle Preprocessing
    ----------------------
    The protocol provides several preprocessing controls:

    - Crop Box Size:
        * Crops particles and volumes to a specified size.
        * Used for high-resolution reconstruction stages.
        * Default behavior preserves original dimensions.

    - Training Box Size:
        * Fourier downsamples cropped particles.
        * Reduces computational cost during neural-network training.
        * Recommended to remain below the consensus FSC resolution limit.

    - Particle Scale Filtering:
        * Removes particles below a minimum scale threshold.
        * Useful for excluding junk or low-quality particles.

    - Particle Count Limitation:
        * Restricts the number of particles used.
        * Ensures compatibility with 3D Flex requirements.
        * Final particle count must be divisible by 1000.

    cryoSPARC Integration
    ---------------------
    The protocol directly interfaces with cryoSPARC jobs through:
        - Job enqueueing.
        - Input connection mapping.
        - Parameter serialization.
        - Execution monitoring.
        - Automatic waiting and synchronization.

    Intermediate cryoSPARC outputs are copied into the protocol workspace
    and converted into compatible metadata formats for Scipion workflows.

    Metadata and Conversion
    -----------------------
    After cryoSPARC execution:
        - cryoSPARC .cs metadata files are converted into RELION STAR files.
        - Particle metadata is reconstructed and linked back to original particles.
        - Alignment and subset information are preserved.
        - Volume headers and sampling rates are corrected automatically.

    Output Generation
    -----------------
    The protocol produces:

    - Output Particles:
        * Filtered and prepared particle set.
        * Preserves original acquisition and alignment information.

    - Output Volume:
        * Consensus map generated during preparation.
        * Sampling rate adjusted to match processed particles.

    Validation Rules
    ----------------
    Before execution, the protocol validates:
        - Presence of CTF information.
        - Availability of 3D alignments.
        - Valid particle counts for 3D Flex training.

    Validation helps prevent incompatible datasets from entering the
    computationally expensive training stages.

    Compute Configuration
    ---------------------
    - Supports cryoSPARC compute lane integration.
    - GPU usage is configurable through compute settings.
    - Optimized for cryoSPARC pipeline execution environments.

    Practical Recommendations
    -------------------------
    - Use high-quality consensus refinements as input volumes.
    - Keep training box sizes moderate to reduce computational cost.
    - Apply scale filtering to remove noisy particles.
    - Ensure particle counts satisfy divisibility constraints.
    - Verify alignment quality before launching 3D Flex preparation.

    Biological Perspective
    ----------------------
    3D Flex preparation is a critical preprocessing stage for studying
    continuous conformational variability in cryo-EM datasets.

    Proper preparation directly impacts:
        - Training stability.
        - Reconstruction quality.
        - Accuracy of flexibility analysis.
        - Interpretation of molecular dynamics and structural heterogeneity.

    Careful particle filtering and preprocessing improve the reliability
    of downstream flexible reconstruction and motion modeling workflows.

    """