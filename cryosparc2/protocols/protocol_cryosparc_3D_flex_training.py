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
from enum import Enum

from pwem import getMatchingFiles
from pwem.objects import ParticleFlex, SetOfParticles, SetOfParticlesFlex
from pwem.protocols import ProtFlexBase
from pyworkflow import BETA
from pyworkflow.protocol import LEVEL_ADVANCED
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        BooleanParam)
from . import ProtCryosparcBase
from ..convert import *
from ..utils import *
from ..constants import *


class outputs(Enum):
    Particles = SetOfParticles


class ProtCryoSparc3DFlexTraining(ProtCryosparcBase, ProtFlexBase):
    """
    Uses a mesh and prepared particles (at a downsampled resolution) to train
    a 3DFlex model. Parameters control the number of latent dimensions,
    size of the model, and training hyperparameters. This job outputs
    checkpoints during training.
    """
    """
    ProtCryoSparc3DFlexTraining — 3D Flex Model Training Protocol

    Overview
    --------
    Trains a cryoSPARC 3D Flex deformation model using prepared particle
    datasets and tetrahedral meshes generated from previous preprocessing
    stages.

    The protocol learns continuous conformational variability directly from
    cryo-EM particle images by combining:
        - A neural deformation model.
        - A canonical density representation.
        - A tetrahedral mesh deformation framework.

    Training produces latent-space representations and checkpoint models
    that can later be used for:
        - Flexible reconstruction.
        - Motion visualization.
        - Conformational trajectory generation.
        - Continuous heterogeneity analysis.

    Main Objectives
    ---------------
    - Learn continuous molecular flexibility.
    - Train latent-space deformation models.
    - Infer particle-specific conformational coordinates.
    - Generate reusable model checkpoints.
    - Capture structural variability beyond discrete classification.

    Inputs and Workflow
    -------------------
    - 3D Flex Data Prepare Protocol:
        * Provides preprocessed particles.
        * Supplies compatible particle metadata.
        * Defines training box sizes and preprocessing conditions.

    - 3D Flex Mesh Prepare Protocol:
        * Provides tetrahedral deformation meshes.
        * Defines structural rigidity constraints.
        * Supplies geometric deformation topology.

    The workflow follows these stages:
        1. Initialize cryoSPARC project environment.
        2. Connect prepared particles and deformation mesh.
        3. Configure latent-space and network parameters.
        4. Launch cryoSPARC 3D Flex training job.
        5. Monitor training execution.
        6. Extract latent coordinates and checkpoints.
        7. Generate particle outputs with associated latent embeddings.

    Latent Space Modeling
    ---------------------
    The protocol learns a continuous latent representation describing
    conformational variability across particles.

    - Number of Latent Dimensions:
        * Controls the complexity of learned motions.
        * Low-dimensional spaces simplify interpretation.
        * Higher-dimensional spaces capture more complex variability.

    Recommended strategy:
        - Start with 2 latent dimensions.
        - Increase dimensionality only if additional motions are observed.

    Learned latent coordinates represent particle-specific conformational
    states within the continuous deformation landscape.

    Neural Network Architecture
    ---------------------------
    The deformation model is implemented as a neural flow generator network.

    Main configurable parameters:
        - Number of Layers:
            Controls model depth and expressiveness.
        - Hidden Units:
            Defines network capacity per layer.

    Larger networks:
        * Capture more complex motions.
        * Increase computational cost.
        * May increase overfitting risk.

    Smaller networks:
        * Train faster.
        * Generalize more robustly.
        * May underrepresent complex flexibility.

    Learning Rate Scheduling
    ------------------------
    Separate learning schedules are used for:
        - Flow deformation parameters.
        - Canonical density map refinement.

    Initial and final learning rates control:
        * Training stability.
        * Convergence speed.
        * Refinement smoothness.

    The protocol progressively anneals:
        - Learning rates.
        - Canonical map resolution.
        - Optimization schedules.

    Rigidity and Deformation Control
    --------------------------------
    Structural smoothness is regulated using rigidity priors.

    - Rigidity Lambda:
        * Controls deformation smoothness.
        * Higher values produce more rigid motions.
        * Lower values allow more flexible local deformations.

    Proper rigidity tuning helps balance:
        - Physical realism.
        - Flexibility representation.
        - Training stability.

    Latent Space Regularization
    ---------------------------
    Several parameters regulate latent-space organization:

    - Noise Injection Standard Deviation:
        * Encourages smooth latent distributions.
        * Prevents unstable latent estimation.
        * Controls latent-space continuity.

    - Latent Centering Strength:
        * Keeps latent coordinates centered near zero.
        * Prevents latent collapse or divergence.

    - Latent Prior Power:
        * Adjusts regularization behavior.
        * Used for advanced latent-space shaping.

    These mechanisms improve:
        - Latent-space interpretability.
        - Smooth conformational transitions.
        - Model robustness.

    Training Schedule
    -----------------
    The protocol performs a default multi-epoch optimization schedule.

    - Extra Epochs:
        * Extend final refinement stages.
        * Useful for difficult or high-resolution datasets.
        * Increase computational cost.

    During training:
        - Resolution gradually increases.
        - Learning rates are annealed.
        - Canonical density maps are progressively refined.

    cryoSPARC Integration
    ---------------------
    The protocol interfaces directly with cryoSPARC through:
        - Job enqueueing.
        - GPU resource allocation.
        - Mesh and particle input connections.
        - Execution monitoring.
        - Automatic synchronization with project workflows.

    Queue-managed execution environments are fully supported.

    Output Generation
    -----------------
    The protocol produces:

    - Trained Model Checkpoints:
        * Saved periodically during training.
        * Reusable for reconstruction and visualization.

    - Latent Coordinate Metadata:
        * Stores learned particle embeddings.
        * Encodes particle-specific conformational states.

    - Flexible Particle Sets:
        * Generated as SetOfParticlesFlex objects.
        * Include associated latent coordinates (Z values).
        * Preserve original particle metadata and CTF information.

    Each output particle receives:
        - Flex metadata.
        - Project and workspace identifiers.
        - Training job references.
        - Learned latent coordinates.

    Downstream Applications
    -----------------------
    Trained models can be used for:
        - Flexible high-resolution refinement.
        - Volume trajectory generation.
        - Continuous motion visualization.
        - Conformational landscape exploration.
        - Structural heterogeneity analysis.

    Latent trajectories can also be generated to visualize continuous
    molecular motions across the learned deformation space.

    Compute Configuration
    ---------------------
    - GPU acceleration is required and strongly recommended.
    - Compatible with cryoSPARC compute lanes.
    - Supports queue-based execution systems.
    - Designed for large-scale neural-network training workflows.

    Practical Recommendations
    -------------------------
    - Start with low latent dimensionality.
    - Use moderate network sizes initially.
    - Monitor latent-space organization during training.
    - Increase rigidity for noisy datasets.
    - Use additional epochs only when necessary.
    - Carefully inspect generated latent trajectories for biological consistency.

    Biological Perspective
    ----------------------
    3D Flex training enables continuous modeling of structural variability
    directly from cryo-EM particle images.

    Unlike discrete classification methods, the protocol captures:
        - Smooth conformational transitions.
        - Coupled domain motions.
        - Continuous molecular rearrangements.
        - Dynamic structural landscapes.

    The learned latent space provides a biologically interpretable framework
    for studying molecular flexibility, functional dynamics, and structural
    heterogeneity in macromolecular complexes.

    """