# **************************************************************************
# *
# *  Authors:     Szu-Chi Chung (phonchi@stat.sinica.edu.tw) 
# *               Yunior C. Fonseca Reyna (cfonseca@cnb.csic.es)
# *
# * SABID Laboratory, Institute of Statistical Science, Academia Sinica
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

import pyworkflow.utils as pwutils
from pyworkflow.object import String
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        LEVEL_ADVANCED, IntParam, Positive,
                                        BooleanParam, EnumParam)

from pwem import ALIGN_PROJ
from pwem.protocols import ProtInitialVolume, ProtClassify3D

from .protocol_base import ProtCryosparcBase
from ..convert import (convertCs2Star, cryosparcToLocation,
                       rowToAlignment)

from ..utils import (addSymmetryParam, addComputeSectionParams,
                     cryosparcValidate, gpusValidate, getSymmetry, enqueueJob,
                     calculateNewSamplingRate, waitForCryosparc,
                     clearIntermediateResults, fixVolume, copyFiles,
                     getOutputPreffix, matchItemRow)
from ..constants import *


class ProtCryoSparcInitialModel(ProtCryosparcBase, ProtInitialVolume,
                                ProtClassify3D):
    """    
    Generate a 3D initial model _de novo_ from 2D particles using
    CryoSparc Stochastic Gradient Descent (SGD) algorithm.
    """

    """
    ProtCryoSparcInitialModel — Ab-Initio 3D Reconstruction Protocol

    Overview
    --------
    Generates one or more initial 3D volumes directly from 2D cryo-EM particles
    using the CryoSPARC stochastic gradient descent (SGD) ab-initio reconstruction
    algorithm. The protocol is designed to estimate low-resolution starting models
    without requiring a pre-existing reference structure.

    This approach is commonly used as the first reconstruction step in a cryo-EM
    workflow before refinement or heterogeneous classification.

    Typical applications include:
        - Generating de novo initial maps from particle datasets.
        - Separating distinct conformations during early reconstruction.
        - Producing multiple candidate starting models for downstream refinement.
        - Initializing heterogeneous reconstruction workflows.

    Inputs and Workflow
    -------------------
    - Input Particles:
        A set of aligned or unaligned cryo-EM particles with associated CTF
        information.

    - The protocol internally:
        * Converts particle metadata into CryoSPARC-compatible format.
        * Launches a CryoSPARC ab-initio reconstruction job.
        * Monitors execution until completion.
        * Converts CryoSPARC outputs into STAR-compatible metadata.
        * Generates Scipion-compatible SetsOfVolumes and 3D classes.

    - Multiple ab-initio classes can be generated simultaneously:
        * Each class is initialized independently.
        * Useful for identifying structural heterogeneity early in processing.

    SGD-Based Ab-Initio Reconstruction
    ----------------------------------
    The reconstruction uses stochastic gradient descent optimization to iteratively
    estimate 3D density maps directly from particle projections.

    Key optimization stages include:
        - Initial low-resolution reconstruction.
        - Progressive resolution annealing.
        - Iterative pose refinement.
        - Batch-based optimization updates.
        - Noise model estimation.

    The workflow gradually increases reconstruction complexity to stabilize
    convergence and reduce reconstruction artefacts.

    Resolution and Frequency Control
    --------------------------------
    Several parameters regulate reconstruction resolution behavior:

    - Initial Resolution:
        Defines the low-frequency starting point for optimization.

    - Maximum Resolution:
        Controls the highest spatial frequency considered during reconstruction.

    - Fourier Radius Step:
        Gradually increases reconstruction frequency coverage during iterations.

    - Initial Fourier Cutoff:
        Applies low-pass filtering to random initial volumes for stability.

    Best practice:
        * Start conservatively with lower resolutions for noisy datasets.
        * Increase maximum resolution only when particle quality supports it.

    Multi-Class Ab-Initio Reconstruction
    ------------------------------------
    The protocol supports simultaneous reconstruction of multiple classes.

    This is useful for:
        - Detecting conformational variability.
        - Separating contaminants or junk particles.
        - Identifying compositional heterogeneity.

    Class similarity annealing parameters help regulate how strongly classes
    influence one another during optimization.

    Recommendations:
        * Use a small number of classes initially.
        * Excessive class numbers may fragment particle populations.

    Noise and Scale Modeling
    ------------------------
    CryoSPARC allows different statistical noise models during reconstruction.

    Available options include:
        - Symmetric noise model (default).
        - White noise model.
        - Coloured noise model.

    Additional experimental options include:
        * Per-micrograph scale correction.
        * Per-particle scale optimization.

    These parameters may improve reconstruction robustness for challenging datasets.

    Minibatch Optimization Strategy
    -------------------------------
    SGD optimization is performed using minibatches of particles.

    Parameters control:
        - Initial minibatch size.
        - Final minibatch size.
        - Automatic minibatch tuning behavior.
        - Transition timing between optimization stages.

    Larger minibatches improve stability but increase computational cost.

    Symmetry and Structural Constraints
    -----------------------------------
    Symmetry can optionally be enforced during reconstruction.

    Supported symmetry groups include:
        - Cyclic (Cn)
        - Dihedral (Dn)
        - Icosahedral (I)
        - Octahedral (O)
        - Tetrahedral (T)

    Practical considerations:
        * C1 symmetry is generally recommended for true ab-initio reconstruction.
        * Incorrect symmetry assignment may introduce severe artefacts.

    Additional structural constraints:
        - Non-negativity enforcement.
        - Real-space centering.
        - Real-space windowing.

    These constraints help stabilize optimization and improve map interpretability.

    Outputs
    -------
    The protocol generates:

    - Output 3D volumes:
        * One reconstructed volume per ab-initio class.
        * Exported as Scipion-compatible volume objects.

    - Output 3D classes:
        * Particle assignments associated with reconstructed volumes.
        * Includes projection alignment information.

    - Converted STAR metadata:
        * CryoSPARC particle metadata transformed into Relion-compatible format.

    - Intermediate reconstruction information:
        * Optional intermediate plots and optimization diagnostics.

    Internal Processing Utilities
    -----------------------------
    The implementation also handles:
        - Automatic file conversion between CryoSPARC and STAR formats.
        - Volume scaling and sampling-rate correction.
        - GPU selection and queue management.
        - Metadata synchronization between particles and classes.
        - Reconstruction job monitoring and error handling.

    Validation and Execution Checks
    -------------------------------
    Before execution, the protocol validates:
        - CryoSPARC environment availability.
        - GPU configuration compatibility.
        - Presence of CTF information in input particles.

    Execution waits until CryoSPARC processing finishes successfully before
    generating outputs.

    Practical Recommendations
    -------------------------
    - Begin with a single class for homogeneous datasets.
    - Use multiple classes when structural variability is suspected.
    - Keep symmetry disabled (C1) unless symmetry is unquestionably known.
    - Use conservative resolution limits for noisy or low-particle datasets.
    - Inspect intermediate reconstructions for instability or collapse.
    - Validate final maps visually before downstream refinement.

    Biological Perspective
    ----------------------
    Ab-initio reconstruction is a critical stage in cryo-EM structure determination
    because it defines the first unbiased estimate of the molecular structure.

    Reliable initial models are essential for:
        * Accurate downstream refinement.
        * Correct particle classification.
        * Detection of conformational heterogeneity.
        * Avoidance of reference bias.

    The quality of the final reconstruction strongly depends on:
        - Particle quality.
        - Dataset homogeneity.
        - Appropriate optimization settings.
        - Correct symmetry usage.

    """