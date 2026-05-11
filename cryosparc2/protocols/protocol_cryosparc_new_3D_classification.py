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


import pyworkflow.utils as pwutils
from pwem.objects import VolumeMask
from pyworkflow.object import String
from pyworkflow.protocol.params import (FloatParam, LEVEL_ADVANCED,
                                        PointerParam, MultiPointerParam,
                                        CsvList, Positive, IntParam,
                                        BooleanParam, EnumParam)

from .protocol_base import ProtCryosparcBase
from ..convert import (convertBinaryVol, convertCs2Star,
                       rowToAlignment, ALIGN_PROJ, cryosparcToLocation)
from ..utils import (addComputeSectionParams, doImportVolumes,
                     get_job_streamlog, calculateNewSamplingRate,
                     cryosparcValidate, gpusValidate, enqueueJob,
                     waitForCryosparc, clearIntermediateResults, fixVolume,
                     copyFiles, getCryosparcVersion, getOutputPreffix, matchItemRow, parse_version)
from ..constants import *


class ProtCryoSparcNew3DClassification(ProtCryosparcBase):
    """
    3D Classification (BETA) is a new job in cryoSPARC v3.3+ to analyze discrete
    heterogeneity in single particle cryo-EM datasets. This job currently
    implements a version of 3D classification without alignment — a
    classification routine that can complement the existing Heterogeneous
    Refinement job in finding new discrete classes of data.
    """

    """
    ProtCryoSparcNew3DClassification — cryoSPARC 3D Classification Protocol

    Overview
    --------
    Performs 3D classification of single-particle cryo-EM datasets using
    cryoSPARC’s classification framework without particle alignment.
    The protocol is designed to identify discrete structural heterogeneity
    by separating particles into multiple 3D classes based on structural
    variability already present in aligned particles.

    Unlike heterogeneous refinement, this protocol does not refine particle
    orientations during classification, making it especially useful for:
        - Discovering subtle conformational states.
        - Separating compositional heterogeneity.
        - Identifying structurally distinct particle populations.
        - Pre-sorting datasets before refinement.

    Inputs and Workflow
    -------------------
    - Input Particles:
        * Requires particles with CTF information and projection alignment.
        * Multiple particle stacks may be concatenated.
        * A minimum dataset size is enforced for stable classification.

    - Initial Volumes:
        * Optional input reference volumes.
        * Required when using “input” initialization mode.
        * The number of input volumes must match the number of classes.

    - Solvent Mask:
        * Defines the global molecular region used during classification.
        * If absent, masks are automatically generated from reconstructed volumes.

    - Focus Mask:
        * Optional secondary mask for focused classification.
        * Helps isolate flexible or biologically relevant regions.

    Classification Strategy
    -----------------------
    The protocol implements iterative expectation-maximization (EM)
    optimization to classify particles into discrete structural groups.

    Main configurable parameters include:
        - Number of classes.
        - Target reconstruction resolution.
        - Number of online EM epochs.
        - Final full-batch EM iterations.
        - Batch size per class.

    The workflow supports scalable GPU execution and integration with
    cryoSPARC job scheduling.

    Initialization Modes
    --------------------
    Three initialization strategies are available:

    - Simple:
        * Generates initial references from random particle subsets.
        * Recommended for general exploratory classification.

    - PCA:
        * Builds initial references using principal component analysis
          of multiple reconstructions.
        * Useful for detecting subtle variability.

    - Input:
        * Uses externally supplied volumes as initial classes.
        * Best when prior structural knowledge is available.

    Proper initialization is critical for stable convergence and
    biologically meaningful classes.

    PCA-Based Initialization
    ------------------------
    PCA initialization introduces additional controls:
        - Number of reconstructions.
        - Number of PCA components.
        - Particles per reconstruction.

    This strategy helps detect continuous or weak variability patterns
    before discrete clustering.

    Masking and Focused Classification
    ----------------------------------
    The protocol supports automatic and user-defined masking.

    Automatic mask generation includes:
        - Threshold-based mask estimation.
        - Near/far expansion distances.
        - Optional inclusion of negative densities.

    Focused classification can improve separation of:
        - Flexible domains.
        - Ligand occupancy states.
        - Small conformational changes.
        - Partial assemblies.

    Classification Controls
    -----------------------
    Several advanced controls influence class behavior:

    - Hard Classification:
        * Forces particles into a single class per iteration.

    - Class Similarity Annealing:
        * Controls similarity constraints between classes.
        * Gradually relaxes similarity during optimization.

    - Effective Sample Size (ESS) Tuning:
        * Stabilizes early classification stages.

    - FSC Split Control:
        * Preserves or regenerates half-set assignments.

    These parameters allow balancing:
        * Class diversity.
        * Stability.
        * Noise robustness.
        * Convergence speed.

    Image and Reconstruction Processing
    -----------------------------------
    Additional preprocessing and reconstruction options include:
        - High-pass filtering.
        - Anisotropic magnification correction.
        - SSD particle caching.
        - Volume compression into zip archives.
        - Intermediate plotting and diagnostics.

    The protocol dynamically adapts parameters depending on
    the installed cryoSPARC version for compatibility.

    Outputs
    -------
    The protocol generates several output objects:

    - 3D Classes:
        * Particle assignments for each structural class.

    - Output Volumes:
        * Representative reconstructed maps for every class.

    - Solvent Mask:
        * Automatically generated or propagated classification mask.

    - Metadata Files:
        * STAR files converted from cryoSPARC outputs.
        * Volume references and class assignments.

    The protocol also establishes relations between:
        - Input particles.
        - Classified particles.
        - Output volumes.
        - Generated masks.

    Internal Workflow
    -----------------
    The execution pipeline includes:
        1. Input conversion and validation.
        2. cryoSPARC project/job initialization.
        3. Parameter translation to cryoSPARC format.
        4. GPU-aware job submission.
        5. Monitoring of cryoSPARC execution.
        6. Conversion of cryoSPARC outputs into Scipion-compatible objects.
        7. Creation of classes, volumes, and metadata structures.

    Intermediate cryoSPARC outputs are copied, converted,
    and reorganized for downstream compatibility.

    Validation and Consistency Checks
    ---------------------------------
    Before execution, the protocol validates:
        - Presence of CTF information.
        - Projection alignment availability.
        - Minimum particle count.
        - GPU compatibility.
        - Consistency between class count and input volumes.
        - Correct initialization mode selection.

    These checks prevent unstable or invalid classifications.

    Practical Recommendations
    -------------------------
    - Use at least several thousand particles for reliable classification.
    - Start with simple initialization for exploratory analysis.
    - Use PCA initialization for subtle conformational variability.
    - Use focused masks for flexible domains or ligand regions.
    - Avoid identical initial volumes when using input mode.
    - Monitor class distributions and convergence behavior carefully.
    - Enable hard classification only for strongly separated states.

    Biological Perspective
    ----------------------
    3D classification is one of the key tools for studying structural
    heterogeneity in cryo-EM datasets.

    Successful classification depends on:
        - Dataset quality.
        - Proper masking strategy.
        - Realistic class number selection.
        - Appropriate initialization.
        - Careful interpretation of reconstructed classes.

    The protocol enables identification of biologically meaningful
    conformational and compositional states that may otherwise remain hidden
    during consensus refinement workflows.
    """
