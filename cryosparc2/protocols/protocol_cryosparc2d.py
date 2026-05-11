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

import pwem.protocols as pwprot
from pwem import ALIGN_2D
from pyworkflow.object import String
from pyworkflow.protocol.params import (PointerParam, FloatParam, IntParam,
                                        BooleanParam, Positive)
import pyworkflow.utils as pwutils

from .protocol_base import ProtCryosparcBase
from ..convert import (rowToAlignment, convertCs2Star, cryosparcToLocation)
from ..utils import (addComputeSectionParams, cryosparcValidate, gpusValidate,
                     enqueueJob, waitForCryosparc, clearIntermediateResults,
                     copyFiles, getOutputPreffix, isCryosparcStandalone)
from ..constants import *


class ProtCryo2D(ProtCryosparcBase, pwprot.ProtClassify2D):
    """ Wrapper to CryoSparc 2D clustering program.
        Classify particles into multiple 2D classes to facilitate stack cleaning
        and removal of junk particles. Also useful as a sanity check to
        investigate particle quality.
    """

    """
    ProtCryo2D — cryoSPARC 2D Classification Protocol

    Overview
    --------
    Performs unsupervised 2D classification of cryo-EM particle images using
    cryoSPARC clustering algorithms.

    The protocol groups particles into multiple 2D classes in order to:
        - Remove junk or low-quality particles.
        - Assess particle dataset quality.
        - Identify preferred orientations.
        - Detect structural heterogeneity.
        - Generate representative class averages.

    2D classification is commonly used as an early validation and cleaning step
    before high-resolution 3D reconstruction workflows.

    Inputs and Workflow
    -------------------
    - Input Particles:
        Set of aligned or unaligned particle images to classify.

    Workflow stages:
        1. Import particles into cryoSPARC.
        2. Configure classification parameters.
        3. Execute cryoSPARC 2D classification.
        4. Retrieve cryoSPARC outputs.
        5. Convert results into Scipion-compatible metadata.
        6. Generate output classes and representatives.

    The protocol automatically manages:
        - cryoSPARC project creation.
        - GPU allocation.
        - Intermediate file handling.
        - Metadata conversion between cryoSPARC and RELION formats.

    2D Classification Strategy
    --------------------------
    The protocol performs iterative expectation-maximization (EM)
    classification to group particles into structurally similar classes.

    Main objectives:
        - Separate meaningful particles from contaminants.
        - Improve dataset homogeneity.
        - Identify rare particle views.
        - Produce interpretable class averages.

    Runtime approximately scales with:
        - Number of particles.
        - Number of requested classes.
        - Internal resolution settings.

    Classification Parameters
    -------------------------
    The protocol exposes several cryoSPARC classification controls:

    - Number of Classes:
        Defines the number of 2D clusters generated.

    - Maximum Resolution:
        Controls the highest spatial frequency used during alignment
        and reconstruction.

    - Initial Classification Uncertainty:
        Regulates how quickly the algorithm converges toward stable classes.

    Higher uncertainty values:
        * Increase class diversity.
        * Preserve alternative views longer.
        * Reduce premature convergence.

    Lower uncertainty values:
        * Produce faster convergence.
        * Often isolate junk particles more aggressively.

    Circular Masking
    ----------------
    Optional circular masking can be applied during classification.

    Features:
        - Removes corner noise.
        - Focuses alignment on central density.
        - Stabilizes class averaging.

    Mask parameters include:
        - Inner diameter.
        - Outer diameter.
        - Smooth transition region.

    Best practices:
        - Use default masking for globular particles.
        - Adjust mask diameter for elongated or flexible particles.
        - Avoid excessively tight masks that truncate signal.

    Re-centering and Alignment Stability
    ------------------------------------
    The protocol supports iterative class recentering to avoid drift.

    Re-centering options:
        - Density threshold definition.
        - Binary or density-weighted center-of-mass calculation.

    Benefits:
        - Maintains centered class averages.
        - Reduces edge artefacts.
        - Improves alignment consistency.

    Filament and Helical Support
    ----------------------------
    Specialized support is included for filamentous or helical particles.

    Optional features:
        - Vertical alignment of class averages.
        - In-plane rotation estimation.

    This is useful for:
        - Helical assemblies.
        - Filament reconstruction workflows.
        - Directionally constrained particles.

    Noise and Regularization Models
    -------------------------------
    Multiple regularization and noise modeling strategies are available.

    Options include:
        - FRC-based regularization.
        - Full FRC regularization.
        - White noise models.
        - Sigma annealing schedules.

    These mechanisms help:
        - Prevent overfitting.
        - Improve stability for low-SNR datasets.
        - Maintain robust classification behavior.

    Iterative Optimization
    ----------------------
    The classification process combines:
        - Online-EM iterations.
        - Final full-dataset refinement passes.

    Adjustable parameters include:
        - Batch size per class.
        - Number of EM iterations.
        - Final refinement iterations.
        - Zeropadding factor.
        - Initial reference scaling.

    These controls allow optimization for:
        - Large datasets.
        - Small particles.
        - Noisy datasets.
        - High-throughput workflows.

    GPU and Compute Management
    --------------------------
    The protocol automatically manages computational resources.

    Features:
        - GPU allocation.
        - Queue-system compatibility.
        - Standalone and cluster execution modes.
        - cryoSPARC lane assignment.

    The implementation dynamically adapts depending on:
        - Queue usage.
        - Number of available GPUs.
        - cryoSPARC deployment mode.

    Output Generation
    -----------------
    After classification, the protocol:

        - Retrieves cryoSPARC outputs.
        - Converts .cs metadata into STAR files.
        - Generates Scipion SetOfClasses2D objects.
        - Associates particles with class assignments.
        - Creates representative class averages.

    Output objects include:
        - Classified particles.
        - 2D class representatives.
        - Alignment transformations.
        - Updated metadata relationships.

    Metadata and Class Reconstruction
    ---------------------------------
    The protocol reconstructs class information by:
        - Parsing cryoSPARC metadata tables.
        - Recovering particle-to-class assignments.
        - Restoring alignment parameters.
        - Scaling class averages when necessary.

    This ensures interoperability between:
        - cryoSPARC outputs.
        - Scipion workflows.
        - RELION-compatible metadata formats.

    Practical Recommendations
    -------------------------
    - Start with 50–100 classes for heterogeneous datasets.
    - Increase iterations for low-SNR particles.
    - Use circular masking for noisy datasets.
    - Enable FRC regularization to reduce overfitting.
    - Inspect class averages visually before downstream processing.
    - Use larger class counts to identify contaminants and rare views.

    For filament datasets:
        - Enable vertical alignment options.
        - Verify orientation consistency manually.

    Biological Perspective
    ----------------------
    2D classification is one of the most important quality-control steps
    in single-particle cryo-EM processing.

    Well-defined class averages indicate:
        - Good particle alignment.
        - Structural consistency.
        - Proper particle picking.
        - Sufficient signal-to-noise ratio.

    Poor or noisy classes may reveal:
        - Ice contamination.
        - Aggregation.
        - Mis-picked particles.
        - Structural flexibility.

    Effective 2D classification significantly improves the reliability
    of downstream 3D reconstruction and refinement workflows.

    """