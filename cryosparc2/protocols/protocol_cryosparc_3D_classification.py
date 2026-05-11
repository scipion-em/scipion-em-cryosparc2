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
from pyworkflow.object import String
from pyworkflow.protocol.params import (FloatParam, LEVEL_ADVANCED,
                                        PointerParam, MultiPointerParam,
                                        CsvList, Positive, IntParam,
                                        BooleanParam, StringParam, EnumParam)

from .protocol_base import ProtCryosparcBase
from ..convert import (convertBinaryVol, convertCs2Star,
                       rowToAlignment, ALIGN_PROJ, cryosparcToLocation)
from ..utils import (addSymmetryParam, addComputeSectionParams, doImportVolumes,
                     get_job_streamlog, calculateNewSamplingRate,
                     cryosparcValidate, gpusValidate, getSymmetry, enqueueJob,
                     waitForCryosparc, clearIntermediateResults, fixVolume,
                     copyFiles, getOutputPreffix, matchItemRow)
from ..constants import *


class ProtCryoSparc3DClassification(ProtCryosparcBase):
    """
    Heterogeneous Refinement simultaneously classifies particles and refines
    structures from n initial structures, usually obtained following an
    Ab-Initio Reconstruction. This facilitates the ability to look for small
    differences between structures which may not be obvious at low resolutions,
    and also to re-classify particles to aid in sorting.
    """

    """
    ProtCryoSparc3DClassification — cryoSPARC 3D Heterogeneous Refinement Protocol

    Overview
    --------
    Performs simultaneous 3D classification and refinement of cryo-EM particles
    using cryoSPARC Heterogeneous Refinement algorithms.

    The protocol classifies particles into multiple structural states while
    refining an independent 3D reconstruction for each class.

    Typical applications include:
        - Structural heterogeneity analysis.
        - Separation of conformational states.
        - Removal of damaged or inconsistent particles.
        - Refinement of multiple structural populations.
        - Exploration of compositional variability.

    The protocol is commonly applied after:
        - Ab-initio reconstruction.
        - Initial homogeneous refinement.
        - Preliminary dataset cleaning.

    Inputs and Workflow
    -------------------
    - Input Particles:
        Set of particles with associated CTF information.

    - Initial Volumes:
        Multiple starting 3D references representing potential structural states.

    Workflow stages:
        1. Import particles into cryoSPARC.
        2. Import multiple reference volumes.
        3. Configure heterogeneous refinement parameters.
        4. Execute cryoSPARC multi-reference refinement.
        5. Retrieve classified particles and refined maps.
        6. Convert outputs into Scipion-compatible objects.

    The same reference volume may be reused multiple times to:
        - Break symmetry in classification.
        - Detect hidden heterogeneity.
        - Encourage alternative particle assignments.

    Multi-Reference Classification Strategy
    ---------------------------------------
    The protocol performs iterative particle assignment and refinement
    against multiple competing 3D references.

    During refinement:
        - Particles are probabilistically assigned to classes.
        - Each class generates an independent reconstruction.
        - Assignments evolve iteratively during optimization.

    This approach enables:
        - Separation of subtle conformational changes.
        - Detection of rare particle populations.
        - Improved structural homogeneity within classes.

    Refinement and Downsampling
    ---------------------------
    Refinement is performed using an internal reconstruction box size.

    Features:
        - Automatic Fourier cropping of particles.
        - GPU memory optimization.
        - Resolution-aware refinement scaling.

    Smaller box sizes:
        * Reduce computational cost.
        * Improve performance.
        * Lower GPU memory usage.

    Larger box sizes:
        * Preserve high-resolution information.
        * Increase computational requirements.

    Symmetry Handling
    -----------------
    Global symmetry can be applied to all refinement classes.

    Supported symmetry groups include:
        - Cyclic (Cn)
        - Dihedral (Dn)
        - Icosahedral (I)
        - Octahedral (O)
        - Tetrahedral (T)

    Proper symmetry selection is critical for:
        - Accurate reconstructions.
        - Improved convergence.
        - Enhanced resolution estimation.

    Classification Modes
    --------------------
    The protocol supports both:
        - Soft classification.
        - Hard classification.

    Soft classification:
        * Particles maintain probabilistic membership across classes.
        * Better for continuous heterogeneity.

    Hard classification:
        * Particles are assigned to a single class per iteration.
        * Produces sharper class separation.

    Hard classification may improve:
        - Discrete state separation.
        - Dataset cleaning.
        - Structural partitioning.

    Optimization and Online-EM Parameters
    -------------------------------------
    The refinement process relies on iterative Online-EM optimization.

    Adjustable optimization parameters include:
        - Learning rates.
        - Learning-rate decay schedules.
        - Batch size per class.
        - Assignment convergence thresholds.
        - Resolution convergence criteria.

    Additional controls allow:
        - Initial random assignments.
        - Final full-dataset refinement iterations.
        - Half-map decay regulation.

    These parameters influence:
        - Convergence speed.
        - Stability of classification.
        - Separation quality between classes.

    Noise Modeling and Regularization
    ---------------------------------
    Multiple noise models are available:

        - Symmetric noise
        - White noise
        - Coloured noise

    Noise handling affects:
        - Reconstruction stability.
        - High-resolution behavior.
        - Overfitting resistance.

    Additional regularization controls include:
        - FSC-based filtering strategies.
        - Shared filtering across classes.
        - Initial low-pass filtering of references.

    These mechanisms improve robustness for:
        - Low-SNR datasets.
        - Small particle populations.
        - Heterogeneous reconstructions.

    GPU and Compute Management
    --------------------------
    The protocol automatically manages cryoSPARC computational resources.

    Features include:
        - GPU selection.
        - Queue-system compatibility.
        - Cluster and standalone execution support.
        - Lane assignment management.

    GPU usage dynamically adapts depending on:
        - Queue availability.
        - cryoSPARC configuration.
        - Execution environment.

    Output Generation
    -----------------
    After refinement, the protocol generates:

        - Classified particle sets.
        - Refined 3D volumes for each class.
        - Particle alignment parameters.
        - Scipion-compatible metadata objects.

    Outputs include:
        - SetOfClasses3D
        - SetOfVolumes
        - Updated particle transformations
        - Representative maps per class

    cryoSPARC outputs are automatically:
        - Retrieved from project folders.
        - Converted into STAR metadata.
        - Linked with Scipion objects.
        - Associated with reconstructed volumes.

    Metadata Reconstruction and Volume Handling
    -------------------------------------------
    The protocol reconstructs class information by:
        - Parsing cryoSPARC metadata tables.
        - Recovering class assignments.
        - Restoring projection alignments.
        - Scaling reconstructed volumes appropriately.

    Additional utilities:
        - Volume fixing and normalization.
        - Sampling rate recalculation.
        - Representative volume assignment.

    This guarantees compatibility between:
        - cryoSPARC reconstructions.
        - Scipion workflows.
        - RELION-style metadata structures.

    Iteration Tracking and Monitoring
    ---------------------------------
    The protocol monitors cryoSPARC stream logs to:
        - Detect the latest refinement iteration.
        - Recover output iteration identifiers.
        - Synchronize output reconstruction files.

    This enables robust extraction of:
        - Final particle assignments.
        - Latest refined maps.
        - Iteration-dependent outputs.

    Validation and Quality Control
    ------------------------------
    Validation checks ensure:
        - cryoSPARC installation availability.
        - GPU compatibility.
        - Presence of particle CTF information.
        - At least two initial reference volumes.

    These constraints are critical because heterogeneous refinement
    fundamentally depends on:
        - Multi-reference competition.
        - Accurate CTF correction.
        - Reliable projection matching.

    Practical Recommendations
    -------------------------
    - Use biologically meaningful initial references whenever possible.
    - Start with low-resolution initial maps to avoid model bias.
    - Use at least two volumes for meaningful heterogeneity analysis.
    - Increase batch size for difficult classification problems.
    - Use soft classification for continuous conformational landscapes.
    - Use hard classification for discrete structural states.
    - Carefully inspect particle distributions across classes.
    - Verify reconstructed volumes visually before downstream refinement.

    For highly heterogeneous datasets:
        - Increase the number of initial random assignment iterations.
        - Use lower initial resolution limits.
        - Monitor convergence carefully.

    Biological Perspective
    ----------------------
    Heterogeneous refinement is one of the most powerful tools for studying
    structural variability in cryo-EM datasets.

    The protocol enables identification of:
        - Conformational flexibility.
        - Compositional variability.
        - Ligand binding states.
        - Dynamic molecular rearrangements.

    Successful classification can reveal biologically meaningful states that
    would otherwise be averaged out in homogeneous refinement.

    Reliable heterogeneous refinement depends on:
        - Appropriate initial references.
        - Balanced classification parameters.
        - Sufficient particle counts.
        - Careful interpretation of reconstructed classes.

    """