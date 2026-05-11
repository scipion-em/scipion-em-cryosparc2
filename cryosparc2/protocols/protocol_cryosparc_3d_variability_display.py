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
import os.path

from pwem.objects import Volume, SetOfVolumes
from pyworkflow import BETA
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        LEVEL_ADVANCED, EnumParam,
                                        Positive, BooleanParam, StringParam)
from pwem.protocols import ProtRefine3D

from . import ProtCryosparcBase
from ..convert import *
from ..utils import *
from ..constants import *


class ImportVolumeOutputs(enum.Enum):
    component = SetOfVolumes


class ProtCryoSparc3DVariabilityDisplay(ProtCryosparcBase, ProtRefine3D):
    """
    Protocol to create various versions of a 3D variability result that can be
    used for display
    """

    """
    ProtCryoSparc3DVariabilityDisplay — 3D Variability Visualization Protocol

    Overview
    --------
    Generates interpretable visualizations and particle groupings from a 
    previously computed CryoSPARC 3D Variability Analysis. The protocol 
    reconstructs volumes along variability components and organizes results 
    into trajectories, clusters, or intermediate conformational states.

    This protocol is intended for:
        - Exploring continuous structural heterogeneity.
        - Visualizing conformational transitions.
        - Grouping particles into structural states.
        - Generating movies or intermediate reconstructions for analysis.

    Inputs and Workflow
    -------------------
    - 3D Variability Analysis Protocol:
        Uses the output particles and variability components generated
        from a previous 3D variability analysis job.

    - The workflow automatically:
        * Launches a CryoSPARC variability display job.
        * Copies generated outputs into the Scipion workspace.
        * Converts CryoSPARC `.cs` metadata into `.star` files.
        * Organizes reconstructed volumes and particle subsets.
        * Creates visualization-ready outputs.

    Output Modes
    ------------
    The protocol supports three visualization strategies:

    - Cluster Mode:
        * Groups particles into discrete conformational clusters.
        * Produces one representative volume per cluster.
        * Generates associated particle subsets.
        * Useful for identifying distinct structural states.

    - Simple Mode:
        * Produces a linear trajectory along each variability component.
        * Generates sequential frames representing smooth motion.
        * Ideal for quick exploratory visualization.

    - Intermediates Mode:
        * Reconstructs multiple intermediate conformations.
        * Captures non-linear transitions more accurately.
        * Can optionally export particle subsets contributing to each frame.

    Variability Component Processing
    --------------------------------
    - Users may restrict visualization to selected components.
    - Components can be filtered using comma-separated indices.
    - Intermediate trajectories are generated independently for each mode.
    - Optional rolling-window averaging improves continuity of transitions.

    Volume Processing and Filtering
    -------------------------------
    Several options are available to refine visualization outputs:

    - Downsampling:
        * Reduces box size for lighter computations and visualization.

    - Cropping:
        * Removes unnecessary peripheral regions after downsampling.

    - Resolution Filtering:
        * Applies low-pass filtering to smooth reconstructed volumes.

    - High-pass Filtering:
        * Removes large-scale low-frequency variability.

    - Handedness Flipping:
        * Allows inversion of output volume handedness.

    Particle Classification and Metadata
    ------------------------------------
    - In cluster mode:
        * Particle subsets are converted into Scipion-compatible classes.
        * Representative volumes are assigned automatically.
        * STAR metadata files are generated for each cluster.

    - In intermediates mode:
        * Optional particle subsets can be exported for each frame.
        * Enables downstream classification or focused refinement.

    Outputs
    -------
    Depending on the selected mode, the protocol generates:

    - Set of reconstructed variability volumes.
    - Clustered 3D classes with representative maps.
    - Intermediate trajectory frames.
    - Particle subsets linked to variability states.
    - Metadata files compatible with RELION and Scipion.

    Practical Recommendations
    -------------------------
    - Use Cluster Mode when discrete conformations are expected.
    - Use Simple Mode for rapid exploratory visualization.
    - Use Intermediates Mode for continuous motions or flexible systems.
    - Apply filtering to reduce noise in highly heterogeneous datasets.
    - Downsample outputs for faster visualization of large volumes.
    - Inspect generated trajectories carefully to avoid overinterpretation.

    Workflow Integration
    --------------------
    - Integrates directly with CryoSPARC 3D Variability Analysis.
    - Converts CryoSPARC outputs into Scipion-compatible objects.
    - Supports downstream classification, visualization, and refinement.
    - Automatically manages particle-to-class assignments.

    Biological Perspective
    ----------------------
    - Structural variability often reflects biologically relevant motions.
    - This protocol helps reveal:
        * Domain rearrangements.
        * Flexible regions.
        * Continuous conformational transitions.
        * Multiple functional states within a dataset.

    - Careful interpretation is essential since:
        * Variability may include noise-driven components.
        * Intermediate reconstructions may not always represent
          physically stable conformations.
        * Clustering results depend strongly on dataset quality
          and parameter selection.

    """


