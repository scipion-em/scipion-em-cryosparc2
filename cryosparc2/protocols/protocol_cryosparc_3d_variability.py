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
from pyworkflow import BETA
from pyworkflow.protocol.params import (PointerParam, FloatParam, Positive,
                                        BooleanParam, EnumParam)
from pwem.objects import Volume
from pwem.protocols import ProtRefine3D

from . import ProtCryosparcBase
from ..convert import *
from ..utils import *
from ..constants import *


class ProtCryoSparc3DVariability(ProtCryosparcBase, ProtRefine3D):
    """
    Protocol to compute the principle modes of variability with a dataset
    of aligned particles
    """


    """
    ProtCryoSparc3DVariability — 3D Variability Analysis Protocol

    Overview
    --------
    Computes the principal modes of structural variability from a set of
    aligned cryo-EM particles using CryoSPARC 3D Variability Analysis.
    The protocol estimates orthogonal variability components that describe
    continuous or discrete conformational changes within the dataset.

    This approach is useful for:
        - Exploring molecular flexibility and heterogeneity.
        - Identifying dominant conformational motions.
        - Separating flexible regions from rigid structural cores.
        - Visualizing continuous transitions between structural states.

    Inputs and Workflow
    -------------------
    - Input Particles:
        Aligned particles containing orientation and CTF information.
        The protocol requires particles with associated CTF models.

    - Input Mask:
        Defines the molecular region included in the variability analysis.
        Proper masking is critical to avoid noise-driven variability.

    Typical workflow:
        1. Load aligned particles from a refinement workflow.
        2. Apply a focused mask around the structure of interest.
        3. Define the number of variability modes to compute.
        4. Run iterative optimization to estimate principal components.
        5. Generate variability volumes and updated particle metadata.

    Variability Analysis Parameters
    -------------------------------
    - Number of Modes (var_K):
        Defines how many orthogonal variability components are solved.
        Larger values may capture more subtle motions but increase
        computational complexity.

    - Number of Iterations:
        Controls optimization convergence.
        Higher values may improve stability for heterogeneous datasets.

    - Symmetry:
        Optional symmetry constraints can be applied during analysis.
        Proper symmetry selection improves robustness and reduces noise.

    Frequency Filtering
    -------------------
    Variability estimation can be stabilized using low-pass and high-pass
    filtering strategies.

    - Filter Resolution:
        Defines the spatial frequency limit of variability estimation.

    - Filter Order:
        Controls the sharpness of the applied filter.

    - High-Pass Order:
        Reduces contributions from low-frequency global variations.

    These filters help emphasize biologically relevant conformational
    changes while suppressing noise or reconstruction artifacts.

    Noise and Orthogonalization Models
    ----------------------------------
    - Gram-Schmidt Orthogonalization:
        Ensures computed variability modes remain orthogonal and independent.

    - White Noise Model:
        Uses simplified white noise assumptions instead of colored noise.
        Depending on the dataset, one model may produce cleaner variability
        components than the other.

    Per-Particle Scaling
    --------------------
    The protocol supports different strategies for handling particle scale
    factors during variability estimation:

        * None:
            All particles use identical scaling.

        * Input:
            Uses scale values inherited from previous refinements.

        * Optimal:
            Computes per-particle optimal scaling dynamically during analysis.

    Dynamic scaling often improves robustness for heterogeneous datasets.

    Regularization and Stability
    ----------------------------
    - Lambda Regularization:
        Stabilizes optimization and helps avoid numerical divergence or
        reconstruction artifacts.

    Increasing lambda values may improve convergence stability in noisy
    datasets, although excessive regularization can oversmooth variability.

    Outputs
    -------
    - Variability Volume:
        Reference volume representing the analyzed structure.

    - Output Particles:
        Updated particle set containing transformed metadata and variability
        information.

    - Converted STAR Metadata:
        CryoSPARC outputs are converted into STAR-compatible files for
        downstream Scipion and Relion workflows.

    Practical Recommendations
    -------------------------
    - Use high-quality aligned particles as input.
    - Apply focused masks around biologically relevant regions.
    - Start with 2–3 variability modes before increasing complexity.
    - Use optimal particle scaling for heterogeneous datasets.
    - Increase lambda if instability or artifacts appear.
    - Inspect resulting variability maps carefully to distinguish
      biological motion from noise-driven components.

    Biological Perspective
    ----------------------
    3D Variability Analysis provides insight into molecular dynamics by
    identifying dominant structural motions directly from cryo-EM data.

    Key elements for meaningful interpretation include:
        * Reliable particle alignments.
        * Appropriate masking strategies.
        * Careful selection of variability modes.
        * Validation of biologically plausible motions.

    The protocol is especially valuable for studying flexible complexes,
    domain rearrangements, and continuous conformational landscapes.
    """

