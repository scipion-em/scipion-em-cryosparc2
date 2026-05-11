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
from pyworkflow.protocol.params import (PointerParam, FloatParam, BooleanParam, FileParam, StringParam)
from . import ProtCryosparcBase
from ..convert import *
from ..utils import *
from ..constants import *


class ProtCryoSparc3DFlexMeshPrepare(ProtCryosparcBase):
    """
    Prepares particles for use in 3DFlex training and reconstruction. At the same
    way, Takes in a consensus (rigid) refinement density map, plus optionally
    a segmentation and generates a tetrahedral mesh for 3DFlex.
    """
    """
    ProtCryoSparc3DFlexMeshPrepare — 3D Flex Mesh Preparation Protocol

    Overview
    --------
    Prepares tetrahedral meshes for cryoSPARC 3D Flex workflows using a
    consensus reconstruction volume and optional segmentation information.

    The protocol generates structural meshes that define how molecular
    deformations are modeled during 3D Flex training and reconstruction.
    These meshes provide the geometric framework used by the deformation
    model to capture continuous conformational variability.

    The protocol can:
        - Generate solvent masks automatically.
        - Use externally provided masks.
        - Create tetrahedral meshes from consensus volumes.
        - Define segmented submeshes for domain-specific flexibility.
        - Configure rigidity weighting across structural regions.

    Inputs and Workflow
    -------------------
    - 3D Flex Data Prepare Protocol:
        * Provides the processed consensus volume.
        * Supplies compatible dimensions and preprocessing metadata.
        * Serves as the structural reference for mesh generation.

    - Optional Input Mask:
        * Binary solvent mask defining molecular regions.
        * Must match the dimensions of the prepared consensus volume.
        * If not provided, the protocol generates a mask automatically.

    The workflow follows these stages:
        1. Initialize cryoSPARC project environment.
        2. Load consensus volume from 3D Flex Data Prepare.
        3. Generate or import solvent mask.
        4. Configure tetrahedral mesh parameters.
        5. Launch cryoSPARC mesh preparation job.
        6. Wait for execution completion.
        7. Export generated mesh visualization files.

    Solvent Mask Preparation
    ------------------------
    When no external mask is provided, the protocol automatically generates
    a solvent mask from the consensus reconstruction.

    The mask generation workflow includes:
        - Low-pass filtering of the input volume.
        - Density thresholding.
        - Morphological dilation.
        - Soft padding expansion.

    Important parameters:
        * Filter Resolution:
            Controls smoothing before thresholding.
        * Threshold Level:
            Defines density cutoff for mask generation.
        * Mask Dilation:
            Expands the mask boundary outward.
        * Soft Padding:
            Smooths mask edges to avoid sharp transitions.

    Best practice:
        - Use conservative threshold values for noisy datasets.
        - Increase padding for flexible peripheral regions.
        - Provide external masks for complex assemblies.

    Tetrahedral Mesh Generation
    ---------------------------
    The protocol creates tetrahedral meshes used by the deformation model.

    - Base Number of Tetrahedral Cells:
        * Defines mesh granularity.
        * Higher values generate finer meshes.
        * Finer meshes improve local flexibility representation but
          increase computational cost.

    Mesh quality directly affects:
        - Deformation smoothness.
        - Motion accuracy.
        - Training stability.
        - Computational efficiency.

    Segmentation and Submesh Definition
    -----------------------------------
    Optional segmentation files can be used to define structural domains.

    Supported segmentation formats:
        - UCSF Chimera Segger (.seg)
        - MRC segmentation maps (.mrc)

    Segment connectivity rules:
        - Connections define how submeshes are fused.
        - Relationships must form a tree structure.
        - Cyclic connections are not allowed.
        - Breadth-first ordering is required.

    This enables:
        - Domain-specific flexibility modeling.
        - Hierarchical motion representation.
        - Controlled deformation propagation.

    Rigid Segment Configuration
    ---------------------------
    Specific regions can be constrained as rigid domains.

    - Rigid Segments:
        * Receive increased rigidity weights.
        * Resist excessive deformation.
        * Help preserve stable structural cores.

    This is particularly useful for:
        - Multi-domain complexes.
        - Rigid-body motions.
        - Preventing unrealistic distortions.

    Rigidity Weighting
    ------------------
    The protocol applies spatially varying rigidity penalties across the mesh.

    - Dense regions:
        * Receive higher rigidity values.
        * Preserve structural consistency.

    - Low-density or empty regions:
        * Receive lower rigidity values.
        * Allow flexible expansion and contraction.

    Optional peripheral stiffening:
        - Stabilizes weak peripheral densities.
        - Helps reduce overfitting in noisy datasets.
        - May oversmooth sharp motion boundaries if overused.

    cryoSPARC Integration
    ---------------------
    The protocol interfaces directly with cryoSPARC through:
        - Job enqueueing.
        - Parameter serialization.
        - Input connection mapping.
        - Execution monitoring.
        - Automatic synchronization with project workflows.

    Intermediate outputs are managed automatically and cleaned after execution.

    Output Generation
    -----------------
    The protocol generates:
        - Tetrahedral mesh representations.
        - PDB mesh visualization files.
        - Mesh geometry required for downstream 3D Flex training.

    Generated PDB mesh files are exported for visualization purposes
    and can be inspected using molecular visualization software.

    Validation Rules
    ----------------
    Before execution, the protocol validates:
        - cryoSPARC environment compatibility.
        - Presence of a valid 3D Flex Data Prepare protocol.
        - Dimensional consistency between masks and consensus volumes.

    Dimension mismatches are prevented to ensure mesh generation stability.

    Compute Configuration
    ---------------------
    - Uses cryoSPARC compute lane integration.
    - GPU acceleration is not required.
    - Designed for distributed cryoSPARC execution environments.

    Practical Recommendations
    -------------------------
    - Use high-quality consensus reconstructions.
    - Start with moderate mesh density values.
    - Use segmentation files for multi-domain systems.
    - Apply rigid constraints to structurally stable regions.
    - Avoid excessively fine meshes for noisy datasets.
    - Carefully tune rigidity weighting to balance flexibility and stability.

    Biological Perspective
    ----------------------
    Mesh preparation is a fundamental step in continuous heterogeneity
    analysis using 3D Flex.

    The generated tetrahedral mesh defines how structural deformations
    propagate across the molecular volume and strongly influences:
        - Motion realism.
        - Domain flexibility interpretation.
        - Structural continuity.
        - Accuracy of learned conformational landscapes.

    Proper mesh design improves the biological interpretability of
    continuous molecular motions reconstructed from cryo-EM datasets.

    """


