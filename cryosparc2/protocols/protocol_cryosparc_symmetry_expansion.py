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


import pyworkflow.utils as pwutils
from pyworkflow.object import String
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        IntParam)

from .protocol_base import ProtCryosparcBase
from ..convert import (convertCs2Star, readSetOfParticles)
from ..utils import (addComputeSectionParams, cryosparcValidate, gpusValidate,
                     enqueueJob, waitForCryosparc, clearIntermediateResults,
                     addSymmetryParam, getSymmetry, copyFiles)


class ProtCryoSparcSymmetryExpansion(ProtCryosparcBase):
    """ Duplicate particles around a point-group symmetry.
    """


    """
    ProtCryoSparcSymmetryExpansion — Symmetry Expansion Protocol

    Overview
    --------
    Expands particle datasets according to a specified symmetry group using 
    cryoSPARC symmetry operations. The protocol duplicates particles and assigns 
    new orientations corresponding to all symmetry-related views.

    Symmetry expansion is commonly used in cryo-EM workflows to:
        - Analyze asymmetric or flexible regions within symmetric complexes.
        - Perform focused classification or local refinement.
        - Increase angular sampling for downstream analysis.
        - Explore symmetry-related conformational variability.

    Inputs and Workflow
    -------------------
    - Input Particles:
        * Requires particles with projection alignment information.
        * Particle orientations are used to generate symmetry-related copies.

    - Symmetry Definition:
        * Supports standard point-group symmetries:
            - Cyclic (Cn)
            - Dihedral (Dn)
            - Icosahedral (I)
            - Octahedral (O)
            - Tetrahedral (T)
        * Examples:
            - C1 (no symmetry)
            - C4
            - D7

    - Helical Symmetry Support:
        * Optional parameters for helical datasets:
            - Helical twist (degrees)
            - Helical rise (Å)
            - Helical symmetry order
        * Values are typically obtained from previous helical refinement jobs.

    Symmetry Expansion Strategy
    ---------------------------
    - Each particle is duplicated according to the specified symmetry operators.
    - New projection orientations are assigned while preserving particle metadata.
    - Particularly useful for:
        * Localized reconstruction approaches.
        * Masked refinement of flexible domains.
        * Symmetry-relaxed structural analysis.

    Parameter Configuration
    -----------------------
    - Symmetry Group:
        * Defines the symmetry operators applied to particles.
        * Strongly affects the number of expanded particles generated.

    - Helical Parameters:
        * Used only for helical reconstruction workflows.
        * Must be consistent with the refinement symmetry definition.

    Compute and Processing Workflow
    -------------------------------
    - Initializes a cryoSPARC project and processing environment.
    - Converts input particles into cryoSPARC-compatible format.
    - Launches a cryoSPARC symmetry expansion job.
    - Monitors execution and waits for job completion.
    - Converts cryoSPARC output files back into STAR format.
    - Reconstructs the expanded particle set with updated alignments.

    Outputs
    -------
    - Expanded particle set containing:
        * Original particle metadata.
        * Symmetry-related orientations.
        * Updated alignment parameters.

    - Output particles preserve:
        * Sampling rate
        * Dimensions
        * Alignment information

    Practical Recommendations
    -------------------------
    - Use symmetry expansion before focused classification or masked refinement.
    - Ensure the selected symmetry matches the biological assembly.
    - Avoid unnecessary expansion for highly heterogeneous datasets.
    - Helical parameters should only be enabled for filamentous structures.
    - Expanded datasets can become very large; monitor storage requirements.

    Biological Perspective
    ----------------------
    - Symmetry expansion enables analysis of local structural variability 
      hidden by global symmetry averaging.
    - Particularly valuable for:
        * Flexible domains
        * Ligand binding studies
        * Partial occupancy analysis
        * Symmetry-breaking events
    - Helps recover biologically relevant asymmetry within otherwise symmetric particles.

    """