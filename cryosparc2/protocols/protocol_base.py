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
import ast
import requests
import logging
logger = logging.getLogger(__name__)

import pwem.protocols as pw
import pyworkflow.object as pwobj
import pyworkflow.utils as pwutils
from pwem.objects import FSC

from ..constants import V3_3_1, excludedFSCValues, fscValues, V4_0_0, V4_1_0, RELIONCOLUMNS
from ..convert import convertBinaryVol, writeSetOfParticles, ImageHandler
from ..utils import (getProjectPath, createEmptyProject,
                     createEmptyWorkSpace, getProjectName,
                     getCryosparcProjectsDir, createProjectContainerDir,
                     doImportParticlesStar, doImportVolumes, killJob, clearJob,
                     get_job_streamlog, getSystemInfo, getJobStatus,
                     STOP_STATUSES, getCryosparcVersion, getProjectInformation,
                     getCryosparcProjectId, _getLicenceFromFile, doImportMicrographs, getCryosparcProjectsList,
                     getCryosparcWorkSpaces, parse_version)


class ProtCryosparcBase(pw.EMProtocol):
    """
    This class contains the common functions for all Cryosparc protocols.
    """

    """
    ProtCryosparcBase — Base cryoSPARC Integration Protocol

    Overview
    --------
    Base class providing the common infrastructure required for integrating
    cryoSPARC workflows within Scipion protocols.

    The class centralizes:
        - cryoSPARC project and workspace initialization.
        - Import/export of cryo-EM data objects.
        - File conversion and scaling utilities.
        - FSC extraction and parsing.
        - cryoSPARC job management and cleanup.
        - Compatibility handling between cryoSPARC versions.

    This class is designed as a reusable foundation for specialized
    cryoSPARC processing protocols.

    Project and Workspace Initialization
    ------------------------------------
    The protocol automatically manages cryoSPARC projects and workspaces.

    Main responsibilities:
        - Detect existing cryoSPARC projects.
        - Create new projects when necessary.
        - Create or reuse cryoSPARC workspaces.
        - Store project identifiers and paths internally.
        - Maintain synchronization between Scipion and cryoSPARC.

    Version-aware behavior is implemented to ensure compatibility across
    different cryoSPARC releases.

    Utility Initialization
    ----------------------
    Internal utility variables are initialized to:
        - Define project directory names.
        - Generate project paths.
        - Create project container directories.
        - Prepare temporary and working paths.

    This guarantees a consistent filesystem structure for all derived protocols.

    Input Conversion and Import Workflow
    ------------------------------------
    The protocol handles automatic import of multiple cryo-EM data types:

        - Particles
        - Volumes
        - Masks
        - Focus masks
        - Micrographs

    Input preparation includes:
        - STAR file generation.
        - Binary volume conversion.
        - cryoSPARC-specific import formatting.
        - File linking and path adaptation.

    The workflow ensures all imported datasets become compatible with
    cryoSPARC internal job execution.

    Particle and File Management
    ----------------------------
    Particle file handling utilities provide:
        - Dynamic filename pattern generation.
        - Path remapping.
        - Input particle retrieval.
        - Pointer access to imported datasets.

    These utilities simplify interoperability between Scipion and cryoSPARC
    filesystem conventions.

    Volume and Mask Import System
    -----------------------------
    The protocol supports importing:
        - Reference volumes.
        - Binary masks.
        - Focus masks.
        - Half maps.

    Features include:
        - Automatic suffix generation depending on cryoSPARC version.
        - Half-map management for gold-standard refinement workflows.
        - Temporary file conversion before import.
        - Internal tracking of generated cryoSPARC jobs.

    The implementation guarantees compatibility with multiple cryoSPARC
    naming conventions introduced across versions.

    Image Scaling and Resolution Handling
    -------------------------------------
    Utilities are provided for scaling cryoSPARC-generated averages
    to match original particle dimensions.

    Scaling workflow:
        - Detect binning differences.
        - Preserve original particle dimensions.
        - Apply spline interpolation or stack scaling.
        - Generate scaled output stacks.

    This is especially useful when cryoSPARC internally downsamples data
    during classification or refinement.

    FSC Extraction and Processing
    -----------------------------
    The protocol provides automated FSC retrieval from cryoSPARC jobs.

    FSC-related features:
        - Download FSC files directly from cryoSPARC services.
        - Parse FSC metadata and resolution curves.
        - Generate Scipion-compatible FSC objects.
        - Support phase-randomized masked map calculations.
        - Compute corrected FSC estimations.

    The implementation supports both:
        - Legacy cryoSPARC web endpoints.
        - Newer API-based communication introduced in later versions.

    FSC Parsing and Data Conversion
    -------------------------------
    FSC files are processed to:
        - Extract frequency and correlation values.
        - Convert raw cryoSPARC outputs into Scipion FSC datasets.
        - Filter unsupported FSC columns.
        - Generate additional corrected FSC curves.

    Parsed FSC datasets are then exposed as protocol outputs.

    Job Monitoring and Iteration Tracking
    ------------------------------------
    Utilities are included to inspect cryoSPARC stream logs and extract:

        - Last refinement iteration.
        - FSC iteration identifiers.
        - Estimated map resolution.
        - Estimated B-factor values.

    This allows downstream protocols to monitor refinement progress
    and retrieve quantitative reconstruction metrics.

    Job Control and Abort Handling
    ------------------------------
    The protocol includes cleanup mechanisms for interrupted executions.

    Abort handling responsibilities:
        - Detect running cryoSPARC jobs.
        - Kill unfinished jobs safely.
        - Clear cryoSPARC job states.
        - Prevent orphan processes.

    This improves workflow robustness and resource management.

    Compatibility and Version Awareness
    -----------------------------------
    The class includes explicit handling for cryoSPARC version differences.

    Examples:
        - Different project identifiers before and after v4.
        - API endpoint changes in FSC retrieval.
        - Version-dependent filename suffixes.
        - Updated workspace and import conventions.

    This ensures stable operation across multiple cryoSPARC environments.

    Extensibility
    -------------
    The class is intended to be subclassed by specialized cryoSPARC protocols.

    Derived protocols can reuse:
        - Import utilities.
        - FSC processing.
        - Job management.
        - File scaling logic.
        - cryoSPARC communication methods.

    Placeholder methods such as:
        - _createModelFile()

    are designed for protocol-specific implementations.

    Biological and Workflow Perspective
    -----------------------------------
    cryoSPARC workflows frequently involve complex interactions between:
        - Particle datasets.
        - Volumetric reconstructions.
        - Masks and focused refinements.
        - FSC-based validation metrics.

    This base protocol standardizes those operations, allowing derived
    protocols to focus on reconstruction algorithms rather than
    infrastructure management.

    By centralizing cryoSPARC interoperability, the class improves:
        - Reproducibility.
        - Workflow consistency.
        - Data traceability.
        - Multi-version compatibility.

    """




