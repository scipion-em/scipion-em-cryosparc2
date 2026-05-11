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

from pyworkflow.protocol.params import (FloatParam, Positive, IntParam,
                                        BooleanParam, EnumParam)
from .protocol_cryosparc_homogeneous_refine import ProtCryoSparc3DHomogeneousRefine


class ProtCryoSparcNewNonUniformRefine3D(ProtCryoSparc3DHomogeneousRefine):
    """ Apply non-uniform refinement to achieve higher resolution and map
    quality, especially for membrane proteins. Non-uniform refinement
    iteratively accounts for regions of a structure that have disordered or
    flexible density causing local loss of resolution. Accounting for these
    regions and dynamically estimating their locations can significantly
    improve resolution in other regions as well as overall map quality by
    impacting the alignment of particles and reducing the tendency for
    refinement algorithms to over-fit disordered regions.
    """

    """
    ProtCryoSparcNewNonUniformRefine3D — Non-Uniform 3D Refinement Protocol

    Overview
    --------
    Performs cryoSPARC Non-Uniform (NU) Refinement to improve 3D reconstruction 
    quality and achievable resolution in cryo-EM datasets. This protocol extends 
    standard homogeneous refinement by dynamically accounting for structural 
    heterogeneity, flexible regions, and locally disordered density.

    Non-uniform refinement is particularly beneficial for:
        - Membrane proteins.
        - Flexible macromolecular complexes.
        - Multi-domain assemblies.
        - Partially disordered structures.
        - Datasets with strong local resolution variation.

    The protocol iteratively adapts regularization across the volume, reducing 
    overfitting in flexible regions while preserving high-resolution information 
    in rigid areas.

    Inputs and Workflow
    -------------------
    This protocol inherits the standard workflow from cryoSPARC Homogeneous 
    Refinement and augments it with advanced non-uniform regularization.

    Main inputs:
        - Input particles with alignment information.
        - Initial 3D reference volume.
        - Optional masks and refinement constraints.
        - Symmetry definitions when applicable.

    Recommended practices:
        * Use a high-quality consensus map as starting reference.
        * Ensure accurate particle alignments before refinement.
        * Apply appropriate masks to isolate meaningful density.
        * Avoid overly aggressive masking of flexible regions.

    Non-Uniform Refinement Strategy
    -------------------------------
    The core goal of non-uniform refinement is to treat different structural 
    regions according to their local signal quality and flexibility.

    Key principles:
        - Rigid regions receive stronger high-resolution refinement.
        - Flexible or poorly resolved regions are regularized adaptively.
        - Alignment stability improves by reducing bias from disordered density.
        - Local overfitting is minimized through cross-validation-driven filtering.

    This adaptive behaviour often improves:
        * Global resolution.
        * Local map interpretability.
        * Particle alignment accuracy.
        * Structural consistency.

    Adaptive Marginalization
    ------------------------
    The protocol supports adaptive marginalization over particle poses and shifts.

    Characteristics:
        - Uses auto-tuned adaptive sampling strategies.
        - Efficiently integrates pose uncertainty.
        - Particularly useful for:
            * Small particles.
            * Noisy datasets.
            * Flexible molecular systems.

    Benefits include:
        - Improved convergence stability.
        - Better refinement robustness.
        - Reduced alignment bias.

    Non-Uniform Regularization
    --------------------------
    Non-uniform regularization dynamically estimates flexible and disordered 
    regions during refinement.

    The refinement process:
        - Detects local variability across the map.
        - Adjusts filtering strength spatially.
        - Preserves signal in stable regions.
        - Suppresses noise amplification in flexible areas.

    This cross-validation-optimal regularization is one of the defining features 
    of cryoSPARC NU refinement.

    Filter Types
    ------------
    The protocol allows different regularization filter models:

        - Butterworth:
            Smooth frequency transitions and stable default behaviour.

        - Rectangular:
            Sharper frequency cutoffs with less smoothing.

        - Gaussian:
            Smooth probabilistic attenuation of frequencies.

    Filter selection influences:
        * Local smoothness.
        * Transition sharpness.
        * Sensitivity to local disorder.

    Advanced Regularization Parameters
    ---------------------------------
    - Filter Order:
        Controls the sharpness of Butterworth filtering.
        Higher values produce steeper transitions.

    - Adaptive Window Factor (AWF):
        Regulates local cross-validation behaviour.

        Lower AWF:
            * Faster transitions between regions.
            * More localized adaptation.

        Higher AWF:
            * More stable local estimation.
            * Smoother regularization behaviour.

    The default settings are generally appropriate for most datasets.

    Outputs
    -------
    The protocol generates:
        - Refined 3D reconstruction.
        - Updated particle alignments.
        - Half maps for FSC validation.
        - Resolution estimation metrics.
        - Improved locally regularized density maps.

    Expected improvements may include:
        * Sharper local features.
        * Reduced overfitting artefacts.
        * Better side-chain visibility.
        * Improved transmembrane density quality.

    Practical Recommendations
    -------------------------
    - Enable non-uniform refinement for heterogeneous datasets.
    - Use adaptive marginalization for small or noisy particles.
    - Start with default filter settings before fine-tuning.
    - Inspect local resolution maps after refinement.
    - Avoid excessive manual regularization unless necessary.
    - Compare FSC curves against homogeneous refinement results.

    Performance Considerations
    --------------------------
    Non-uniform refinement is computationally more demanding than standard 
    homogeneous refinement because of:
        - Adaptive local regularization.
        - Dynamic cross-validation.
        - Spatially varying filtering operations.

    GPU acceleration is strongly recommended for efficient execution.

    Biological Perspective
    ----------------------
    Non-uniform refinement is highly valuable for studying biologically relevant 
    structural variability.

    Typical applications include:
        - Flexible membrane transporters.
        - Dynamic ribonucleoprotein assemblies.
        - Multi-state enzymatic complexes.
        - Conformationally heterogeneous particles.
        - Regions with partial occupancy or mobility.

    By accounting for local disorder during refinement, the protocol improves 
    the interpretability of structurally important regions while maintaining 
    robust global reconstruction quality.
    """