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
        Performs non-uniform 3D refinement in cryoSPARC to improve map quality
        and local resolution by dynamically accounting for flexible or
        disordered regions within cryo-EM reconstructions.

        AI Generated:

        Non-Uniform Refinement (ProtCryoSparcNewNonUniformRefine3D) — User Manual
            Overview

            The Non-Uniform Refinement protocol applies advanced cryoSPARC
            refinement strategies to improve the quality, interpretability,
            and local resolution of cryo-EM reconstructions. Its primary goal
            is to handle structural heterogeneity and local flexibility more
            effectively than conventional homogeneous refinement approaches.

            In practical cryo-EM workflows, this protocol is especially
            valuable for membrane proteins, flexible macromolecular
            assemblies, multi-domain complexes, and particles containing
            regions with variable conformational stability. Traditional
            refinement methods often overfit highly flexible densities or
            distribute alignment errors unevenly across the map. Non-uniform
            refinement addresses these limitations by dynamically estimating
            local structural variability and adapting the regularization
            strategy during refinement.

            From a biological perspective, this approach frequently produces
            sharper densities, improved secondary-structure definition, and
            enhanced interpretability of flexible regions that would otherwise
            remain poorly resolved in standard refinements.

            Inputs and General Workflow

            The protocol extends the standard homogeneous refinement workflow
            and therefore assumes that particles already possess reliable
            alignment parameters obtained from previous refinement or
            reconstruction stages. The refinement uses an input reference map
            together with particle images to iteratively optimize particle
            orientations, shifts, and reconstructed densities.

            Because this protocol inherits the behavior of homogeneous
            refinement, users typically apply it after consensus refinement,
            ab-initio reconstruction, or heterogeneous refinement steps.
            Non-uniform refinement is often considered a late-stage
            optimization procedure intended to maximize final map quality.

            A good biological practice is to begin with a reasonably stable
            and globally consistent reconstruction before applying
            non-uniform refinement. Extremely heterogeneous or poorly aligned
            datasets may require additional classification or cleaning before
            meaningful local improvements can be achieved.

            Non-Uniform Regularization

            The central feature of this protocol is the non-uniform
            regularization strategy. During refinement, the algorithm
            dynamically estimates which regions of the map are flexible,
            disordered, or weakly supported by the data. These regions are
            treated differently from highly ordered areas to reduce the risk
            of overfitting and improve overall reconstruction stability.

            Biologically, this is particularly important for systems with
            continuous flexibility, peripheral domains, membrane-associated
            regions, or compositional heterogeneity. Flexible regions often
            degrade the alignment quality of otherwise stable structural
            cores. By accounting for these local variations, the protocol can
            improve not only the flexible densities themselves but also the
            quality of globally ordered regions.

            The refinement process continuously adapts its regularization
            behavior during optimization, producing reconstructions that are
            often cleaner and more interpretable than those obtained through
            classical homogeneous refinement alone.

            Adaptive Marginalization

            The protocol includes adaptive marginalization over poses and
            shifts. This strategy efficiently integrates uncertainty in
            particle alignment parameters by using adaptive sampling during
            refinement.

            In biological datasets containing small particles, weak signal,
            or substantial conformational variability, adaptive
            marginalization can improve alignment robustness and reduce the
            tendency to converge toward unstable local minima. This often
            results in improved density continuity and more reliable local
            structural features.

            For many cryo-EM datasets, enabling adaptive marginalization is
            beneficial and generally recommended unless computational
            efficiency becomes a limiting factor.

            Filter Types and Regularization Control

            The protocol provides several filter models used during
            cross-validation-based regularization, including Butterworth,
            rectangular, and Gaussian filters. These filters influence how
            local frequency information is treated during refinement.

            The Butterworth filter is typically the default and most broadly
            applicable option because it provides smooth transitions between
            frequency regions while maintaining stable refinement behavior.
            Gaussian filtering may produce smoother regularization behavior in
            certain datasets, whereas rectangular filtering applies sharper
            cutoffs that may be useful in specialized situations.

            Additional parameters such as filter order and Adaptive Window
            Factor (AWF) allow advanced users to fine-tune the balance
            between local sensitivity and refinement stability. Lower AWF
            values tend to favor faster local transitions between flexible
            and rigid regions, while higher values emphasize more stable
            cross-validation behavior.

            In most biological workflows, the default parameter values are
            appropriate and should only be modified when optimization
            artifacts or unusually heterogeneous datasets require additional
            tuning.

            Outputs and Their Interpretation

            The protocol generates a refined 3D reconstruction with improved
            local regularization and optimized particle alignment parameters.
            Compared with conventional homogeneous refinement, the resulting
            maps frequently display sharper density features, improved local
            contrast, and reduced noise in flexible regions.

            From a biological standpoint, improvements are often most visible
            in peripheral domains, membrane interfaces, flexible loops, or
            regions affected by conformational variability. Enhanced local
            detail may facilitate model building, ligand interpretation, or
            structural analysis of dynamic regions.

            Resolution estimates and FSC-based validation remain important
            indicators of reconstruction quality, but visual inspection is
            equally critical. Users should verify that improvements represent
            biologically meaningful signal rather than artificially sharpened
            noise.

            Practical Recommendations

            In most cryo-EM refinement workflows, non-uniform refinement is
            best applied after achieving a stable consensus reconstruction.
            Using high-quality particles and accurate initial alignments
            significantly improves the effectiveness of the protocol.

            For membrane proteins and flexible complexes, enabling
            non-uniform refinement and adaptive marginalization generally
            produces substantial gains in local map quality. Default filter
            settings are usually sufficient, and unnecessary parameter tuning
            should be avoided unless specific refinement problems are
            observed.

            Users working with highly dynamic assemblies should carefully
            evaluate whether unresolved variability originates from
            continuous flexibility or from the presence of discrete structural
            states. In cases of strong compositional heterogeneity,
            classification approaches may still be required before refinement.

            Final Perspective

            Non-uniform refinement represents one of the most important
            advances in modern cryo-EM image processing because it explicitly
            accounts for local structural variability during optimization.
            Rather than treating the entire reconstruction as uniformly
            ordered, the protocol adapts refinement behavior according to the
            confidence and flexibility of individual regions.

            For biological users, this often translates into improved map
            interpretability, better-resolved structural features, and more
            reliable downstream analysis. Successful application depends on
            starting from a stable consensus reconstruction and understanding
            how local flexibility influences both alignment quality and map
            resolution across the structure.
        """
    _label = '3D non-uniform refinement'
    _className = "nonuniform_refine_new"
    ewsParamsName = []

    def _defineParams(self, form):
        ProtCryoSparc3DHomogeneousRefine._defineParams(self, form)

        # ------------[Non-uniform Refinement]-----------------

        form.addSection(label='Advanced Refinement')

        form.addParam('refine_do_marg', BooleanParam, default=True,
                      label="Adaptive Marginalization",
                      help='Efficiently marginalize over poses and shifts '
                           'using an auto-tuning adaptive sampling strategy. '
                           'Can improve results on small molecules.')

        form.addParam('refine_nu_enable', BooleanParam, default=True,
                      label="Non-uniform refine enable",
                      help='Enable cross-validation-optimal non-uniform '
                           'regularization during refinement.')

        form.addParam('refine_nu_filtertype', EnumParam,
                      choices=['butterworth', 'rect', 'gaussian'],
                      default=0,
                      label="Non-uniform filter type",
                      help='butterworth, rect, or gaussian')

        form.addParam('refine_nu_order', IntParam,
                      default=8,
                      validator=[Positive],
                      label="Non-uniform filter order",
                      help='Order of the butterworth filter used for '
                           'cross-validation-optimal regularization. Default t'
                           'o 8, probably no need to change this.')

        form.addParam('refine_nu_awf', FloatParam,
                      default=3,
                      validator=[Positive],
                      label="Non-uniform AWF",
                      help='Adaptive Window Factor for cross-validation-optimal '
                           'regularization. Trade off between fast transitions '
                           'between regions (AWF should be lower) and more '
                           'accurate local cross-validation test (AWF should '
                           'be higher). Default of 3 is good, can try as low '
                           'as 1.5 ')


    def _defineParamsName(self):
        """ Define a list with all protocol parameters names"""
        ProtCryoSparc3DHomogeneousRefine._defineParamsName(self)
        self._paramsName += ['refine_do_marg', 'refine_nu_enable',
                             'refine_nu_filtertype', 'refine_nu_order',
                             'refine_nu_awf']