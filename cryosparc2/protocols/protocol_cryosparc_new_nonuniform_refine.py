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
    """
    Applies non-uniform refinement to improve cryo-EM map quality and
    resolution by accounting for flexible or heterogeneous regions within
    macromolecular structures. This refinement strategy is especially useful
    for membrane proteins, dynamic assemblies, and complexes containing local
    disorder, where traditional homogeneous refinement approaches may lose
    resolution or overfit flexible density.

    AI Generated:

    Non-Uniform Refinement (ProtCryoSparcNewNonUniformRefine3D) — User Manual
        Overview

        The Non-Uniform Refinement protocol performs advanced single-particle
        cryo-EM refinement designed to improve reconstruction quality in cases
        where structural flexibility, local disorder, or compositional
        variability limit the performance of standard refinement procedures.
        The protocol dynamically models differences in local signal quality
        across the structure, allowing stable regions to refine to higher
        resolution while reducing the negative influence of flexible regions.

        In practical biological workflows, this protocol is commonly used after
        an initial consensus refinement has already produced a reasonable map.
        It is particularly effective for membrane proteins, ribonucleoprotein
        assemblies, multi-domain complexes, and particles containing flexible
        peripheral regions. By adapting regularization locally across the map,
        the refinement process can significantly improve interpretability, map
        continuity, and downstream atomic modeling.

        Biological Motivation

        In many cryo-EM datasets, different regions of the same particle do not
        behave uniformly during imaging and reconstruction. Some domains may be
        highly ordered and rigid, while others remain mobile, partially occupied,
        or conformationally heterogeneous. Conventional refinement approaches
        often treat the entire structure uniformly, which can lead to reduced
        alignment accuracy and excessive smoothing or overfitting.

        Non-uniform refinement addresses this challenge by estimating how local
        regions contribute differently to the reconstruction. Flexible regions
        receive stronger regularization while stable regions retain high-resolution
        information. From a biological perspective, this often results in sharper
        transmembrane helices, clearer secondary structure features, improved side
        chain visibility, and better preservation of structurally stable cores.

        Inputs and Refinement Context

        The protocol typically requires aligned particles together with an initial
        three-dimensional reference map. In most workflows, users begin with a
        homogeneous or consensus refinement and subsequently apply non-uniform
        refinement as a higher-precision optimization step.

        Best results are generally obtained when the starting reconstruction is
        already reasonably well aligned and free of severe heterogeneity. If the
        dataset contains multiple conformational states, substantial compositional
        variability, or major classification uncertainty, additional classification
        or heterogeneity analysis may still be necessary before refinement.

        Adaptive Marginalization

        One important refinement feature is adaptive marginalization, which improves
        the treatment of alignment uncertainty during particle refinement. Instead
        of relying on a single rigid alignment estimate, the protocol evaluates
        multiple possible alignments in an adaptive manner, helping reduce alignment
        bias and improving robustness for difficult datasets.

        This strategy can be particularly beneficial for smaller complexes, weakly
        scattering particles, or datasets with low signal-to-noise ratio. In these
        situations, alignment ambiguity can significantly degrade reconstruction
        quality, and adaptive marginalization often improves convergence and final
        resolution.

        Non-Uniform Regularization

        The central component of the protocol is non-uniform regularization. This
        procedure dynamically estimates which regions of the map are well ordered
        and which are locally disordered or poorly resolved. Stable regions are
        refined more aggressively, whereas flexible regions are constrained to
        avoid over-interpretation of noise.

        Biologically, this approach is highly valuable for systems containing
        flexible loops, mobile domains, detergent belts, lipid-associated regions,
        or partially occupied conformations. Instead of forcing the entire map to
        behave uniformly, the protocol adapts locally to the underlying structural
        variability.

        Filter Models and Spatial Regularization

        The protocol allows different filtering behaviors that determine how local
        regularization transitions between ordered and disordered regions. Smooth
        filtering approaches are often preferred for highly continuous structures,
        while sharper transitions may be useful when boundaries between rigid and
        flexible regions are well defined.

        In most biological applications, the default filtering behavior performs
        well and should only be adjusted by experienced users investigating difficult
        datasets. Excessive modification of regularization behavior may produce
        unstable reconstructions or artificially sharpened features.

        Adaptive Window Factor

        The adaptive window factor controls the balance between local sensitivity
        and refinement stability during non-uniform regularization. Lower values
        increase responsiveness to rapid local transitions, whereas higher values
        provide smoother and more conservative estimation.

        For highly flexible systems, moderate adjustment of this parameter may
        improve map interpretability. However, extreme settings can amplify noise
        or suppress meaningful structural detail. In routine cryo-EM refinement,
        the default configuration is generally appropriate.

        Outputs and Interpretation

        The protocol produces a refined three-dimensional reconstruction together
        with updated particle alignment information and refinement statistics.
        Compared to standard homogeneous refinement, the resulting maps frequently
        display improved local contrast, better-defined secondary structure, and
        enhanced interpretability in rigid regions.

        Users should still interpret flexible or weak-density regions cautiously.
        Although non-uniform refinement improves local regularization, it does not
        eliminate intrinsic conformational variability or compensate for missing
        information within the experimental data.

        Practical Recommendations

        In most cryo-EM workflows, non-uniform refinement is best applied after an
        initial consensus refinement has converged successfully. For membrane
        proteins and structurally heterogeneous complexes, it is often beneficial
        to perform this refinement routinely before downstream model building or
        local refinement procedures.

        If refinement becomes unstable or produces unrealistic sharpening, users
        should verify particle quality, masking strategy, symmetry assignment, and
        upstream classification results. Strong biological heterogeneity is often
        better addressed through focused classification or multibody analysis before
        additional refinement cycles.

        Final Perspective

        Non-uniform refinement has become one of the most important refinement
        strategies in modern cryo-EM because it directly addresses the uneven local
        behavior present in many biological macromolecules. By adapting refinement
        strength according to local structural order, the protocol improves both
        global resolution and biological interpretability while reducing the risk
        of overfitting flexible density regions.
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