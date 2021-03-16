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

from pyworkflow import NEW
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
    _label = '3D non-uniform refinement'
    _className = "nonuniform_refine_new"
    _devStatus = NEW

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
                      label="Non-uniform filter order",
                      help='Adaptive Window Factor for cross-validation-optimal '
                           'regularization. Trade off between fast transitions '
                           'between regions (AWF should be lower) and more '
                           'accurate local cross-validation test (AWF should '
                           'be higher). Default of 3 is good, can try as low '
                           'as 1.5 ')

    def _defineParamsName(self):
        """ Define a list with all protocol parameters names"""
        ProtCryoSparc3DHomogeneousRefine._defineParamsName(self)
        self._paramsName += ['refine_do_marg', 'refine_do_marg',
                             'refine_nu_filtertype', 'refine_nu_order',
                             'refine_nu_awf']