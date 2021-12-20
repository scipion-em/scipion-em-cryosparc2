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

from pyworkflow.object import String
from pyworkflow.protocol.params import (FloatParam, LEVEL_ADVANCED, Positive,
                                        IntParam, BooleanParam, EnumParam,
                                        StringParam)

from .protocol_cryorefine import ProtCryoSparcRefine3D
from ..utils import (getSymmetry, enqueueJob, waitForCryosparc,
                     clearIntermediateResults)
from ..constants import *


class ProtCryoSparcNonUniformRefine3D(ProtCryoSparcRefine3D):
    """ Apply non-uniform refinement to achieve higher resolution and map
    quality, especially for membrane proteins. Non-uniform refinement
    iteratively accounts for regions of a structure that have disordered or
    flexible density causing local loss of resolution. Accounting for these
    regions and dynamically estimating their locations can significantly
    improve resolution in other regions as well as overall map quality by
    impacting the alignment of particles and reducing the tendency for
    refinement algorithms to over-fit disordered regions.
    """
    _label = '3D non-uniform refinement(Legacy)'
    _className = "nonuniform_refine"

    def _defineParams(self, form):

        ProtCryoSparcRefine3D._defineParams(self, form)

        # ------------[Non-uniform Refinement 2]-----------------

        form.addSection(label='Advanced Refinement')

        form.addParam('refine_locproc_start_res', FloatParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=8.0,
                      validators=[Positive],
                      label="Local processing start resolution (A)",
                      help='The resolution user sets to start local processing.')

        form.addParam('refine_early_stopping', IntParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=1,
                      validators=[Positive],
                      label="Difference between radwn",
                      help='Difference in radius wave number between refinement '
                           'iterations for early stopping criteria')

        form.addParam('locres_use_mask', BooleanParam,
                      default=True,
                      label="Only compute at voxels within mask",
                      help='Whether or not to compute local resolution only '
                           'inside a mask around the structure. If true, the '
                           'mask is either the user-selected mask given as '
                           'input to this job, or else the refinement mask used '
                           'in the refinement result selected as input for this '
                           'job')

        form.addParam('locres_step_size', IntParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=2,
                      validators=[Positive],
                      label="Step size",
                      help='Sub-sampling step size - the local resolution is '
                           'only computed every this many voxels in each '
                           'dimension. This causes the output of local '
                           'resolution to be a smaller box with larger voxels')

        form.addParam('locres_stop_subsampling', IntParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=5,
                      validators=[Positive],
                      label="Resolution to stop subsampling",
                      help='Stop sub-sampling if the GS-FSC falls below this '
                           'resolution threshold.')

        form.addParam('locres_fsc_thresh', FloatParam,
                      default=0.5,
                      validators=[Positive],
                      label="FSC threshold",
                      help='Threshold for local FSC for measuring local '
                           'resolution value. 0.5 and 0.143 can both be used')

        form.addParam('locres_awf', IntParam,
                      default=6,
                      validators=[Positive],
                      label="Adaptive Window Factor",
                      help='Multiple of the global FSC resolution to use as '
                           'the box size. Dependent on FSC threshold. Default '
                           'of 6 gives reasonable box sizes in most cases for '
                           'FSC threshold of 0.5. >35 is recommended for FSC '
                           'threshold of 0.143.')

        form.addParam('locres_zeropad_factor', FloatParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=2.0,
                      validators=[Positive],
                      label="Zeropadding Factor",
                      help='Zeropadding factor applied to the selected box '
                           'width to improve resolution accuracy')

        form.addParam('locres_cap', BooleanParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=False,
                      label="Cap local resolution estimate",
                      help='Cap local resolution estimate below global '
                           'resolution estimate')

        form.addParam('locres_compute_facility', EnumParam,
                      choices=['GPU', 'CPU'],
                      expertLevel=LEVEL_ADVANCED,
                      default=0,
                      label="Use GPU or CPU for computation",
                      help='he computation facility to use for local processing. '
                           'Options are \'CPU\' and \'GPU\'. Leave as default')

        form.addParam('locres_anneal_factor', FloatParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=0.0,
                      label="Annealing factor",
                      help='Annealing factor to anneal towards global FSC. 0 '
                           'means no annealing and 1 locks to global FSC.')

        form.addParam('locres_fsc_weighting', BooleanParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=False,
                      label="Enable FSC weighting",
                      help='Enable FSC weighting to consider over-fitting '
                           'across mask boundary.')

        form.addParam('localfilter_kernel_type', StringParam,
                      default='lanczos',
                      label="Filter type to apply",
                      help="The filter type to apply on the refined structure. "
                           "Options are 'lanczos' and 'gaussian'.")

        form.addParam('cross_validation', BooleanParam,
                      default=True,
                      label="Enable local cross validation",
                      help='Use local cross validation to estimate local '
                           'resolution during refinement iterations. Leave as '
                           'default')

        form.addParam('temporal_locres', BooleanParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=False,
                      label="Enable temporal locres estimate",
                      help='Take temporal local resolution estimate into '
                           'account. Leave as default')

        form.addParam('use_phenix_sharpen', BooleanParam,
                      default=False,
                      label="Use Phenix's style of sharpening",
                      help="Use Phenix's style of sharpening when outputting "
                           "sharpened maps. Leave as default")

        form.addParam('use_monores', BooleanParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=False,
                      label="Enable MonoRes")

        form.addParam('extra_sharpen', BooleanParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=False,
                      label="Sharpen before local processing")

    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        ProtCryoSparcRefine3D._insertAllSteps(self)

    def _defineParamsName(self):
        """ Define a list with all protocol parameters names"""
        ProtCryoSparcRefine3D._defineParamsName(self)
        self._paramsName += ['refine_locproc_start_res',
                             'refine_early_stopping',
                             'locres_use_mask',
                             'locres_step_size',
                             'locres_stop_subsampling',
                             'locres_fsc_thresh',
                             'locres_awf',
                             'locres_zeropad_factor',
                             'locres_cap',
                             'locres_compute_facility',
                             'locres_anneal_factor',
                             'locres_fsc_weighting',
                             'localfilter_kernel_type',
                             'cross_validation',
                             'temporal_locres',
                             'use_phenix_sharpen',
                             'use_monores',
                             'extra_sharpen',
                             'compute_use_ssd']

    def doRunRefine(self):
        """
        :return:
        """
        if self.mask.get() is not None:
            input_group_connect = {"particles": self.particles.get(),
                                   "volume": self.volume.get(),
                                   "mask": self.mask.get()}
        else:
            input_group_connect = {"particles": self.particles.get(),
                                   "volume": self.volume.get()}

        params = {}

        for paramName in self._paramsName:
            if (paramName != 'refine_symmetry' and
                    paramName != 'refine_noise_model' and
                    paramName != 'locres_compute_facility' and
                    paramName != 'refine_mask' and
                    paramName != 'refine_N'):
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif (paramName == 'refine_N' and
                  int(self.getAttributeValue(paramName)) > 0):
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif paramName == 'refine_symmetry':
                symetryValue = getSymmetry(self.symmetryGroup.get(),
                                           self.symmetryOrder.get())

                params[str(paramName)] = symetryValue
            elif paramName == 'refine_noise_model':
                params[str(paramName)] = str(
                    NOISE_MODEL_CHOICES[self.refine_noise_model.get()])
            elif paramName == 'refine_mask':
                params[str(paramName)] = str(
                    REFINE_MASK_CHOICES[self.refine_mask.get()])
            elif paramName == 'locres_compute_facility':
                params[str(paramName)] = str(COMPUTE_FACILITY_CHOICES[self.locres_compute_facility.get()])

        # Determinate the GPUs to use (in dependence of
        # the cryosparc version)
        try:
            gpusToUse = self.getGpuList()
        except Exception:
            gpusToUse = False

        runRefineJob = enqueueJob(self._className, self.projectName.get(),
                                    self.workSpaceName.get(),
                                    str(params).replace('\'', '"'),
                                    str(input_group_connect).replace('\'', '"'),
                                    self.lane, gpusToUse)

        self.runRefine = String(runRefineJob.get())
        self.currenJob.set(runRefineJob.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(), self.runRefine.get(),
                         "An error occurred in the Refinement process. "
                         "Please, go to cryosPARC software for more "
                         "details.")
        clearIntermediateResults(self.projectName.get(), self.runRefine.get())





