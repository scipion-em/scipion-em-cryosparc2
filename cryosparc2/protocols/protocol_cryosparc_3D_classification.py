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

from pyworkflow.protocol.params import (FloatParam, LEVEL_ADVANCED,
                                        PointerParam)

from . import ProtCryosparcBase
from ..convert import writeSetOfParticles
from ..utils import *
from ..constants import *


class ProtCryoSparc3DClassification(ProtCryosparcBase):
    """
    Heterogeneous Refinement simultaneously classifies particles and refines
    structures from n initial structures, usually obtained following an
    Ab-Initio Reconstruction. This facilitates the ability to look for small
    differences between structures which may not be obvious at low resolutions,
    and also to re-classify particles to aid in sorting.
    """
    _label = '3D Classification'

    def _initialize(self):
        self._defineFileNames()

    def _defineFileNames(self):
        """ Centralize how files are called. """
        myDict = {
            'input_particles': self._getTmpPath('input_particles.star'),
            'out_particles': self._getExtraPath('output_particle.star'),
            'stream_log': self._getPath() + '/stream.log'
        }
        self._updateFilenamesDict(myDict)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      label="Input particles", important=True,
                      validators=[Positive],
                      help='Select the input images from the project.')
        form.addParam('refVolume', PointerParam,
                      pointerClass='Volume',
                      important=True,
                      label="Initial volumes",
                      help='Multiple initial volumes to refine and classify. '
                           'The same input volume can be connected multiple '
                           'times.')
        form.addParallelSection(threads=1, mpi=1)

        # --------------[Heterogeneous Refinement]---------------------------
        form.addSection(label='Heterogeneous Refinement')
        form.addParam('multirefine_N', IntParam, default=128,
                      label="Refinement box size (Voxels)",
                      help='Box size of each volume during refinement. '
                           'Particles will automatically be downsampled '
                           '(Fourier cropped) to this box size on the fly. '
                           'Keep this as small as possible to limit GPU '
                           'memory usage.')

        addSymmetryParam(form)

        form.addParam('multirefine_sharp_bfactor', IntParam, default=-100,
                      label="Plotting bfactor",
                      help='B-Factor to apply to the structures before '
                           'plotting, to enhance medium/high resolution '
                           'detail in plots. Outputs are not affected by this '
                           'parameter.')

        form.addParam('multirefine_force_hard_class', BooleanParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=False,
                      label="Force hard classification",
                      help='Force hard classification, so that each particle '
                           'is only assigned to one class at every iteration, '
                           'rather than having partial assignment to all '
                           'classes.')

        form.addParam('multirefine_batch_size_per_class', IntParam, default=1000,
                      expertLevel=LEVEL_ADVANCED,
                      label="Batch size per class",
                      help='Number of images per class used in each batch (one '
                           'batch per iteration of online-EM). Larger values '
                           'slow the algorithm down but can provide better '
                           'classification results.')

        form.addParam('multirefine_update_rule', StringParam,
                      default="online_em",
                      expertLevel=LEVEL_ADVANCED,
                      label="Optimization method",
                      help='Optimization method to use.')

        form.addParam('multirefine_online_em_lr_rand', FloatParam,
                      default=0.2,
                      expertLevel=LEVEL_ADVANCED,
                      label="O-EM learning rate during randomization",
                      help='Not recommended to change.')

        form.addParam('multirefine_online_em_lr_init', FloatParam,
                      default=0.1,
                      expertLevel=LEVEL_ADVANCED,
                      label="O-EM learning rate init",
                      help='Not recommended to change.')

        form.addParam('multirefine_online_em_lr_hl', IntParam,
                      default=50,
                      expertLevel=LEVEL_ADVANCED,
                      label="O-EM learning rate halflife (iters)",
                      help='Not recommended to change.')

        form.addParam('multirefine_halfmap_decay', FloatParam,
                      default=0.9,
                      expertLevel=LEVEL_ADVANCED,
                      label="Halfmap decay constant",
                      help='Not recommended to change.')

        form.addParam('multirefine_res_init', FloatParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=20,
                      label="Initial resolution (A)",
                      help='Initial low-pass resolution applied to input '
                           'volumes before classification or reconstruction.')

        form.addParam('multirefine_bp_res_factor', FloatParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=1.5,
                      label="Backprojection resolution factor",
                      help='Backproject at this multiple of the best resolution '
                           'amongst classes. Not recommended to change.')

        form.addParam('multirefine_use_max_fsc', BooleanParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=True,
                      label="Use max FSC over classes for filtering",
                      help='Use the maximum FSC across classes to filter all '
                           'classes. This prevents smaller classes from being '
                           'over-filtered during reconstruction.')

        form.addParam('multirefine_assignment_conv_eps', FloatParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=0.05,
                      label="Assignment convergence criteria",
                      help='Fraction of the batch that is allowed to have '
                           'changed classes in the past iteration to be '
                           'considered converged.')

        form.addParam('multirefine_assignment_conv_eps', IntParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=1,
                      label="Resolution convergence criteria",
                      help='Maximum change in resolution (in Fourier shells) '
                           'between iterations to be considered converged.')

        form.addParam('multirefine_num_rand_assign_iters', IntParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=5,
                      label="Number of initial random assignment iterations",
                      help='Number of iterations to perform initially with '
                           'random assignments to break symmetry of multiple '
                           'identical initial references.')

        form.addParam('multirefine_num_final_full_iters', IntParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=2,
                      label="Number of final full iterations",
                      help='Number of final full iterations through the entire '
                           'dataset. Generally 2 iterations is enough.')

        form.addParam('multirefine_noise_model', EnumParam,
                      expertLevel=LEVEL_ADVANCED,
                      choices=['symmetric', 'white', 'coloured'],
                      default=0,
                      label='Noise model',
                      help='Noise model to use. Valid options are white, '
                           'coloured or symmetric. Symmetric is the default, '
                           'meaning coloured with radial symmetry. '
                           'Not recommended to change.')

        form.addParam('multirefine_noise_init_sigmascale', IntParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=3,
                      label="Noise initial sigma-scale",
                      help='Scale factor initially applied to the base noise '
                           'estimate. Not recommended to change.')

        # --------------[Compute settings]---------------------------
        form.addSection(label="Compute settings")
        addComputeSectionParams(form, allowMultipleGPUs=False)

    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._defineFileNames()
        self._defineParamsName()
        self._initializeCryosparcProject()
        self._insertFunctionStep("convertInputStep")
        self._insertFunctionStep('processStep')
        # self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions ------------------------------
    def processStep(self):
        self.vol = self.importVolume.get() + '.imported_volume.map'
        print(pwutils.yellowStr("3D Classification started..."), flush=True)
        self.do3DClasification()

    def _defineParamsName(self):
        """ Define a list with all protocol parameters names"""
        self._paramsName = ['multirefine_N',
                            'multirefine_symmetry',
                            'multirefine_sharp_bfactor',
                            'multirefine_force_hard_class',
                            'multirefine_batch_size_per_class',
                            'multirefine_update_rule',
                            'multirefine_online_em_lr_rand',
                            'multirefine_online_em_lr_init',
                            'multirefine_online_em_lr_hl',
                            'multirefine_halfmap_decay',
                            'multirefine_res_init',
                            'multirefine_bp_res_factor',
                            'multirefine_use_max_fsc',
                            'multirefine_assignment_conv_eps',
                            'multirefine_num_rand_assign_iters',
                            'multirefine_num_final_full_iters',
                            'multirefine_noise_model',
                            'multirefine_noise_init_sigmascale',
                            'compute_use_ssd']
        self.lane = str(self.getAttributeValue('compute_lane'))

    def do3DClasification(self):
        """
        """
        className = "hetero_refine"
        input_group_conect = {"particles": str(self.par),
                              "volume": str(self.vol)}
        # {'particles' : 'JXX.imported_particles' }
        params = {}

        for paramName in self._paramsName:
            if (paramName != 'multirefine_symmetry' and
                    paramName != 'multirefine_noise_model'):
                params[str(paramName)] = str(self.getAttributeValue(paramName))

            elif paramName == 'refine_symmetry':
                symetryValue = getSymmetry(self.symmetryGroup.get(),
                                           self.symmetryOrder.get())

                params[str(paramName)] = symetryValue
            elif paramName == 'multirefine_noise_model':
                params[str(paramName)] = str(NOISE_MODEL_CHOICES[self.multirefine_noise_model.get()])

        # Determinate the GPUs to use (in dependence of
        # the cryosparc version)
        try:
            gpusToUse = self.getGpuList()
        except Exception:
            gpusToUse = False

        self.run3dClassification = enqueueJob(className, self.projectName.get(),
                                    self.workSpaceName.get(),
                                    str(params).replace('\'', '"'),
                                    str(input_group_conect).replace('\'',
                                                                    '"'),
                                    self.lane, gpusToUse)

        self.currenJob.set(self.run3dClassification.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(), self.run3dClassification.get(),
                         "An error occurred in the 3D Classification process. "
                         "Please, go to cryosPARC software for more "
                         "details.")