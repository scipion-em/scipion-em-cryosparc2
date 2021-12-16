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

from pkg_resources import parse_version

from pwem import SCIPION_SYM_NAME

import pyworkflow.utils as pwutils
from pyworkflow import BETA
from pyworkflow.object import String
from pyworkflow.protocol.params import (FloatParam, Positive, IntParam,
                                        BooleanParam, EnumParam, PointerParam)

from .protocol_cryorefine import ProtCryoSparcRefine3D
from ..utils import (getSymmetry, enqueueJob, waitForCryosparc,
                     clearIntermediateResults, addComputeSectionParams,
                     cryosparcValidate, gpusValidate, getCryosparcVersion)
from ..constants import *


class ProtCryoSparcHelicalRefine3D(ProtCryoSparcRefine3D):
    """ Reconstruct and refine a homogeneous helical assembly, with or without
    imposition and refinement of symmetry parameters. Helical Refinement (BETA)
    uses an algorithm that is conceptually similar to Egelman's Iterative
    Helical Real Space Reconstruction (IHRSR) algorithm, while incorporating
    the same maximum likelihood framework, accelerated branch-and-bound
    alignment algorithm, and optional Non-Uniform regularization as used in
    other cryoSPARC refinement jobs.
    """
    _label = '3D helical refinement'
    _devStatus = BETA
    _fscColumns = 4
    _protCompatibility = [V3_0_0, V3_1_0, V3_2_0, V3_3_0, V3_3_1]
    _className = "helix_refine"

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      label="Input particles", important=True,
                      validators=[Positive],
                      help='Particle stacks to use. Multiple stacks will '
                           'be concatenated.')
        form.addParam('referenceVolume', PointerParam, pointerClass='Volume',
                      default=None,
                      allowsNull=True,
                      label="Initial volume",
                      help='Initial volume to use for helical refinement.')
        form.addParam('refMask', PointerParam, pointerClass='VolumeMask',
                      default=None,
                      label='Mask to be applied to this map(Optional)',
                      allowsNull=True,
                      help='Initial volume mask raw data.')

        form.addSection(label='Helical Refinement')

        form.addParam('refine_init_twist', FloatParam,
                      default=None,
                      allowsNull=True,
                      label="Helical twist estimate (degrees)",
                      help='Angular distance between adjacent subunits, in the '
                           'range (-180,180). Positive and negative values '
                           'correspond to right and left-handed helices, '
                           'respectively.')

        form.addParam('refine_init_shift', FloatParam,
                      default=None,
                      allowsNull=True,
                      label="Helical rise estimate (A)",
                      help='Positive non-zero translation distance (along '
                           'helical axis) between adjacent subunits.')

        form.addParam('refine_limit_shifts', BooleanParam,
                      default=True,
                      label="Limit shifts along the helical axis",
                      help='Limit alignment shifts along the meridian to +/- '
                           '0.5 helical rises, once all particles are seen. '
                           'For fairly rigid helices, can improve resolution'
                           ' by ensuring symmetry averaging occurs over only '
                           'central asymmetric units. Not recommended for '
                           'helices with significant flexibility.')

        form.addParam('refine_hsym_order', FloatParam,
                      default=None,
                      allowsNull=True,
                      label="Maximum symmetry order to apply during reconstruction",
                      help='The maximum amount of helical symmetry to impose '
                           'during reconstruction; i.e., the maximum number of '
                           '(twist, rise) pairs to backproject each particle '
                           'image with. For particles picked outside of the '
                           'filament tracer or template picker, this should be '
                           'set to the distance between extracted boxes, '
                           'divided by the helical rise. If left as None, '
                           'will be calculated based on the inter-box distance '
                           '(if available). Set to 1 for no helical Fourier'
                           ' space symmetrization.')

        form.addParam('refine_sym_enforce_r', FloatParam,
                      default=8,
                      validators=[Positive],
                      label="Resolution to begin real-space symmetrization",
                      help='At what resolution (A) to begin enforce symmetry '
                           'in real-space prior to alignment; may improve'
                           ' particle alignments. Set to 0 to disable '
                           'real-space symmetry enforcement. If symmetry '
                           'parameters are being searched, it\'s recommended '
                           'to set this to a fairly high resolution '
                           '(e.g. ~5 A).')

        form.addParam('symmetryGroup', EnumParam,
                      choices=[CS_SYM_NAME[SYM_CYCLIC] +
                               " (" + SCIPION_SYM_NAME[SYM_CYCLIC] + ")",
                               CS_SYM_NAME[SYM_DIHEDRAL_Y] +
                               " (" + SCIPION_SYM_NAME[SYM_DIHEDRAL_Y] + ")"],
                      default=SYM_CYCLIC,
                      label="Symmetry",
                      help="Symmetry String (C, D). E.g. C1, D7, C4, etc. "
                           "Only cyclic and dihedral symmetries are supported."
                      )

        form.addParam('symmetryOrder', IntParam, default=1,
                      condition='symmetryGroup==%d or symmetryGroup==%d' %
                                (SYM_DIHEDRAL_Y - 1, SYM_CYCLIC),
                      label='Point group symmetry',
                      validators=[Positive],
                      help='Order of symmetry.')

        form.addSection(label='Non-Uniform Refinement')

        form.addParam('nu_refine', BooleanParam,
                      default=False,
                      label="Use Non-Uniform Refinement?",
                      help='Use Non-Uniform regularization during refinement '
                           'to achieve higher resolution and map quality.')

        form.addSection(label='Initial Model')

        form.addParam('refine_res_init', FloatParam,
                      default=20,
                      validators=[Positive],
                      label="Initial lowpass resolution (A)",
                      help='Lowpass filter resolution applied to input '
                           'structure. Values between 15 and 35 Angstroms '
                           'may produce best results.')

        form.addParam('refine_initmodel_numimages', FloatParam,
                      default=5000,
                      validators=[Positive],
                      label="Number of images for initial density generation",
                      help='Number of images used in generating the initial '
                           'density.')

        form.addParam('use_cylindrical_model', BooleanParam,
                      default=False,
                      label="Generate a cylindrical initial model?",
                      help='Whether or not to generate a cylindrical initial model')

        form.addParam('filament_outer_diameter', FloatParam,
                      default=None,
                      allowsNull=True,
                      label="Filament Outer Diameter (Angstrom)",
                      help='Approximate outer diameter of the filament in Angstroms')

        form.addParam('filament_inner_diameter', FloatParam,
                      default=0,
                      label="Filament Inner Diameter (Angstrom)",
                      help='Approximate inner diameter of the filament in Angstroms')

        form.addParam('filament_far_dist_A', FloatParam,
                      default=6,
                      label="Far distance (Angstrom)",
                      help='Distance over which the model is padded, with voxel values fading to 0.')

        form.addSection(label='Refinement')

        form.addParam('refine_res_align_max', FloatParam,
                      default=None,
                      allowsNull=True,
                      label="Maximum align resolution (A)",
                      help='Manual override for maximum resolution that is '
                           'used for alignment. This value is normally '
                           'set by the GS-FSC')

        form.addParam('refine_res_gsfsc_split', FloatParam,
                      default=20,
                      validators=[Positive],
                      label="GSFSC split resolution (A)",
                      help='Resolution beyond which two GS-FSC halves are '
                           'independent')

        form.addParam('refine_mask', EnumParam,
                      choices=["dynamic", "static"],
                      default=0,
                      label="Mask",
                      help='Type of masking to use. Either "dynamic", or '
                           '"static".')

        form.addParam('refine_dynamic_mask_thresh_factor', FloatParam,
                      default=0.3,
                      validators=[Positive],
                      label="Dynamic mask threshold (0-1)",
                      help='Level set threshold for selecting regions that are '
                           'included in the dynamic mask.')

        # --------------[Compute settings]---------------------------
        form.addSection(label="Compute settings")
        addComputeSectionParams(form, allowMultipleGPUs=False)

    def _insertAllSteps(self):
        ProtCryoSparcRefine3D._insertAllSteps(self)

    def processStep(self):
        print(pwutils.yellowStr("Refinement started..."), flush=True)
        self.doRunRefine()

    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        validateMsgs = cryosparcValidate()
        if not validateMsgs:
            csVersion = getCryosparcVersion()
            if [version for version in self._protCompatibility
                if parse_version(version) >= parse_version(csVersion)]:
                validateMsgs = gpusValidate(self.getGpuList(), checkSingleGPU=True)
                if not validateMsgs:
                    if self.referenceVolume.get() is None and not self.use_cylindrical_model.get():
                        validateMsgs.append("Cannot generate initial model "
                                            "without in-plane rotation "
                                            "information. Please input an "
                                            "initial model from a previous ab-initio or "
                                            "refinement protocol, or activate the "
                                            "'Generate a cylindrical initial "
                                            "model?' parameter")

                    if self.use_cylindrical_model.get() and self.filament_outer_diameter.get() is None:
                        validateMsgs.append("Must set the filament outer diameter to use a cylindrical model")

            else:
                validateMsgs.append("The protocol is not compatible with the "
                                    "cryoSPARC version %s" % csVersion)
        return validateMsgs

    def _defineParamsName(self):
        """ Define a list with all protocol parameters names"""
        self._paramsName = ['refine_init_twist',
                            'refine_init_shift',
                            'refine_hsym_order',
                            'refine_limit_shifts',
                            'refine_sym_enforce_r',
                            'refine_pg_symmetry',
                            'nu_refine',
                            'refine_res_init',
                            'refine_initmodel_numimages',
                            'use_cylindrical_model',
                            'refine_res_align_max',
                            'refine_res_gsfsc_split',
                            'refine_mask',
                            'refine_dynamic_mask_thresh_factor',
                            'filament_outer_diameter',
                            'filament_inner_diameter',
                            'filament_far_dist_A',
                            'compute_use_ssd']
        self.lane = str(self.getAttributeValue('compute_lane'))

    def doRunRefine(self):
        input_group_connect = {"particles": self.particles.get()}
        if self.volume.get() is not None:
            input_group_connect["volume"] = self.volume.get()
        if self.mask.get() is not None:
            input_group_connect["mask"] = self.mask.get()
        params = {}

        for paramName in self._paramsName:
            if (paramName != 'refine_hsym_order' and
                    paramName != 'refine_pg_symmetry' and
                    paramName != "refine_init_twist" and
                    paramName != "refine_init_shift" and
                    paramName != 'filament_outer_diameter' and
                    paramName != 'refine_res_align_max'):
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif paramName == 'refine_pg_symmetry':
                symetryValue = getSymmetry(self.symmetryGroup.get(),
                                           self.symmetryOrder.get())
                params[str(paramName)] = symetryValue
            elif self.getAttributeValue(paramName) is not None and float(self.getAttributeValue(paramName)) > 0:
                params[str(paramName)] = str(self.getAttributeValue(paramName))

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