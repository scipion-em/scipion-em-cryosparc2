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

import pwem.objects as pwobj
import pwem.protocols as pwprot
import pyworkflow.utils as pwutils
from pyworkflow import NEW
from pyworkflow.protocol.params import *

from .protocol_base import ProtCryosparcBase
from ..convert import (defineArgs, convertCs2Star, createItemMatrix,
                       setCryosparcAttributes)
from ..utils import (addSymmetryParam, addComputeSectionParams,
                     calculateNewSamplingRate,
                     cryosparcValidate, gpusValidate, getSymmetry,
                     waitForCryosparc, clearIntermediateResults, enqueueJob,
                     getCryosparcVersion, fixVolume, copyFiles)
from ..constants import (md, NOISE_MODEL_CHOICES, REFINE_MASK_CHOICES, V3_0_0,
                         V3_1_0, V3_2_0, V3_3_0, V3_3_1, REFINE_FILTER_TYPE)


class ProtCryoSparc3DHomogeneousRefine(ProtCryosparcBase, pwprot.ProtRefine3D):
    """ Protocol to refine a 3D map using cryosparc.
        Rapidly refine a single homogeneous structure to high-resolution and
        validate using the gold-standard FSC.
    """
    _label = '3D homogeneous refinement'
    _devStatus = NEW
    _fscColumns = 6
    _className = "homo_refine_new"
    _protCompatibility = [V3_0_0, V3_1_0, V3_2_0, V3_3_0, V3_3_1]

    # --------------------------- DEFINE param functions ----------------------
    def _defineFileNames(self):
        """ Centralize how files are called within the protocol. """
        myDict = {
                  'input_particles': self._getTmpPath('input_particles.star'),
                  'out_particles': self._getPath() + '/output_particle.star',
                  'stream_log': self._getPath()+'/stream.log'
                  }
        self._updateFilenamesDict(myDict)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      label="Input particles", important=True,
                      validators=[Positive],
                      help='Select the input images from the project.')
        form.addParam('referenceVolume', PointerParam, pointerClass='Volume',
                      important=True,
                      label="Input volume",
                      help='Initial reference 3D map, it should have the same '
                           'dimensions and the same pixel size as your input '
                           'particles.')
        form.addParam('refMask', PointerParam, pointerClass='VolumeMask',
                      default=None,
                      label='Mask to be applied to this map(Optional)',
                      allowsNull=True,
                      help='A volume mask containing a (soft) mask with '
                           'the same dimensions as the reference(s), '
                           'and values between 0 and 1, with 1 being 100% '
                           'protein and 0 being 100% solvent. The '
                           'reconstructed reference map will be multiplied '
                           'by this mask. If no mask is given, a soft '
                           'spherical mask based on the <radius> of the '
                           'mask for the experimental images will be '
                           'applied.')

        # --------------[Homogeneous Refinement]---------------------------
        form.addSection(label='Homogeneous Refinement')
        addSymmetryParam(form, help="Symmetry String (C, D, I, O, T). E.g. C1, "
                                    "D7, C4, etc")

        form.addParam('refine_symmetry_do_align', BooleanParam, default=True,
                      label="Do symmetry alignment",
                      help='Align the input structure to the symmetry axes')

        form.addParam('refine_do_init_scale_est', BooleanParam, default=True,
                      label="Re-estimate greyscale level of input reference")

        form.addParam('refine_num_final_iterations', IntParam, default=0,
                      expertLevel=LEVEL_ADVANCED,
                      label="Number of extra final passes",
                      help='Number of extra passes through the data to do '
                           'after the GS-FSC resolution has stopped improving')

        form.addParam('refine_res_init', IntParam, default=30,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label="Initial lowpass resolution (A)",
                      help='Applied to input structure')

        form.addParam('refine_res_gsfsc_split', IntParam, default=20,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label="GSFSC split resolution (A)",
                      help='Resolution beyond which two GS-FSC halves are '
                           'independent')

        form.addParam('refine_highpass_res', IntParam, default=None,
                      allowsNull=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Highpass resolution (A)")

        form.addParam('refine_clip', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Enforce non-negativity",
                      help='Clip negative density. Probably should be false')

        form.addParam('refine_window', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Skip interpolant premult",
                      help='Softly window the structure in real space with a '
                           'spherical window. Should be true')

        form.addParam('refine_skip_premult', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Window structure in real space",
                      help='Leave this as true')

        form.addParam('refine_ignore_dc', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Ignore DC component",
                      help='Ignore the DC component of images. Should be true')

        form.addParam('refine_batchsize_init', IntParam, default=0,
                      expertLevel=LEVEL_ADVANCED,
                      label="Initial batchsize",
                      help='Number of images used in the initial iteration. '
                           'Set to zero to autotune')

        form.addParam('refine_batchsize_epsilon', FloatParam, default=0.001,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label="Batchsize epsilon",
                      help='Controls batch size when autotuning batchsizes. '
                           'Set closer to zero for larger batches')

        form.addParam('refine_batchsize_snrfactor', FloatParam, default=40.0,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label="Batchsize snrfactor",
                      help='Specifies the desired improvement in SNR from the '
                           'images when autotuning batchsizes. Directly '
                           'multiplies the number of images in the batch')

        form.addParam('refine_scale_min', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Minimize over per-particle scale")

        form.addParam('refine_scale_start_iter', IntParam, default=0,
                      expertLevel=LEVEL_ADVANCED,
                      label="Scale min/use start iter",
                      help='Iteration to start minimizing over per-particle scale')

        form.addParam('refine_noise_model', EnumParam,
                      choices=['symmetric', 'white', 'coloured'],
                      default=0,
                      label="Noise model:",
                      help='Noise model to be used. Valid options are white, '
                           'coloured or symmetric. Symmetric is the default, '
                           'meaning coloured with radial symmetry')

        form.addParam('refine_noise_priorw', IntParam, default=50,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label="Noise priorw",
                      help='Weight of the prior for estimating noise (units of '
                           '# of images)')

        form.addParam('refine_noise_initw', IntParam, default=200,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label="Noise initw",
                      help='Weight of the initial noise estimate (units of # '
                           'of images)')

        form.addParam('refine_noise_init_sigmascale', IntParam, default=3,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label="Noise initial sigma-scale",
                      help='Scale factor initially applied to the base noise '
                           'estimate')

        form.addParam('refine_mask', EnumParam,
                      choices=['dynamic', 'static', 'null'],
                      default=0,
                      expertLevel=LEVEL_ADVANCED,
                      label="Mask:",
                      help='Type of masking to use. Either "dynamic", '
                           '"static", or "null"')

        form.addParam('refine_dynamic_mask_thresh_factor', FloatParam,
                      default=0.2,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label="Dynamic mask threshold (0-1)",
                      help='Level set threshold for selecting regions that are '
                           'included in the dynamic mask. Probably don\'t need '
                           'to change this')

        form.addParam('refine_dynamic_mask_near_ang', FloatParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=6.0,
                      validators=[Positive],
                      label="Dynamic mask near (A)",
                      help='Controls extent to which mask is expanded. At the '
                           'near distance, the mask value is 1.0 (in A)')

        form.addParam('refine_dynamic_mask_far_ang', FloatParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=14.0,
                      validators=[Positive],
                      label="Dynamic mask far (A)",
                      help='Controls extent to which mask is expanded. At the '
                           'far distance the mask value becomes 0.0 (in A)')

        form.addParam('refine_dynamic_mask_start_res', IntParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=12,
                      validators=[Positive],
                      label="Dynamic mask start resolution (A)",
                      help='Map resolution at which to start dynamic masking '
                           '(in A)')

        form.addParam('refine_dynamic_mask_use_abs', BooleanParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=False,
                      label="Dynamic mask use absolute value",
                      help='Include negative regions if they are more negative '
                           'than the threshold')

        form.addSection(label='Defocus Refinement')
        form.addParam('refine_defocus_refine', BooleanParam,
                      default=True,
                      label="Optimize per-particle defocus",
                      help='Minimize over per-particle defocus at each '
                           'iteration of refinement. The optimal defocus will'
                           ' be used for backprojection as well, and will be '
                           'written out at each iteration. Defocus refinement '
                           'will start only once refinement with current '
                           'defocus values converges. Beware that with '
                           'small/disordered proteins, defocus refinement '
                           'may actually make resolutions worse.')

        form.addParam('crl_num_plots', IntParam,
                      default=3,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label="Num. particles to plot",
                      help='Number of particles to make plots for. After '
                           'this many, stop plotting to save time.')

        form.addParam('crl_min_res_A', FloatParam,
                      default=20,
                      condition='refine_defocus_refine==True',
                      validators=[Positive],
                      label="Minimum Fit Res (A)",
                      help='The minimum resolution to use during refinement '
                           'of image aberrations.')

        form.addParam('crl_df_range', FloatParam,
                      default=2000,
                      condition='refine_defocus_refine==True',
                      validators=[Positive],
                      label="Defocus Search Range (A +/-)",
                      help='Defocus search range in Angstroms, searching '
                           'both above and below the input defocus by this '
                           'amount')

        form.addSection(label='Global CTF Refinement')

        form.addParam('refine_ctf_global_refine', BooleanParam,
                      default=True,
                      label="Optimize per-group CTF params",
                      help='Optimize the per-exposure-group CTF parameters '
                           '(for higher-order aberrations) at each iteration '
                           'of refinement. The optimal CTF will be used for '
                           'backprojection as well, and will be written out '
                           'at each iteration. CTF refinement will start only '
                           'once refinement with current CTF values converges. '
                           'Beware that with small/disordered proteins, CTF'
                           ' refinement may actually make resolutions worse.')

        form.addParam('crg_num_plots', IntParam,
                      default=3,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label="Num. groups to plot",
                      help='Number of exposure groups to make plots for. '
                           'After this many, stop plotting to save time.')

        form.addParam('crg_min_res_A', FloatParam,
                      default=10,
                      condition="refine_ctf_global_refine == True",
                      validators=[Positive],
                      label="Minimum Fit Res (A)",
                      help='The minimum resolution to use during refinement '
                           'of image aberrations.')

        form.addParam('crg_do_tilt', BooleanParam,
                      condition="refine_ctf_global_refine == True",
                      default=True,
                      label="Fit Tilt",
                      help='Whether to fit beam tilt.')

        form.addParam('crg_do_trefoil', BooleanParam,
                      condition="refine_ctf_global_refine == True",
                      default=True,
                      label="Fit Trefoil",
                      help='Whether to fit beam trefoil.')

        form.addParam('crg_do_spherical', BooleanParam,
                      condition="refine_ctf_global_refine == True",
                      default=True,
                      label="Fit Spherical Aberration",
                      help='Whether to fit spherical aberration.')

        form.addParam('crg_do_tetrafoil', BooleanParam,
                      condition="refine_ctf_global_refine == True",
                      default=True,
                      label="Fit Tetrafoil",
                      help='Whether to fit beam tetrafoil.')

        # --------------[Compute settings]---------------------------
        form.addSection(label="Compute settings")
        addComputeSectionParams(form, allowMultipleGPUs=True)

    # --------------------------- INSERT steps functions -----------------------

    def _insertAllSteps(self):
        self._defineFileNames()
        self._defineParamsName()
        self._initializeCryosparcProject()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.processStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions ------------------------------
    def _getInputVolume(self):
        return self.referenceVolume.get()

    def processStep(self):
        print(pwutils.yellowStr("Refinement started..."), flush=True)
        self.doRunRefine()

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to Relion file
        """
        self._initializeUtilsVariables()
        idd, itera = self.findLastIteration(self.runRefine.get())
        csOutputFolder = os.path.join(self.projectPath, self.projectName.get(),
                                      self.runRefine.get())
        csOutputPattern = "cryosparc_%s_%s_%s" % (self.projectName.get(),
                                                  self.runRefine.get(),
                                                  itera)
        csParticlesName = csOutputPattern + "_particles.cs"

        fnVolName = csOutputPattern + "_volume_map.mrc"
        half1Name = csOutputPattern + "_volume_map_half_A.mrc"
        half2Name = csOutputPattern + "_volume_map_half_B.mrc"

        # Copy the CS output volume and half to extra folder
        copyFiles(csOutputFolder, self._getExtraPath(), files=[csParticlesName,
                                                               fnVolName,
                                                               half1Name,
                                                               half2Name])

        csFile = os.path.join(self._getExtraPath(), csParticlesName)

        outputStarFn = self._getFileName('out_particles')
        argsList = [csFile, outputStarFn]

        parser = defineArgs()
        args = parser.parse_args(argsList)
        convertCs2Star(args)

        fnVol = os.path.join(self._getExtraPath(), fnVolName)
        half1 = os.path.join(self._getExtraPath(), half1Name)
        half2 = os.path.join(self._getExtraPath(), half2Name)

        imgSet = self._getInputParticles()
        vol = pwobj.Volume()
        fixVolume([fnVol, half1, half2])
        vol.setFileName(fnVol)
        vol.setSamplingRate(calculateNewSamplingRate(vol.getDim(),
                                                     imgSet.getSamplingRate(),
                                                     imgSet.getDim()))
        vol.setHalfMaps([half1, half2])

        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(imgSet)
        self._fillDataFromIter(outImgSet)

        self._defineOutputs(outputVolume=vol)
        self._defineSourceRelation(self.inputParticles.get(), vol)
        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(self.inputParticles.get(), outImgSet)
        self.createFSC(idd, imgSet, vol)

    def _validate(self):
        validateMsgs = cryosparcValidate()
        if not validateMsgs:
            csVersion = getCryosparcVersion()
            if [version for version in self._protCompatibility
                if parse_version(version) >= parse_version(csVersion)]:
                validateMsgs = gpusValidate(self.getGpuList())
            else:
                validateMsgs.append("The protocol is not compatible with the "
                                    "cryoSPARC version %s" % csVersion)
        return validateMsgs

    def _summary(self):
        summary = []
        if (not hasattr(self, 'outputVolume') or
                not hasattr(self, 'outputParticles') or
                not hasattr(self, 'outputFSC')):
            summary.append("Output objects not ready yet.")
        else:
            summary.append("Input Particles: %s" %
                           self.getObjectTag('inputParticles'))
            summary.append("Input Volume: %s" %
                           self.getObjectTag('referenceVolume'))
            summary.append("Input Mask: %s" %
                           self.getObjectTag('refMask'))
            summary.append("Symmetry: %s" %
                           getSymmetry(self.symmetryGroup.get(),
                                       self.symmetryOrder.get())
                           )
            summary.append("------------------------------------------")
            summary.append("Output particles %s" %
                           self.getObjectTag('outputParticles'))
            summary.append("Output volume %s" %
                           self.getObjectTag('outputVolume'))

            if self.hasAttribute('mapResolution'):
                summary.append(
                    "\nMap Resolution: %s" % self.mapResolution.get())
            if self.hasAttribute('estBFactor'):
                summary.append(
                    '\nEstimated Bfactor: %s' % self.estBFactor.get())
        return summary

    # -------------------------- UTILS functions ------------------------------

    def _fillDataFromIter(self, imgSet):
        outImgsFn = 'particles@' + self._getFileName('out_particles')
        imgSet.setAlignmentProj()
        imgSet.copyItems(self._getInputParticles(),
                         updateItemCallback=self._createItemMatrix,
                         itemDataIterator=md.iterRows(outImgsFn,
                                                      sortByLabel=md.RLN_IMAGE_ID))

    def _createItemMatrix(self, particle, row):
        createItemMatrix(particle, row, align=pwobj.ALIGN_PROJ)
        setCryosparcAttributes(particle, row,
                               md.RLN_PARTICLE_RANDOM_SUBSET)

    def _defineParamsName(self):
        """ Define a list with all protocol parameters names"""
        self._paramsName = ['refine_symmetry',
                            'refine_symmetry_do_align',
                            'refine_do_init_scale_est',
                            'refine_highpass_res',
                            'refine_num_final_iterations',
                            'refine_res_init',
                            'refine_res_gsfsc_split',
                            'refine_clip',
                            'refine_window', 'refine_skip_premult',
                            'refine_ignore_dc',
                            'refine_batchsize_init',
                            'refine_batchsize_snrfactor',
                            'refine_batchsize_epsilon',
                            'refine_scale_min',
                            'refine_scale_start_iter',
                            'refine_noise_model', 'refine_noise_priorw',
                            'refine_noise_initw',
                            'refine_noise_init_sigmascale',
                            'refine_mask',
                            'refine_dynamic_mask_thresh_factor',
                            'refine_dynamic_mask_near_ang',
                            'refine_dynamic_mask_far_ang',
                            'refine_dynamic_mask_start_res',
                            'refine_dynamic_mask_use_abs',
                            'refine_defocus_refine', 'crl_num_plots',
                            'crl_min_res_A', 'crl_df_range',
                            'refine_ctf_global_refine',
                            'crg_num_plots', 'crg_min_res_A', 'crg_do_tilt',
                            'crg_do_trefoil', 'crg_do_spherical',
                            'crg_do_tetrafoil',
                            'compute_use_ssd']
        self.lane = str(self.getAttributeValue('compute_lane'))

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
                    paramName != 'refine_mask' and
                    paramName != 'refine_highpass_res' and
                    paramName != 'refine_nu_filtertype'):
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif (paramName == 'refine_highpass_res' and self.getAttributeValue(paramName) is not None and
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
            elif paramName == 'refine_nu_filtertype':
                params[str(paramName)] = str(
                    REFINE_FILTER_TYPE[self.refine_nu_filtertype.get()])

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

        self.runRefine = pwobj.String(runRefineJob.get())
        self.currenJob.set(self.runRefine.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(), self.runRefine.get(),
                         "An error occurred in the Refinement process. "
                         "Please, go to cryosPARC software for more "
                         "details.")
        clearIntermediateResults(self.projectName.get(), self.runRefine.get())



