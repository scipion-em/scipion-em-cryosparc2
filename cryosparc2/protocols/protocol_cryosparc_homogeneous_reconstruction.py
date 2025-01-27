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
import ast
import os

import emtable
from pkg_resources import parse_version

from pwem import ALIGN_PROJ

import pyworkflow.utils as pwutils
from pyworkflow import NEW
from pyworkflow.object import String
from pyworkflow.protocol.params import (PointerParam, LEVEL_ADVANCED, IntParam,
                                        BooleanParam, EnumParam, FloatParam)
from pwem.objects import Volume

from .protocol_base import ProtCryosparcBase
from ..convert import (convertCs2Star, createItemMatrix,
                       setCryosparcAttributes)
from ..utils import (addComputeSectionParams, calculateNewSamplingRate,
                     cryosparcValidate, gpusValidate, enqueueJob,
                     waitForCryosparc, clearIntermediateResults, fixVolume,
                     copyFiles, addSymmetryParam, getSymmetry,
                     getCryosparcVersion, get_job_streamlog, getOutputPreffix)
from ..constants import *


class ProtCryoSparcHomogeneousReconstruct(ProtCryosparcBase):
    """ Create a 3D reconstruction from input particles that already have alignments in 3D.
    """
    _label = 'homogeneous reconstruction'
    _className = "homo_reconstruct"
    _devStatus = NEW
    _fscColumns = 6
    _protCompatibility = [V3_3_0, V3_3_1, V3_3_2, V4_0_0, V4_0_1, V4_0_2,
                          V4_0_3, V4_1_0, V4_1_1, V4_1_2, V4_2_0, V4_2_1,
                          V4_3_1, V4_4_0, V4_4_1, V4_5_1, V4_5_3, V4_6_0, V4_6_1, V4_6_2]
    ewsParamsName = []

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
                      help='Select the experimental particles.')
        form.addParam('refMask', PointerParam, pointerClass='VolumeMask',
                      label='Mask used for computing FSC',
                      allowsNull=True,
                      help="The mask is used for computing"
                           "the FSC between the reconstructed half-maps")

        form.addSection(label='Homogeneous Reconstruction')

        form.addParam('refine_N_cmp', IntParam, default=None,
                      allowsNull=True,
                      label="Reconstruction box size (voxels)",
                      help='The volume size to use for refinement. Default: raw particle size')

        addSymmetryParam(form, help="Symmetry String (C, D, I, O, T). E.g. C1, "
                                    "D7, C4, etc")

        form.addParam('refine_helical_twist', IntParam, default=None,
                      allowsNull=True,
                      label="Helical twist (degrees)",
                      help='Helical twist to apply in degrees')

        form.addParam('refine_helical_shift', IntParam, default=None,
                      allowsNull=True,
                      label="Helical rise (Angstroms)",
                      help='Helical rise to apply in Angstroms')

        form.addParam('refine_hsym_order', IntParam, default=None,
                      allowsNull=True,
                      label="Helical symmetry order to apply",
                      help='The maximum amount of helical symmetry to impose '
                           'during reconstruction; i.e., the maximum number '
                           'of (twist, rise) pairs to backproject each '
                           'particle image with. For particles picked outside '
                           'of the filament tracer or template picker, this '
                           'should be set to the distance between extracted '
                           'boxes, divided by the helical rise. If left as '
                           'None, will be calculated based on the inter-box '
                           'distance (if available). Set to 1 for no helical '
                           'symmetry during reconstruction.')

        form.addParam('refine_gs_resplit', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Force re-do GS split",
                      help='Force re-splitting the particles into two random '
                           'gold-standard halves. If this is not set, split is '
                           'preserved from input alignments')

        form.addParam('refine_do_flip', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Flip the reconstruction hand",
                      help='Flip the hand of the reconstruction, and correct '
                           'the alignments. If flipping the hand and applying '
                           'helical symmetry, ensure that the sign of the '
                           'helical twist is inverted to match the hand of '
                           'the new symmetry.')

        form.addParam('refine_ignore_tilt', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Ignore the tilt",
                      help='Ignore the tilt')


        form.addParam('refine_ignore_trefoil', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Ignore trefoil",
                      help='Ignore trefoil')

        form.addParam('refine_ignore_tetra', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Ignore tetra",
                      help='Ignore tetra')

        form.addParam('refine_ignore_anisomag', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Ignore anisomag",
                      help='Ignore the anisomag')

        csVersion = getCryosparcVersion()
        if parse_version(csVersion) >= parse_version(V3_3_1):

            form.addParam('refine_do_ews_correct', BooleanParam, default=False,
                          expertLevel=LEVEL_ADVANCED,
                          label="Do EWS correction",
                          help='Whether or not to correct for the curvature of the '
                               'Ewald Sphere.')

            form.addParam('refine_ews_zsign', EnumParam,
                          choices=['positive', 'negative'],
                          default=0,
                          expertLevel=LEVEL_ADVANCED,
                          label="EWS curvature sign",
                          help='Whether to use positive or negative curvature in '
                               'Ewald Sphere correction.')

            form.addParam('refine_ews_simple', EnumParam,
                          choices=['simple', 'iterative'],
                          default=0,
                          expertLevel=LEVEL_ADVANCED,
                          label="EWS correction method",
                          help='Whether to use the simple insertion method, or to '
                               'use an iterative optimization method, for Ewald '
                               'Sphere correction.')

            form.addParam('refine_fsc_mask_opt', BooleanParam, default=False,
                          expertLevel=LEVEL_ADVANCED,
                          label="Optimize FSC mask",
                          help='Whether or not to optimize the mask used for '
                               'calculating FSCs')

            form.addParam('refine_override_filter', BooleanParam, default=False,
                          expertLevel=LEVEL_ADVANCED,
                          label="Override FSC Filtering",
                          help='Whether to override the FSC-filtering and use '
                               'manual filtering and sharpening (set below) '
                               'instead. FSC is still computed, just not used '
                               'for filtering.')

            form.addParam('refine_override_filter_res', FloatParam, default=None,
                          allowsNull=True,
                          expertLevel=LEVEL_ADVANCED,
                          label="Override Filtering Resolution",
                          help='Override filter corner resolution (A)')

            form.addParam('refine_override_filter_order', IntParam, default=8,
                          expertLevel=LEVEL_ADVANCED,
                          label="Override Filtering Order",
                          help='Override filter order (Butterworth)')

            form.addParam('refine_override_filter_bfactor', FloatParam, default=0,
                          expertLevel=LEVEL_ADVANCED,
                          label="Override Filtering Bfactor",
                          help='Override sharpening B-factor. Negative to sharpen')

            self.ewsParamsName = ['refine_do_ews_correct',
                                  'refine_ews_zsign', 'refine_fsc_mask_opt',
                                  'refine_override_filter',
                                  'refine_override_filter_res',
                                  'refine_ews_simple',
                                  'refine_override_filter_order',
                                   'refine_override_filter_bfactor']

        form.addParam('recon_do_expand_mask', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Wide mask final output",
                      help='Apply a wide mask to the final map, to reduce file '
                           'size after compression.')

        form.addParam('intermediate_plots', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Show plots from intermediate steps",
                      help='Hide plots from intermediate steps to speed up '
                           'processing.')

        form.addParam('refine_compute_batch_size', IntParam, default=None,
                      expertLevel=LEVEL_ADVANCED,
                      allowsNull=True,
                      label="GPU batch size of images",
                      help='Batch size of images to process at a time on the '
                           'GPU. If you run out of GPU memory, try setting '
                           'this to a small number to override the '
                           'auto-detect procedure.')

        # --------------[Compute settings]---------------------------
        form.addSection(label="Compute settings")
        addComputeSectionParams(form, allowMultipleGPUs=False)

    # --------------------------- INSERT steps functions -----------------------

    def _insertAllSteps(self):
        self._defineFileNames()
        self._defineParamsName()
        self._initializeCryosparcProject()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.processStep)
        self._insertFunctionStep(self.createOutputStep)

    def processStep(self):
        self.info(pwutils.yellowStr("Homogeneous Reconstruction started..."))
        self.doHomogeneousReconstruction()

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to Relion file
        """
        self._initializeUtilsVariables()
        idd = self.findLastIteration(self.runHomogeneousReconstruction.get())

        csOutputFolder = os.path.join(self.projectDir.get(),
                                       self.runHomogeneousReconstruction.get())

        csOutputPattern = "%s%s" % (getOutputPreffix(self.projectName.get()),
                                    self.runHomogeneousReconstruction.get())
        csParticlesName = csOutputPattern + "_particles.cs"
        fnVolName = csOutputPattern + "_volume_map.mrc"
        half1Name = csOutputPattern + "_volume_map_half_A.mrc"
        half2Name = csOutputPattern + "_volume_map_half_B.mrc"

        # Copy the CS output to extra folder
        copyFiles(csOutputFolder, self._getExtraPath(), files=[csParticlesName,
                                                               fnVolName,
                                                               half1Name,
                                                               half2Name])

        csFile = os.path.join(self._getExtraPath(), csParticlesName)
        outputStarFn = self._getFileName('out_particles')
        argsList = [csFile, outputStarFn]
        convertCs2Star(argsList)

        fnVol = os.path.join(self._getExtraPath(), fnVolName)
        half1 = os.path.join(self._getExtraPath(), half1Name)
        half2 = os.path.join(self._getExtraPath(), half2Name)
        imgSet = self._getInputParticles()
        vol = Volume()
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
        self._defineSourceRelation(self.inputParticles, vol)
        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(self.inputParticles, outImgSet)
        self.createFSC(idd, imgSet, vol)

    def _fillDataFromIter(self, imgSet):
        outImgsFn = 'particles@' + self._getFileName('out_particles')
        imgSet.setAlignmentProj()
        imgSet.copyItems(self._getInputParticles(),
                         updateItemCallback=self._createItemMatrix,
                         itemDataIterator=emtable.Table.iterRows(outImgsFn))

    def _createItemMatrix(self, particle, row):
        createItemMatrix(particle, row, align=ALIGN_PROJ)
        setCryosparcAttributes(particle, row,
                               RELIONCOLUMNS.rlnRandomSubset.value)

    def findLastIteration(self, jobName):
        get_job_streamlog(self.projectName.get(),
                          jobName,
                          self._getFileName('stream_log'))

        # Get the metadata information from stream.log
        with open(self._getFileName('stream_log')) as f:
            data = f.readlines()

        x = ast.literal_eval(data[0])

        # Find the ID of last iteration and the map resolution
        for y in x:
            if 'text' in y:
                z = str(y['text'])

                if z.startswith('FSC, after mask auto-tightening'):
                    idd = y['imgfiles'][2]['fileid']

        return idd

    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        """ Should be overwritten in subclasses to
               return summary message for NORMAL EXECUTION.
               """
        validateMsgs = cryosparcValidate()
        if not validateMsgs:
            validateMsgs = gpusValidate(self.getGpuList(), checkSingleGPU=True)
            if not validateMsgs:
                particles = self._getInputParticles()
                if not particles.hasCTF():
                    validateMsgs.append("The Particles has not associated a "
                                        "CTF model")
                    if not validateMsgs and not particles.hasAlignment3D():
                        validateMsgs.append("The Particles has not a 3D "
                                            "alignment")

        return validateMsgs

    def _summary(self):
        summary = []
        if (not hasattr(self, 'outputVolume') or
                not hasattr(self, 'outputParticles')):
            summary.append("Output objects not ready yet.")
        else:
            summary.append("Input Particles: %s" %
                           self.getObjectTag('inputParticles'))
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

    def _defineParamsName(self):
        """ Define a list with all protocol parameters names"""
        self._paramsName = ['refine_N_cmp', 'refine_gs_resplit',
                            'refine_do_flip', 'refine_ignore_tilt',
                            'refine_ignore_trefoil', 'refine_ignore_tetra',
                            'refine_ignore_anisomag', 'refine_fsc_mask_opt',
                            'recon_do_expand_mask',
                            'intermediate_plots', 'refine_symmetry',
                            'refine_compute_batch_size', 'refine_helical_twist',
                            'refine_helical_shift', 'refine_hsym_order',
                            'compute_use_ssd'] + self.ewsParamsName

        self.lane = str(self.getAttributeValue('compute_lane'))

    def doHomogeneousReconstruction(self):
        input_group_connect = {"particles": self.particles.get()}
        if self.mask.get() is not None:
            input_group_connect["mask"] = self.mask.get()

        params = {}

        for paramName in self._paramsName:
            if (paramName != 'refine_N_cmp' and
                    paramName != 'refine_compute_batch_size' and
                    paramName != 'refine_symmetry' and
                    paramName != 'refine_helical_twist' and
                    paramName != 'refine_helical_shift' and
                    paramName != 'refine_hsym_order' and
                    paramName != 'refine_fsc_mask_opt' and
                    paramName != 'refine_ews_zsign' and
                    paramName != 'refine_ews_simple' and
                    paramName != 'refine_override_filter_res'):
                params[str(paramName)] = str(self.getAttributeValue(paramName))

            elif paramName == 'refine_fsc_mask_opt':
                params[str(paramName)] = str("True")
            elif paramName == 'refine_N_cmp' and self.refine_N_cmp.get() is not None:
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif paramName == 'refine_compute_batch_size' and self.refine_compute_batch_size.get() is not None:
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif paramName == 'refine_symmetry':
                symetryValue = getSymmetry(self.symmetryGroup.get(),
                                           self.symmetryOrder.get())
                params[str(paramName)] = symetryValue
            elif paramName == 'refine_helical_twist' and self.refine_helical_twist.get() is not None:
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif paramName == 'refine_helical_shift' and self.refine_helical_shift.get() is not None:
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif paramName == 'refine_hsym_order' and self.refine_hsym_order.get() is not None:
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif paramName == 'refine_ews_zsign':
                params[str(paramName)] = str(EWS_CURVATURE_SIGN[self.refine_ews_zsign.get()])
            elif paramName == 'refine_ews_simple':
                params[str(paramName)] = str(EWS_CORRECTION_METHOD[self.refine_ews_simple.get()])
            elif paramName == 'refine_override_filter_res' and self.refine_override_filter_res.get() is not None:
                params[str(paramName)] = str(self.getAttributeValue(paramName))

        # Determinate the GPUs to use (in dependence of
        # the cryosparc version)
        try:
            gpusToUse = self.getGpuList()
        except Exception:
            gpusToUse = False

        runHomogeneousReconstructionJob = enqueueJob(self._className, self.projectName.get(),
                                                       self.workSpaceName.get(),
                                                       str(params).replace('\'', '"'),
                                                       str(input_group_connect).replace('\'', '"'),
                                                       self.lane, gpusToUse)
        self.runHomogeneousReconstruction = String(runHomogeneousReconstructionJob.get())
        self.currenJob.set(runHomogeneousReconstructionJob.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(), self.runHomogeneousReconstruction.get(),
                         "An error occurred in the homogeneous reconstruction process. "
                         "Please, go to cryoSPARC software for more "
                         "details.", self)
        clearIntermediateResults(self.projectName.get(), self.runHomogeneousReconstruction.get())








