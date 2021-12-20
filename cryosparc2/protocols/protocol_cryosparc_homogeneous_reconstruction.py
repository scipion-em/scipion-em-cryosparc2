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

from pkg_resources import parse_version

from pwem import ALIGN_PROJ

import pyworkflow.utils as pwutils
from pyworkflow.object import String
from pyworkflow.protocol.params import (PointerParam, LEVEL_ADVANCED, IntParam,
                                        BooleanParam)
from pwem.objects import Volume

from .protocol_base import ProtCryosparcBase
from ..convert import (defineArgs, convertCs2Star, createItemMatrix,
                       setCryosparcAttributes)
from ..utils import (addComputeSectionParams, calculateNewSamplingRate,
                     cryosparcValidate, gpusValidate, enqueueJob,
                     waitForCryosparc, clearIntermediateResults, fixVolume,
                     copyFiles, addSymmetryParam, getSymmetry,
                     getCryosparcVersion, get_job_streamlog)
from ..constants import *


class ProtCryoSparcHomogeneousReconstruct(ProtCryosparcBase):
    """ Create a 3D reconstruction from input particles that already have alignments in 3D.
    """
    _label = 'homogeneous reconstruction'
    _className = "homo_reconstruct"
    _fscColumns = 6
    _protCompatibility = [V3_0_0, V3_1_0, V3_2_0, V3_3_0, V3_3_1]

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
                      label='Mask to be applied to this map',
                      important=True,
                      allowsNull=False,
                      help="Provide a soft mask where the protein density "
                           "you wish to subtract from the experimental "
                           "particles is white (1) and the rest of the "
                           "protein and the solvent is black (0). "
                           "That is: *the mask should INCLUDE the part of the "
                           "volume that you wish to SUBTRACT.*")

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
        print(pwutils.yellowStr("Homogeneous Reconstruction started..."), flush=True)
        self.doHomogeneousReconstruction()

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to Relion file
        """
        self._initializeUtilsVariables()
        idd = self.findLastIteration(self.runHomogeneousReconstruction.get())

        csOutputFolder = os.path.join(self.projectPath, self.projectName.get(),
                                       self.runHomogeneousReconstruction.get())

        csOutputPattern = "cryosparc_%s_%s" % (self.projectName.get(),
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
        parser = defineArgs()
        args = parser.parse_args(argsList)
        convertCs2Star(args)

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
        self._defineSourceRelation(self.inputParticles.get(), vol)
        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(self.inputParticles.get(), outImgSet)
        self.createFSC(idd, imgSet, vol)

    def _fillDataFromIter(self, imgSet):
        outImgsFn = 'particles@' + self._getFileName('out_particles')
        imgSet.setAlignmentProj()
        imgSet.copyItems(self._getInputParticles(),
                         updateItemCallback=self._createItemMatrix,
                         itemDataIterator=md.iterRows(outImgsFn,
                                                      sortByLabel=md.RLN_IMAGE_ID))

    def _createItemMatrix(self, particle, row):
        createItemMatrix(particle, row, align=ALIGN_PROJ)
        setCryosparcAttributes(particle, row,
                               md.RLN_PARTICLE_RANDOM_SUBSET)

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
            csVersion = getCryosparcVersion()
            if [version for version in self._protCompatibility
                if parse_version(version) >= parse_version(csVersion)]:
                validateMsgs = gpusValidate(self.getGpuList(),
                                            checkSingleGPU=True)
                if not validateMsgs:
                    particles = self._getInputParticles()
                    if not particles.hasCTF():
                        validateMsgs.append("The Particles has not associated a "
                                            "CTF model")
                        if not validateMsgs and not particles.hasAlignment3D():
                            validateMsgs.append("The Particles has not a 3D "
                                                "alignment")
            else:
                validateMsgs.append("The protocol is not compatible with the "
                                    "cryoSPARC version %s" % csVersion)

        return validateMsgs

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
                            'compute_use_ssd']
        self.lane = str(self.getAttributeValue('compute_lane'))

    def doHomogeneousReconstruction(self):
        input_group_connect = {"particles": self.particles.get(),
                               "mask": self.mask.get()}

        params = {}

        for paramName in self._paramsName:
            if (paramName != 'refine_N_cmp' and
                    paramName != 'refine_compute_batch_size' and
                    paramName != 'refine_symmetry' and
                    paramName != 'refine_helical_twist' and
                    paramName != 'refine_helical_shift' and
                    paramName != 'refine_hsym_order' and
                    paramName != 'refine_fsc_mask_opt'):
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
                         "Please, go to cryosPARC software for more "
                         "details.")
        clearIntermediateResults(self.projectName.get(), self.runHomogeneousReconstruction.get())








