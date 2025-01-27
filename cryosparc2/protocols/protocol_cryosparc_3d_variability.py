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
from pyworkflow import BETA
from pyworkflow.protocol.params import (PointerParam, FloatParam, Positive,
                                        BooleanParam, EnumParam)
from pwem.objects import Volume
from pwem.protocols import ProtRefine3D

from . import ProtCryosparcBase
from ..convert import *
from ..utils import *
from ..constants import *


class ProtCryoSparc3DVariability(ProtCryosparcBase, ProtRefine3D):
    """
    Protocol to compute the principle modes of variability with a dataset
    of aligned particles
    """
    _label = '3D variability Analysis '
    _devStatus = BETA

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
        form.addParam('refMask', PointerParam, pointerClass='VolumeMask',
                      default=None,
                      label='Input Mask',
                      allowsNull=False,
                      help='Mask raw data')

        form.addParallelSection(threads=1, mpi=1)

        # --------------[3D Variability]---------------------------

        form.addSection(label='3D Variability')
        form.addParam('var_K', IntParam, default=3,
                      label="Number of modes to solve",
                      help='The number of orthogonal principle modes (i.e. '
                           'eigenvectors of the 3D covariance) to solve.')

        form.addParam('var_num_iterations', IntParam, default=20,
                      validators=[Positive],
                      label="Number of iterations",
                      help='Number of iterations')

        # form.addParam('var_N', FloatParam, default=-1,
        #               expertLevel=LEVEL_ADVANCED,
        #               label="Refinement box size (Voxels)",
        #               help='The volume size to use for refinement. If this is '
        #                    '-1, use the full image size. Otherwise images '
        #                    'are automatically downsampled')

        addSymmetryParam(form)

        # form.addParam('var_num_particles', IntParam, default=None,
        #               label="Only use this many particles")

        form.addParam('var_filter_res', FloatParam, default=None,
                      validator=[Positive],
                      label="Filter resolution (A)",
                      help='Resolution at which results are filtered')

        form.addParam('var_filter_order', FloatParam, default=1.5,
                      validators=[Positive],
                      label="Filter order",
                      help='Order of filter')

        # form.addParam('var_highpass_res', StringParam, default='',
        #               expertLevel=LEVEL_ADVANCED,
        #               label="Highpass resolution (A)",
        #               help='Resolution below which variability is ignored')

        form.addParam('var_highpass_order', FloatParam, default=8,
                      validators=[Positive],
                      label="Highpass order",
                      help='Order of filter')

        form.addParam('var_use_gramschmidt', BooleanParam, default=True,
                      label="Use Gram-Schmidt",
                      help='Order of filter')

        form.addParam('var_use_white_noise', BooleanParam, default=False,
                      label="Use white noise model",
                      help='Use a white noise model (default) instead of a '
                           'colored noise model. One or the other may work '
                           'better depending on the dataset.')

        form.addParam('var_use_scales', EnumParam,
                      choices=['none', 'input', 'optimal'],
                      default=2,
                      label="Per-particle scale",
                      help='How to treat per-particle scale factors. No scales '
                           'means all particles have scale 1.0 (useful for very '
                           'small particles or strange cases). Input scales '
                           'means to use the input scale factors (which may '
                           'have come from another refinement). Optimal means '
                           'to compute per-particle optimal scales on the fly '
                           'during 3D variability.')

        form.addParam('var_lambda', FloatParam, default=0.01,
                      validators=[Positive],
                      label="Lambda",
                      help='Stabilizing coefficient - try larger values if '
                           'optimization diverges and creates artefacts. '
                           'In v2.12 this was changed to a normalized '
                           'fractional value, where 0.01 should work for '
                           'almost all datasets.')
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

    # --------------------------- STEPS functions ------------------------------
    def processStep(self):
        self.info(pwutils.yellowStr("3D Variability started..."))
        self.doRun3DVariability()

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to Relion file
        """
        self._initializeUtilsVariables()
        csOutputFolder = os.path.join(self.projectDir.get(),
                                      self.run3DVariability.get())
        csParticlesName = (getOutputPreffix(self.projectName.get()) +
                           self.run3DVariability.get() + "_particles.cs")

        csFile = os.path.join(csOutputFolder, csParticlesName)

        # Copy the particles to scipion output folder
        os.system("cp -r " + csFile + " " + self._getExtraPath())
        csFile = os.path.join(self._getExtraPath(), csParticlesName)

        outputStarFn = self._getFileName('out_particles')
        argsList = [csFile, outputStarFn]

        convertCs2Star(argsList)

        fnVolName = (getOutputPreffix(self.projectName.get()) +
                     self.run3DVariability.get() + "_map.mrc")

        fnVol = os.path.join(csOutputFolder, fnVolName)

        # Copy the volumes to extra folder
        os.system("cp -r " + fnVol + " " + self._getExtraPath())
        fnVol = os.path.join(self._getExtraPath(), fnVolName)

        imgSet = self._getInputParticles()
        vol = Volume()
        vol.setFileName(fnVol)
        vol.setSamplingRate(calculateNewSamplingRate(vol.getDim(),
                                                     imgSet.getSamplingRate(),
                                                     imgSet.getDim()))
        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(imgSet)
        self._fillDataFromIter(outImgSet)

        self._defineOutputs(outputVolume=vol)
        self._defineSourceRelation(self.inputParticles.get(), vol)
        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(self.inputParticles.get(), outImgSet)

    #  ----------------------------UTILS functions --------------------------

    def cleanTmp(self):
        """ Delete all files and subdirectories under Tmp folder. For this
        protocol we need to keep the tmp folder content"""
        pass

    def _validate(self):
        validateMsgs = cryosparcValidate()
        if not validateMsgs:
            validateMsgs = gpusValidate(self.getGpuList(), checkSingleGPU=True)
            if not validateMsgs:
                particles = self._getInputParticles()
                if not particles.hasCTF():
                    validateMsgs.append("The Particles has not associated a "
                                        "CTF model")
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
                           self.getObjectTag('refVolume'))
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
        return summary

    # -------------------------- UTILS functions ------------------------------

    def _fillDataFromIter(self, imgSet):
        outImgsFn = 'particles@' + self._getFileName('out_particles')
        imgSet.setAlignmentProj()
        imgSet.copyItems(self._getInputParticles(),
                         updateItemCallback=self._createItemMatrix,
                         itemDataIterator=emtable.Table.iterRows(
                             fileName=outImgsFn))

    def _createItemMatrix(self, particle, row):
        createItemMatrix(particle, row, align=ALIGN_PROJ)
        setCryosparcAttributes(particle, row,
                               RELIONCOLUMNS.rlnRandomSubset.value)

    def _defineParamsName(self):
        """ Define a list with all protocol parameters names"""
        self._paramsName = ['var_K',
                            'var_filter_res',
                            'var_filter_order',
                            'var_highpass_order',
                            'var_use_gramschmidt',
                            'var_use_white_noise',
                            'var_use_scales',
                            'var_num_iterations',
                            'var_lambda',
                            'compute_use_ssd',
                            'var_symmetry']
        self.lane = str(self.getAttributeValue('compute_lane'))

    def doRun3DVariability(self):
        """
        :return:
        """
        className = "var_3D"
        input_group_conect = {"particles": str(self.particles),
                              "mask": str(self.mask)}
        params = {}

        for paramName in self._paramsName:
            if (paramName != 'var_symmetry' and
                    paramName != 'var_use_scales'):
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif paramName == 'var_symmetry':
                symetryValue = getSymmetry(self.symmetryGroup.get(),
                                           self.symmetryOrder.get())
                params[str(paramName)] = symetryValue
            elif paramName == 'var_use_scales':
                params[str(paramName)] = str(VAR_USE_SCALES[self.var_use_scales.get()])

        # Determinate the GPUs to use (in dependence of
        # the cryosparc version)
        try:
            gpusToUse = self.getGpuList()
        except Exception:
            gpusToUse = False

        self.run3DVariability = enqueueJob(className, self.projectName.get(),
                                           self.workSpaceName.get(),
                                           str(params).replace('\'', '"'),
                                           str(input_group_conect).replace('\'',
                                                                           '"'),
                                           self.lane, gpusToUse)

        self.currenJob.set(self.run3DVariability.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(), self.run3DVariability.get(),
                         "An error occurred in the 3D Variability process. "
                         "Please, go to cryosPARC software for more "
                         "details.", self)



