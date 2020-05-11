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
from pwem.protocols import ProtOperateParticles
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        LEVEL_ADVANCED)

from . import ProtCryosparcBase
from ..convert import *
from ..utils import *
from ..constants import *


class ProtCryoSparcSubtract(ProtCryosparcBase, ProtOperateParticles):
    """ Signal subtraction protocol of cryoSPARC.
        Subtract projections of a masked volume from particles.
        """
    _label = 'subtract projection'

    def _initialize(self):
        self._createFilenameTemplates()

    def _createFilenameTemplates(self):
        """ Centralize how files are called. """
        myDict = {
            'input_particles': self._getTmpPath('input_particles.star'),
            'out_particles': self._getExtraPath('output_particle.star')
        }
        self._updateFilenamesDict(myDict)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      label="Input particles", important=True,
                      help='Select the experimental particles.')
        form.addParam('refVolume', PointerParam, pointerClass='Volume',
                      label="Input map to be projected",
                      important=True,
                      help='Provide the input volume that will be used to '
                           'calculate projections, which will be subtracted '
                           'from the experimental particles. Make sure this '
                           'map was calculated by RELION from the same '
                           'particles as above, and preferably with those '
                           'orientations, as it is crucial that the absolute '
                           'greyscale is the same as in the experimental '
                           'particles.')
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

        form.addParallelSection(threads=1, mpi=1)

        # -----------[Particles Subtraction]------------------------
        form.addSection(label="Particle Subtraction")

        # form.addParam('n_particles', IntParam, default=-1,
        #               label='Number of particles to subtract:',
        #               help='Leave as -1 to process all particles.')

        form.addParam('inner_radius', FloatParam, default=0.85,
                      validators=[Positive],
                      label='Inner radius of reference window',
                      help='Inner radius of the windowing applied to the '
                           'particles used to generate the input structure.')

        form.addParam('outer_radius', FloatParam, default=0.99,
                      validators=[Positive],
                      label='Outer radius of reference window',
                      help='Outer radius of the windowing applied to the '
                           'particles used to generate the input structure')

        form.addParam('use_premult', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Use premultiplier for scaling",
                      help="Premultiplier to scale initial real space mode")

        form.addParam('use_halfmaps', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Use halfmaps for Gold-Standard Subtraction",
                      help="Subtract each halfmap from particles used to "
                           "generate it. Disabling this parameter breaks the "
                           "assumptions of gold-standard FSC calculation.")

        form.addParam('lpf_volume', FloatParam, default=0.0,
                      label='Low-pass Filter Input Structure (A)',
                      help='Apply a lowpass filter to the specified reoslution '
                           'in Angstroms to the input volume before subtraction.'
                           'Leave 0.0 to ignore a lowpass filter')

        form.addParam('mask_threshold', FloatParam, default=0.0,
                      label='Mask threshold',
                      help='The threshold of binarization of the mask. Must be '
                           'set to dilate or pad. Leave 0.0 to skip mask processing.')

        form.addParam('mask_fill_holes', BooleanParam, default=False,
                      label="Fill holes",
                      help="Fill the holes in the binarized mask")

        # form.addParam('mask_dilator', FloatParam, default=0.0,
        #               label='Mask Dilation Radius',
        #               help='The radius of the spherical mask used to dilate '
        #                    'the mask. Leave 0 to skip dilation.')
        #
        # form.addParam('mask_pad', FloatParam, default=0.0,
        #               label='Cosine padding width',
        #               help='The width of the cosine-padded region at the edge'
        #                    'of the mask. Leave 0 to skip padding.')

        # --------------[Compute settings]---------------------------
        form.addSection(label="Compute settings")
        addComputeSectionParams(form, allowMultipleGPUs=False)

    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._defineParamsName()
        self._initializeCryosparcProject()
        self._insertFunctionStep("convertInputStep")
        self._insertFunctionStep('processStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions ------------------------------
    def processStep(self):
        self.vol = self.importVolume.get() + '.imported_volume.map'
        self.mask = self.importMask.get() + '.imported_mask.mask'

        print(pwutils.yellowStr("Particles Subtraction started..."), flush=True)
        self.doPartStract()

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to Relion file
        """
        self._initializeUtilsVariables()
        outputStarFn = self._getFileName('out_particles')

        # Create the output folder
        os.system("cp -r " + self.projectPath + "/" + self.projectName.get() +
                  '/' + self.runPartStract.get() + " " + self._getExtraPath())

        csFileName = "subtracted_particles.cs"
        csFile = os.path.join(self._getExtraPath(), self.runPartStract.get(),
                              csFileName)

        argsList = [csFile, outputStarFn]

        parser = defineArgs()
        args = parser.parse_args(argsList)
        convertCs2Star(args)

        imgSet = self._getInputParticles()

        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(imgSet)
        self._fillDataFromIter(outImgSet)

        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(imgSet, outImgSet)

    def _fillDataFromIter(self, imgSet):

        outImgsFn = self._getFileName('out_particles')
        readSetOfParticles(outImgsFn, imgSet,
                           postprocessImageRow=self._updateItem,
                           alignType=ALIGN_PROJ)

    def _updateItem(self, item, row):
        newFn = row.getValue(md.RLN_IMAGE_NAME)
        index, file = cryosparcToLocation(newFn)
        item.setLocation((index, self._getExtraPath(file)))
        item.setSamplingRate(calculateNewSamplingRate(item.getDim(),
                                                      self._getInputParticles().getSamplingRate(),
                                                      self._getInputParticles().getDim()))

    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        """ Should be overwritten in subclasses to
            return summary message for NORMAL EXECUTION.
        """
        validateMsgs = cryosparcValidate()
        if not validateMsgs:
            validateMsgs = gpusValidate(self.getGpuList(), checkSingleGPU=True)
            if not validateMsgs:
                self._validateDim(self._getInputParticles(),
                                  self.refVolume.get(),
                                  validateMsgs, 'Input particles',
                                  'Input volume')
        return validateMsgs

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputParticles'):
            summary.append("Output Particles not ready yet.")
        else:
            summary.append("Input Particles: %s" %
                           self.getObjectTag('inputParticles'))
            summary.append("Reference Volume: %s" %
                           self.getObjectTag('refVolume'))
            summary.append("Reference Mask: %s" %
                           self.getObjectTag('refMask'))

            summary.append("Inner radius of the window: %s" %
                           str(self.inner_radius.get()))

            summary.append("Outer radius of the window: %s" %
                           str(self.outer_radius.get()))

            summary.append("--------------------------------------------------")
            summary.append("Output particles %s" %
                           self.getObjectTag('outputParticles'))
        return summary

    # ---------------Utils Functions-------------------------------------------

    def _defineParamsName(self):
        """ Define a list with all protocol parameters names"""
        self._paramsName = ['inner_radius',
                            'outer_radius',
                            'use_premult',
                            'use_halfmaps',
                            'lpf_volume',
                            'mask_threshold',
                            'mask_fill_holes',
                            'compute_use_ssd']
        self.lane = str(self.getAttributeValue('compute_lane'))

    def doPartStract(self):
        """
        :return:
        """
        className = "particle_subtract"
        input_group_conect = {"particles": str(self.par),
                              "volume": str(self.vol),
                              "mask": str(self.mask)}
        # {'particles' : 'JXX.imported_particles' }
        params = {}

        for paramName in self._paramsName:
            if paramName != 'lpf_volume' and paramName != 'mask_threshold':
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif paramName == 'lpf_volume':
                if float(self.getAttributeValue(paramName)) > 0:
                    params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif paramName == 'mask_threshold':
                if int(self.getAttributeValue(paramName)) > 0:
                    params[str(paramName)] = str(self.getAttributeValue(paramName))

        # Determinate the GPUs to use (in dependence of
        # the cryosparc version)
        try:
            gpusToUse = self.getGpuList()
        except Exception:
            gpusToUse = False

        self.runPartStract = enqueueJob(className, self.projectName.get(),
                                  self.workSpaceName.get(),
                                  str(params).replace('\'', '"'),
                                  str(input_group_conect).replace('\'', '"'),
                                  self.lane, gpusToUse)

        self.currenJob.set(self.runPartStract.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(), self.runPartStract.get(),
                         "An error occurred in the particles subtraction process. "
                         "Please, go to cryosPARC software for more "
                         "details.")
