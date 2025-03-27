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
from pwem.objects import Volume
from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam
from . import ProtCryosparcBase
from ..convert import *
from ..utils import *
from ..constants import *


class ProtCryoSparc3DFlexDataPrepare(ProtCryosparcBase):
    """
    Prepares particles for use in 3DFlex training and reconstruction. At the same
    way, takes in a consensus (rigid) refinement density map, plus optionally
    a segmentation and generates a tetrahedral mesh for 3DFlex.
    """
    _label = '3D flex data prepare'
    _devStatus = BETA
    _protCompatibility = [V4_1_0, V4_1_1, V4_1_2, V4_2_0, V4_2_1, V4_3_1,
                          V4_4_0, V4_4_1, V4_5_1, V4_5_3, V4_6_0, V4_6_1, V4_6_2]

    # --------------------------- DEFINE param functions ----------------------
    def _defineFileNames(self):
        """ Centralize how files are called within the protocol. """
        myDict = {
            'input_particles': self._getTmpPath('input_particles.star'),
            'out_particles': self._getPath() + '/output_particle.star',
            'stream_log': self._getPath() + '/stream.log'
        }
        self._updateFilenamesDict(myDict)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      label="Input particles",
                      important=True,
                      help='Particle stacks to use.')
        form.addParam('refVolume', PointerParam, pointerClass='Volume',
                      label='Initial volume',
                      important=True,
                      help='Initial volume raw data')

        form.addSection(label='Data Prepare')
        form.addParam('box_size_pix', IntParam, default=None,
                      allowsNull=True,
                      label="Crop box size (pix)",
                      help="Crop input particles and volume to this box size. "
                           "This is the box size that will be used for high "
                           "resolution reconstruction with 3D Flex. Particles "
                           "are cropped to this box size first, and then "
                           "downsampled to the Training Box Size for training "
                           "time. Default (None) means to keep the original "
                           "box size of the particles.")

        form.addParam('bin_size_pix', IntParam, default=128,
                      allowsNull=True,
                      label="Training box size (pix)",
                      help="Downsample cropped particles (via Fourier cropping)"
                           " to this box size for training the 3D Flex model. "
                           "This should be chosen to limit 3D Flex training to "
                           "a resolution below the gold-standard FSC resolution "
                           "of the consensus reconstruction, in order to ensure "
                           "high resolution reconstructions can be validated. "
                           "Box sizes over 256 may become prohibitively slow")

        form.addParam('alpha_min', IntParam, default=None,
                      allowsNull=True,
                      label="Min. scale to keep",
                      help="Only keep particles with scale factor above this "
                           "value. Useful for discarding particles that might "
                           "be junk.")

        form.addParam('keep_num_particles', IntParam, default=None,
                      allowsNull=True,
                      label="Num. particles to use",
                      help="Only keep the first X particles. The final number"
                           " of particles used during 3D Flex training "
                           "and reconstruction must be divisible by 1000. "
                           "If this is None (default) then the number of input "
                           "particles will be rounded down to the nearest 1000.")

        # --------------[Compute settings]---------------------------
        form.addSection(label="Compute settings")
        addComputeSectionParams(form, allowMultipleGPUs=False, needGPU=False)

    def _insertAllSteps(self):
        self._defineFileNames()
        self._defineParamsPrepareName()
        self._initializeCryosparcProject()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.dataPrepareStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions ------------------------------
    def dataPrepareStep(self):
        self.info(pwutils.yellowStr("3D Flex Data Preparation started..."))
        self.doRun3DFlexDataPrepare()

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to Relion file
        """
        self.info(pwutils.yellowStr("Creating the output..."))
        self._initializeUtilsVariables()
        outputStarFn = self._getFileName('out_particles')
        self.info(pwutils.yellowStr("outputStarFn: %s " % outputStarFn))
        csOutputFolder = os.path.join(self.projectDir.get(), self.run3DFlexDataPrepJob.get())
        self.info("csOutputFolder: %s " % csOutputFolder)
        # csFileName = "subtracted_particles.cs"
        csFileName = "%s_passthrough_particles.cs" % self.run3DFlexDataPrepJob.get()
        self.info("csFileName: %s " % csOutputFolder)
        # Create the output folder
        copyFiles(csOutputFolder,  os.path.join(self._getExtraPath(), self.run3DFlexDataPrepJob.get()))
        self.info("copyFolder: src-> dst %s %s" % (csOutputFolder, os.path.join(self._getExtraPath(), self.run3DFlexDataPrepJob.get())))
        csFile = os.path.join(self._getExtraPath(), self.run3DFlexDataPrepJob.get(), csFileName)
        self.info("csFile (metadata): %s " % csFile)
        argsList = [csFile, outputStarFn]
        self.info("starFile: %s " % outputStarFn)
        convertCs2Star(argsList)
        self.info("Convert to star done.")
        self.info("Creating the particles output")
        imgSet = self._getInputParticles()
        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(imgSet)
        outImgSet.setSamplingRate(imgSet.getSamplingRate())
        imgSetDict = {}
        for img in imgSet.iterItems():
            imgSetDict[str(img.getIndex()) + '@' + os.path.basename(img.getFileName())] = img.clone()
        self._fillDataFromIter(outImgSet, imgSetDict)
        self.info("Creating the consensus volume  output")

        csMapName = "%s_map.mrc" % self.run3DFlexDataPrepJob.get()
        fnVol = os.path.join(self._getExtraPath(), self.run3DFlexDataPrepJob.get(), csMapName)
        vol = Volume()
        fixVolume(fnVol)
        vol.setFileName(fnVol)
        vol.setSamplingRate(calculateNewSamplingRate(vol.getDim(),
                                                     imgSet.getSamplingRate(),
                                                     imgSet.getDim()))

        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(imgSet, outImgSet)
        self._defineOutputs(outputVolume=vol)

    # ------------------------- Utils methods ----------------------------------

    def _fillDataFromIter(self, outImgSet, imgSetDict):
        filename = 'particles@' + self._getFileName('out_particles')
        for imgRow in emtable.Table.iterRows(filename):
            fileName = imgRow.get(RELIONCOLUMNS.rlnImageName.value)
            splitFileName = fileName.split('@')
            key = str(int(splitFileName[0])) + '@' + '_'.join(splitFileName[1].split('_')[1:])
            if key in imgSetDict:
                img = imgSetDict[key]
                outImgSet.append(img)

    def _createItemMatrix(self, particle, row):
        createItemMatrix(particle, row, align=ALIGN_PROJ)
        setCryosparcAttributes(particle, row, RELIONCOLUMNS.rlnRandomSubset.value)

    def _validate(self):
        validateMsgs = cryosparcValidate()
        if not validateMsgs:
            particles = self._getInputParticles()
            if not particles.hasCTF():
                validateMsgs.append(
                    "The Particles has not associated a CTF model")
                if not particles.hasAlignment3D():
                    validateMsgs.append("The Particles has not aligned")
                    if self.keep_num_particles.get() % 1000 != 0:
                        validateMsgs.append("The final number of particles "
                                            "used during 3D Flex training and "
                                            "reconstruction must be divisible by 1000")

        return validateMsgs

    def _defineParamsPrepareName(self):
        """ Define a list with 3D Flex Prepare Data parameters names"""
        self._paramsPrepareName = ['box_size_pix', 'bin_size_pix', 'alpha_min',
                                   'keep_num_particles']
        self.lane = str(self.getAttributeValue('compute_lane'))

    def doRun3DFlexDataPrepare(self):
        self._className = "flex_prep"
        input_group_connect = {"particles": self.particles.get(),
                               "volume": self.volume.get()}
        params = {}

        for paramName in self._paramsPrepareName:
            if self.getAttributeValue(paramName) is not None:
                params[str(paramName)] = str(self.getAttributeValue(paramName))

        run3DFlexDataPrepJob = enqueueJob(self._className, self.projectName.get(),
                                  self.workSpaceName.get(),
                                  str(params).replace('\'', '"'),
                                  str(input_group_connect).replace('\'', '"'),
                                  self.lane, False)

        self.run3DFlexDataPrepJob = String(run3DFlexDataPrepJob.get())
        self.currenJob.set(self.run3DFlexDataPrepJob.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(), self.run3DFlexDataPrepJob.get(),
                         "An error occurred in the 3D Flex Data Preparation process. "
                         "Please, go to cryoSPARC software for more "
                         "details.", self)
        clearIntermediateResults(self.projectName.get(), self.run3DFlexDataPrepJob.get())








