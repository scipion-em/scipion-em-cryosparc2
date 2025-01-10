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
import os.path
import zipfile

from pwem import getMatchingFiles
from pwem.protocols import ProtFlexBase
from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam
from . import ProtCryosparcBase
from ..convert import *
from ..utils import *
from ..constants import *

class ProtCryoSparc3DFlexGenerator(ProtCryosparcBase, ProtFlexBase):
    """
    Takes in a checkpoint from training and generates volume series from it, to show what the model is learning
    about the motion of the particle. This job can be run while training is ongoing to see progress along the way.
    This job can also optionally take in a high-resolution density map (e.g., from 3D Flex Reconstruction) and
    will upscale the deformation model and apply deformations to the high resolution map.
    """
    _label = '3D flex generator'
    _devStatus = BETA
    _protCompatibility = [V4_1_0, V4_1_1, V4_1_2, V4_2_0, V4_2_1, V4_3_1, V4_4_0, V4_4_1, V4_5_1,
                          V4_5_3, V4_6_0, V4_6_1, V4_6_2]

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('flexTraining', PointerParam,
                      pointerClass='ProtCryoSparc3DFlexTraining',
                      label="3D flex training protocol",
                      important=True,
                      help='3DFlex model')

        form.addParam('flex_gen_num_pts', IntParam, default=41,
                      label='Nomber of frames per series',
                      help="Number of volumes to generate per series (unless custom latent components connected). "
                            "This should be an odd number to ensure that the zero position along each latent "
                            "dimension is sampled at the middle of the series")

        # --------------[Compute settings]---------------------------
        form.addSection(label="Compute settings")
        addComputeSectionParams(form, allowMultipleGPUs=False, needGPU=True)

    def _insertAllSteps(self):
        self._defineParamsName()
        self._initializeCryosparcProject()
        self._insertFunctionStep(self.generatorStep)
        self._insertFunctionStep(self.createOutputStep)

    def generatorStep(self):
        self.info(pwutils.yellowStr("3D Flex Generator started..."))
        self.doRun3DFlexGenerator()

    def createOutputStep(self):
        self._initializeUtilsVariables()
        self.info(pwutils.yellowStr("Creating the output..."))

        csOutputFolder = os.path.join(self.projectDir.get(),
                                      self.run3DGeneratorJob.get())

        pattern = csOutputFolder + '/%s_series_*.zip' % self.run3DGeneratorJob.get()
        csSeries = getMatchingFiles(pattern, True)
        fileNameList = [os.path.basename(fileName) for fileName in csSeries]
        copyFiles(csOutputFolder, self._getExtraPath(), files=fileNameList)
        outputFolder = self._getExtraPath()
        for fileName in fileNameList:
            filePath = os.path.join(outputFolder, fileName)
            with zipfile.ZipFile(filePath, 'r') as fileZip:
                extractFolder = os.path.join(outputFolder, os.path.splitext(fileName)[0])
                os.mkdir(extractFolder)
                fileZip.extractall(extractFolder)
            os.remove(filePath)

    def _validate(self):
        validateMsgs = cryosparcValidate()
        if not validateMsgs:
            flexTrainingProt = self.flexTraining.get()
            flexTrainingJob = getJob(flexTrainingProt.projectName.get(), flexTrainingProt.run3DFlexTrainJob.get())
            checkPoints = eval(flexTrainingJob[1])['output_result_groups']
            for output in checkPoints:
                if output['name'] == 'flex_model':
                    numItems = output['num_items']
                    if not numItems:
                        validateMsgs.append("There is no data to show yet")
                    break
        return validateMsgs

    def _defineParamsName(self):
        """ Define a list with 3D Flex Training parameters names"""
        self._paramsName = ['flex_gen_num_pts']
        self.lane = str(self.getAttributeValue('compute_lane'))

    def doRun3DFlexGenerator(self):
        self._className = "flex_generate"

        try:
            gpusToUse = self.getGpuList()
        except Exception:
            gpusToUse = False

        protocolTraining = self.flexTraining.get().run3DFlexTrainJob.get()
        input_group_connect = {"flex_model": "%s.flex_model" % protocolTraining}
        params = {}

        for paramName in self._paramsName:
            if self.getAttributeValue(paramName) is not None:
                params[str(paramName)] = str(self.getAttributeValue(paramName))

        run3DGeneratorJob = enqueueJob(self._className,
                                   self.projectName.get(),
                                   self.workSpaceName.get(),
                                   str(params).replace('\'', '"'),
                                   str(input_group_connect).replace('\'','"'),
                                   self.lane, gpusToUse)

        self.run3DGeneratorJob = String(run3DGeneratorJob.get())
        self.currenJob.set(self.run3DGeneratorJob.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(),
                         self.run3DGeneratorJob.get(),
                         "An error occurred in the 3D Flex Generator process. "
                         "Please, go to cryoSPARC software for more "
                         "details.", self)
        clearIntermediateResults(self.projectName.get(),
                                 self.run3DGeneratorJob.get())
