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

import os

from pwem.objects import SetOfParticles

import pyworkflow.utils as pwutils
from pyworkflow import NEW
from pyworkflow.object import String
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        IntParam)

from .protocol_base import ProtCryosparcBase
from ..convert import (defineArgs, convertCs2Star, readSetOfParticles)
from ..utils import (addComputeSectionParams, cryosparcValidate, gpusValidate,
                     enqueueJob, waitForCryosparc, clearIntermediateResults,
                     addSymmetryParam, getSymmetry, copyFiles)


class ProtCryoSparcSymmetryExpansion(ProtCryosparcBase):
    """ Duplicate particles around a point-group symmetry.
    """
    _label = 'symmetry expansion'
    _className = "sym_expand"
    _devStatus = NEW

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

        addSymmetryParam(form, help="Symmetry String (C, D, I, O, T). E.g. C1,"
                                    "D7, C4, etc. Particles will be "
                                    "symmetry-expanded at this symmetry.")

        form.addParam('sym_twist_deg', FloatParam, default=None,
                      allowsNull=True,
                      label='Helical twist (degrees)',
                      help='Helical twist for symmetry expansion. This can be '
                           'found in the final iteration of the source Helical '
                           'Refinement job streamlog.')

        form.addParam('sym_rise_A', FloatParam, default=None,
                      allowsNull=True,
                      label='Helical rise (A)',
                      help='Helical rise for symmetry expansion. This can be '
                           'found in the final iteration of the source '
                           'Helical Refinement job streamlog.')

        form.addParam('sym_num_rises', IntParam, default=None,
                      allowsNull=True,
                      label='Helical symmetry order (integer)',
                      help='Helical symmetry order for symmetry expansion. '
                           'This can be found in the final iteration of the '
                           'source Helical Refinement job streamlog.')

        # --------------[Compute settings]---------------------------
        form.addSection(label="Compute settings")
        addComputeSectionParams(form, allowMultipleGPUs=False)

    # --------------------------- INSERT steps functions ----------------------

    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._defineParamsName()
        self._initializeCryosparcProject()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.processStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions ------------------------------

    def processStep(self):
        print(pwutils.yellowStr("Symmetry Expansion started..."), flush=True)
        self.doSymmetryExpansion()

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to Relion file
        """
        self._initializeUtilsVariables()
        outputStarFn = self._getFileName('out_particles')
        csFileName = "particles_expanded.cs"

        # Copy the CS output expanded particles to extra folder
        copyFiles(os.path.join(self.projectPath, self.projectName.get(),
                               self.runSymExp.get()),
                  self._getExtraPath(), files=[csFileName])

        csFile = os.path.join(self._getExtraPath(), csFileName)

        argsList = [csFile, outputStarFn]

        parser = defineArgs()
        args = parser.parse_args(argsList)
        convertCs2Star(args)
        imgSet = self._getInputParticles()
        self.setFilePattern(imgSet.getFirstItem().getFileName())
        outImgSet = SetOfParticles.create(self._getExtraPath())
        outImgSet.copyInfo(imgSet)
        readSetOfParticles('particles@' + outputStarFn, outImgSet,
                           alignType=imgSet.getAlignment(),
                           postprocessImageRow=self.updateParticlePath)

        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(imgSet, outImgSet)

    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        """ Should be overwritten in subclasses to
            return summary message for NORMAL EXECUTION.
        """
        validateMsgs = cryosparcValidate()
        if not validateMsgs:
            validateMsgs = gpusValidate(self.getGpuList(),
                                        checkSingleGPU=True)
        return validateMsgs

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputParticles'):
            summary.append("Output Particles not ready yet.")
        else:
            summary.append("Input Particles: %s" %
                           self.getObjectTag('inputParticles'))
            summary.append(
                "--------------------------------------------------")
            summary.append("Output particles %s" %
                           self.getObjectTag('outputParticles'))
        return summary

    # ---------------Utils Functions-------------------------------------------

    def _defineParamsName(self):
        """ Define a list with all protocol parameters names"""
        self._paramsName = ['sym_symmetry',
                            'sym_twist_deg',
                            'sym_rise_A',
                            'sym_num_rises',
                            'compute_use_ssd']
        self.lane = str(self.getAttributeValue('compute_lane'))

    def doSymmetryExpansion(self):
        """
        Launch a symmetry expansion job
        """
        input_group_connect = {"particles": self.particles.get()}
        params = {}

        for paramName in self._paramsName:
            if paramName == 'sym_symmetry':
                symetryValue = getSymmetry(self.symmetryGroup.get(),
                                           self.symmetryOrder.get())
                params[str(paramName)] = symetryValue
            elif paramName == 'sym_num_rises':
                if self.getAttributeValue(paramName) is not None and int(self.getAttributeValue(paramName)) > 0:
                    params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif self.getAttributeValue(paramName) is not None and float(self.getAttributeValue(paramName)) > 0:
                params[str(paramName)] = str(self.getAttributeValue(paramName))

        # Determinate the GPUs to use (in dependence of
        # the cryosparc version)
        try:
            gpusToUse = self.getGpuList()
        except Exception:
            gpusToUse = False

        runSymExpJob = enqueueJob(self._className, self.projectName.get(),
                                    self.workSpaceName.get(),
                                    str(params).replace('\'', '"'),
                                    str(input_group_connect).replace('\'', '"'),
                                    self.lane, gpusToUse)

        self.runSymExp = String(runSymExpJob.get())
        self.currenJob.set(self.runSymExp.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(), self.runSymExp.get(),
                         "An error occurred in the particles subtraction process. "
                         "Please, go to cryosPARC software for more "
                         "details.")
        clearIntermediateResults(self.projectName.get(), self.runSymExp.get())

