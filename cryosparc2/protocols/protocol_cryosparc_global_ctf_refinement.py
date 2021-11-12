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

from pwem import ALIGN_PROJ
import pwem.protocols as pwprot

import pyworkflow.utils as pwutils
from pyworkflow import BETA
from pyworkflow.object import String
from pyworkflow.protocol.params import (PointerParam, FloatParam, IntParam,
                                        LEVEL_ADVANCED, Positive, BooleanParam)

from .protocol_base import ProtCryosparcBase
from ..convert import (defineArgs, convertCs2Star, createItemMatrix,
                       setCryosparcAttributes)
from ..utils import (addComputeSectionParams, cryosparcValidate, gpusValidate,
                     enqueueJob, waitForCryosparc, clearIntermediateResults,
                     copyFiles)
import pwem.emlib.metadata as md


class ProtCryoSparcGlobalCtfRefinement(ProtCryosparcBase, pwprot.ProtParticles):
    """
    Wrapper protocol for the Cryosparc's per-particle Global CTF refinement.
    Performs per-exposure-group CTF parameter refinement of higher-order
    aberrations, against a given 3D reference
    """
    _label = 'global ctf refinement'
    _devStatus = BETA
    _className = "ctf_refine_global"

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
                      help='Provide a set of particles for global '
                           'CTF refinement.')
        form.addParam('refVolume', PointerParam, pointerClass='Volume',
                      important=True,
                      label="Input volume",
                      help='Provide a reference volume for global '
                           'CTF refinement.')
        form.addParam('refMask', PointerParam, pointerClass='VolumeMask',
                      label='Mask to be applied to this map',
                      important=True,
                      help="Provide a soft mask. if mask is present, use that, "
                           "otherwise use mask_refine if present, otherwise "
                           "fail")

        # -----------[Global CTF Refinement]------------------------
        form.addSection(label="Global CTF Refinement")
        form.addParam('crg_num_iters', IntParam, default=1,
                      validators=[Positive],
                      label='Number of iterations',
                      help='Number of times refinement of the various CTF '
                           'parameters is done. Using 2 or more iterations '
                           'allows changes in one parameter (eg. anisomag) to '
                           'affect the estimation of other parameters '
                           '(eg. tetrafoil).')
        form.addParam('crg_num_plots', IntParam, default=50,
                      validators=[Positive],
                      label='Num. groups to plot',
                      help='Number of exposure groups to make plots for. '
                           'After this many, stop plotting to save time.')

        form.addParam('crg_plot_binfactor', IntParam, default=1,
                      validators=[Positive],
                      expertLevel=LEVEL_ADVANCED,
                      label='Binning to apply to plots',
                      help='Binning makes it easier to see tilt/trefoil/tetrafoil '
                           'etc in the data plots, but does not change the '
                           'results')

        form.addParam('crg_min_res_A', FloatParam, default=10,
                      validators=[Positive],
                      label='Minimum Fit Res (A)',
                      help='The minimum resolution to use during refinement of '
                           'image aberrations.')

        form.addParam('crg_do_tilt', BooleanParam, default=True,
                      label="Fit Tilt",
                      help="Whether to fit beam tilt.")

        form.addParam('crg_do_trefoil', BooleanParam, default=True,
                      label="Fit Trefoil",
                      help="Whether to fit beam trefoil.")

        form.addParam('crg_do_spherical', BooleanParam, default=True,
                      label="Fit Spherical Aberration",
                      help="Whether to fit spherical aberration.")

        form.addParam('crg_do_tetrafoil', BooleanParam, default=True,
                      label="Fit Tetrafoil",
                      help="Whether to fit beam tetrafoil.")

        # --------------[Compute settings]---------------------------
        form.addSection(label="Compute settings")
        addComputeSectionParams(form, allowMultipleGPUs=False)

    # --------------------------- INSERT steps functions -----------------------

    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._defineParamsName()
        self._initializeCryosparcProject()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.processStep)
        self._insertFunctionStep(self.createOutputStep)

    # -------------------------- UTILS functions ------------------------------

    def _defineParamsName(self):
        """ Define a list with all protocol parameters names"""
        self._paramsName = ['crg_num_iters',
                            'crg_num_plots',
                            'crg_plot_binfactor',
                            'crg_min_res_A',
                            'crg_do_tilt',
                            'crg_do_trefoil',
                            'crg_do_spherical',
                            'crg_do_tetrafoil',
                            'compute_use_ssd']
        self.lane = str(self.getAttributeValue('compute_lane'))

    # --------------------------- STEPS functions ------------------------------
    def processStep(self):
        print(pwutils.yellowStr("Ctf Refinement started..."), flush=True)
        self.doGlobalCtfRefinement()

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to Relion file
        """
        self._initializeUtilsVariables()
        outputStarFn = self._getFileName('out_particles')
        csOutputFolder = os.path.join(self.projectPath, self.projectName.get(),
                                      self.runGlobalCtfRefinement.get())
        csFileName = "particles.cs"

        # Copy the CS output particles to extra folder
        copyFiles(csOutputFolder, self._getExtraPath(), files=[csFileName])

        csFile = os.path.join(self._getExtraPath(), csFileName)

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
                                  self._getInputVolume(),
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
            summary.append("Reference Mask: %s" %
                           self.getObjectTag('refMask'))
            summary.append("Number of Iterations: %s" %
                           self.crg_num_iters.get())
            summary.append("--------------------------------------------------")
            summary.append("Output particles %s" %
                           self.getObjectTag('outputParticles'))
        return summary

    def doGlobalCtfRefinement(self):
        """
        :return:
        """
        input_group_connect = {"particles": self.particles.get(),
                               "volume": self.volume.get(),
                               "mask": self.mask.get()}
        input_result_connect = None
        if self._getInputVolume().hasHalfMaps():
            input_result_connect = {"volume.0.map_half_A": self.importVolumeHalfA.get(),
                                    "volume.0.map_half_B": self.importVolumeHalfB.get()}
        params = {}

        for paramName in self._paramsName:
            params[str(paramName)] = str(self.getAttributeValue(paramName))

        # Determinate the GPUs to use (in dependence of
        # the cryosparc version)
        try:
            gpusToUse = self.getGpuList()
        except Exception:
            gpusToUse = False

        runGlobalCtfRefinementJob = enqueueJob(self._className, self.projectName.get(),
                                                 self.workSpaceName.get(),
                                                 str(params).replace('\'', '"'),
                                                 str(input_group_connect).replace('\'', '"'),
                                                 self.lane, gpusToUse,
                                                 result_connect=input_result_connect)

        self.runGlobalCtfRefinement = String(runGlobalCtfRefinementJob)
        self.currenJob.set(runGlobalCtfRefinementJob.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(), self.runGlobalCtfRefinement.get(),
                         "An error occurred in the particles subtraction process. "
                         "Please, go to cryosPARC software for more "
                         "details.")