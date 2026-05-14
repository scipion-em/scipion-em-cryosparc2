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


import pyworkflow.utils as pwutils
from pyworkflow.object import String
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        IntParam)

from .protocol_base import ProtCryosparcBase
from ..convert import (convertCs2Star, readSetOfParticles)
from ..utils import (addComputeSectionParams, cryosparcValidate, gpusValidate,
                     enqueueJob, waitForCryosparc, clearIntermediateResults,
                     addSymmetryParam, getSymmetry, copyFiles)


class ProtCryoSparcSymmetryExpansion(ProtCryosparcBase):
    """ Duplicate particles around a point-group symmetry.
    """

    """
        CryoSPARC Symmetry Expansion Protocol — User Manual

        Overview

        The ProtCryoSparcSymmetryExpansion protocol is a Scipion wrapper
        designed to execute cryoSPARC symmetry expansion jobs within an
        integrated cryo-electron microscopy workflow. Its primary purpose is
        to duplicate particle orientations according to a specified symmetry
        group so that each particle can be represented multiple times under
        all symmetry-related orientations. In cryo-EM structural analysis,
        symmetry expansion is particularly important for studying local
        conformational variability, asymmetric structural features, and
        flexible regions inside otherwise symmetric macromolecular complexes.

        The protocol inherits from ProtCryosparcBase and acts as a bridge
        between Scipion and cryoSPARC. Internally, it prepares the particle
        metadata, defines the symmetry parameters, launches the cryoSPARC
        symmetry expansion job, monitors execution, and imports the expanded
        particle dataset back into the Scipion project for downstream
        processing and analysis.

        Inputs and Workflow

        The protocol requires an input set of particles containing projection
        alignment information. These particles typically originate from a
        previous refinement or reconstruction step in which particle
        orientations have already been estimated. The presence of alignment
        parameters is essential because symmetry expansion generates new
        symmetry-related orientations from the existing alignment solutions.

        The workflow begins with initialization of filename templates and
        protocol parameters. The protocol then prepares the cryoSPARC project
        environment, converts the input metadata, executes the symmetry
        expansion job, and finally imports the generated expanded particle
        dataset into Scipion.

        During execution, the protocol creates temporary STAR and cryoSPARC
        metadata files that allow compatibility between both software
        environments. Once the cryoSPARC job finishes, the generated
        `particles_expanded.cs` file is converted into STAR format and used
        to reconstruct the expanded Scipion particle set.

        Symmetry Definition

        The central parameter of the protocol is the symmetry definition.
        Users may specify standard point-group symmetries such as cyclic,
        dihedral, tetrahedral, octahedral, or icosahedral symmetry. Examples
        include C1, C4, D7, and other commonly used cryo-EM symmetry groups.

        Biologically, symmetry expansion allows each particle to contribute
        multiple symmetry-equivalent orientations. This is especially useful
        when investigating local asymmetry inside globally symmetric
        assemblies. For example, ligand binding, conformational flexibility,
        or partial occupancy may only affect a subset of symmetry-related
        regions. Symmetry expansion enables these regions to be analyzed
        independently during focused classification or local refinement.

        The protocol also supports helical symmetry parameters including
        helical twist, helical rise, and helical symmetry order. These values
        are commonly obtained from previous helical refinement jobs and are
        necessary when processing filamentous or helical biological systems.
        In these cases, symmetry expansion generates multiple symmetry-related
        representations along the helical lattice, improving the analysis of
        local structural variability within repetitive filament assemblies.

        Parameter Management and Validation

        Internally, the protocol stores all symmetry-related parameters inside
        a dedicated parameter list that is later transmitted to cryoSPARC.
        Parameters are only included when biologically meaningful values are
        provided. For example, helical rise and twist values are ignored if
        they are undefined or non-positive.

        Before execution, the protocol validates cryoSPARC accessibility and
        GPU configuration to ensure compatibility with the current execution
        environment. GPU validation is especially important because cryoSPARC
        processing workflows rely heavily on GPU acceleration for efficient
        large-scale cryo-EM computation.

        GPU Management and cryoSPARC Integration

        During execution, the protocol dynamically determines whether GPU
        resources should be assigned directly or managed through a queue
        system. This behavior allows compatibility with both standalone
        workstations and distributed computing infrastructures commonly used
        in cryo-EM facilities.

        The actual symmetry expansion job is launched through the cryoSPARC
        scheduling interface using `enqueueJob`. Once submitted, the protocol
        continuously monitors execution status using `waitForCryosparc`,
        ensuring synchronization between Scipion and cryoSPARC execution
        states. Intermediate cryoSPARC files are removed after completion in
        order to reduce unnecessary storage consumption.

        Outputs and Biological Interpretation

        After successful execution, the protocol generates an expanded set of
        particles in which every original particle has been duplicated across
        all symmetry-related orientations defined by the selected symmetry
        group. The expanded particles preserve the original acquisition and
        alignment metadata while introducing the additional transformed
        orientations generated during expansion.

        The resulting particle set is imported back into Scipion as an output
        `SetOfParticles` object and maintains the same dimensionality,
        alignment type, and sampling rate as the original dataset. A
        transformation relationship is also established between the original
        particles and the expanded particles so that downstream protocols can
        track their correspondence.

        From a biological perspective, symmetry expansion is not intended to
        improve global resolution directly. Instead, it provides a framework
        for studying local heterogeneity within symmetric assemblies. This is
        particularly valuable when analyzing flexible domains, asymmetric
        ligand binding, partial occupancy events, or localized conformational
        transitions that would otherwise be averaged out during standard
        symmetric reconstruction.

        In practical cryo-EM workflows, symmetry-expanded particles are
        commonly used for focused classification, local refinement, masked
        variability analysis, and detailed investigation of asymmetric
        structural features embedded within highly symmetric complexes.
        """
    _label = 'symmetry expansion'
    _className = "sym_expand"

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
        self.info(pwutils.yellowStr("Symmetry Expansion started..."))
        self.doSymmetryExpansion()

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to Relion file
        """
        self._initializeUtilsVariables()
        outputStarFn = self._getFileName('out_particles')
        csOutputFolder = os.path.join(self.projectDir.get(),
                                      self.runSymExp.get())
        csFileName = "particles_expanded.cs"

        # Copy the CS output expanded particles to extra folder
        copyFiles(csOutputFolder, self._getExtraPath(), files=[csFileName])

        csFile = os.path.join(self._getExtraPath(), csFileName)

        argsList = [csFile, outputStarFn]

        convertCs2Star(argsList)
        imgSet = self._getInputParticles()
        self.setFilePattern(imgSet.getFirstItem().getFileName())
        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(imgSet)
        outImgSet.setDim(imgSet.getDim())
        self._fillDataFromIter(outImgSet)

        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(imgSet, outImgSet)

    def _fillDataFromIter(self, imgSet):
        outImgsFn = 'particles@' + self._getFileName('out_particles')
        readSetOfParticles(outImgsFn, imgSet,
                           postprocessImageRow=self.updateParticlePath,
                           alignType=imgSet.getAlignment(),
                           samplingRate=imgSet.getSamplingRate())

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
            if not self.useQueueForSteps() and not self.useQueue():  # not using queue system
                gpusToUse = self.getGpuList()
            else:  # using queue system
                gpusToUse = False
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
                         "Please, go to cryoSPARC software for more "
                         "details.", self)
        clearIntermediateResults(self.projectName.get(), self.runSymExp.get())

