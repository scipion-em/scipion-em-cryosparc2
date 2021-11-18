# **************************************************************************
# *
# *  Authors:     Szu-Chi Chung (phonchi@stat.sinica.edu.tw) 
# *               Yunior C. Fonseca Reyna (cfonseca@cnb.csic.es)
# *
# * SABID Laboratory, Institute of Statistical Science, Academia Sinica
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
                                        LEVEL_ADVANCED, IntParam, Positive,
                                        BooleanParam, EnumParam)

from pwem import ALIGN_PROJ
from pwem.protocols import ProtInitialVolume, ProtClassify3D

from .protocol_base import ProtCryosparcBase
from ..convert import (defineArgs, convertCs2Star, cryosparcToLocation,
                       rowToAlignment)

from ..utils import (addSymmetryParam, addComputeSectionParams,
                     cryosparcValidate, gpusValidate, getSymmetry, enqueueJob,
                     calculateNewSamplingRate, waitForCryosparc,
                     clearIntermediateResults, fixVolume, copyFiles)
from ..constants import *


class ProtCryoSparcInitialModel(ProtCryosparcBase, ProtInitialVolume,
                                ProtClassify3D):
    """    
    Generate a 3D initial model _de novo_ from 2D particles using
    CryoSparc Stochastic Gradient Descent (SGD) algorithm.
    """
    _label = 'initial model'
    _className = "homo_abinit"
    # --------------------------- DEFINE param functions ----------------------

    def _defineFileNames(self):
        """ Centralize how files are called within the protocol. """
        myDict = {
                  'input_particles': self._getTmpPath('input_particles.star'),
                  'out_particles': self._getExtraPath() + '/output_particle.star',
                  'out_class': self._getExtraPath() + '/output_class.star'
                  }
        self._updateFilenamesDict(myDict)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      label="Input particles", important=True,
                      help='Select the input images from the project.')

        # --------------[Ab-Initio reconstruction]---------------------------

        form.addSection(label='Ab-Initio reconstruction')

        form.addParam('abinit_K', IntParam, default=1,
                      validators=[Positive],
                      label='Number of Ab-Initio classes:',
                      help='The number of classes. Each class will be randomly '
                           'initialized independently, unless an initial '
                           'structure was provided, in which case each class '
                           'will be a random variant of the initial structure')

        form.addParam('abinit_max_res', FloatParam, default=12.0,
                      validators=[Positive],
                      label='Maximum resolution (Angstroms):',
                      help='Maximum frequency to consider')

        form.addParam('abinit_init_res', FloatParam, default=35.0,
                      validators=[Positive],
                      label='Initial resolution (Angstroms):',
                      help='Starting frequency to consider')

        form.addParam('abinit_num_init_iters', IntParam, default=200,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label='Number of initial iterations:',
                      help='Number of initial iterations before annealing starts')

        form.addParam('abinit_num_final_iters', IntParam, default=300,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label='Number of final iterations:',
                      help='Number of final iterations after annealing ends')

        form.addParam('abinit_radwn_step', FloatParam, default=0.04,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label='Fourier radius step:',
                      help='Increase in Fourier radius a each iteration')

        form.addParam('abinit_window', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label='Window structures in real space:',
                      help='Softly window the reconstructions in real space at '
                           'each iteration')

        form.addParam('abinit_center', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label='Center structures in real space:',
                      help='Center the reconstructions in real space at each '
                           'iteration')

        form.addParam('abinit_scale_mg_correct', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Correct for per-micrograph optimal scales:',
                      help='(Experimental) Estimate and compute optimal scales '
                           'per micrograph')

        form.addParam('abinit_scale_compute', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Compute per-image optimal scales:',
                      help='(Experimental) Estimate and compute optimal scales '
                           'per image')

        form.addParam('abinit_mom', FloatParam, default=0,
                      expertLevel=LEVEL_ADVANCED,
                      label='SGD Momentum:',
                      help='Momentum for stochastic gradient descent')

        form.addParam('abinit_sparsity', FloatParam, default=0,
                      expertLevel=LEVEL_ADVANCED,
                      label='Sparsity prior:',
                      help='')

        form.addParam('abinit_minisize_init', IntParam, default=90,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label='Initial minibatch size:',
                      help='Number of images per minibatch at the beginning. '
                           'Set to zero to autotune')

        form.addParam('abinit_minisize', IntParam, default=300,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label='Final minibatch size:',
                      help='Final number of images per minibatch. Set to zero '
                           'to autotune')

        form.addParam('abinit_minisize_epsilon', FloatParam, default=0.05,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label='Abinit minisize epsilon:',
                      help='Parameter that controls batch size when autotuning '
                           'minibatch size. Set closer to zero for larger '
                           'batches')

        form.addParam('abinit_minisize_minp', FloatParam, default=0.01,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label='Abinit minisize minp:',
                      help='Parameter that controls how the batch size adjusts '
                           'to low probability classes when autotuning '
                           'minibatch sizes')

        form.addParam('abinit_minisize_num_init_iters', IntParam, default=300,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label='Initial minibatch size num iters:',
                      help='When to switch to final number of images per '
                           'minibatch')

        form.addParam('abinit_noise_model', EnumParam,
                      choices=['symmetric', 'white', 'coloured'],
                      default=0,
                      label='Noise model:',
                      help='Noise model to use. Valid options are white, '
                           'coloured or symmetric. Symmetric is the default, '
                           'meaning coloured with radial symmetry')

        form.addParam('abinit_noise_priorw', IntParam, default=50,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label='Noise priorw:',
                      help='Weight of the prior for estimating noise '
                           '(units of # of images)')

        form.addParam('abinit_noise_initw', IntParam, default=5000,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label='Noise initw:',
                      help='Weight of the initial noise estimate '
                           '(units of # of images)')

        form.addParam('abinit_class_anneal_beta', FloatParam, default=0.1,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label='Class similarity:',
                      help='Expected similarity of structures from different '
                           'classes. A number between 0 and 1. 0 means classes '
                           'are independent, 1 means classes are very similar)')

        form.addParam('abinit_class_anneal_start', IntParam, default=300,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label='Class similarity anneal start iter:',
                      help='Start point for annealing the similarity factor')

        form.addParam('abinit_class_anneal_end', IntParam, default=350,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label='Class similarity anneal end iter:',
                      help='Finish point for annealing the similarity factor')

        form.addParam('abinit_target_initial_ess_fraction', FloatParam,
                      default=0.011,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label='Target 3D ESS Fraction:',
                      help='Fraction of poses at the first iteration that '
                           'should have significant probability (used for '
                           'auto-tuning initial noise sigma-scale)')

        addSymmetryParam(form, help="Symmetry enforced (C, D, I, O, T). Eg. "
                                    "C1, D7, C4 etc. Enforcing symmetry above "
                                    "C1 is not recommended for ab-initio "
                                    "reconstruction")

        form.addParam('abinit_r_grid', FloatParam, default=25,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label='Abinit_r_grid:')

        form.addParam('abinit_high_lr_duration', FloatParam, default=100,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label='Initial learning rate duration:',
                      help='How long to apply the initial learning rate')

        form.addParam('abinit_high_lr', FloatParam, default=0.4,
                      validators=[Positive],
                      label='Initial learning rate:',
                      help='Learning rate (step size) used at the start of '
                           'optimization to help make rapid progress')

        form.addParam('abinit_nonneg', BooleanParam, default=True,
                      label='Enforce non-negativity:',
                      help='Enforce non-negativity of structures in real '
                           'space during optimization. Non-negativity is '
                           'recommended for ab-initio reconstruction')

        form.addParam('abinit_ignore_dc', BooleanParam, default=True,
                      label='Ignore DC component:',
                      help='Ignore the DC component of images. Should be true')

        form.addParam('abinit_init_radwn_cutoff', IntParam, default=7,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label='Initial structure lowpass (Fourier radius):',
                      help='Lowpass filter cutoff in Fourier radius for '
                           'initial random structures')

        form.addParam('abinit_search_start_iter', IntParam, default=200,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label='Abinit_search_start_iter:')

        form.addParam('abinit_use_engine', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label='Use fast codepaths:')

        form.addParam('intermediate_plots', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label='Show plots from intermediate steps:')

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
        print(pwutils.yellowStr("Ab Initial Model Generation Started..."), flush=True)
        self.doRunAbinit()

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to Relion file
        """
        print(pwutils.yellowStr("Creating the output..."), flush=True)
        self._initializeUtilsVariables()
        csOutputFolder = os.path.join(self.projectPath, self.projectName.get(),
                                      self.runAbinit.get())
        csFileName = "cryosparc_%s_%s_final_particles.cs" % (self.projectName.get(),
                                                             self.runAbinit.get())

        outputFolder = os.path.join(self._getExtraPath(), self.runAbinit.get())

        # Copy the CS output to extra folder
        copyFiles(csOutputFolder, outputFolder)

        csFile = os.path.join(outputFolder, csFileName)
        outputClassFn = self._getFileName('out_particles')
        argsList = [csFile, outputClassFn]
        parser = defineArgs()
        args = parser.parse_args(argsList)
        convertCs2Star(args)

        # Create model files for 3D classification
        self._createModelFile()

        imgSet = self._getInputParticles()
        classes3D = self._createSetOfClasses3D(imgSet)
        self._fillClassesFromIter(classes3D, self._getFileName('out_particles'))

        self._defineOutputs(outputClasses=classes3D)
        self._defineSourceRelation(self.inputParticles.get(), classes3D)

        # create a SetOfVolumes and define its relations
        volumes = self._createSetOfVolumes()
        vol = None

        for class3D in classes3D:
            vol = class3D.getRepresentative()
            vol.setObjId(class3D.getObjId())
            volumes.append(vol)

        volumes.setSamplingRate(vol.getSamplingRate())

        self._defineOutputs(outputVolumes=volumes)
        self._defineSourceRelation(self.inputParticles.get(), volumes)

    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        validateMsgs = cryosparcValidate()
        if not validateMsgs:
            validateMsgs = gpusValidate(self.getGpuList(), checkSingleGPU=True)
        return validateMsgs

    def _summary(self):
        summary = []
        if (not hasattr(self, 'outputVolumes') or
                not hasattr(self, 'outputClasses')):
            summary.append("Output objects not ready yet.")
        else:
            summary.append("Input Particles: %s" %
                           self.getObjectTag('inputParticles'))
            summary.append("Number of Ab-Initio classes: %s" %
                           str(self.abinit_K.get()))

            summary.append("Symmetry: %s" %
                           getSymmetry(self.symmetryGroup.get(),
                                       self.symmetryOrder.get()))
            summary.append("------------------------------------------")
            summary.append("Output volume %s" %
                           self.getObjectTag('outputVolumes'))
            summary.append("Output classes %s" %
                           self.getObjectTag('outputClasses'))

        return summary

    # --------------------------- UTILS functions ---------------------------
    def _loadClassesInfo(self, filename):
        """ Read some information about the produced CryoSparc Classes
        from the star file.
        """
        self._classesInfo = {}  # store classes info, indexed by class id

        modelStar = md.MetaData(filename)

        for classNumber, row in enumerate(md.iterRows(modelStar)):
            index, fn = cryosparcToLocation(row.getValue('rlnReferenceImage'))
            # Store info indexed by id, we need to store the row.clone() since
            # the same reference is used for iteration
            scaledFile = self._getScaledAveragesFile(fn, force=True)
            self._classesInfo[classNumber+1] = (index, scaledFile, row.clone())

    def _fillClassesFromIter(self, clsSet, filename):
        """ Create the SetOfClasses3D """
        outImgsFn = 'images@' + self._getFileName('out_class')
        self._loadClassesInfo(outImgsFn)
        clsSet.classifyItems(updateItemCallback=self._updateParticle,
                             updateClassCallback=self._updateClass,
                             itemDataIterator=md.iterRows('particles@'+ filename,
                                                          sortByLabel=md.RLN_IMAGE_ID))

    def _updateParticle(self, item, row):
        item.setClassId(row.getValue(md.RLN_PARTICLE_CLASS))
        item.setTransform(rowToAlignment(row, ALIGN_PROJ))

    def _updateClass(self, item):
        classId = item.getObjId()
        if classId in self._classesInfo:
            index, fn, row = self._classesInfo[classId]
            fixVolume(fn)
            item.setAlignmentProj()
            vol = item.getRepresentative()
            vol.setLocation(index, fn)
            vol.setSamplingRate(calculateNewSamplingRate(vol.getDim(),
                                                         self._getInputParticles().getSamplingRate(),
                                                         self._getInputParticles().getDim()))

    def _createModelFile(self):
        with open(self._getFileName('out_class'), 'w') as output_file:
            output_file.write('\n')
            output_file.write('data_images')
            output_file.write('\n\n')
            output_file.write('loop_')
            output_file.write('\n')
            output_file.write('_rlnReferenceImage')
            output_file.write('\n')
            for i in range(int(self.abinit_K.get())):
                row = ("%s/%s/cryosparc_%s_%s_class_%02d_final_volume.mrc\n"
                       % (self._getExtraPath(), self.runAbinit.get(),
                          self.projectName.get(), self.runAbinit.get(), i))
                output_file.write(row)

    def _defineParamsName(self):
        """ Define a list with all protocol parameters names"""
        self._paramsName = ['abinit_K', 'abinit_max_res', 'abinit_init_res',
                            'abinit_num_init_iters', 'abinit_num_final_iters',
                            'abinit_radwn_step', 'abinit_window',
                            'abinit_center',
                            'abinit_scale_mg_correct', 'abinit_scale_compute',
                            'abinit_mom', 'abinit_sparsity',
                            'abinit_minisize_init',
                            'abinit_minisize', 'abinit_minisize_epsilon',
                            'abinit_minisize_minp',
                            'abinit_minisize_num_init_iters',
                            'abinit_noise_model', 'abinit_noise_priorw',
                            'abinit_noise_initw', 'abinit_class_anneal_beta',
                            'abinit_class_anneal_start',
                            'abinit_class_anneal_end',
                            'abinit_target_initial_ess_fraction',
                            'abinit_symmetry',
                            'abinit_r_grid', 'abinit_high_lr_duration',
                            'abinit_high_lr',
                            'abinit_nonneg', 'abinit_ignore_dc',
                            'abinit_init_radwn_cutoff',
                            'abinit_search_start_iter', 'abinit_use_engine',
                            'intermediate_plots', 'compute_use_ssd']
        self.lane = str(self.getAttributeValue('compute_lane'))

    def doRunAbinit(self):
        """self._program + "  \'do_run_abinit(\"" + self.projectName +
        "\", \"" + self.workSpaceName + "\", \"\'" + self._user + "\'\", \""
        + self.par + "\",\"\'1\'\")\'")
        """
        input_group_connect = {"particles": self.particles.get()}
        params = {}

        if self.hasExpert():
            for paramName in self._paramsName:
                if paramName != 'abinit_symmetry' and paramName != 'abinit_noise_model':
                    params[str(paramName)] = str(self.getAttributeValue(paramName))
                elif paramName == 'abinit_symmetry':
                    symetryValue = getSymmetry(self.symmetryGroup.get(),
                                               self.symmetryOrder.get())
                    params[str(paramName)] = symetryValue
                elif paramName == 'abinit_noise_model':
                    params[str(paramName)] = str(NOISE_MODEL_CHOICES[self.abinit_noise_model.get()])

        # Determinate the GPUs to use (in dependence of
        # the cryosparc version)
        try:
            gpusToUse = self.getGpuList()
        except Exception:
            gpusToUse = False

        runAbinitJob = enqueueJob(self._className, self.projectName.get(),
                                    self.workSpaceName.get(),
                                    str(params).replace('\'', '"'),
                                    str(input_group_connect).replace('\'', '"'),
                                    self.lane,
                                    gpusToUse)

        self.runAbinit = String(runAbinitJob.get())
        self.currenJob.set(self.runAbinit.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(), self.runAbinit.get(),
                         "An error occurred in the initial volume process. "
                         "Please, go to cryosPARC software for more "
                         "details.")
        clearIntermediateResults(self.projectName.get(), self.runAbinit.get(), wait=7)
