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
import pwem.protocols as pwprot
from pwem import ALIGN_2D
from pyworkflow.object import String
from pyworkflow.protocol.params import (PointerParam, FloatParam, IntParam,
                                        BooleanParam, Positive)
import pyworkflow.utils as pwutils

from .protocol_base import ProtCryosparcBase
from ..convert import (rowToAlignment, defineArgs, convertCs2Star,
                       cryosparcToLocation)
from ..utils import (addComputeSectionParams, cryosparcValidate, gpusValidate,
                     enqueueJob, waitForCryosparc, clearIntermediateResults,
                     copyFiles)
from ..constants import *


class ProtCryo2D(ProtCryosparcBase, pwprot.ProtClassify2D):
    """ Wrapper to CryoSparc 2D clustering program.
        Classify particles into multiple 2D classes to facilitate stack cleaning
        and removal of junk particles. Also useful as a sanity check to
        investigate particle quality.
    """
    _label = '2d classification'
    IS_2D = True
    _className = "class_2D"
    
    def __init__(self, **args):
        pwprot.ProtClassify2D.__init__(self, **args)
        if self.numberOfMpi.get() < 2:
            self.numberOfMpi.set(2)
    
    def _defineFileNames(self):
        """ Centralize how files are called within the protocol. """
        myDict = {
                  'input_particles': self._getTmpPath('input_particles.star'),
                  'out_particles': self._getExtraPath() + '/output_particle.star',
                  'out_class': self._getExtraPath() + '/output_class.star',
                  'out_class_m2': self._getExtraPath() + '/output_class_m2.star'
                  }
        self._updateFilenamesDict(myDict)
    
    # --------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      label="Input particles", important=True,
                      help='Select the input images from the project.')

        # ----------- [2D Classification] --------------------------------

        form.addSection(label="2D Classification")
        form.addParam('numberOfClasses', IntParam, default=50,
                      validators=[Positive],
                      label='Number of classes:',
                      help='The number of 2D classes into which to sort the '
                           'dataset. Runtime is approximately linear in the '
                           'number of classes. Often, sorting the same dataset '
                           'into different numbers of classes can be helpful in '
                           'finding junk particles or rare views.')

        form.addParam('maximunResolution', IntParam, default=6,
                      validators=[Positive],
                      label='Maximum resolution (A)',
                      help='The maximum resolution in Angstroms to consider when '
                           'aligning and reconstructing 2D classes. This setting '
                           'controls the box size that is used internally, and '
                           'higher resolutions can slow down processing.')
        form.addParam('initialClassification', FloatParam, default=2.0,
                      label='Initial classification uncertainty factor',
                      validators=[Positive],
                      help='This factor (a number greater than 1) controls the '
                           'initial search for 2D references. A value of 1.0 '
                           'indicates that the search should quickly become '
                           'certain about classes and assignments, leading to '
                           'the algorithm finding more "junk" classes. A value '
                           'larger than 1.0 (usually between 2 and 10) causes '
                           'the algorithm to remain uncertain about classes and '
                           'assignments for more iterations, resulting in more '
                           'diversity of "good" classes.')

        form.addParam('useCircular2D', BooleanParam, default=True,
                      label='Use circular mask on 2D classes?',
                      help='Whether or not to apply a circular window to the 2D '
                           'classes during classification. This ensures that '
                           'each 2D class has no density outside the circular '
                           'window. By default, the window is a circle that only '
                           'masks out the corners of the 2D classes.')

        form.addParam('class2D_window_inner_A', FloatParam, default=None,
                      label='Circular mask diameter (A)',
                      help='The inner diameter (in Angstroms) of the window '
                            'that is applied to 2D classes during '
                            'classification. If None, the window only masks out '
                            'the corners of each 2D class.',
                      allowsNull=True,
                      condition='useCircular2D==True')

        form.addParam('class2D_window_outer_A', FloatParam, default=None,
                      label='Circular mask diameter outer (A)',
                      help='The outer diameter (in Angstroms) of the window. '
                           'If None, outer diameter is 20 percent larger than '
                           'inner diameter. The window mask transitions '
                           'smoothly between inner and outer diameters.',
                      allowsNull=True,
                      condition='useCircular2D==True')

        form.addParam('reCenter2D', BooleanParam, default=True,
                      label='Re-center 2D classes',
                      help='Whether or not to re-center 2D class references at '
                           'every iteration to avoid drift of density away from '
                           'the center of the box. This option is often '
                           'important to keep classes centered and avoid '
                           'artefacts near the edges of the box.')

        form.addParam('reCenterMask', FloatParam, default=0.2,
                      validators=[Positive],
                      label='Re-center mask threshold',
                      help='2D classes are recentered by computing the '
                           'center-of-mass (COM) of pixels that are above this '
                           'threshold value. The threshold is relative to the '
                           'maximum density value in the reference, so 0.2 means '
                           'pixels with greater than 20%% of the maximum density.')

        form.addParam('reCenterMaskBinary', BooleanParam, default=False,
                      label='Re-center mask binary',
                      help='If True, compute the COM for re-centering by equally '
                           'weighting every pixel that was above the threshold. '
                           'If False, weight every pixel by its greyscale '
                           'density value.')

        form.addParam('forceMaxover', BooleanParam, default=True,
                      label='Force Max over poses/shifts',
                      help='If True, maximize over poses and shifts when '
                           'aligning particles to references. If False, '
                           'marginalize over poses and shifts to account for '
                           'alignment uncertainty. This is generally not '
                           'necessary, but can provide better results with very '
                           'small or low SNR particles.')

        form.addParam('ctfFlipPhases', BooleanParam, default=False,
                      label='CTF flip phases only',
                      help='Treat the CTF by flipping phases only, rather that '
                           'correctly accounting for amplitude and phase. Not '
                           'recommended.')

        form.addParam('numberFinalIterator', IntParam, default=1,
                      validators=[Positive],
                      label='Number of final full iterations',
                      help='The number of final full passes through the dataset '
                           'at the end of classification. Usually only one full '
                           'pass is needed.')

        form.addParam('numberOnlineEMIterator', IntParam, default=20,
                      validators=[Positive],
                      label='Number of online-EM iterations',
                      help='The total number of iterations of online-EM to '
                           'perform. Typically 20 is enough, but for small or '
                           'low SNR particles, or when classifying subsets that '
                           'have few distinct views, a larger number like 40 '
                           'can help.')

        form.addParam('batchSizeClass', IntParam, default=100,
                      validators=[Positive],
                      label='Batchsize per class',
                      help='The number of particles per class to use during each '
                           'iteration of online-EM. For small or low SNR '
                           'particles, this can be increased to 200.')

        form.addParam('initialScale2D', IntParam, default=1,
                      validators=[Positive],
                      label='2D initial scale',
                      help='Initial scale of random starting references. Not '
                           'recommended to change.')
        form.addParam('zeropadFactor', IntParam, default=2,
                      validators=[Positive],
                      label='2D zeropad factor',
                      help='Zeropadding factor. For very large box particles, '
                           'this can be reduced to speed up computation and '
                           'reduce memory requirements.')

        form.addParam('useFRCRegularized', BooleanParam, default=True,
                      label='Use FRC based regularizer',
                      help='Use an FRC based regularizer to avoid overfitting '
                           'during classification.')

        form.addParam('useFullFRC', BooleanParam, default=True,
                      label='Use full FRC')

        form.addParam('iterationToStartAnneling', IntParam, default=2,
                      validators=[Positive],
                      label='Iteration to start annealing sigma',
                      help='Iteration at which noise model should be annealed. '
                           'Not recommended to change.')

        form.addParam('iterationToStartAnneal', IntParam, default=15,
                      validators=[Positive],
                      label='Number of iteration to anneal sigma',
                      help='Number of iterations over which to anneal noise '
                           'model. Not recommended to change.')

        form.addParam('useWhiteNoiseModel', BooleanParam, default=False,
                      label='Use white noise model',
                      help='Force the use of a white noise model.')

        # ----------- [Compute settings] --------------------------------
        form.addSection(label="Compute settings")
        addComputeSectionParams(form)

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
        """
        Classify particles into multiples 2D classes
        """
        print(pwutils.yellowStr("2D Classifications Started..."), flush=True)
        self.doRunClass2D()

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to Relion file
        """
        print(pwutils.yellowStr("Creating the output..."), flush=True)
        self._initializeUtilsVariables()

        csOutputFolder = os.path.join(self.projectPath, self.projectName.get(),
                                      self.runClass2D.get())
        _numberOfIterSuffix = self._getNumberOfIterSuffix()
        csOutputPattern = "cryosparc_%s_%s%s" % (self.projectName.get(),
                                                 self.runClass2D.get(),
                                                 _numberOfIterSuffix)
        csParticlesName = csOutputPattern + "_particles.cs"
        csClassAveragesName = csOutputPattern + "_class_averages.cs"
        mrcFileName = csOutputPattern + "_class_averages.mrc"

        # Copy the CS output to extra folder
        copyFiles(csOutputFolder, self._getExtraPath(), files=[csParticlesName,
                                                               csClassAveragesName,
                                                               mrcFileName])

        csPartFile = os.path.join(self._getExtraPath(), csParticlesName)
        outputStarFn = self._getFileName('out_particles')
        argsList = [csPartFile, outputStarFn]

        parser = defineArgs()
        args = parser.parse_args(argsList)
        convertCs2Star(args)

        csClassAverageFile = os.path.join(self._getExtraPath(), csClassAveragesName)
        outputClassFn = self._getFileName('out_class')
        argsList = [csClassAverageFile, outputClassFn]

        parser = defineArgs()
        args = parser.parse_args(argsList)
        convertCs2Star(args)

        self._createModelFile()
        self._loadClassesInfo(self._getFileName('out_class_m2'))
        inputParticles = self._getInputParticles()
        classes2DSet = self._createSetOfClasses2D(inputParticles)
        self._fillClassesFromLevel(classes2DSet)
  
        self._defineOutputs(outputClasses=classes2DSet)
        self._defineSourceRelation(self.inputParticles.get(), classes2DSet)

    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        validateMsgs = cryosparcValidate()
        if not validateMsgs:
            validateMsgs = gpusValidate(self.getGpuList())
        return validateMsgs

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputClasses'):
            summary.append("Output classes not ready yet.")
        else:
            summary.append("Input Particles: %s" %
                           self.getObjectTag('inputParticles'))
            summary.append("Classified into *%d* classes." %
                           self.numberOfClasses.get())
            summary.append("Output set: %s" %
                           self.getObjectTag('outputClasses'))

        return summary

    def _methods(self):
        methods = "We classified input particles %s (%d items) " % (
            self.getObjectTag('inputParticles'),
            self._getInputParticles().getSize())
        methods += "into %d classes using CryoSparc " % self.numberOfClasses.get()
        return [methods]
    
    # --------------------------- UTILS functions ------------------------------
    def _loadClassesInfo(self, filename):
        """ Read some information about the produced 2D classes
        from the metadata file.
        """
        self._classesInfo = {}  # store classes info, indexed by class id

        # mdClasses = md.MetaData(filename)
        mdClasses = md.MetaData('%s@%s' % ('particles', filename))

        for classNumber, row in enumerate(md.iterRows(mdClasses)):
            index, fn = cryosparcToLocation(row.getValue('rlnImageName'))
            # Store info indexed by id, we need to store the row.clone() since
            # the same reference is used for iteration
            scaledFile = self._getScaledAveragesFile(fn)
            self._classesInfo[classNumber + 1] = (index, scaledFile, row.clone())
        self._numClass = index

    def _fillClassesFromLevel(self, clsSet):
        """ Create the SetOfClasses2D from a given iteration. """

        # the particle with orientation parameters (all_parameters)
        xmpMd = 'particles@' + self._getFileName("out_particles")

        clsSet.classifyItems(updateItemCallback=self._updateParticle,
                             updateClassCallback=self._updateClass,
                             itemDataIterator=md.iterRows(xmpMd,
                                                          sortByLabel=md.MDL_ITEM_ID))  # relion style

    def _updateParticle(self, item, row):
        item.setClassId(row.getValue(md.RLN_PARTICLE_CLASS))
        item.setTransform(rowToAlignment(row, ALIGN_2D))
        
    def _updateClass(self, class2D):
        classId = class2D.getObjId()
        if classId in self._classesInfo:
            index, fn, row = self._classesInfo[classId]
            class2D.setAlignment2D()
            class2Drep = class2D.getRepresentative()
            class2Drep.setLocation(index, fn)
            class2Drep.setSamplingRate(class2D.getSamplingRate())

    def _createModelFile(self):
        with open(self._getFileName('out_class'), 'r') as input_file, \
                open(self._getFileName('out_class_m2'), 'w') as output_file:
            for line in input_file:
                if "@" in line:
                    row = "%s@%s/%s"
                    classNumber = line.split('@')[0]
                    image = line.split('/')[1]
                    output_file.write(row % (classNumber, self._getExtraPath(), image))
                else:
                    output_file.write(line)

    def _getNumberOfIterSuffix(self):
        _numberOfIter = (self.numberOnlineEMIterator.get() +
                         self.numberFinalIterator.get() - 1)
        _numberOfIterSuffix = "_00%s" % str(self.numberOnlineEMIterator.get())
        if _numberOfIter > 9:
            _numberOfIterSuffix = "_0%s" % str(_numberOfIter)
        if _numberOfIter > 99:
            _numberOfIterSuffix = "_%s" % str(_numberOfIter)
        return _numberOfIterSuffix

    def _defineParamsName(self):
        """ Define a list with all protocol parameters names"""
        self.lane = str(self.getAttributeValue('compute_lane'))

    def assignParamValue(self):
        params = {"class2D_K": str(self.numberOfClasses.get()),
                  "class2D_max_res": str(self.maximunResolution.get()),
                  "class2D_sigma_init_factor": str(self.initialClassification.get()),
                  "class2D_window": str(self.useCircular2D.get()),
                  "class2D_recenter": str(self.reCenter2D.get()),
                  "class2D_recenter_thresh": str(self.reCenterMask.get()),
                  "class2D_recenter_binary": str(self.reCenterMaskBinary.get()),
                  "class2D_force_max": str(self.forceMaxover.get()),
                  "class2D_ctf_phase_flip_only": str(self.ctfFlipPhases.get()),
                  "class2D_num_full_iter": str(self.numberFinalIterator.get()),
                  "class2D_num_full_iter_batch": str(self.numberOnlineEMIterator.get()),
                  "class2D_num_full_iter_batchsize_per_class": str(self.batchSizeClass.get()),
                  "class2D_init_scale": str(self.initialScale2D.get()),
                  "class2D_zp_factor": str(self.zeropadFactor.get()),
                  "class2D_use_frc_reg": str(self.useFRCRegularized.get()),
                  "class2D_use_frc_reg_full": str(self.useFullFRC.get()),
                  "class2D_sigma_init_iter": str(self.iterationToStartAnneling.get()),
                  "class2D_sigma_num_anneal_iters": str(self.iterationToStartAnneal.get()),
                  "class2D_sigma_use_white": str(self.useWhiteNoiseModel.get()),
                  "intermediate_plots": str('False'),
                  "compute_use_ssd": str(self.compute_use_ssd.get())}
        if self.class2D_window_inner_A.get() is not None:
                params["class2D_window_inner_A"] = str(self.class2D_window_inner_A.get())
        if self.class2D_window_outer_A.get() is not None:
                params["class2D_window_outer_A"] = str(self.class2D_window_outer_A.get())
        return params

    def doRunClass2D(self):
        """
        do_run_class_2D:  do_job(job_type, puid='P1', wuid='W1',
                                 uuid='devuser', params={},
                                 input_group_connects={})
        returns: the new uid of the job that was created
        """
        # {'particles' : 'JXX.imported_particles' }
        input_group_connect = {"particles": self.particles.get()}

        # Determinate the GPUs or the number of GPUs to use (in dependence of
        # the cryosparc version)
        try:
            gpusToUse = self.getGpuList()
            numberGPU = len(gpusToUse)
        except Exception:
            gpusToUse = False
            numberGPU = 1

        params = self.assignParamValue()
        params["compute_num_gpus"] = str(numberGPU)
        runClass2DJob = enqueueJob(self._className, self.projectName.get(),
                                     self.workSpaceName.get(),
                                     str(params).replace('\'', '"'),
                                     str(input_group_connect).replace('\'', '"'),
                                     self.lane, gpusToUse)

        self.runClass2D = String(runClass2DJob.get())
        self.currenJob.set(runClass2DJob.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(), self.runClass2D.get(),
                         "An error occurred in the 2D classification process. "
                         "Please, go to cryosPARC software for more "
                         "details.")
        clearIntermediateResults(self.projectName.get(), self.runClass2D.get())

