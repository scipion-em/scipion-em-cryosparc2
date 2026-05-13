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
import emtable

import pwem.protocols as pwprot
from pwem import ALIGN_2D
from pyworkflow.object import String
from pyworkflow.protocol.params import (PointerParam, FloatParam, IntParam,
                                        BooleanParam, Positive)
import pyworkflow.utils as pwutils

from .protocol_base import ProtCryosparcBase
from ..convert import (rowToAlignment, convertCs2Star, cryosparcToLocation)
from ..utils import (addComputeSectionParams, cryosparcValidate, gpusValidate,
                     enqueueJob, waitForCryosparc, clearIntermediateResults,
                     copyFiles, getOutputPreffix, isCryosparcStandalone)
from ..constants import *


class ProtCryo2D(ProtCryosparcBase, pwprot.ProtClassify2D):
    """
    Wrapper to CryoSparc 2D clustering program.
    Classify particles into multiple 2D classes to facilitate stack cleaning
    and removal of junk particles. Also useful as a sanity check to
    investigate particle quality.

    AI Generated:

    CryoSPARC 2D Classification (ProtCryo2D) — User Manual
        Overview

        The CryoSPARC 2D Classification protocol organizes particle images into
        groups of structurally similar 2D projections. The main biological goal
        of this procedure is to separate meaningful particle views from noise,
        contaminants, damaged particles, aggregation artifacts, or incorrectly
        picked regions. In most cryo-EM workflows, 2D classification is one of
        the earliest and most important quality-control stages because it allows
        the user to evaluate whether the dataset contains interpretable signal
        before proceeding toward three-dimensional reconstruction.

        For biological users, this protocol is especially valuable for assessing
        sample quality, particle integrity, preferred orientation problems, and
        biochemical heterogeneity. Well-defined class averages usually indicate
        that the sample contains reproducible structural information, while poor
        or noisy classes often reveal ice contamination, broken particles, or
        inaccurate particle picking.

        Inputs and General Workflow

        The protocol requires a set of aligned or unaligned particle images as
        input. These particles are iteratively grouped according to similarity
        in appearance, generating representative class averages that summarize
        common structural views present in the dataset. During the process,
        particles are aligned in-plane and compared against evolving references
        until stable classes emerge.

        In routine cryo-EM processing, users commonly perform several rounds of
        2D classification. Early rounds are often permissive and intended to
        remove obvious junk particles, while later rounds may focus on refining
        the quality of the remaining dataset. Repeating classification with
        different numbers of classes can help reveal rare orientations or small
        structural populations that might otherwise remain hidden.

        Number of Classes and Biological Interpretation

        The number of requested classes strongly influences the outcome and
        interpretation of the analysis. Using relatively few classes tends to
        produce broad averages that capture dominant structural features but may
        merge distinct conformations or orientations together. Increasing the
        number of classes provides finer separation and can reveal structural
        variability, rare views, or minor populations.

        From a biological perspective, there is no universally optimal number
        of classes. Small and homogeneous datasets may require only a modest
        number of classes, whereas large or structurally heterogeneous datasets
        often benefit from many more. Excessively large numbers of classes,
        however, may fragment the signal and produce unstable or noisy averages.

        Resolution and Computational Considerations

        The protocol allows users to define the effective resolution range used
        during classification. Lower-resolution settings generally increase
        robustness and speed, making them suitable for noisy datasets or early
        exploratory processing. Higher-resolution settings can reveal finer
        structural details but require greater computational resources and may
        become unstable when particle quality is limited.

        In biological practice, starting with moderate resolution settings is
        usually recommended. Once clear particle populations have been isolated,
        more demanding classifications can be performed if additional structural
        detail is needed.

        Circular Masking and Particle Centering

        Circular masking is an important feature because it restricts the region
        of the image contributing to classification. In most datasets, applying
        a circular mask improves robustness by suppressing noisy image corners
        and solvent background. The mask diameter should approximately match the
        particle size while avoiding unnecessary surrounding noise.

        For elongated or filamentous assemblies, careful adjustment of masking
        becomes especially important. Overly restrictive masks may remove real
        structural information, whereas masks that are too large can allow noise
        to dominate classification.

        Re-centering options help maintain particles and class averages aligned
        around the image center. This is biologically important because drifting
        class averages can produce blurred reconstructions or distorted views.
        Proper centering is particularly beneficial for asymmetric complexes,
        flexible assemblies, and datasets with broad orientation variability.

        Classification Uncertainty and Structural Diversity

        The protocol includes parameters controlling the degree of uncertainty
        maintained during early classification iterations. Biologically, this
        affects how aggressively particles are separated into classes.

        Lower uncertainty tends to produce rapid convergence and can strongly
        isolate junk particles, but may reduce structural diversity by forcing
        particles prematurely into specific classes. Higher uncertainty allows
        the algorithm to remain more flexible for longer periods, which can help
        preserve subtle conformational differences or weak particle populations.

        For heterogeneous biological systems such as multi-domain proteins,
        flexible complexes, or membrane assemblies, maintaining moderate initial
        uncertainty often improves the recovery of meaningful structural states.

        Filament and Helical Data

        The protocol includes support for filamentous and helical specimens.
        These datasets present unique challenges because particles frequently
        appear in continuous orientations and may exhibit directional ambiguity.
        Vertical alignment of filament classes helps standardize the appearance
        of class averages and facilitates interpretation of repeating helical
        features.

        Biological users working with cytoskeletal assemblies, amyloid fibrils,
        or filament-forming proteins may benefit substantially from enabling
        filament-oriented alignment options during classification.

        Noise Modeling and Regularization

        Noise handling and regularization are central to obtaining reliable 2D
        class averages. The protocol incorporates statistical regularization
        methods intended to reduce overfitting and improve stability during
        iterative refinement. These approaches help ensure that class averages
        represent reproducible structural signal rather than amplification of
        random noise.

        In practical biological applications, the default regularization
        settings are often sufficient. More advanced users may adjust noise
        behavior when working with unusually small particles, very low signal-
        to-noise datasets, or highly heterogeneous samples.

        Outputs and Their Interpretation

        The protocol produces a set of 2D class averages together with particle
        assignments for each class. These outputs allow users to visually assess
        dataset quality and decide which particles should be retained for
        downstream refinement.

        Good classes typically display recognizable structural features, clear
        secondary structure patterns at high quality, and reproducible particle
        orientations. Poor classes often contain diffuse density, inconsistent
        shapes, contamination, ice artifacts, or overlapping particles.

        Biologically, the retained particles define the quality ceiling of the
        remainder of the cryo-EM workflow. Careful inspection and selection of
        classes at this stage significantly influence the success of later 3D
        reconstruction and refinement.

        Practical Recommendations

        In most biological workflows, it is advisable to begin with a moderate
        number of classes and default alignment settings. After inspecting the
        resulting averages, users can iteratively refine the dataset by removing
        poor classes and rerunning classification on the cleaned particles.

        Datasets with severe heterogeneity or contamination often benefit from
        multiple sequential rounds of classification. Conversely, very clean and
        homogeneous samples may require only minimal cleaning before proceeding
        to ab initio reconstruction or high-resolution refinement.

        When class averages appear unstable or poorly centered, adjusting mask
        sizes, enabling re-centering, or increasing the number of iterations can
        substantially improve interpretability.

        Final Perspective

        For most cryo-EM practitioners, 2D classification is far more than a
        computational sorting procedure. It is a biologically meaningful quality
        assessment stage that determines whether the dataset contains coherent
        structural information suitable for downstream analysis. Careful tuning
        of class diversity, masking strategy, and alignment behavior allows the
        protocol to adapt to a wide range of biological specimens, from rigid
        globular proteins to highly flexible or filamentous assemblies.
    """
    _label = '2D classification'
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
                      allowsPointers=True,
                      condition='useCircular2D==True')

        form.addParam('class2D_window_outer_A', FloatParam, default=None,
                      label='Circular mask diameter outer (A)',
                      help='The outer diameter (in Angstroms) of the window. '
                           'If None, outer diameter is 20 percent larger than '
                           'inner diameter. The window mask transitions '
                           'smoothly between inner and outer diameters.',
                      allowsNull=True,
                      allowsPointers=True,
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
        form.addParam('class2D_estimate_in_plane_pose', BooleanParam, default=False,
                      label='Align filament classes vertically',
                      help='Set to True only if the particle images are of filamentous/helical assemblies - this will'
                           ' align all class averages vertically in the second last iteration, enabling estimation '
                           'of in-plane rotation. Note that this will not attempt to estimate the relative polarity '
                           'of class averages.')
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
        self.info(pwutils.yellowStr("2D Classifications Started..."))
        self.doRunClass2D()

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to Relion file
        """
        self.info(pwutils.yellowStr("Creating the output..."))
        self._initializeUtilsVariables()

        csOutputFolder = os.path.join(self.projectDir.get(),
                                      self.runClass2D.get())
        _numberOfIterSuffix = self._getNumberOfIterSuffix()

        csOutputPattern = "%s%s%s" % (getOutputPreffix(self.projectName.get()),
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
        convertCs2Star(argsList)

        csClassAverageFile = os.path.join(self._getExtraPath(),
                                          csClassAveragesName)
        outputClassFn = self._getFileName('out_class')
        argsList = [csClassAverageFile, outputClassFn]

        convertCs2Star(argsList)

        self._createModelFile()
        self._loadClassesInfo(self._getFileName('out_class_m2'))
        # Use the pointer with extended (indirect)
        classes2DSet = self._createSetOfClasses2D(self.inputParticles)
        self._fillClassesFromLevel(classes2DSet)

        self._defineOutputs(outputClasses=classes2DSet)
        self._defineSourceRelation(self.inputParticles, classes2DSet)

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

        mdFileName = '%s@%s' % ('particles', filename)
        table = emtable.Table(fileName=filename)

        for classNumber, row in enumerate(table.iterRows(mdFileName)):
            index, fn = cryosparcToLocation(
                row.get(RELIONCOLUMNS.rlnImageName.value))

            # Store info indexed by id, we need to store the row.clone() since
            # the same reference is used for iteration
            scaledFile = self._getScaledAveragesFile(fn)
            self._classesInfo[classNumber + 1] = (index, scaledFile, row)
        self._numClass = index

    def _fillClassesFromLevel(self, clsSet):
        """ Create the SetOfClasses2D from a given iteration. """

        # the particle with orientation parameters (all_parameters)
        xmpMd = 'particles@' + self._getFileName("out_particles")

        clsSet.classifyItems(updateItemCallback=self._updateParticle,
                             updateClassCallback=self._updateClass,
                             itemDataIterator=emtable.Table.iterRows(
                                 xmpMd),
                             raiseOnNextFailure=False,
                             cancelNextWhenAppendIsFalse=True)  # relion style

    def _updateParticle(self, item, row):
        item.setClassId(row.get(RELIONCOLUMNS.rlnClassNumber.value))
        samplingRate = item.getSamplingRate()
        item.setTransform(rowToAlignment(row, ALIGN_2D, samplingRate))

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
                    output_file.write(
                        row % (classNumber, self._getExtraPath(), image))
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
                  "class2D_sigma_init_factor": str(
                      self.initialClassification.get()),
                  "class2D_window": str(self.useCircular2D.get()),
                  "class2D_recenter": str(self.reCenter2D.get()),
                  "class2D_recenter_thresh": str(self.reCenterMask.get()),
                  "class2D_recenter_binary": str(self.reCenterMaskBinary.get()),
                  "class2D_estimate_in_plane_pose": str(self.class2D_estimate_in_plane_pose.get()),
                  "class2D_force_max": str(self.forceMaxover.get()),
                  "class2D_ctf_phase_flip_only": str(self.ctfFlipPhases.get()),
                  "class2D_num_full_iter": str(self.numberFinalIterator.get()),
                  "class2D_num_full_iter_batch": str(
                      self.numberOnlineEMIterator.get()),
                  "class2D_num_full_iter_batchsize_per_class": str(
                      self.batchSizeClass.get()),
                  "class2D_init_scale": str(self.initialScale2D.get()),
                  "class2D_zp_factor": str(self.zeropadFactor.get()),
                  "class2D_use_frc_reg": str(self.useFRCRegularized.get()),
                  "class2D_use_frc_reg_full": str(self.useFullFRC.get()),
                  "class2D_sigma_init_iter": str(
                      self.iterationToStartAnneling.get()),
                  "class2D_sigma_num_anneal_iters": str(
                      self.iterationToStartAnneal.get()),
                  "class2D_sigma_use_white": str(self.useWhiteNoiseModel.get()),
                  "intermediate_plots": str('False'),
                  "compute_use_ssd": str(self.compute_use_ssd.get())}
        if self.class2D_window_inner_A.get() is not None:
            params["class2D_window_inner_A"] = str(
                self.class2D_window_inner_A.get())
        if self.class2D_window_outer_A.get() is not None:
            params["class2D_window_outer_A"] = str(
                self.class2D_window_outer_A.get())
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
            if not self.useQueueForSteps() and not self.useQueue(): # not using queue system
                gpusToUse = self.getGpuList()
                numberGPU = len(gpusToUse)
            else: # using queue system
                gpusToUse = False
                numberGPU = 1
        except Exception:
            gpusToUse = False
            numberGPU = 1
        params = self.assignParamValue()
        if not isCryosparcStandalone():  # Cluster case
            gpusToUse = False
            numberGPU = self.compute_num_gpus.get()

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
                         "Please, go to cryoSPARC software for more "
                         "details.", self)
        clearIntermediateResults(self.projectName.get(), self.runClass2D.get())

