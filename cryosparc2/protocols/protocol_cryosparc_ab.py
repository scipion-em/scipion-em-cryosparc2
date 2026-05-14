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

import pyworkflow.utils as pwutils
from pyworkflow.object import String
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        LEVEL_ADVANCED, IntParam, Positive,
                                        BooleanParam, EnumParam)

from pwem import ALIGN_PROJ
from pwem.protocols import ProtInitialVolume, ProtClassify3D

from .protocol_base import ProtCryosparcBase
from ..convert import (convertCs2Star, cryosparcToLocation,
                       rowToAlignment)

from ..utils import (addSymmetryParam, addComputeSectionParams,
                     cryosparcValidate, gpusValidate, getSymmetry, enqueueJob,
                     calculateNewSamplingRate, waitForCryosparc,
                     clearIntermediateResults, fixVolume, copyFiles,
                     getOutputPreffix, matchItemRow)
from ..constants import *


class ProtCryoSparcInitialModel(ProtCryosparcBase, ProtInitialVolume,
                                ProtClassify3D):
    """
    Generate a 3D initial model _de novo_ from 2D particles using
    CryoSparc Stochastic Gradient Descent (SGD) algorithm.
    """

    """
        Generates ab initio 3D initial models from 2D particle images using
        the CryoSPARC stochastic gradient descent (SGD) reconstruction algorithm.
        The protocol is designed to estimate one or multiple initial volumes
        directly from particle data without requiring a prior structural model.

        AI Generated:

        Initial Model Reconstruction (ProtCryoSparcInitialModel) — User Manual
            Overview

            The Initial Model protocol reconstructs one or several 3D ab initio
            volumes from 2D cryo-EM particles using the CryoSPARC homogeneous
            ab initio reconstruction algorithm. Its primary objective is to
            generate an initial structural model suitable for downstream
            refinement, heterogeneous classification, or structural exploration
            in cryo-electron microscopy workflows.

            From a biological perspective, this protocol is typically used
            during the earliest stages of structural determination, when no
            reliable reference structure is available. Instead of relying on
            template-based approaches, the method reconstructs structures
            directly from experimental particle images through iterative
            stochastic optimization. This makes it especially valuable for
            novel macromolecular complexes, unknown conformational states,
            or datasets with significant structural variability.

            Inputs and General Workflow

            The protocol requires a set of particles with associated CTF
            information. These particles represent the experimental projections
            used to estimate the initial 3D reconstruction. Proper preprocessing,
            including particle extraction and CTF estimation, is essential for
            obtaining biologically meaningful results.

            The reconstruction process begins by initializing one or several
            random 3D volumes. During optimization, CryoSPARC iteratively
            estimates particle orientations, class assignments, and volume
            updates using stochastic gradient descent. Over successive iterations,
            the structures progressively evolve from low-resolution random
            densities toward biologically interpretable maps.

            The protocol supports generation of multiple ab initio classes.
            This functionality is particularly important for heterogeneous
            datasets containing different conformations, compositional states,
            or partially damaged particles. Each class evolves independently,
            allowing the algorithm to separate distinct structural populations.

            Number of Classes and Structural Heterogeneity

            One of the most biologically relevant parameters is the number
            of ab initio classes. Using a single class is appropriate for
            homogeneous datasets where all particles are expected to belong
            to the same structural state. In contrast, heterogeneous datasets
            often benefit from multiple classes because the protocol can
            separate different conformations or assemblies during reconstruction.

            Increasing the number of classes may improve structural separation,
            but it also increases computational complexity and the risk of
            generating noisy or poorly populated classes. In practical cryo-EM
            workflows, users often start with a small number of classes and
            iteratively refine the strategy based on the biological quality
            of the resulting maps.

            Resolution Progression and Optimization Strategy

            The reconstruction follows a progressive optimization scheme.
            The initial resolution parameter controls the starting frequency
            information used during early iterations. Beginning with low-resolution
            information stabilizes orientation estimation and avoids overfitting
            noise. As optimization progresses, the maximum resolution parameter
            determines the highest spatial frequency incorporated into the model.

            The protocol separates optimization into initial and final iteration
            phases. Early iterations focus on global convergence and orientation
            stabilization, while later iterations progressively refine structural
            details. The Fourier radius step parameter controls how rapidly
            higher-resolution information is introduced during optimization.

            From a biological standpoint, gradual resolution annealing improves
            robustness, especially for difficult datasets with low signal-to-noise
            ratios or strong preferred orientations.

            SGD Optimization and Learning Parameters

            The reconstruction is driven by stochastic gradient descent.
            Several parameters regulate optimization stability and convergence,
            including learning rate, learning rate duration, momentum,
            minibatch size, and sparsity regularization.

            The learning rate controls the magnitude of volume updates during
            optimization. Higher values accelerate convergence but may introduce
            instability, whereas lower values improve stability at the expense
            of runtime. The protocol initially applies an elevated learning rate
            to accelerate exploration of the solution space before transitioning
            to more stable refinement phases.

            Momentum can improve convergence by smoothing stochastic updates,
            although excessive momentum may destabilize difficult datasets.
            Sparsity regularization may help suppress noisy regions in some
            reconstructions, particularly for small particles or weak datasets.

            Minibatch sizes determine how many particles contribute to each
            optimization step. Small minibatches introduce stochastic variability
            that may help escape local minima, whereas larger minibatches provide
            smoother and more stable optimization. The protocol also supports
            automatic minibatch tuning for adaptive optimization behavior.

            Noise Modeling and Stability

            Cryo-EM particle images contain substantial experimental noise,
            making noise estimation a central component of ab initio reconstruction.
            The protocol supports symmetric, white, and coloured noise models.

            The symmetric model is generally recommended for most biological
            datasets because it approximates coloured noise while assuming
            radial symmetry in Fourier space. White noise is computationally
            simpler but biologically less realistic for most cryo-EM experiments.
            The coloured noise model may improve reconstruction quality in
            datasets with complex background characteristics.

            Additional parameters regulate the prior and initialization weights
            of the noise estimation procedure. These parameters influence
            optimization stability during early iterations, particularly when
            dealing with low particle counts or noisy experimental conditions.

            Symmetry and Structural Constraints

            The protocol supports symmetry enforcement during reconstruction,
            including cyclic, dihedral, tetrahedral, octahedral, and icosahedral
            symmetries. Applying symmetry can dramatically improve reconstruction
            quality when the biological complex genuinely possesses the specified
            symmetry.

            However, enforcing incorrect symmetry may introduce severe structural
            artifacts or mask biologically relevant asymmetry. For this reason,
            symmetry above C1 is generally discouraged during exploratory
            ab initio reconstruction unless strong prior biological evidence
            exists.

            The protocol also supports real-space centering and windowing of
            reconstructed volumes. These operations improve numerical stability
            and help maintain compact, centered density distributions during
            optimization.

            Non-negativity constraints can additionally be enforced in real
            space. This is biologically meaningful because electron density
            values are expected to remain physically positive in most regions.

            Particle Alignment and Classification

            During optimization, particles are iteratively assigned orientations
            and, when multiple classes are used, class memberships. The protocol
            automatically stores these alignments and reconstructs representative
            volumes for each class.

            Internally, particle alignment information is converted into
            compatible metadata formats for downstream integration with
            Scipion and RELION-style workflows. Representative volumes are
            associated with their corresponding classes, preserving particle-to-class
            relationships throughout the processing chain.

            Outputs and Their Interpretation

            After completion, the protocol produces a set of reconstructed
            3D volumes together with their associated particle classifications.
            Each output class contains a representative volume and the aligned
            particles assigned to that structural state.

            Biologically, these volumes often represent candidate conformational
            states, compositional assemblies, or distinct structural populations.
            Interpretation should therefore consider both map quality and
            particle occupancy.

            The protocol also exports intermediate metadata files and converted
            STAR-format outputs to facilitate interoperability with external
            cryo-EM software ecosystems.

            Practical Recommendations

            In most routine cryo-EM workflows, it is advisable to begin with
            conservative parameters and a limited number of classes. Excessive
            heterogeneity or aggressive optimization settings can destabilize
            convergence and produce noisy maps.

            For difficult datasets, gradual resolution progression, moderate
            minibatch sizes, and symmetric noise modeling generally improve
            robustness. Enabling non-negativity and real-space centering is
            also recommended for most biological applications.

            Symmetry should only be enforced when strongly supported by prior
            structural knowledge. Incorrect symmetry assumptions are among the
            most common causes of biologically misleading reconstructions.

            Visual inspection of intermediate and final volumes remains essential.
            Successful reconstructions should exhibit coherent structural features,
            stable particle assignments, and biologically interpretable density
            distributions.

            Final Perspective

            Ab initio reconstruction represents one of the most critical stages
            in single-particle cryo-EM analysis because it establishes the
            structural foundation for all subsequent refinement and interpretation.
            Reliable initial models facilitate downstream classification,
            high-resolution refinement, and biological discovery.

            Careful parameter selection, realistic handling of structural
            heterogeneity, and biologically informed interpretation of the
            resulting maps are essential for obtaining meaningful structural
            insights from experimental cryo-EM datasets.
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
        self.info(pwutils.yellowStr("Ab Initial Model Generation Started..."))
        self.doRunAbinit()

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to Relion file
        """
        self.info(pwutils.yellowStr("Creating the output..."))
        self._initializeUtilsVariables()
        csOutputFolder = os.path.join(self.projectDir.get(),
                                      self.runAbinit.get())
        csFileName = "%s%s_final_particles.cs" % (getOutputPreffix(self.projectName.get()),
                                                  self.runAbinit.get())

        outputFolder = os.path.join(self._getExtraPath(), self.runAbinit.get())

        # Copy the CS output to extra folder
        copyFiles(csOutputFolder, outputFolder)

        csFile = os.path.join(outputFolder, csFileName)
        outputClassFn = self._getFileName('out_particles')
        argsList = [csFile, outputClassFn]
        convertCs2Star(argsList)

        # Create model files for 3D classification
        self._createModelFile()

        imgSet = self._getInputParticlesPointer()
        classes3D = self._createSetOfClasses3D(imgSet)
        self._fillClassesFromIter(classes3D, self._getFileName('out_particles'))

        self._defineOutputs(outputClasses=classes3D)
        self._defineSourceRelation(imgSet, classes3D)

        # create a SetOfVolumes and define its relations
        volumes = self._createSetOfVolumes()
        vol = None

        for class3D in classes3D:
            vol = class3D.getRepresentative()
            vol.setObjId(class3D.getObjId())
            volumes.append(vol)

        volumes.setSamplingRate(vol.getSamplingRate())

        self._defineOutputs(outputVolumes=volumes)
        self._defineSourceRelation(imgSet, volumes)

    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        validateMsgs = cryosparcValidate()
        if not validateMsgs:
            validateMsgs = gpusValidate(self.getGpuList(), checkSingleGPU=True)
            if not validateMsgs:
                particles = self._getInputParticles()
                if not particles.hasCTF():
                    validateMsgs.append(
                        "The Particles has not associated a "
                        "CTF model")
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

        table = emtable.Table(fileName=filename)

        for classNumber, row in enumerate(table.iterRows(filename)):
            index, fn = cryosparcToLocation(row.get(RELIONCOLUMNS.rlnReferenceImage.value))
            # Store info indexed by id, we need to store the row.clone() since
            # the same reference is used for iteration
            scaledFile = self._getScaledAveragesFile(fn, force=True)
            self._classesInfo[classNumber+1] = (index, scaledFile, row)

    def _fillClassesFromIter(self, clsSet, filename):
        """ Create the SetOfClasses3D """
        outImgsFn = 'particles@' + filename
        self._loadClassesInfo(self._getFileName('out_class'))
        clsSet.classifyItems(updateItemCallback=self._updateParticle,
                             updateClassCallback=self._updateClass,
                             itemDataIterator=emtable.Table.iterRows(outImgsFn),
                             raiseOnNextFailure=False,
                             cancelNextWhenAppendIsFalse=True)

    def _updateParticle(self, item, row):
        if row is not None and matchItemRow(item, row):
            if row.hasColumn(RELIONCOLUMNS.rlnClassNumber.value):
                item.setClassId(row.get(RELIONCOLUMNS.rlnClassNumber.value))
            else:
                item.setClassId(1)
            samplingRate = item.getSamplingRate()
            item.setTransform(rowToAlignment(row, ALIGN_PROJ, samplingRate))
        else:
            item._appendItem = False

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
                row = ("%s/%s/%s%s_class_%02d_final_volume.mrc\n"
                       % (self._getExtraPath(), self.runAbinit.get(),
                          getOutputPreffix(self.projectName.get()), self.runAbinit.get(), i))
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
            if not self.useQueueForSteps() and not self.useQueue():  # not using queue system
                gpusToUse = self.getGpuList()
            else:  # using queue system
                gpusToUse = False
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
                         "Please, go to cryoSPARC software for more "
                         "details.", self)
        clearIntermediateResults(self.projectName.get(), self.runAbinit.get(), wait=7)
