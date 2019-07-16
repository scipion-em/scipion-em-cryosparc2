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
import commands
import pyworkflow.em.metadata as md
from pyworkflow.em import ALIGN_PROJ
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        LabelParam, IntParam, Positive,
                                        EnumParam, StringParam,
                                        BooleanParam, PathParam,
                                        LEVEL_ADVANCED)
from pyworkflow.em.data import String
from pyworkflow.em.protocol import ProtInitialVolume, ProtClassify3D
from cryosparc2.convert import *
from cryosparc2.utils import *
from cryosparc2.constants import *

relionPlugin = pwutils.importFromPlugin("relion.convert", doRaise=True)


class ProtCryoSparcInitialModel(ProtInitialVolume, ProtClassify3D):
    """    
    Generate a 3D initial model _de novo_ from 2D particles using
    CryoSparc Stochastic Gradient Descent (SGD) algorithm.
    """
    _label = 'initial model'
    # --------------------------- DEFINE param functions ----------------------
    def _defineFileNames(self):
        """ Centralize how files are called within the protocol. """
        myDict = {
                  'input_particles': self._getPath('input_particles.star'),
                  'out_particles': self._getPath() + '/output_particle.star',
                  'out_class': self._getPath() + '/output_class.star'
                  }
        self._updateFilenamesDict(myDict)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      label="Input particles", important=True,
                      help='Select the input images from the project.')
        form.addParallelSection(threads=1, mpi=1)

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

        addSymmetryParam(form)

        # form.addParam('abinit_symmetry', StringParam, default="C1",
        #               label='Symmetry:',
        #               help='Symmetry enforced (C, D, I, O, T). Eg. C1, D7, C4 '
        #                    'etc. Enforcing symmetry above C1 is not '
        #                    'recommended for ab-initio reconstruction')

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

        form.addSection(label='Compute settings')

        form.addParam('compute_use_ssd', BooleanParam, default=True,
                      label='Cache particle images on SSD:')

    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._defineFileNames()
        self._defineParamsName()
        self._initializeCryosparcProject()
        self._insertFunctionStep("convertInputStep")
        self._insertFunctionStep('processStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions ------------------------------
    def convertInputStep(self):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file. 
        """
        print(pwutils.greenStr("Importing Particles..."))
        imgSet = self._getInputParticles()
        relionPlugin.writeSetOfParticles(imgSet,
                                         self._getFileName('input_particles'),
                                         outputDir=self._getExtraPath(),
                                         fillMagnification=True,
                                         fillRandomSubset=True)
        self._importParticles()
        while getJobStatus(self.projectName.get(), self.importedParticles.get()) not in STOP_STATUSES:
            waitJob(self.projectName.get(), self.importedParticles.get())

    def processStep(self):

        print(pwutils.greenStr("Ab Initial Model Generation Started..."))
        self.runAbinit = String(self.doRunAbinit()[-1].split()[-1])
        self.currenJob.set(self.runAbinit)
        self._store(self)

        while getJobStatus(self.projectName.get(), self.runAbinit.get()) not in STOP_STATUSES:
            waitJob(self.projectName.get(), self.runAbinit.get())

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to Relion file
        """
        self._initializeUtilsVariables()
        print (pwutils.greenStr("Creating the output..."))

        csFile = os.path.join(self.projectPath, self.projectName.get(),
                              self.runAbinit.get(), ("cryosparc_" +
                                                     self.projectName.get() +
                                                     "_" + self.runAbinit.get() +
                                                     "_final_particles.cs"))

        outputClassFn = self._getFileName('out_particles')
        argsList = [csFile, outputClassFn]

        parser = defineArgs()
        args = parser.parse_args(argsList)
        convertCs2Star(args)

        # Link the folder on SSD to scipion directory
        os.system("ln -s " + self.projectPath + "/" + self.projectName.get() +
                  '/' + self.runAbinit.get() + " " + self._getExtraPath())

        # Create model files for 3D classiffication
        with open(self._getFileName('out_class'), 'w') as output_file:
            output_file.write('\n')
            output_file.write('data_images')
            output_file.write('\n\n')
            output_file.write('loop_')
            output_file.write('\n')
            output_file.write('_rlnReferenceImage')
            output_file.write('\n')
            for i in range(int(self.abinit_K.get())):
                output_file.write("%02d"%(i+1)+"@"+self._getExtraPath() + "/" + self.runAbinit.get() + "/cryosparc_" +\
                self.projectName.get() + "_"+self.runAbinit.get() + "_class_%02d"%i +
                                  "_final_volume.mrc\n")

        imgSet = self._getInputParticles()
        classes3D = self._createSetOfClasses3D(imgSet)
        self._fillClassesFromIter(classes3D, self._getFileName('out_particles'))

        self._defineOutputs(outputClasses=classes3D)
        self._defineSourceRelation(self.inputParticles.get(), classes3D)

        # create a SetOfVolumes and define its relations
        volumes = self._createSetOfVolumes()
        volumes.setSamplingRate(imgSet.getSamplingRate())

        for class3D in classes3D:
            vol = class3D.getRepresentative()
            vol.setObjId(class3D.getObjId())
            volumes.append(vol)

        self._defineOutputs(outputVolumes=volumes)
        self._defineSourceRelation(self.inputParticles.get(), volumes)

    def setAborted(self):
        """ Set the status to aborted and updated the endTime. """
        ProtInitialVolume.setAborted(self)
        killJob(str(self.projectName.get()), str(self.currenJob.get()))

    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        validateMsgs = cryosparcExist()
        if not validateMsgs:
            validateMsgs = isCryosparcRunning()
            if not validateMsgs:
                particles = self._getInputParticles()
                if not particles.hasCTF():
                    validateMsgs.append("The Particles has not associated a "
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

    def _getInputParticles(self):
        return self.inputParticles.get()

    def _loadClassesInfo(self, filename):
        """ Read some information about the produced CryoSparc Classes
        from the star file.
        """
        self._classesInfo = {}  # store classes info, indexed by class id

        modelStar = md.MetaData(filename)

        for classNumber, row in enumerate(md.iterRows(modelStar)):
            index, fn = relionPlugin.relionToLocation(row.getValue('rlnReferenceImage'))
            # Store info indexed by id, we need to store the row.clone() since
            # the same reference is used for iteration
            self._classesInfo[classNumber+1] = (index, fn, row.clone())

    def _fillClassesFromIter(self, clsSet, filename):
        """ Create the SetOfClasses3D """
        self._loadClassesInfo(self._getFileName('out_class'))
        outImgsFn = self._getFileName('out_class')
        clsSet.classifyItems(updateItemCallback=self._updateParticle,
                             updateClassCallback=self._updateClass,
                             itemDataIterator=md.iterRows(filename,
                                                          sortByLabel=md.RLN_IMAGE_ID))

    def _updateParticle(self, item, row):
        item.setClassId(row.getValue(md.RLN_PARTICLE_CLASS))
        item.setTransform(relionPlugin.rowToAlignment(row, ALIGN_PROJ))

    def _updateClass(self, item):
        classId = item.getObjId()
        if classId in self._classesInfo:
            index, fn, row = self._classesInfo[classId]
            fn += ":mrc"
            item.setAlignmentProj()
            item.getRepresentative().setLocation(index, fn)

    def _initializeUtilsVariables(self):
        """
        Initialize all utils cryoSPARC variables
        """
        self._program = getCryosparcProgram()
        self._user = getCryosparcUser()
        self._ssd = getCryosparcSSD()

        # Create a cryoSPARC project dir
        self.projectDirName = suffix + self.getProject().getShortName()
        self.projectPath = pwutils.join(self._ssd, self.projectDirName)
        self.projectDir = createProjectDir(self.projectPath)

    def _initializeCryosparcProject(self):
        """
        Initialize the cryoSPARC project and workspace
        """
        self._initializeUtilsVariables()
        # create empty project or load an exists one
        folderPaths = getProjectPath(self.projectPath)
        if not folderPaths:
            self.a = createEmptyProject(self.projectPath, self.projectDirName)
            self.projectName = self.a[-1].split()[-1]
        else:
            self.projectName = str(folderPaths[0])

        self.projectName = String(self.projectName)
        self._store(self)

        # create empty workspace
        self.b = createEmptyWorkSpace(self.projectName, self.getRunName(),
                                      self.getObjComment())
        self.workSpaceName = String(self.b[-1].split()[-1])
        self._store(self)

    def _importParticles(self):

        print("Importing Particles")

        # import_particles_star
        self.c = self.doImportParticlesStar()

        self.importedParticles = String(self.c[-1].split()[-1])
        self._store(self)

        self.currenJob = String(self.importedParticles.get())
        self._store(self)

        self.par = String(self.importedParticles.get() + '.imported_particles')

    def doImportParticlesStar(self):
        """
        do_import_particles_star(puid, wuid, uuid, abs_star_path,
                                 abs_blob_path=None, psize_A=None)
        returns the new uid of the job that was created
        """
        cmd = """ 'do_import_particles_star("%s","%s", "%s", "%s", "%s", "%s")'"""
        import_particles_cmd = (self._program + cmd % (
            self.projectName.get(), self.workSpaceName.get(),
            self._user,
            os.path.join(os.getcwd(),
            self._getFileName('input_particles')),
            os.path.join(os.getcwd(),
            self._getExtraPath()),
            str(self._getInputParticles().getSamplingRate())
        ))
        print(pwutils.greenStr(import_particles_cmd))
        return commands.getstatusoutput(import_particles_cmd)

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

    def doRunAbinit(self):
        """self._program + "  \'do_run_abinit(\"" + self.projectName +
        "\", \"" + self.workSpaceName + "\", \"\'" + self._user + "\'\", \""
        + self.par + "\",\"\'1\'\")\'")
        """
        className = "homo_abinit"
        input_group_conect = {"particles": str(self.par)}
        # {'particles' : 'JXX.imported_particles' }
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

        return doJob(className, self.projectName.get(), self.workSpaceName.get(),
                     str(params).replace('\'', '"'),
                     str(input_group_conect).replace('\'', '"'))

