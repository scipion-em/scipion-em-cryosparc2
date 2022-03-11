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

import emtable
from pkg_resources import parse_version

import pyworkflow.utils as pwutils
from pwem.objects import VolumeMask
from pyworkflow import BETA
from pyworkflow.object import String
from pyworkflow.protocol.params import (FloatParam, LEVEL_ADVANCED,
                                        PointerParam, MultiPointerParam,
                                        CsvList, Positive, IntParam,
                                        BooleanParam, EnumParam)

from .protocol_base import ProtCryosparcBase
from ..convert import (convertBinaryVol, defineArgs, convertCs2Star,
                       rowToAlignment, ALIGN_PROJ, cryosparcToLocation)
from ..utils import (addComputeSectionParams, doImportVolumes,
                     get_job_streamlog, calculateNewSamplingRate,
                     cryosparcValidate, gpusValidate, enqueueJob,
                     waitForCryosparc, clearIntermediateResults, fixVolume,
                     copyFiles, getCryosparcVersion)
from ..constants import *


class ProtCryoSparcNew3DClassification(ProtCryosparcBase):
    """
    3D Classification (BETA) is a new job in cryoSPARC v3.3+ to analyze discrete
    heterogeneity in single particle cryo-EM datasets. This job currently
    implements a version of 3D classification without alignment — a
    classification routine that can complement the existing Heterogeneous
    Refinement job in finding new discrete classes of data.
    """
    _label = '3D Classification'
    _className = "class_3D"
    _devStatus = BETA
    _protCompatibility = [V3_3_1]

    def _initialize(self):
        self._defineFileNames()

    def _defineFileNames(self):
        """ Centralize how files are called. """
        myDict = {
            'input_particles': self._getTmpPath('input_particles.star'),
            'out_particles': self._getExtraPath('output_particle.star'),
            'stream_log': self._getPath() + '/stream.log',
            'out_class': self._getExtraPath() + '/output_class.star'
        }
        self._updateFilenamesDict(myDict)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      label="Input particles", important=True,
                      help='Particle stacks to use. Multiple stacks will be '
                           'concatenated.')
        form.addParam('refVolumes', MultiPointerParam,
                      pointerClass='Volume',
                      allowsNull=True,
                      label="Initial volumes",
                      help='Multiple initial volume raw data')
        form.addParam('refMask', PointerParam, pointerClass='VolumeMask',
                      label='Static mask',
                      allowsNull=True,
                      help="Mask for focussed classification. If none supplied,"
                           " a mask will be created from initial structures "
                           "(bootstrapped or input).")

        # --------------[3D Classification]---------------------------
        form.addSection(label='3D Classification (without alignment)')

        form.addParam('class3D_N_K', IntParam, default=2,
                      label="Number of classes",
                      validators=[Positive],
                      help='Number of classes. This value must be grater than 2')

        form.addParam('class3D_target_res', IntParam, default=4,
                      label="Target resolution (A)",
                      validators=[Positive],
                      help='Set based off highest resolution at which variation '
                           'occurs. This parameter determines the reconstruction '
                           'box size during classification (which also accounts '
                           'for input particle extent and optimal sizes for FFT ops).')

        form.addParam('class3D_num_oem_epochs', IntParam, default=5,
                      label="Number of O-EM epochs",
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      help='Number of `epochs` (traversals through the entire '
                           'dataset) to perform during online expectation '
                           'maximization (O-EM).')

        form.addParam('class3D_num_oem_epochs', IntParam, default=5,
                      label="Number of O-EM epochs",
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      help='Number of `epochs` (traversals through the entire '
                           'dataset) to perform during online expectation '
                           'maximization (O-EM).')

        form.addParam('class3D_num_final_full_iters', IntParam, default=2,
                      label="Number of final full iterations",
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      help='Number of final `full batch` iterations through '
                           'the entire dataset. Generally 2 iterations is enough.')

        form.addParam('class3D_batch_size_per_class', IntParam, default=1000,
                      label="Batch size per class",
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      help='Number of images per class used in each batch (one '
                           'batch per iteration of online-EM). Larger values '
                           'slow the algorithm down but can provide better'
                           ' classification results.')

        form.addParam('class3D_init_res', IntParam, default=20,
                      label="Init structure lowpass resolution (A)",
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      help='Init structure lowpass resolution (A)')

        form.addParam('class3D_init_mode', EnumParam, default=0,
                      label='Initialization mode',
                      choices=['simple', 'PCA', 'input'],
                      help="simple mode: reconstruct initial volumes from "
                           "randomly selected particle subsets. PCA mode: "
                           "reconstruct initial volumes from cluster-average "
                           "maps computed after applying PCA to a large set of "
                           "reconstructions. input mode: use input volumes "
                           "(Warning: volumes should *not* be identical!)")

        form.addParam('class3D_N_ik', IntParam, default=1000,
                      label="PCA/simple: particles per reconstruction",
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      help='Number of particles used to compute each density '
                           'initialization. Recommended settings -- '
                           'simple mode: 1000-5000; PCA mode: 100-500.')

        form.addParam('class3D_PCA_num_reconstructions', IntParam, default=100,
                      label="PCA: number of reconstructions",
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      help='Number of reconstructions which will be used in '
                           'PCA.')

        form.addParam('class3D_PCA_num_components', IntParam, default=5,
                      label="PCA: number of components",
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      help='PCA: number of components')

        form.addParam('class3D_filter_hp_res', IntParam, default=None,
                      allowsNull=True,
                      label="Highpass resolution (A)",
                      expertLevel=LEVEL_ADVANCED,
                      help='Resolution below which variability is ignored')

        form.addParam('class3D_filter_hp_order', IntParam, default=8,
                      label="Highpass order",
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      help='Order of filter')

        form.addParam('class3D_mask_thresh_factor', FloatParam, default=0.2,
                      label="Auto mask threshold (0-1)",
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      help='(If no input mask is provided) Level set threshold '
                           'for selecting regions that are included in the mask '
                           'made from initial structures.')

        form.addParam('class3D_mask_near_ang', FloatParam, default=10.0,
                      label="Auto mask near (A)",
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      help='(If no input mask is provided) Controls extent to'
                           ' which mask is expanded. At the near distance, '
                           'the mask value is 1.0 (in A)')

        form.addParam('class3D_mask_far_ang', FloatParam, default=20.0,
                      label="Auto mask far (A)",
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      help='(If no input mask is provided) Controls extent to '
                           'which mask is expanded. At the far distance the '
                           'mask value becomes 0.0 (in A)')

        form.addParam('class3D_mask_use_abs', BooleanParam, default=False,
                      label="Auto mask use absolute value",
                      expertLevel=LEVEL_ADVANCED,
                      help='(If no input mask is provided) Include negative '
                           'regions if they are more negative than the '
                           'threshold')

        form.addParam('class3D_class_anneal_beta', FloatParam, default=0.5,
                      label="Class similarity",
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      help='Expected similarity of structures from different '
                           'classes. A number between 0 and 1. 0 means classes '
                           'are independent, 1 means classes are very similar')

        form.addParam('class3D_class_anneal_beta_tune', BooleanParam,
                      default=True,
                      label="Auto-tune initial class similarity",
                      expertLevel=LEVEL_ADVANCED,
                      help='Tune initial class similarity based on a target '
                           'class ESS (Effective Sample Size). An ESS close to '
                           'the number of classes at the outset helps to '
                           'ensure the classification does not prematurely '
                           'converge.')

        form.addParam('class3D_class_anneal_beta_tune_essf', FloatParam,
                      default=0.75,
                      label="Target class ESS fraction",
                      expertLevel=LEVEL_ADVANCED,
                      help='Target ESS fraction (of total classes) for '
                           'similarity tuning. A number between 0 and 1.')

        form.addParam('class3D_class_anneal_end', IntParam,
                      default=10,
                      validators=[Positive],
                      label="Class similarity anneal end iter",
                      expertLevel=LEVEL_ADVANCED,
                      help='Class similarity will be annealed to 0 throughout '
                           'the optimization. This is the end iteration when '
                           'class similarity should equal 0.')

        form.addParam('class3D_bootstrap_noise', BooleanParam,
                      default=True,
                      label="Estimate noise based on initial maps",
                      expertLevel=LEVEL_ADVANCED,
                      help='Estimate the initial noise distribution from '
                           'initial maps. If false, noise will be initially '
                           'set to a constant value.')

        form.addParam('class3D_use_anisomag', BooleanParam,
                      default=False,
                      label="Correct Anisotropic Magnification",
                      expertLevel=LEVEL_ADVANCED,
                      help='Whether or not to correct for anisotropic '
                           'magnification. The anisomag parameters must '
                           'already have been fit but a global CTF refinement '
                           'protocol or within a previous homogeneous '
                           'refinement.')

        form.addParam('class3D_zip_volumes', BooleanParam,
                      default=True,
                      label="Compress volumes into series",
                      expertLevel=LEVEL_ADVANCED,
                      help='Collect and compress reconstructed volumes into a'
                           ' single zip file output. Useful when running '
                           'classification jobs with a large number of classes.')

        form.addParam('class3D_plot_iters', IntParam,
                      default=10,
                      validators=[Positive],
                      label="Output plots every __ iters",
                      expertLevel=LEVEL_ADVANCED,
                      help='Plot at this cadence (units of O-EM iterations -- '
                           'not epochs).')

        form.addParam('class3D_full_iter_out', BooleanParam,
                      default=False,
                      label="Output data after every full iter",
                      expertLevel=LEVEL_ADVANCED,
                      help='Turn on to output volumes/particles after every '
                           'full iteration.')

        form.addParam('class3D_dist_plots', BooleanParam,
                      default=True,
                      label="Show viewing orientation distribution plots",
                      expertLevel=LEVEL_ADVANCED,
                      help='Whether or not to produce and show viewing '
                           'direction distribution. Turn off to reduce '
                           'runtime and decrease the number of plots shown in '
                           'the stream log.')

        form.addParam('class3D_num_particles', IntParam,
                      default=None,
                      allowsNull=True,
                      label="Number of particles to classify",
                      expertLevel=LEVEL_ADVANCED,
                      help='None means use the entire input particle stack -- '
                           'otherwise use the first N particles.')

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
        if self.refVolumes is not None:
            self.volumes = [vol + self.outputVolumeSuffix for vol in
                            self.importVolumes]
        else:
            self.volumes = None
        print(pwutils.yellowStr("3D Classification started..."), flush=True)
        self.doNew3DClasification()

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to Relion file
        """
        self._initializeUtilsVariables()
        print(pwutils.yellowStr("Creating the output..."), flush=True)

        csOutputFolder = os.path.join(self.projectPath, self.projectName.get(),
                                      self.run3dClassification.get())
        itera = self.findLastIteration(self.run3dClassification.get())

        csParticlesName = "cryosparc_%s_%s_00%s_particles.cs" % (
                                                 self.projectName.get(),
                                                 self.run3dClassification.get(),
                                                 itera)
        # Copy the CS output particles to extra folder
        copyFiles(csOutputFolder, self._getExtraPath(), files=[csParticlesName])

        csFile = os.path.join(self._getExtraPath(), csParticlesName)

        outputStarFn = self._getFileName('out_particles')
        argsList = [csFile, outputStarFn]

        parser = defineArgs()
        args = parser.parse_args(argsList)
        convertCs2Star(args)

        self._createModelFile(csOutputFolder, itera)

        # Create 3D classes
        imgSet = self._getInputParticles()
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
        self._defineSourceRelation(self.inputParticles.get(), volumes)

        # Create a 3D mask
        volMask = VolumeMask()
        maskFileName = "cryosparc_%s_%s__mask.mrc" % (
            self.projectName.get(),
            self.run3dClassification.get())
        copyFiles(csOutputFolder, self._getExtraPath(), files=[maskFileName])
        maskFilePath = os.path.join(self._getExtraPath(), maskFileName)
        volMask.setFileName(maskFilePath)
        sr = self._getInputParticles().getSamplingRate()
        volMask.setSamplingRate(sr)

        self._defineOutputs(outputMask=volMask)
        self._defineSourceRelation(self.inputParticles.get(), self.outputMask)

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
            self._classesInfo[classNumber + 1] = (index, scaledFile, row)

    def _fillClassesFromIter(self, clsSet, filename):
        """ Create the SetOfClasses3D """
        xmpMd = 'micrographs@' + filename
        self._loadClassesInfo(self._getFileName('out_class'))
        clsSet.classifyItems(updateItemCallback=self._updateParticle,
                             updateClassCallback=self._updateClass,
                             itemDataIterator=emtable.Table.iterRows(xmpMd))

    def _updateParticle(self, item, row):
        if row.get(RELIONCOLUMNS.rlnAnglePsi.value):
            item.setClassId(row.get(RELIONCOLUMNS.rlnClassNumber.value))
            item.setTransform(rowToAlignment(row, ALIGN_PROJ))
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

    def _createModelFile(self, csOutputFolder, itera):
        # Create model files for 3D classification
        with open(self._getFileName('out_class'), 'w') as output_file:
            output_file.write('\n')
            output_file.write('data_images')
            output_file.write('\n\n')
            output_file.write('loop_')
            output_file.write('\n')
            output_file.write('_rlnReferenceImage')
            output_file.write('\n')
            numOfClass = self.class3D_N_K.get()
            for i in range(numOfClass):
                csVolName = ("cryosparc_%s_%s_class_%02d_00%s_volume.mrc" %
                             (self.projectName.get(),
                              self.run3dClassification.get(), i, itera))

                copyFiles(csOutputFolder, self._getExtraPath(),
                          files=[csVolName])

                row = ("%s/cryosparc_%s_%s_class_%02d_00%s_volume.mrc\n" %
                       (self._getExtraPath(), self.projectName.get(),
                        self.run3dClassification.get(), i, itera))
                output_file.write(row)

    def findLastIteration(self, jobName):
        import ast
        get_job_streamlog(self.projectName.get(),
                          jobName,
                          self._getFileName('stream_log'))

        # Get the metadata information from stream.log
        with open(self._getFileName('stream_log')) as f:
            data = f.readlines()

        x = ast.literal_eval(data[0])

        # Find the ID of last iteration and the map resolution
        for y in x:
            if 'text' in y:
                z = str(y['text'])
                if z.startswith('Batch Class Distribution (Iteration:'):
                    it = z.split(': ')[1].split(')')[0]
                    if int(it) > 99:
                        itera = it
                    else:
                        itera = '0' + it if int(it) > 10 else '00' + it

        return itera

    # --------------------------- INFO functions -------------------------------

    def _validate(self):
        validateMsgs = cryosparcValidate()
        if not validateMsgs:
            csVersion = getCryosparcVersion()
            if [version for version in self._protCompatibility
                if parse_version(version) >= parse_version(csVersion)]:
                validateMsgs = gpusValidate(self.getGpuList(),
                                            checkSingleGPU=True)
                if not validateMsgs:
                    particles = self._getInputParticles()
                    if not particles.hasCTF():
                        validateMsgs.append("The Particles has not associated a "
                                            "CTF model")
                    if not validateMsgs and not particles.hasAlignmentProj():
                        validateMsgs.append("The Particles has not "
                                            "alignment")

                    inputVolumes = self._getInputVolume()
                    if not validateMsgs:
                        if inputVolumes is not None and inputVolumes:
                            if self.class3D_init_mode.get() != 2:
                                validateMsgs.append('Input volumes detected, please set initialization mode to `input` or clear volume inputs.')
                            if len(inputVolumes) != self.class3D_N_K.get():
                                validateMsgs.append('No. of input volumes must equal no. of classes')
                        else:
                            if self.class3D_init_mode.get() != 0:
                                validateMsgs.append('Please connect input volumes or change initialization mode to `simple`')
                            elif len(particles) < 1000:
                                validateMsgs.append('Not Enough Images! The set of particles must contain at least 1000 images')
                        if not validateMsgs:
                            if self.class3D_N_K.get() < 2:
                                validateMsgs.append('The number of classes must be grater than 2')

            else:
                validateMsgs.append("The protocol is not compatible with the "
                                    "cryoSPARC version %s" % csVersion)

        return validateMsgs

    # --------------------------- UTILS functions ------------------------------

    def _getInputVolume(self):
        if self.hasAttribute('refVolumes') and self.refVolumes is not None:
            self.vols = []
            for volume in self.refVolumes:
                vol = volume.get()
                self.vols.append(vol)
            return self.vols
        return None

    def _importVolume(self):
        self.importVolumes = CsvList()
        self._initializeVolumeSuffix()
        for vol in self.vols:
            self.vol_fn = os.path.join(os.getcwd(),
                                       convertBinaryVol(
                                           vol,
                                           self._getTmpPath()))
            self.importVolume = doImportVolumes(self, self.vol_fn, vol,
                                                'map',
                                                'Importing volume...')
            self.importVolumes.append(self.importVolume.get())
            self.currenJob.set(self.importVolume.get())

    def _defineParamsName(self):
        self._paramsName = ['class3D_N_K', 'class3D_target_res',
                            'class3D_num_oem_epochs',
                            'class3D_num_oem_epochs',
                            'class3D_num_final_full_iters',
                            'class3D_batch_size_per_class', 'class3D_init_res',
                            'class3D_init_mode', 'class3D_N_ik',
                            'class3D_PCA_num_reconstructions',
                            'class3D_PCA_num_components',
                            'class3D_filter_hp_res', 'class3D_filter_hp_order',
                            'class3D_mask_thresh_factor',
                            'class3D_mask_near_ang',
                            'class3D_mask_far_ang',
                            'class3D_mask_use_abs',
                            'class3D_class_anneal_beta',
                            'class3D_class_anneal_beta_tune',
                            'class3D_class_anneal_beta_tune_essf',
                            'class3D_class_anneal_end',
                            'class3D_use_anisomag',
                            'class3D_zip_volumes',
                            'class3D_plot_iters',
                            'class3D_full_iter_out',
                            'class3D_dist_plots',
                            'class3D_num_particles',
                            'compute_use_ssd']

        self.lane = str(self.getAttributeValue('compute_lane'))


    def _summary(self):
        summary = []
        if (not hasattr(self, 'outputVolumes') or
                not hasattr(self, 'outputClasses')):
            summary.append("Output objects not ready yet.")
        else:
            summary.append("Input Particles: %s" %
                           self.getObjectTag('inputParticles'))
            summary.append("Initial volumes: %s" %
                           self.getObjectTag('refVolumes'))

            summary.append("------------------------------------------")
            summary.append("Output volumes %s" %
                           self.getObjectTag('outputVolumes'))
            summary.append("Output classes %s" %
                           self.getObjectTag('outputClasses'))

        return summary

    def doNew3DClasification(self):
        """
        """
        group_connect = {}
        input_group_connect = {"particles": self.particles.get()}
        if self.volumes is not None:
            group_connect["volume"] = self.volumes

        if self.mask.get() is not None:
            group_connect["mask"] = [self.mask]
        params = {}


        for paramName in self._paramsName:
            if (paramName != 'class3D_filter_hp_res' and
                    paramName != 'class3D_num_particles' and paramName != 'class3D_init_mode'):
                params[str(paramName)] = str(self.getAttributeValue(paramName))

            elif paramName == 'class3D_init_mode':
                params[str(paramName)] = str(CLASS_3D_INIT_MODE[self.class3D_init_mode.get()])

            elif paramName == 'class3D_filter_hp_res' and self.class3D_filter_hp_res.get() is not None:
                params[str(paramName)] = str(self.class3D_filter_hp_res.get())

            elif paramName == 'class3D_num_particles' and self.class3D_num_particles.get() is not None:
                params[str(paramName)] = str(self.class3D_num_particles.get())

        # Determinate the GPUs to use (in dependence of
        # the cryosparc version)
        try:
            gpusToUse = self.getGpuList()
        except Exception:
            gpusToUse = False

        run3dClassificationJob = enqueueJob(self._className,
                                              self.projectName.get(),
                                              self.workSpaceName.get(),
                                              str(params).replace('\'', '"'),
                                              str(input_group_connect).replace('\'', '"'),
                                              self.lane, gpusToUse,
                                              group_connect=group_connect)

        self.run3dClassification = String(run3dClassificationJob.get())
        self.currenJob.set(run3dClassificationJob.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(), self.run3dClassification.get(),
                         "An error occurred in the 3D Classification process. "
                         "Please, go to cryosPARC software for more "
                         "details.")
        clearIntermediateResults(self.projectName.get(), self.run3dClassification.get())
