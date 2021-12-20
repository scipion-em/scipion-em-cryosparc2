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
from pyworkflow import BETA
from pyworkflow.object import String
from pyworkflow.protocol.params import (FloatParam, LEVEL_ADVANCED,
                                        PointerParam, MultiPointerParam,
                                        CsvList, Positive, IntParam,
                                        BooleanParam, StringParam, EnumParam)

from .protocol_base import ProtCryosparcBase
from ..convert import (convertBinaryVol, defineArgs, convertCs2Star,
                       rowToAlignment, ALIGN_PROJ, cryosparcToLocation)
from ..utils import (addSymmetryParam, addComputeSectionParams, doImportVolumes,
                     get_job_streamlog, calculateNewSamplingRate,
                     cryosparcValidate, gpusValidate, getSymmetry, enqueueJob,
                     waitForCryosparc, clearIntermediateResults, fixVolume,
                     copyFiles)
from ..constants import *


class ProtCryoSparc3DClassification(ProtCryosparcBase):
    """
    Heterogeneous Refinement simultaneously classifies particles and refines
    structures from n initial structures, usually obtained following an
    Ab-Initio Reconstruction. This facilitates the ability to look for small
    differences between structures which may not be obvious at low resolutions,
    and also to re-classify particles to aid in sorting.
    """
    _label = '3D classification'
    _className = "hetero_refine"
    _devStatus = BETA

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
                      validators=[Positive],
                      help='Select the input images from the project.')
        form.addParam('refVolumes', MultiPointerParam,
                      pointerClass='Volume',
                      important=True,
                      label="Initial volumes",
                      help='Multiple initial volumes to refine and classify. '
                           'The same input volume can be connected multiple '
                           'times.')

        # --------------[Heterogeneous Refinement]---------------------------
        form.addSection(label='Heterogeneous Refinement')
        form.addParam('multirefine_N', IntParam, default=128,
                      label="Refinement box size (Voxels)",
                      help='Box size of each volume during refinement. '
                           'Particles will automatically be downsampled '
                           '(Fourier cropped) to this box size on the fly. '
                           'Keep this as small as possible to limit GPU '
                           'memory usage.')

        addSymmetryParam(form, help="Symmetry String (C, D, I, O, T). E.g. "
                                    "C1, D7, C4, etc. Symmetry is applied to "
                                    "all classes.")

        form.addParam('multirefine_sharp_bfactor', IntParam, default=-100,
                      label="Plotting bfactor",
                      help='B-Factor to apply to the structures before '
                           'plotting, to enhance medium/high resolution '
                           'detail in plots. Outputs are not affected by this '
                           'parameter.')

        form.addParam('multirefine_force_hard_class', BooleanParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=False,
                      label="Force hard classification",
                      help='Force hard classification, so that each particle '
                           'is only assigned to one class at every iteration, '
                           'rather than having partial assignment to all '
                           'classes.')

        form.addParam('multirefine_batch_size_per_class', IntParam, default=1000,
                      expertLevel=LEVEL_ADVANCED,
                      label="Batch size per class",
                      help='Number of images per class used in each batch (one '
                           'batch per iteration of online-EM). Larger values '
                           'slow the algorithm down but can provide better '
                           'classification results.')

        form.addParam('multirefine_update_rule', StringParam,
                      default="online_em",
                      expertLevel=LEVEL_ADVANCED,
                      label="Optimization method",
                      help='Optimization method to use.')

        form.addParam('multirefine_online_em_lr_rand', FloatParam,
                      default=0.2,
                      expertLevel=LEVEL_ADVANCED,
                      label="O-EM learning rate during randomization",
                      help='Not recommended to change.')

        form.addParam('multirefine_online_em_lr_init', FloatParam,
                      default=0.1,
                      expertLevel=LEVEL_ADVANCED,
                      label="O-EM learning rate init",
                      help='Not recommended to change.')

        form.addParam('multirefine_online_em_lr_hl', IntParam,
                      default=50,
                      expertLevel=LEVEL_ADVANCED,
                      label="O-EM learning rate halflife (iters)",
                      help='Not recommended to change.')

        form.addParam('multirefine_halfmap_decay', FloatParam,
                      default=0.9,
                      expertLevel=LEVEL_ADVANCED,
                      label="Halfmap decay constant",
                      help='Not recommended to change.')

        form.addParam('multirefine_res_init', FloatParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=20,
                      label="Initial resolution (A)",
                      help='Initial low-pass resolution applied to input '
                           'volumes before classification or reconstruction.')

        form.addParam('multirefine_bp_res_factor', FloatParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=1.5,
                      label="Backprojection resolution factor",
                      help='Backproject at this multiple of the best resolution '
                           'amongst classes. Not recommended to change.')

        form.addParam('multirefine_use_max_fsc', BooleanParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=True,
                      label="Use max FSC over classes for filtering",
                      help='Use the maximum FSC across classes to filter all '
                           'classes. This prevents smaller classes from being '
                           'over-filtered during reconstruction.')

        form.addParam('multirefine_assignment_conv_eps', FloatParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=0.05,
                      label="Assignment convergence criteria",
                      help='Fraction of the batch that is allowed to have '
                           'changed classes in the past iteration to be '
                           'considered converged.')

        form.addParam('multirefine_assignment_conv_eps', IntParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=1,
                      label="Resolution convergence criteria",
                      help='Maximum change in resolution (in Fourier shells) '
                           'between iterations to be considered converged.')

        form.addParam('multirefine_num_rand_assign_iters', IntParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=5,
                      label="Number of initial random assignment iterations",
                      help='Number of iterations to perform initially with '
                           'random assignments to break symmetry of multiple '
                           'identical initial references.')

        form.addParam('multirefine_num_final_full_iters', IntParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=2,
                      label="Number of final full iterations",
                      help='Number of final full iterations through the entire '
                           'dataset. Generally 2 iterations is enough.')

        form.addParam('multirefine_noise_model', EnumParam,
                      expertLevel=LEVEL_ADVANCED,
                      choices=['symmetric', 'white', 'coloured'],
                      default=0,
                      label='Noise model',
                      help='Noise model to use. Valid options are white, '
                           'coloured or symmetric. Symmetric is the default, '
                           'meaning coloured with radial symmetry. '
                           'Not recommended to change.')

        form.addParam('multirefine_noise_init_sigmascale', IntParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=3,
                      label="Noise initial sigma-scale",
                      help='Scale factor initially applied to the base noise '
                           'estimate. Not recommended to change.')

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
    def _getInputVolume(self):
        if self.hasAttribute('refVolumes'):
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
            self.importVolume = doImportVolumes(self, self.vol_fn, vol, 'map',
                                                'Importing volume...')
            self.importVolumes.append(self.importVolume.get())
            self.currenJob.set(self.importVolume.get())

    def processStep(self):
        self.volumes = [vol + self.outputVolumeSuffix for vol in self.importVolumes]
        print(pwutils.yellowStr("3D Classification started..."), flush=True)
        self.do3DClasification()

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to Relion file
        """
        self._initializeUtilsVariables()
        print(pwutils.yellowStr("Creating the output..."), flush=True)

        csOutputFolder = os.path.join(self.projectPath, self.projectName.get(),
                              self.run3dClassification.get())
        itera = self.findLastIteration(self.run3dClassification.get())

        csParticlesName = "cryosparc_%s_%s_000%s_particles.cs" % (self.projectName.get(),
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

    # --------------------------- UTILS functions ---------------------------
    def _loadClassesInfo(self, filename):
        """ Read some information about the produced CryoSparc Classes
        from the star file.
        """
        self._classesInfo = {}  # store classes info, indexed by class id

        modelStar = md.MetaData(filename)

        for classNumber, row in enumerate(md.iterRows(modelStar)):
            index, fn = cryosparcToLocation(
                row.getValue('rlnReferenceImage'))
            # Store info indexed by id, we need to store the row.clone() since
            # the same reference is used for iteration
            scaledFile = self._getScaledAveragesFile(fn, force=True)
            self._classesInfo[classNumber + 1] = (index, scaledFile, row.clone())

    def _fillClassesFromIter(self, clsSet, filename):
        """ Create the SetOfClasses3D """
        xmpMd = 'micrographs@' + filename
        self._loadClassesInfo(self._getFileName('out_class'))
        clsSet.classifyItems(updateItemCallback=self._updateParticle,
                             updateClassCallback=self._updateClass,
                             itemDataIterator=md.iterRows(xmpMd,
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
            numOfClass = len(self.importVolumes)
            for i in range(numOfClass):
                csVolName = ("cryosparc_%s_%s_class_%02d_000%s_volume.mrc" %
                             (self.projectName.get(),
                              self.run3dClassification.get(), i, itera))

                copyFiles(csOutputFolder, self._getExtraPath(), files=[csVolName])

                row = ("%s/cryosparc_%s_%s_class_%02d_000%s_volume.mrc\n" %
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
                if z.startswith('Done iteration'):
                    itera = z.split(' ')[2]

        return itera

    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        validateMsgs = cryosparcValidate()
        if not validateMsgs:
            validateMsgs = gpusValidate(self.getGpuList(),
                                        checkSingleGPU=True)
            if not validateMsgs:
                volumes = self._getInputVolume()
                if volumes is not None and len(volumes) < 2:
                    validateMsgs.append("The number of initial volumes must "
                                        "be equal or greater than 2")
        return validateMsgs

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
            summary.append("Symmetry: %s" %
                           getSymmetry(self.symmetryGroup.get(),
                                       self.symmetryOrder.get()))
            summary.append("------------------------------------------")
            summary.append("Output volumes %s" %
                           self.getObjectTag('outputVolumes'))
            summary.append("Output classes %s" %
                           self.getObjectTag('outputClasses'))

        return summary

    # --------------------------- UTILS functions ---------------------------

    def _defineParamsName(self):
        """ Define a list with all protocol parameters names"""
        self._paramsName = ['multirefine_N',
                            'multirefine_symmetry',
                            'multirefine_sharp_bfactor',
                            'multirefine_force_hard_class',
                            'multirefine_batch_size_per_class',
                            'multirefine_update_rule',
                            'multirefine_online_em_lr_rand',
                            'multirefine_online_em_lr_init',
                            'multirefine_online_em_lr_hl',
                            'multirefine_halfmap_decay',
                            'multirefine_res_init',
                            'multirefine_bp_res_factor',
                            'multirefine_use_max_fsc',
                            'multirefine_assignment_conv_eps',
                            'multirefine_num_rand_assign_iters',
                            'multirefine_num_final_full_iters',
                            'multirefine_noise_model',
                            'multirefine_noise_init_sigmascale',
                            'intermediate_plots',
                            'distribution_plots',
                            'compute_use_ssd']
        self.lane = str(self.getAttributeValue('compute_lane'))

    def do3DClasification(self):
        """
        """
        input_group_connect = {"particles": self.particles.get()}
        group_connect = {"volume": self.volumes}
        params = {}

        for paramName in self._paramsName:
            if (paramName != 'multirefine_symmetry' and
                    paramName != 'multirefine_noise_model' and
                    paramName != 'intermediate_plots' and
                    paramName != 'distribution_plots'):
                params[str(paramName)] = str(self.getAttributeValue(paramName))

            elif paramName == 'multirefine_symmetry':
                symetryValue = getSymmetry(self.symmetryGroup.get(),
                                           self.symmetryOrder.get())

                params[str(paramName)] = symetryValue
            elif paramName == 'multirefine_noise_model':
                params[str(paramName)] = str(NOISE_MODEL_CHOICES[self.multirefine_noise_model.get()])

            elif paramName == 'intermediate_plots' or paramName == 'distribution_plots':
                params[str(paramName)] = str("False")

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