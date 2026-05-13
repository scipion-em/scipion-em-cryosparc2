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

import pyworkflow.utils as pwutils
from pyworkflow.object import String
from pyworkflow.protocol.params import (FloatParam, LEVEL_ADVANCED,
                                        PointerParam, MultiPointerParam,
                                        CsvList, Positive, IntParam,
                                        BooleanParam, StringParam, EnumParam)

from .protocol_base import ProtCryosparcBase
from ..convert import (convertBinaryVol, convertCs2Star,
                       rowToAlignment, ALIGN_PROJ, cryosparcToLocation)
from ..utils import (addSymmetryParam, addComputeSectionParams, doImportVolumes,
                     get_job_streamlog, calculateNewSamplingRate,
                     cryosparcValidate, gpusValidate, getSymmetry, enqueueJob,
                     waitForCryosparc, clearIntermediateResults, fixVolume,
                     copyFiles, getOutputPreffix, matchItemRow)
from ..constants import *


class ProtCryoSparc3DClassification(ProtCryosparcBase):
    """
    Performs heterogeneous 3D refinement and particle classification in cryo-EM datasets using multiple initial reference volumes. The protocol is designed to separate structurally distinct particle populations while simultaneously refining the corresponding 3D reconstructions, enabling the identification of conformational variability, compositional heterogeneity, or differences in structural quality within a dataset.

    AI Generated:

    3D Heterogeneous Refinement (ProtCryoSparc3DClassification) — User Manual
        Overview

        The 3D Heterogeneous Refinement protocol classifies particles into
        multiple structural groups while refining an independent 3D map for
        each group. This strategy is particularly important in cryo-EM studies
        where the sample contains several conformational states, partially
        assembled complexes, ligand-bound and ligand-free populations, or
        damaged and low-quality particles mixed together.

        Rather than treating all particles as belonging to a single structure,
        the protocol continuously evaluates similarities between particles and
        multiple reference maps. As refinement progresses, particles are
        redistributed among the classes and the maps improve iteratively. This
        approach allows subtle biological differences to emerge even when they
        are not initially obvious at low resolution.

        In practical cryo-EM workflows, heterogeneous refinement is commonly
        applied after Ab-Initio Reconstruction or after obtaining preliminary
        consensus refinements. The protocol is often one of the key stages for
        understanding structural dynamics and biochemical variability.

        Inputs and Biological Context

        The protocol requires a particle dataset together with at least two
        initial reference volumes. These references define the starting point
        for the classification process and ideally represent distinct
        conformational or compositional states. In many workflows, the initial
        maps are generated from ab-initio reconstruction methods or previous
        classification steps.

        Biological interpretation strongly depends on the quality and diversity
        of the starting references. If all references are too similar, the
        classification may fail to separate meaningful structural states.
        Conversely, if references are unrealistically different or contain
        artifacts, particles may be forced into biologically incorrect classes.

        Input particles should already contain accurate imaging information and
        reliable alignment parameters. Proper preprocessing, including motion
        correction, CTF estimation, and particle polishing when appropriate,
        generally improves the stability and interpretability of the refinement.

        Classification and Refinement Strategy

        The protocol performs classification and refinement simultaneously.
        During each refinement cycle, particles are probabilistically assigned
        to the reference map that best explains their signal. At the same time,
        each map is reconstructed and improved using the particles associated
        with that class.

        This iterative approach is especially powerful for detecting continuous
        or discrete conformational variability. For example, flexible domains,
        binding events, rotational rearrangements, or assembly intermediates
        may emerge naturally as distinct structural classes.

        In many biological systems, some particle populations are much smaller
        than others. The protocol is capable of preserving minority classes,
        although very small populations may still become unstable or noisy if
        insufficient particle numbers are available.

        Symmetry Considerations

        Symmetry can be applied during refinement to improve signal quality and
        reconstruction stability. For highly symmetric particles such as viral
        capsids or oligomeric assemblies, symmetry enforcement often produces
        significant resolution improvements.

        However, symmetry should only be imposed when biologically justified.
        Incorrect symmetry can obscure meaningful asymmetry, hide flexible
        regions, or merge distinct conformational states into artificial
        averages. For complexes with partial symmetry breaking or asymmetric
        ligand binding, lower symmetry or fully asymmetric refinement may be
        more appropriate.

        Hard and Soft Classification

        The protocol supports both probabilistic and strict particle assignment
        strategies. In probabilistic classification, particles may contribute
        partially to multiple classes during refinement. This approach is often
        beneficial when structural differences are subtle or when conformational
        transitions are continuous.

        Hard classification forces each particle into a single class at every
        iteration. This strategy may improve class separation in some cases,
        particularly when the expected structural states are highly distinct.
        However, it can also introduce instability when transitions between
        states are gradual.

        From a biological perspective, soft classification is generally more
        tolerant of structural continua, whereas hard classification emphasizes
        discrete state separation.

        Resolution and Box Size Considerations

        The refinement box size determines the working resolution and memory
        requirements of the protocol. Smaller box sizes reduce computational
        cost and are often sufficient during early exploratory analyses.
        Larger box sizes preserve higher-frequency information but require
        substantially more GPU memory and processing time.

        Initial low-pass filtering is commonly used to stabilize early
        refinement iterations. Starting from lower resolution information helps
        prevent overfitting and encourages robust convergence, especially when
        initial references are noisy or uncertain.

        In practical workflows, users often begin with moderate resolutions and
        later continue refinement using higher-resolution settings once stable
        classes have emerged.

        Noise Modeling and Optimization

        The protocol includes several optimization and noise modeling strategies
        that influence refinement behavior. Most standard cryo-EM workflows can
        rely on the default parameters, which are designed to provide stable
        convergence across a wide range of datasets.

        Advanced optimization settings become relevant mainly in difficult
        cases, such as highly heterogeneous samples, extremely noisy datasets,
        or very large particle collections. Careful tuning may improve class
        separation, but excessive parameter manipulation can also destabilize
        refinement or introduce overfitting.

        In general biological practice, conservative parameter choices are
        recommended unless there is a clear experimental reason to modify the
        optimization behavior.

        Outputs and Their Interpretation

        The protocol produces a set of refined 3D classes together with their
        associated particle assignments. Each class corresponds to a refined
        volume representing one structural population within the dataset.

        The resulting maps can reveal distinct conformations, binding states,
        assembly intermediates, or quality differences among particles.
        Particle assignments may also be used for downstream refinement,
        focused classification, variability analysis, or atomic modeling.

        Interpretation should always consider the particle distribution across
        classes. Large classes generally produce more stable reconstructions,
        while very small classes may contain either rare biological states or
        poorly aligned particles. Visual inspection and biological consistency
        remain essential for distinguishing meaningful heterogeneity from noise.

        Practical Recommendations

        For most cryo-EM projects, it is advisable to begin with diverse but
        biologically plausible initial references. If classification collapses
        into nearly identical classes, stronger structural diversity among the
        starting maps may be needed. Conversely, unrealistic starting maps may
        bias the refinement toward incorrect solutions.

        Moderate box sizes and default optimization settings are usually
        sufficient for initial exploration. Once meaningful classes appear,
        selected subsets can be refined further at higher resolution using more
        specialized refinement strategies.

        Biological users should carefully inspect not only the final resolution
        values but also the structural interpretability of each class. Flexible
        regions, domain motions, and ligand densities often provide more
        meaningful insight than global resolution alone.

        Final Perspective

        Heterogeneous refinement is one of the most biologically informative
        stages in modern cryo-EM analysis because it transforms structural
        variability into interpretable three-dimensional states. Successful
        application depends on thoughtful selection of initial references,
        careful interpretation of class distributions, and biological awareness
        of the expected conformational landscape.

        When used appropriately, the protocol enables researchers to move
        beyond consensus averaging and directly investigate the dynamic and
        heterogeneous nature of macromolecular systems.
    """
    _label = '3D Heterogeneous Refinement'
    _className = "hetero_refine"

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
        self.info(pwutils.yellowStr("3D Heterogeneous Refinement started..."))
        self.do3DClasification()

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to Relion file
        """
        self._initializeUtilsVariables()
        self.info(pwutils.yellowStr("Creating the output..."))

        csOutputFolder = os.path.join(self.projectDir.get(),
                              self.run3dClassification.get())
        itera = self.findLastIteration(self.run3dClassification.get())

        csParticlesName = "%s%s_000%s_particles.cs" % (getOutputPreffix(self.projectName.get()),
                                                       self.run3dClassification.get(),
                                                       itera)
        # Copy the CS output particles to extra folder
        copyFiles(csOutputFolder, self._getExtraPath(), files=[csParticlesName])

        csFile = os.path.join(self._getExtraPath(), csParticlesName)

        outputStarFn = self._getFileName('out_particles')
        argsList = [csFile, outputStarFn]

        convertCs2Star(argsList)

        self._createModelFile(csOutputFolder, itera)

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
                             itemDataIterator=emtable.Table.iterRows(xmpMd),
                             raiseOnNextFailure=False,
                             cancelNextWhenAppendIsFalse=True)

    def _updateParticle(self, item, row):
        if matchItemRow(item, row):
            item.setClassId(row.get(RELIONCOLUMNS.rlnClassNumber.value))
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
                csVolName = ("%s%s_class_%02d_000%s_volume.mrc" %
                             (getOutputPreffix(self.projectName.get()),
                              self.run3dClassification.get(), i, itera))

                copyFiles(csOutputFolder, self._getExtraPath(), files=[csVolName])

                row = ("%s/%s%s_class_%02d_000%s_volume.mrc\n" %
                       (self._getExtraPath(), getOutputPreffix(self.projectName.get()),
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
                particles = self._getInputParticles()
                if not particles.hasCTF():
                    validateMsgs.append(
                        "The Particles has not associated a "
                        "CTF model")
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
            if not self.useQueueForSteps() and not self.useQueue():  # not using queue system
                gpusToUse = self.getGpuList()
            else:  # using queue system
                gpusToUse = False
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
                         "Please, go to cryoSPARC software for more "
                         "details.", self)
        clearIntermediateResults(self.projectName.get(), self.run3dClassification.get())