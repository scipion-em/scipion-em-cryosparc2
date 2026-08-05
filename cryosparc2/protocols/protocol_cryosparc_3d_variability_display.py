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
import os.path

from pwem.objects import Volume, SetOfVolumes
from pyworkflow import BETA
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        LEVEL_ADVANCED, EnumParam,
                                        Positive, BooleanParam, StringParam)
from pwem.protocols import ProtRefine3D

from . import ProtCryosparcBase
from ..convert import *
from ..utils import *
from ..constants import *


class ImportVolumeOutputs(enum.Enum):
    component = SetOfVolumes


class ProtCryoSparc3DVariabilityDisplay(ProtCryosparcBase, ProtRefine3D):
    """
    Protocol to create various versions of a 3D variability result that can be
    used for display.

    AI Generated:

    3D Variability Display (ProtCryoSparc3DVariabilityDisplay) — User Manual
        Overview

        The 3D Variability Display protocol generates interpretable visual
        representations from previously computed 3D variability analysis
        results. Its primary purpose is to help biological users explore and
        visualize continuous structural heterogeneity present within cryo-EM
        datasets. Instead of focusing on a single static reconstruction, the
        protocol produces multiple related volumes that describe how molecular
        structures change along variability dimensions identified during
        analysis.

        In practical cryo-EM workflows, this protocol is especially useful for
        studying conformational flexibility, domain motions, ligand-induced
        rearrangements, and compositional heterogeneity. It allows researchers
        to transform abstract variability components into visually meaningful
        structural transitions that can be interpreted biologically and
        communicated more effectively.

        Inputs and General Workflow

        The protocol requires the output of a previous 3D variability analysis.
        This input contains both the reconstructed reference volume and the
        variability information describing structural changes across the
        particle population. The protocol then generates representative outputs
        that summarize these motions using several visualization strategies.

        Depending on the selected mode, the workflow can produce linear motion
        trajectories, discrete structural clusters, or intermediate
        reconstructions that sample conformational landscapes more densely.
        These outputs can subsequently be inspected in molecular visualization
        software, compared structurally, or used for downstream biological
        interpretation.

        Biological Interpretation of Variability

        Structural variability in cryo-EM datasets often reflects meaningful
        biological processes such as substrate binding, enzymatic cycling,
        maturation, assembly dynamics, or flexible domain rearrangements. This
        protocol helps convert statistical variability information into a set
        of interpretable structural states.

        Users should keep in mind that variability components do not always
        correspond to perfectly discrete conformations. In many cases, they
        represent continuous transitions between related molecular states.
        Consequently, the generated outputs should be interpreted as a visual
        approximation of structural motion rather than a strict classification
        into independent biological species.

        Output Modes

        The protocol provides several visualization strategies adapted to
        different biological questions.

        The simple mode generates a sequence of volumes arranged along each
        variability dimension. This creates a movie-like representation of
        structural motion and is often the best starting point for exploratory
        analysis. It is particularly useful for identifying large-scale domain
        movements, hinge motions, or gradual conformational transitions.

        The cluster mode separates particles into groups based on their
        positions within variability space and reconstructs representative
        volumes for each cluster. This mode is valuable when the biological
        system contains partially discrete conformations or compositional
        states. The resulting classes can help distinguish alternative
        assemblies, ligand occupancies, or structural substates.

        The intermediates mode reconstructs multiple volumes distributed along
        the variability trajectory. Compared to the simple mode, this approach
        can better represent non-linear structural transitions and more subtle
        conformational landscapes. It is especially useful for systems with
        continuous flexibility or gradual structural rearrangements.

        Number of Frames or Clusters

        The number of generated frames or clusters controls how finely the
        conformational landscape is sampled. Lower values produce simpler and
        easier-to-interpret outputs, while larger values may reveal more subtle
        transitions but increase computational complexity and interpretation
        difficulty.

        In biological practice, moderate values are often sufficient to capture
        dominant motions. Excessively large numbers of frames may introduce
        redundancy or emphasize noise rather than meaningful structural
        variation. Users are encouraged to balance interpretability with
        structural detail.

        Component Selection

        The protocol allows users to restrict visualization to selected
        variability components. This becomes important when only a subset of
        components contains biologically relevant motion. Focusing on the most
        interpretable dimensions often simplifies downstream analysis and
        improves clarity during visualization.

        Different variability components may capture distinct biological
        phenomena. For example, one component may describe domain breathing
        motions while another reflects compositional changes or local
        flexibility. Careful inspection of each component is therefore
        recommended before selecting those most suitable for interpretation.

        Filtering and Resolution Considerations

        Optional low-pass and high-pass filtering controls help optimize the
        visual appearance of reconstructed volumes. Low-pass filtering is often
        beneficial when emphasizing large conformational changes or reducing
        noise in flexible regions. High-pass filtering may help reveal finer
        structural details but should be applied cautiously to avoid
        over-interpreting noisy features.

        Resolution filtering does not increase structural information content.
        Instead, it modifies the balance between signal clarity and detail
        visibility. Biological conclusions should therefore always be supported
        by consistent structural evidence across the dataset.

        Downsampling and Cropping

        The protocol supports downsampling and cropping of output volumes. These
        options are useful for reducing storage requirements, accelerating
        visualization, or focusing attention on specific structural regions.

        Downsampling is particularly advantageous during exploratory analysis
        of very large complexes or when rapid inspection of multiple states is
        required. Cropping can help isolate flexible domains or remove large
        solvent regions that may distract from biologically important motions.

        Cluster Visualization and Particle Subsets

        In cluster mode, the protocol can generate visual plots describing the
        distribution of particles within variability space. These plots help
        users understand whether conformational states are continuous, partially
        overlapping, or clearly separated into distinct populations.

        In intermediates mode, the protocol can optionally produce particle
        subsets associated with individual frames. This functionality is useful
        when users wish to inspect the particles contributing to particular
        conformational states or perform additional downstream processing on
        selected structural regions of the variability landscape.

        Handedness and Reconstruction Options

        The protocol includes the possibility of flipping the handedness of
        output volumes. This option only affects the displayed maps and should
        be used carefully, particularly in workflows where handedness has
        already been validated experimentally.

        Reconstruction steps may also be skipped in certain exploratory
        scenarios. This enables rapid evaluation of clustering behavior or
        variability partitioning before committing computational resources to
        full reconstruction procedures.

        Outputs and Their Interpretation

        Depending on the selected mode, the protocol produces sets of volumes,
        structural trajectories, clustered reconstructions, or particle
        classifications. These outputs are intended to facilitate visual and
        biological interpretation of conformational variability.

        The resulting volumes preserve the sampling properties and spatial
        reference frame of the original variability analysis, allowing direct
        comparison between generated states. In cluster-oriented workflows,
        representative classes can also be associated with corresponding
        particle subsets for additional refinement or validation.

        Practical Recommendations

        For most biological investigations, it is advisable to begin with the
        simple mode to identify dominant motions and assess whether variability
        is continuous or discrete. Once the main structural transitions are
        understood, cluster or intermediates modes can provide more detailed
        structural interpretation.

        When studying highly flexible complexes, intermediate reconstructions
        often provide the most informative representation of conformational
        landscapes. Conversely, systems containing distinct structural states
        are frequently better analyzed using cluster mode.

        Users should visually inspect generated trajectories carefully and
        compare them with known biochemical or structural information. Apparent
        motions may occasionally reflect noise, alignment uncertainty, or local
        reconstruction artifacts rather than biologically meaningful dynamics.

        Final Perspective

        The visualization of structural variability has become a central aspect
        of modern cryo-EM analysis because many biological systems function
        through dynamic conformational changes rather than static structures.
        This protocol provides a practical framework for transforming
        variability analysis into interpretable structural representations that
        can reveal molecular motions, functional transitions, and heterogeneous
        states in biologically meaningful ways.
    """
    _label = '3D variability Display'
    _devStatus = BETA
    partClassDict = NestedDict(depth=2)

    def _initialize(self):
        self._defineFileNames()

    def _defineFileNames(self):
        """ Centralize how files are called. """
        myDict = {
            'out_class': self._getExtraPath() + '/output_class.star'
        }
        self._updateFilenamesDict(myDict)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('input3DVariablityAnalisysProt', PointerParam,
                      label="3D variability analysis protocol",
                      pointerClass='ProtCryoSparc3DVariability',
                      help='Particle stacks to use.')
        form.addParallelSection(threads=1, mpi=1)

        # --------------[3D Variability Display]---------------------------
        form.addSection(label='3D Variability Output')

        form.addParam('var_output_mode', EnumParam,
                      choices=['cluster', 'simple', 'intermediates'],
                      default=0,
                      label="Output mode",
                      help='simple mode: output a simple linear "movie" of '
                           'volumes along each dimension. Number of frames is '
                           'the next param. '
                           'cluster mode: fit clusters to reaction coordinates '
                           'and output volumes and particles from each cluster. '
                           'Number of clusters is the next param. '
                           'intermediates mode: - reconstruct multiple '
                           'intermediate volumes along each variability '
                           'dimension, useful for better visualizing '
                           'non-linear changes. Number of intermediate frames '
                           'is the next param.')
        form.addParam('var_num_frames', IntParam, default=20,
                      validators=[Positive],
                      label="Number of frames/clusters",
                      help='Number of clusters for cluster mode, or number of '
                           'frames for simple and intermediates mode.')

        form.addParam('var_range_percentile', FloatParam, default=3,
                      validators=[Positive],
                      label="Min/Max Range Percentile (%)",
                      help='Percentile for computing the min and max for '
                           'each component.')

        form.addParam('var_N', IntParam, default=None,
                      allowsNull=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Downsample to box size",
                      help='Downsample the output volumes to this size.')

        form.addParam('var_M', IntParam, default=None,
                      expertLevel=LEVEL_ADVANCED,
                      allowsNull=True,
                      label="Crop to size (after downsample)",
                      help='Crop the output volumes to this size after '
                           'downsampling.')

        # form.addParam('var_num_particles', IntParam, default=None,
        #               validators=[Positive],
        #               label="Only use this many particles",
        #               help='Only use this many particles for reconstructions.')

        form.addParam('var_component_filter', StringParam, default='0',
                      expertLevel=LEVEL_ADVANCED,
                      label="Only use these components",
                      help='Select a subset of 3DVA components to use. Provide a comma-separated, zero-based list '
                           '(e.g., `0,3,4`). Leave as None to use all components. WARNING: This will renumber '
                           'components upon output.')

        form.addParam('var_filter_res', FloatParam, default=None,
                      allowsNull=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Filter resolution (A)",
                      help='Resolution at which results are filtered')

        form.addParam('var_filter_order', FloatParam, default=1.5,
                      validators=[Positive],
                      label="Filter order",
                      help='Order of filter')

        form.addParam('var_highpass_res', FloatParam, default=None,
                      allowsNull=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Highpass resolution (A)",
                      help='Resolution at which results are filtered')

        form.addParam('var_highpass_order', FloatParam, default=8,
                      validators=[Positive],
                      expertLevel=LEVEL_ADVANCED,
                      label="Highpass order",
                      help='Order of filter')

        form.addParam('var_vol_flip', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Flip hand of outputs",
                      help='Flip the hand out output volumes (including series). '
                           'Note that this does not flip output '
                           'alignments - only the output volumes.')

        form.addParam('var_skip_reconstruction', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      condition='var_output_mode == 0 or var_output_mode == 2',
                      label="Skip reconstruction",
                      help='Skip reconstructions in cluster and '
                           'intermediate mode - helpful to quickly inspect '
                           'whether a clustering or division of particles will'
                           ' be helpful.')

        form.addParam('var_cluster_3D_plots', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      condition='var_output_mode == 0',
                      label="Cluster mode: 3D plots",
                      help='Make 3D instead of 2D plots to show clusters.')

        form.addParam('var_intermediate_output_frame_particles', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      condition='var_output_mode==2',
                      label="Intermediates: output particle subsets",
                      help='Output particle stacks that contribute to each intermediate frame along one component.')

        form.addParam('var_intermediate_output_component', IntParam, default=0,
                      condition='var_output_mode==2 and var_intermediate_output_frame_particles==True',
                      expertLevel=LEVEL_ADVANCED,
                      label="Intermediates: output particle subsets",
                      help='Select a single component/series for which particle subsets will be output. Input must be '
                           'a zero-based integer in the range [0, K-1], where K is the number of components that '
                           'exist after the (optional) filtering.')

        form.addParam('var_intermediate_width', IntParam, default=2,
                      expertLevel=LEVEL_ADVANCED,
                      condition='var_output_mode==2',
                      label="Intermediates: window (frames)",
                      help=('Size of rolling window in each bin (in number of frames). Set to a '
                            'positive integer to define a triangular weighting function for each volume. '
                            'Set to 0 for a tophat weighting function. Set to -1 to create bins with equal '
                            'number of particles. '))

        form.addSection(label="Compute settings")
        addComputeSectionParams(form, allowMultipleGPUs=False)

    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._defineFileNames()
        self._defineParamsName()
        self._initializeCryosparcProject()
        self._insertFunctionStep(self.processStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions ------------------------------
    def processStep(self):
        self.info(pwutils.yellowStr("3D Variability started..."))
        self.doRun3DVariabilityDisplay()

    def generateStartFiles(self):
        self.info(pwutils.yellowStr("Copying files from CS to Scipion folder..."))
        csOutputFolder = os.path.join(self.projectDir.get(),
                                      self.run3DVariabilityDisplay.get())
        if self.var_output_mode.get() == 0:  # Cluster mode

            # Copy the CS output to extra folder
            outputFolder = os.path.join(self._getExtraPath(), 'outputs')
            pwutils.cleanPath(outputFolder)
            copyFiles(csOutputFolder, outputFolder)

            # Creating the folder where .star files will be created
            starFilesPath = os.path.join(outputFolder, 'clustersFiles')
            os.makedirs(starFilesPath)

            # Generating the .star file por all clusters
            allClusterSeries = '%s%s_series_all_clusters.cs' % (getOutputPreffix(self.projectName.get()),
                                                             self.run3DVariabilityDisplay.get())
            allClusterSeriesStar = '%s%s_series_all_clusters.star' % (getOutputPreffix(self.projectName.get()),
                                                             self.run3DVariabilityDisplay.get())

            filePath = os.path.join(outputFolder, allClusterSeries)
            outputStarFn = os.path.join(starFilesPath, allClusterSeriesStar)
            argsList = [filePath, outputStarFn]
            convertCs2Star(argsList)

            # Creating the .star cluster per cluster
            numOfClusters = self.var_num_frames.get()
            self.info(pwutils.yellowStr("Generating .star files for all clusters"))
            for i in range(numOfClusters):
                self.info(pwutils.yellowStr("Processing cluster %d ..." % i))

                clusterParticlesPattern = '%s%s_cluster_%03d_particles.cs' % (getOutputPreffix(self.projectName.get()),
                                                                              self.run3DVariabilityDisplay.get(), i)


                passthroughParticlesPattern = '%s%s_passthrough_particles_cluster_%d.cs' % (getOutputPreffix(self.projectName.get()),
                                                                                            self.run3DVariabilityDisplay.get(), i)
                patterns = [clusterParticlesPattern, passthroughParticlesPattern]
                outputs = ['output_cluster_particle%03d.star' % i,
                           'output_passthrough_particle%03d.star' % i]
                # Generating .star files
                for j in range(len(patterns)):
                    filePath = os.path.join(outputFolder, patterns[j])
                    outputStarFn = os.path.join(starFilesPath, outputs[j])
                    argsList = [filePath, outputStarFn]
                    convertCs2Star(argsList)
                self._partClassDict(os.path.join(starFilesPath, outputs[0]), i)
        else:
            # Copy the CS output to extra folder
            outputFolder = self._getExtraPath()
            pwutils.cleanPath(outputFolder)
            copyFiles(csOutputFolder, outputFolder)
            # Creating the .star cluster particles
            numOfComponets = int(self.input3DVariablityAnalisysProt.get().var_K)

            for i in range(numOfComponets):
                self.info(pwutils.yellowStr("Extracting component %d ..." % i))
                # Creating the folder where clusters will be unzipped
                componetFilesPath = self._getExtraPath('component%03d' % i)
                os.makedirs(componetFilesPath)
                componentPattern = "%s%s_component_%03d.zip" % (getOutputPreffix(self.projectName.get()),
                                                                self.run3DVariabilityDisplay.get(), i)
                cmd = "unzip %s -d %s" % (self._getExtraPath(componentPattern), componetFilesPath)
                os.system(cmd)

            if self.var_output_mode.get() == 2 and self.var_intermediate_output_frame_particles.get():
                clusterNumber = self.var_num_frames.get()
                component = self.var_intermediate_output_component.get()
                componetFilesPath = self._getExtraPath('component%03d' % component)

                self._create2DModelFile(componetFilesPath)

                for i in range(clusterNumber):
                    componentCSParticlesPattern = '%s%s_particles_series_%d_frame_%d.cs' % (getOutputPreffix(self.projectName.get()),
                                                                                          self.run3DVariabilityDisplay.get(),
                                                                                          component, i)
                    componentStarPattern = '%s%s_particles_series_%d_frame_%d.star' % (getOutputPreffix(self.projectName.get()),
                                                                                       self.run3DVariabilityDisplay.get(),
                                                                                       component, i)
                    outputStarPath = os.path.join(componetFilesPath, componentStarPattern)

                    filePath = os.path.join(outputFolder, componentCSParticlesPattern)
                    argsList = [filePath, outputStarPath]
                    convertCs2Star(argsList)

                    self._partClassDict(outputStarPath, i)

    def createOutputStep(self):
        self._initializeUtilsVariables()
        self.generateStartFiles()
        self.info(pwutils.yellowStr("Creating the output..."))

        if self.var_output_mode.get() == 0:
            outputFolder = os.path.join(self._getExtraPath(), 'outputs')
            self._create3DModelFile(outputFolder)

            # Create 3D classes
            imgSet = self.input3DVariablityAnalisysProt.get().outputParticles
            classes3D = self._createSetOfClasses3D(imgSet)
            self._fillClassesFromIter(classes3D)

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
        else:
            numOfComponets = int(self.input3DVariablityAnalisysProt.get().var_K)
            referenceVol = self.input3DVariablityAnalisysProt.get().outputVolume
            for i in range(numOfComponets):
                componetFilesPath = self._getExtraPath('component%03d' % i)
                # create a SetOfVolumes and define its relations
                volSet = self._createSetOfVolumes(suffix='%03d' % i)
                numOfClusters = self.var_num_frames.get()
                for j in range(numOfClusters):
                    framePattern = '%s%s_component_%03d_frame_%03d.mrc' % (getOutputPreffix(self.projectName.get()),
                                                                           self.run3DVariabilityDisplay.get(), i, j)
                    localFileName = os.path.join(componetFilesPath, framePattern)
                    vol = Volume()
                    vol.setFileName(localFileName)
                    vol.setSamplingRate(referenceVol.getSamplingRate())
                    vol.setOrigin(referenceVol.getOrigin())
                    volSet.append(vol)

                volSet.setSamplingRate(referenceVol.getSamplingRate())
                self._defineOutputs(**{ImportVolumeOutputs.component.name + '%03d' % i: volSet})

            if self.var_output_mode.get() == 2 and self.var_intermediate_output_frame_particles.get():
                imgSet = self.input3DVariablityAnalisysProt.get().outputParticles
                classes2DSet = self._createSetOfClasses2D(imgSet)
                self._fillClassesFromIter(classes2DSet)
                component = self.var_intermediate_output_component.get()

                self._defineOutputs(**{"SetOfClass2D_" + ImportVolumeOutputs.component.name + '%03d' % component: classes2DSet})
                self._defineSourceRelation(imgSet, classes2DSet)

    # --------------------------- UTILS functions ------------------------------
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
            self._classesInfo[classNumber + 1] = (index, fn, row)

    def _fillClassesFromIter(self, clsSet):
        """ Create the SetOfClasses3D """
        self._loadClassesInfo(self._getFileName('out_class'))
        clsSet.classifyItems(updateItemCallback=self._updateParticle,
                             updateClassCallback=self._updateClass,
                             itemDataIterator=None)

    def _updateParticle(self, item, row):
        classNumber = self.partClassDict.search([os.path.splitext(os.path.basename(item.getFileName()))[0],
                                                 item.getIndex()])
        if not classNumber:
            setattr(item, '_appendItem', False)
        else:
            item.setClassId(classNumber)

    def _updateClass(self, item):
        classId = item.getObjId()
        imgSet = self.input3DVariablityAnalisysProt.get().outputParticles
        if classId in self._classesInfo:
            index, fn, row = self._classesInfo[classId]
            fixVolume(fn)
            item.setAlignmentProj()
            vol = item.getRepresentative()
            vol.setLocation(index, fn)
            vol.setSamplingRate(calculateNewSamplingRate(vol.getDim(),
                                                         imgSet.getSamplingRate(),
                                                         imgSet.getDim()))

    def _create3DModelFile(self, volumesPath):
        # Create model files for 3D classification
        with open(self._getFileName('out_class'), 'w') as output_file:
            output_file.write('\n')
            output_file.write('data_images')
            output_file.write('\n\n')
            output_file.write('loop_')
            output_file.write('\n')
            output_file.write('_rlnReferenceImage')
            output_file.write('\n')
            numOfClass = self.var_num_frames.get()
            for i in range(numOfClass):
                csVolName = '%s%s_cluster_%03d.mrc' % (getOutputPreffix(self.projectName.get()),
                                                       self.run3DVariabilityDisplay.get(), i)
                row = ("%s/%s\n" % (volumesPath, csVolName))
                output_file.write(row)

    def _create2DModelFile(self, filePath):
        with open(self._getFileName('out_class'), 'w') as output_file:
            output_file.write('\n')
            output_file.write('data_particles')
            output_file.write('\n\n')
            output_file.write('loop_')
            output_file.write('\n')
            output_file.write('_rlnReferenceImage')
            output_file.write('\n')
            component = self.var_intermediate_output_component.get()
            for i in range(self.var_num_frames.get()):
                csVolName = '%s%s_component_%03d_frame_%03d.mrc' % (getOutputPreffix(self.projectName.get()),
                                                                    self.run3DVariabilityDisplay.get(), component, i)

                row = ("%05d@%s/%s\n" % (i+1, filePath, csVolName))
                output_file.write(row)

    def _partClassDict(self, fileName, classNumber):
        with open(os.path.abspath(fileName), 'r') as _file:
            rawLine = _file.readline()

            # Find data_particle block
            while rawLine:
                if rawLine.startswith('data_particles'):
                    break
                rawLine = _file.readline()

            # Find metadatada labels
            while rawLine:
                if rawLine.startswith('loop_'):
                    break
                rawLine = _file.readline()

            rawLine = _file.readline()
            while rawLine.startswith('_'):
                rawLine = _file.readline()

            # Taking the particles ids
            while rawLine:
                line = rawLine.strip().split()[0].split('@')
                partId = int(line[0])
                imageName = '_'.join(line[1].split('_')[1:]).split('.')[0]
                self.partClassDict.insert([imageName, partId], classNumber + 1)
                rawLine = _file.readline()

    def _defineParamsName(self):
        """ Define a list with all protocol parameters names"""
        self._paramsName = ['var_output_mode',
                            'var_num_frames',
                            'var_range_percentile',
                            'var_N',
                            'var_M',
                            'var_filter_res',
                            'var_filter_order',
                            'var_highpass_res',
                            'var_highpass_order',
                            'var_vol_flip',
                            'var_skip_reconstruction',
                            'var_cluster_3D_plots',
                            'var_intermediate_width',
                            'var_component_filter',
                            'var_intermediate_output_frame_particles',
                            'var_intermediate_output_component',
                            'compute_use_ssd']
        self.lane = str(self.getAttributeValue('compute_lane'))

    def _validate(self):
        validateMsgs = cryosparcValidate()
        if not validateMsgs:
            validateMsgs = gpusValidate(self.getGpuList(), checkSingleGPU=True)
        return validateMsgs

    def doRun3DVariabilityDisplay(self):
        """
        :return:
        """
        className = "var_3D_disp"
        varAnalysisJob = str(self.input3DVariablityAnalisysProt.get().run3DVariability)
        input_group_conect = {"particles": str('%s.particles' %varAnalysisJob) ,
                              "volume": str('%s.volume'%varAnalysisJob)}
        params = {}

        for paramName in self._paramsName:
            if (paramName != 'var_output_mode' and paramName != 'var_N' and
                    paramName != 'var_M' and paramName != 'var_filter_res' and
                    paramName != 'var_highpass_res' and paramName != 'var_intermediate_width' and
                    paramName != 'var_component_filter'):
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif paramName == 'var_output_mode':
                params[str(paramName)] = str(VAR_OUTPUT_MODE[self.var_output_mode.get()])
            elif paramName == 'var_N' and self.var_N.get() is not None:
                params[str(paramName)] = str(self.var_N.get())
            elif paramName == 'var_M' and self.var_M.get() is not None:
                params[str(paramName)] = str(self.var_M.get())
            elif paramName == 'var_filter_res' and self.var_filter_res.get() is not None:
                params[str(paramName)] = str(self.var_filter_res.get())
            elif paramName == 'var_highpass_res' and self.var_highpass_res.get() is not None:
                params[str(paramName)] = str(self.var_highpass_res.get())
            elif paramName == 'var_intermediate_width':
                if self.var_intermediate_width.get() is not None or self.var_intermediate_width.get() > 0:
                    params[str(paramName)] = str(self.var_intermediate_width.get())
                else:
                    params[str(paramName)] = '0'

            # elif paramName == 'var_component_filter':
            #     params[str(paramName)] = str(self.var_component_filter.get())

        # Determinate the GPUs to use (in dependence of the cryosparc version)
        try:
            if not self.useQueueForSteps() and not self.useQueue():  # not using queue system
                gpusToUse = self.getGpuList()
            else:  # using queue system
                gpusToUse = False
        except Exception:
            gpusToUse = False

        self.run3DVariabilityDisplay = enqueueJob(className,
                                                  self.projectName.get(),
                                                  self.workSpaceName.get(),
                                                  str(params).replace('\'',
                                                                      '"'),
                                                  str(input_group_conect).replace(
                                                      '\'',
                                                      '"'),
                                                  self.lane, gpusToUse)

        self.currenJob.set(self.run3DVariabilityDisplay.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(),
                         self.run3DVariabilityDisplay.get(),
                         "An error occurred in the 3D Variability process. "
                         "Please, go to cryosPARC software for more "
                         "details.", self)
