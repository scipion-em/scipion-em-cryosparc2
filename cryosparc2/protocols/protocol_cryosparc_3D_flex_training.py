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
from enum import Enum

from pwem import getMatchingFiles
from pwem.objects import ParticleFlex, SetOfParticles, SetOfParticlesFlex
from pwem.protocols import ProtFlexBase
from pyworkflow import BETA
from pyworkflow.protocol import LEVEL_ADVANCED
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        BooleanParam)
from . import ProtCryosparcBase
from ..convert import *
from ..utils import *
from ..constants import *


class outputs(Enum):
    Particles = SetOfParticles


class ProtCryoSparc3DFlexTraining(ProtCryosparcBase, ProtFlexBase):
    """
    Uses a mesh and prepared particles (at a downsampled resolution) to train
    a 3DFlex model. Parameters control the number of latent dimensions,
    size of the model, and training hyperparameters. This job outputs
    checkpoints during training.
    """
    """
        Trains a CryoSPARC 3DFlex model using previously prepared particles and a
        tetrahedral mesh representation of the structure. The protocol learns the
        continuous conformational variability of macromolecular complexes by modeling
        structural deformations in a low-dimensional latent space. During training,
        the protocol generates checkpoints and latent coordinates that describe the
        motion and flexibility of particles across different conformational states.

        AI Generated:

        3D Flex Training (ProtCryoSparc3DFlexTraining) — User Manual
            Overview

            The 3D Flex Training protocol performs the core learning stage of the
            CryoSPARC 3DFlex workflow. Using downsampled particles together with a
            tetrahedral deformation mesh, the protocol trains a neural deformation
            model capable of representing continuous molecular motion directly from
            cryo-EM particle images.

            Unlike traditional discrete classification approaches, this protocol
            models conformational variability as a continuous latent space. Each
            particle receives latent coordinates describing its structural state,
            allowing biological motions to be explored as smooth transitions rather
            than isolated classes. This makes the protocol particularly useful for
            studying flexible proteins, dynamic assemblies, multi-domain complexes,
            and molecular machines exhibiting continuous motion.

            Inputs and General Workflow

            The protocol requires two main inputs. The first is a 3D Flex Data
            Prepare protocol containing the downsampled particles and associated
            metadata. The second is a 3D Flex Mesh Prepare protocol containing the
            tetrahedral mesh used to regularize and constrain structural deformations.

            During execution, the protocol connects both particle and mesh data to
            the CryoSPARC training job and launches the flexible refinement model.
            The training process optimizes both the canonical density map and the
            deformation field simultaneously while learning latent coordinates for
            each particle.

            The workflow is designed to capture biologically meaningful continuous
            motion while maintaining physically realistic deformations through mesh
            regularization and rigidity constraints.

            Latent Space Representation

            One of the central concepts of the protocol is the latent space, defined
            by the number of latent dimensions parameter. This latent space represents
            the continuous conformational landscape explored by the particles.

            Lower-dimensional latent spaces are generally easier to interpret and are
            often sufficient for simple motions such as hinge bending or domain
            rotations. In many biological applications, starting with two latent
            dimensions provides a good balance between interpretability and flexibility.
            More complex systems containing multiple coupled motions may require higher
            dimensional latent representations.

            From a biological perspective, latent coordinates can later be analyzed
            to identify conformational trajectories, structural transitions, or
            distinct motion pathways within heterogeneous datasets.

            Neural Network Architecture

            The protocol allows control over the architecture of the deformation
            network through parameters such as the number of layers and hidden units.
            Larger neural networks can model more complex motions and deformation
            patterns but also increase computational cost and the risk of overfitting.

            For most routine cryo-EM datasets, the default architecture provides a
            good compromise between flexibility and stability. Increasing network size
            may become useful when studying highly dynamic complexes with large-scale
            structural rearrangements.

            Excessively large models, however, may start fitting noise instead of
            biologically meaningful motion, especially in datasets with limited signal
            or poor particle quality.

            Learning Rates and Optimization

            Separate learning rates are defined for the deformation model and the
            canonical density map. This separation allows independent control over
            how rapidly structural deformations and density features evolve during
            optimization.

            Higher learning rates can accelerate convergence but may lead to unstable
            training or unrealistic deformations. Lower learning rates generally
            improve stability at the expense of longer training times.

            The protocol progressively anneals learning rates during training while
            gradually increasing the effective resolution of the canonical density.
            This strategy improves robustness during early optimization and refines
            structural details during later stages.

            Rigidity Regularization

            Rigidity regularization is one of the most biologically important aspects
            of the protocol. The rigidity lambda parameter controls how smoothly the
            tetrahedral mesh deforms during training.

            Higher rigidity values favor smooth, coordinated motions and reduce the
            likelihood of unrealistic local distortions. Lower rigidity values allow
            more flexible and localized deformations but may increase the risk of
            overfitting noise or generating physically implausible motions.

            In practice, relatively rigid settings are often preferred for small
            particles, low signal-to-noise datasets, or systems with subtle motions.
            More flexible settings may be beneficial when studying highly dynamic
            assemblies or large conformational rearrangements.

            Latent Space Regularization

            The protocol includes several parameters dedicated to stabilizing and
            organizing the latent space representation. Noise injection during latent
            inference helps smooth the conformational landscape and prevents unstable
            latent clustering.

            The latent centering parameters constrain latent coordinates to remain
            distributed around the center of latent space. Proper balancing of these
            parameters improves the continuity and interpretability of conformational
            trajectories.

            If latent coordinates accumulate excessively near the boundaries of the
            estimation range, increasing the latent centering strength can stabilize
            training. Conversely, overly concentrated latent distributions may indicate
            excessive regularization.

            Training Outputs

            After training, the protocol produces several important outputs. The most
            significant output is a particle set containing latent coordinates for
            every particle. These latent values encode the conformational state of
            each observation within the learned flexibility landscape.

            The protocol also stores CryoSPARC checkpoint files that preserve the
            trained model and allow downstream protocols such as 3D Flex Generator
            or 3D Flex Reconstruction to continue processing.

            Internally, the protocol extracts latent coordinates from CryoSPARC output
            files and associates them with Scipion particle objects, creating a
            flexible particle dataset compatible with downstream flexibility analysis.

            Biological Interpretation

            Biologically, the trained latent space can reveal continuous molecular
            motions that are often hidden in conventional discrete classifications.
            Smooth trajectories through latent space may correspond to domain motions,
            ligand-induced rearrangements, breathing motions, or assembly transitions.

            Interpretation should nevertheless be performed carefully. Not every latent
            dimension necessarily corresponds to a unique biological motion, and some
            variability may still reflect noise, preferred orientations, or imperfect
            reconstruction conditions.

            Visual inspection of generated trajectories and reconstructed maps remains
            essential to validate that observed motions are structurally and
            biologically meaningful.

            Practical Recommendations

            For most cryo-EM workflows, beginning with two latent dimensions and
            default training parameters provides a robust starting point. If the
            resulting latent space appears too simple or fails to capture observed
            heterogeneity, additional latent dimensions can be introduced gradually.

            When training becomes unstable or generates unrealistic motions, increasing
            rigidity regularization or reducing network complexity often improves
            robustness. Conversely, highly dynamic systems may benefit from slightly
            lower rigidity constraints and richer latent representations.

            Careful inspection of latent distributions, generated trajectories, and
            reconstructed density maps is critical before drawing biological conclusions.

            Final Perspective

            The 3D Flex Training protocol represents the central learning stage of the
            CryoSPARC flexibility framework. Rather than separating particles into
            discrete structural classes, it models structural heterogeneity as a
            continuous conformational landscape.

            For biological users, this provides a powerful framework to investigate
            molecular dynamics directly from cryo-EM data, enabling visualization and
            interpretation of motions that are often inaccessible through conventional
            refinement approaches.
        """
    _label = '3D flex training'
    _devStatus = BETA
    _protCompatibility = [V4_1_0, V4_1_1, V4_1_2, V4_2_0, V4_2_1, V4_3_1, V4_4_0, V4_4_1, V4_5_1,
                          V4_5_3, V4_6_0, V4_6_1, V4_6_2, V4_7_0, V4_7_1]
    _possibleOutputs = outputs

    # --------------------------- DEFINE param functions ----------------------
    def _defineFileNames(self):
        """ Centralize how files are called within the protocol. """
        myDict = {
            'input_particles': self._getTmpPath('input_particles.star'),
            'out_particles': self._getPath() + '/output_particle.star',
            'stream_log': self._getPath() + '/stream.log'
        }
        self._updateFilenamesDict(myDict)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('input3DFlexDataPrepareProt', PointerParam,
                      pointerClass='ProtCryoSparc3DFlexDataPrepare',
                      label="3D flex data prepare protocol",
                      important=True,
                      help='Particle stacks to use.')

        form.addParam('input3DMeshFlexPrepareProt', PointerParam,
                      pointerClass='ProtCryoSparc3DFlexMeshPrepare',
                      label="3D flex mesh prepare protocol",
                      important=True,
                      help='3d Flex mesh.')

        form.addParam('flex_K', IntParam, default=2,
                      label="Number of latent dims",
                      help="Number of latent dimensions in the flex refine "
                           "model. See guide for more details. Typically, "
                           "start with 2 and increase if the data appears to "
                           "have more modes of motion present.")

        form.addParam('flex_num_layers', IntParam, default=6,
                      expertLevel=LEVEL_ADVANCED,
                      label="Number of layers",
                      help="Number of layers in the flow generator network. "
                           "See guide for more details. Larger networks can "
                           "be more expressive but can increase training time "
                           "and propensity for overfitting.")

        form.addParam('flex_hidden_units', IntParam, default=64,
                      label="Number of hidden units",
                      help="Number of hidden units per layer in the flow "
                           "generator network. See guide for more details. "
                           "Larger networks can be more expressive but can "
                           "increase training time and propensity for "
                           "overfitting.")

        form.addParam('flex_lr_flex_init', FloatParam, default=0.005,
                      label="Flow learning rate initial",
                      help="Initial learning rate for flow generator params")

        form.addParam('flex_lr_flex_final', FloatParam, default=0.0005,
                      label="Flow learning rate final",
                      help="Final learning rate for flow generator params")

        form.addParam('flex_lr_density_init', FloatParam, default=0.01,
                      label="Density learning rate initial",
                      help="Initial learning rate for canonical density map")

        form.addParam('flex_lr_density_final', FloatParam, default=0.01,
                      label="Density learning rate final",
                      help="Final learning rate for canonical density map")

        form.addParam('flex_sv_lam', FloatParam, default=2.0,
                      label="Rigidity (lambda)",
                      help="Rigidity prior strength. This modulates the "
                           "rigidity of all tetra elements based on the "
                           "element rigidity weighting. Increasing this value "
                           "ensures motions are relatively more smooth.")

        form.addParam('flex_extra_epochs', IntParam, default=0,
                      label="Number of extra epochs",
                      help="Number of extra additional final epochs of "
                           "training to run. By default, 16 epochs are done "
                           "and during this time, the training schedule anneals "
                           "the learning rates as well as the canonical map"
                           " resolution. The resolution goes up to 80% of the "
                           "training box size Nyquist.")

        form.addParam('flex_latent_samp_std', FloatParam, default=0.15,
                      label="Noise injection stdev.",
                      help="Standard deviation of noise injected during latent "
                           "inference. Latent coordinates typically range "
                           "between (-1.5, 1.5). Larger values of this "
                           "parameter introduce more noise during estimation, "
                           "forcing the deformation model to be smoother over "
                           "the latent space. Smaller values allow the latent "
                           "coordinates to be estimated and retained with more "
                           "precision, but can sometimes lead to poorly "
                           "structured latent spaces. See guide for more details.")

        form.addParam('flex_latent_prior_lam', IntParam, default=20,
                      label="Latent centering strength",
                      help="Strength of prior that keeps latent coordinates "
                           "centered and distributed around (0,0) in the "
                           "latent space. This typically needs to be tuned "
                           "for every dataset, but has relatively little "
                           "effect on results. If you notice many latent coordinates "
                           "reaching the end of the (-1.5, 1.5) estimation "
                           "range, this value should be increased. If latent "
                           "coordinates are very concentrated around (0,0) "
                           "then this value should be decreased.")

        form.addParam('flex_latent_prior_pow', FloatParam, default=4.0,
                      label="Latent centering pow",
                      expertLevel=LEVEL_ADVANCED)

        form.addParam('flex_latent_ext_init', BooleanParam, default=True,
                      label="Initialize latents from input",
                      expertLevel=LEVEL_ADVANCED)

        """
        job.param_add('flex_train', "flex_latent_ext_init_idxs",base_value=None,   title="Initialize latents input indices",  param_type="string", desc="Comma separated list of (zero-based) indices for which input components to use for initializing latent coordinates. This list should be the same length as the latent dimension specified.")

        job.param_add('flex_train', "flex_force_restart", base_value=False,  title="Restart training",  param_type="boolean", desc="Force restart of training even if a model with trained checkpoint is connected.", hidden=True)
        """


        # --------------[Compute settings]---------------------------
        form.addSection(label="Compute settings")
        addComputeSectionParams(form, allowMultipleGPUs=False, needGPU=True)

    def _insertAllSteps(self):
        self._defineFileNames()
        self._defineParamsName()
        self._initializeCryosparcProject()
        self._insertFunctionStep(self.doRun3DFlexTraining)
        self._insertFunctionStep(self.createOutputStep)

    def trainingStep(self):
        self.info(pwutils.yellowStr("3D Flex Training started..."))
        self.doRun3DFlexTraining()

    def createOutputStep(self):
        self._initializeUtilsVariables()
        self.info(pwutils.yellowStr("Creating the output..."))

        csOutputFolder = os.path.join(self.projectDir.get(),
                                      self.run3DFlexTrainJob.get())

        pattern = csOutputFolder + '/*latents*'
        csParticlesName = os.path.basename(getMatchingFiles(pattern, True)[-1])

        pattern = csOutputFolder + '/*train_checkpoint*.tar'
        trainModelRar = os.path.basename(getMatchingFiles(pattern, True)[-1])

        pattern = csOutputFolder + '/*train_checkpoint*.cs'
        trainModelCs = os.path.basename(getMatchingFiles(pattern, True)[-1])


        # Copy the CS output particles to extra folder
        copyFiles(csOutputFolder, self._getExtraPath(), files=[csParticlesName, trainModelRar, trainModelCs])
        csPartFile = os.path.join(self._getExtraPath(), csParticlesName)

        # Taking the zvalues from the .cs file using numpy
        arr = np.load(csPartFile)
        zValues = [[arr[i][j] for j in range(2, len(arr[i]), 2)] for i in range(len(arr))]

        inputSet = self.input3DFlexDataPrepareProt.get()._getInputParticles()
        outImgSet = SetOfParticlesFlex.create(self._getPath(), suffix='', progName=CRYOSPARCFLEX)

        outImgSet.copyInfo(inputSet)
        outImgSet.setHasCTF(inputSet.hasCTF())
        outImgSet.getFlexInfo().setProgName(CRYOSPARCFLEX)
        outImgSet.getFlexInfo().setAttr('projectId', str(self.projectName.get()))
        outImgSet.getFlexInfo().setAttr('workSpaceId', str(self.workSpaceName.get()))
        outImgSet.getFlexInfo().setAttr('trainJobId', str(self.run3DFlexTrainJob.get()))
        outImgSet.getFlexInfo().setAttr('projectPath', self.projectDir.get())

        for particle, zValue in zip(inputSet, zValues):
            outParticle = ParticleFlex(progName=CRYOSPARCFLEX)
            outParticle.copyInfo(particle)
            outParticle.getFlexInfo().setProgName(CRYOSPARCFLEX)

            outParticle.setZFlex(list(zValue))

            outImgSet.append(outParticle)

        self._defineOutputs(**{outputs.Particles.name: outImgSet})
        self._defineSourceRelation(inputSet, outImgSet)

        # This is an example to create a latent trajectory in order to launch the flex generator job
        # arr = np.stack([zValues[20], zValues[21]], axis=0)
        # latentTrajectoryJob = customLatentTrajectory(arr,
        #                                              str(self.projectName.get()),
        #                                              str(self.workSpaceName.get()),
        #                                              str(self.run3DFlexTrainJob.get()))
        #
        # flexGeneratorJob = runFlexGeneratorJob(str(self.run3DFlexTrainJob.get()),
        #                                        latentTrajectoryJob,
        #                                        str(self.projectName.get()),
        #                                        str(self.workSpaceName.get()))



    def _defineParamsName(self):
        """ Define a list with 3D Flex Training parameters names"""
        self._paramsName = ['flex_K', 'flex_num_layers', 'flex_num_layers',
                        'flex_hidden_units', 'flex_lr_flex_init',
                        'flex_lr_flex_final', 'flex_lr_density_init',
                        'flex_lr_density_final', 'flex_sv_lam',
                        'flex_extra_epochs', 'flex_latent_samp_std',
                        'flex_latent_prior_lam', 'flex_latent_prior_pow',
                        'flex_latent_ext_init']
        self.lane = str(self.getAttributeValue('compute_lane'))

    def doRun3DFlexTraining(self):
        self._className = "flex_train"

        try:
            if not self.useQueueForSteps() and not self.useQueue():  # not using queue system
                gpusToUse = self.getGpuList()
            else:  # using queue system
                gpusToUse = False
        except Exception:
            gpusToUse = False

        protocolPrepare = self.input3DFlexDataPrepareProt.get()
        protocolMesh = self.input3DMeshFlexPrepareProt.get()
        varDataPrepJobParticles = str(protocolPrepare.run3DFlexDataPrepJob)
        ## varDataMeshJob = str(protocolMesh.run3DFlexMeshPrepJob)
        varDataMeshJob = str(protocolMesh.run3DFlexMeshPrep)
        input_group_connect = {"particles": "%s.particles" % varDataPrepJobParticles,
                               "flex_mesh": "%s.flex_mesh" % varDataMeshJob}
        params = {}

        for paramName in self._paramsName:
            if self.getAttributeValue(paramName) is not None:
                params[str(paramName)] = str(self.getAttributeValue(paramName))

        run3DTrainJob = enqueueJob(self._className,
                                   self.projectName.get(),
                                   self.workSpaceName.get(),
                                   str(params).replace('\'', '"'),
                                   str(input_group_connect).replace('\'','"'),
                                   self.lane, gpusToUse)

        self.run3DFlexTrainJob = String(run3DTrainJob.get())
        self.currenJob.set(self.run3DFlexTrainJob.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(),
                         self.run3DFlexTrainJob.get(),
                         "An error occurred in the 3D Flex Training process. "
                         "Please, go to cryoSPARC software for more "
                         "details.", self)
        clearIntermediateResults(self.projectName.get(),
                                 self.run3DFlexTrainJob.get())
