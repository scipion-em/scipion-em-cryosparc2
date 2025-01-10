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
    _label = '3D flex training'
    _devStatus = BETA
    _protCompatibility = [V4_1_0, V4_1_1, V4_1_2, V4_2_0, V4_2_1, V4_3_1, V4_4_0, V4_4_1, V4_5_1,
                          V4_5_3, V4_6_0, V4_6_1, V4_6_2]
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
            gpusToUse = self.getGpuList()
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
