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
import zipfile

from pwem import getMatchingFiles
from pwem.protocols import ProtFlexBase
from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam
from . import ProtCryosparcBase
from ..convert import *
from ..utils import *
from ..constants import *


class ProtCryoSparc3DFlexGenerator(ProtCryosparcBase, ProtFlexBase):
    """
    Generates continuous series of cryo-EM volumes from a trained 3DFlex model
    in order to visualize conformational variability and structural motions
    learned during flexible refinement.

    AI Generated:

    3D Flex Generator (ProtCryoSparc3DFlexGenerator) — User Manual
        Overview

        The 3D Flex Generator protocol produces animated or sequential
        representations of structural variability learned by a previously
        trained 3DFlex model. Its primary objective is to transform the
        abstract latent conformational space into interpretable density maps
        that reveal how the particle changes shape across continuous motions.

        In cryo-EM studies, many biologically important systems do not exist
        as a single rigid structure but instead fluctuate between multiple
        conformations. Traditional classification methods often separate only
        a limited number of discrete states, whereas 3DFlex attempts to model
        smooth and continuous structural transitions. The generator protocol
        provides a direct way to visualize these transitions and explore the
        dynamic behavior captured during training.

        Biological Interpretation of Generated Motion

        The generated volume series represents hypothetical structural
        trajectories inferred from experimental particle images. These
        trajectories can correspond to biologically meaningful motions such as
        domain rotations, hinge bending, channel opening, ligand-induced
        rearrangements, breathing motions, or coordinated subunit movements.

        Rather than producing unrelated independent maps, the protocol creates
        a continuous progression through conformational space. Consecutive
        frames therefore represent gradual structural changes instead of abrupt
        transitions between disconnected classes.

        For biological users, this approach can provide insight into molecular
        mechanisms that are difficult to observe using static reconstructions
        alone. Dynamic processes such as substrate transport, conformational
        activation, assembly transitions, or allosteric coupling can often
        become more interpretable when visualized as smooth structural motion.

        Inputs and General Workflow

        The protocol requires a previously trained 3DFlex model generated from
        the training stage of the workflow. The quality and interpretability
        of the generated series therefore depend strongly on the quality of
        the earlier training and mesh preparation steps.

        During execution, the protocol samples points across the learned latent
        conformational space and converts those positions into corresponding
        density maps. The resulting maps are organized as ordered series that
        describe motion along one or more learned flexibility dimensions.

        The protocol can also operate while model training is still in progress.
        This allows users to monitor the evolution of the deformation model and
        evaluate whether biologically meaningful motions are emerging before
        the final training stage is completed.

        Frame Sampling and Motion Resolution

        One of the most important user-controlled settings is the number of
        frames generated within each series. This parameter determines how
        finely the conformational trajectory is sampled.

        Smaller numbers of frames produce shorter and more simplified motion
        trajectories that are easier to inspect quickly. Larger numbers of
        frames generate smoother transitions and more continuous animations,
        which can improve interpretation of subtle structural rearrangements.

        In biological practice, moderate frame counts are often sufficient for
        identifying major conformational trends, whereas highly flexible systems
        with complex motions may benefit from denser sampling to capture gradual
        transitions more accurately.

        Because the latent space is continuous, the generated trajectories do
        not necessarily correspond to experimentally isolated biochemical states.
        Instead, they represent smooth interpolations within the learned
        deformation landscape. Biological interpretation should therefore focus
        on consistent structural trends rather than assuming that every frame
        corresponds to a stable physical intermediate.

        High-Resolution Flexible Visualization

        The protocol can optionally apply learned deformations to higher
        resolution density maps obtained from downstream reconstruction stages.
        This capability is especially valuable because it combines the dynamic
        information learned during flexible modeling with the improved structural
        detail available from high-resolution refinement.

        In practical cryo-EM workflows, this allows users to visualize flexible
        motions while preserving sharper local density features. Domain
        rearrangements, secondary structure displacements, and coordinated
        structural transitions may therefore become easier to interpret.

        However, biological caution remains important. Generated maps represent
        model-based deformations and should not automatically be interpreted as
        experimentally independent reconstructions. Validation against known
        biochemical, structural, or functional information is strongly advised.

        Outputs and Their Interpretation

        The protocol produces ordered series of reconstructed density maps
        corresponding to trajectories through latent conformational space.
        These outputs can be visualized as movies, morphing animations, or
        sequential structural states within standard cryo-EM visualization
        software.

        Generated trajectories are often most informative when examined
        interactively. Users should inspect whether motions appear physically
        plausible, structurally continuous, and biologically meaningful.
        Unrealistic distortions, disconnected density, or abrupt transitions
        may indicate insufficient training, overfitting, noisy data, or poorly
        constrained flexibility models.

        In many cases, multiple independent trajectories may reveal different
        modes of motion. These may correspond to separate biological processes,
        coupled domain rearrangements, or hierarchical conformational changes.

        Practical Recommendations

        For most biological applications, it is advisable to begin by generating
        moderate-length trajectories and visually inspecting the resulting
        motions before performing extensive interpretation. Large conformational
        changes are often easier to identify early in the analysis process.

        If the generated motions appear unstable or physically implausible,
        improving earlier stages of the workflow may help. Better masking,
        more appropriate rigidity settings, additional training, or refined
        mesh preparation can substantially improve the interpretability of
        generated trajectories.

        Flexible systems containing multiple articulated domains frequently
        produce especially informative visualizations because the learned latent
        space captures coordinated inter-domain movement. Smaller or highly
        rigid particles may instead show only subtle local variability.

        Final Perspective

        The 3D Flex Generator protocol transforms latent conformational models
        into visually interpretable structural trajectories that help bridge
        the gap between static cryo-EM reconstructions and dynamic biological
        processes. By revealing continuous molecular motion directly from
        experimental data, it provides an important framework for studying
        structural heterogeneity, conformational landscapes, and functional
        dynamics in cryo-EM research.
    """
    _label = '3D flex generator'
    _devStatus = BETA

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('flexTraining', PointerParam,
                      pointerClass='ProtCryoSparc3DFlexTraining',
                      label="3D flex training protocol",
                      important=True,
                      help='3DFlex model')

        form.addParam('flex_gen_num_pts', IntParam, default=41,
                      label='Nomber of frames per series',
                      help="Number of volumes to generate per series (unless custom latent components connected). "
                            "This should be an odd number to ensure that the zero position along each latent "
                            "dimension is sampled at the middle of the series")

        # --------------[Compute settings]---------------------------
        form.addSection(label="Compute settings")
        addComputeSectionParams(form, allowMultipleGPUs=False, needGPU=True)

    def _insertAllSteps(self):
        self._defineParamsName()
        self._initializeCryosparcProject()
        self._insertFunctionStep(self.generatorStep)
        self._insertFunctionStep(self.createOutputStep)

    def generatorStep(self):
        self.info(pwutils.yellowStr("3D Flex Generator started..."))
        self.doRun3DFlexGenerator()

    def createOutputStep(self):
        self._initializeUtilsVariables()
        self.info(pwutils.yellowStr("Creating the output..."))

        csOutputFolder = os.path.join(self.projectDir.get(),
                                      self.run3DGeneratorJob.get())

        pattern = csOutputFolder + '/%s_series_*.zip' % self.run3DGeneratorJob.get()
        csSeries = getMatchingFiles(pattern, True)
        fileNameList = [os.path.basename(fileName) for fileName in csSeries]
        copyFiles(csOutputFolder, self._getExtraPath(), files=fileNameList)
        outputFolder = self._getExtraPath()
        for fileName in fileNameList:
            filePath = os.path.join(outputFolder, fileName)
            with zipfile.ZipFile(filePath, 'r') as fileZip:
                extractFolder = os.path.join(outputFolder, os.path.splitext(fileName)[0])
                os.mkdir(extractFolder)
                fileZip.extractall(extractFolder)
            os.remove(filePath)

    def _validate(self):
        validateMsgs = cryosparcValidate()
        if not validateMsgs:
            flexTrainingProt = self.flexTraining.get()
            flexTrainingJob = getJob(flexTrainingProt.projectName.get(), flexTrainingProt.run3DFlexTrainJob.get())
            checkPoints = eval(flexTrainingJob[1])['output_result_groups']
            for output in checkPoints:
                if output['name'] == 'flex_model':
                    numItems = output['num_items']
                    if not numItems:
                        validateMsgs.append("There is no data to show yet")
                    break
        return validateMsgs

    def _defineParamsName(self):
        """ Define a list with 3D Flex Training parameters names"""
        self._paramsName = ['flex_gen_num_pts']
        self.lane = str(self.getAttributeValue('compute_lane'))

    def doRun3DFlexGenerator(self):
        self._className = "flex_generate"

        try:
            if not self.useQueueForSteps() and not self.useQueue():  # not using queue system
                gpusToUse = self.getGpuList()
            else:  # using queue system
                gpusToUse = False
        except Exception:
            gpusToUse = False

        protocolTraining = self.flexTraining.get().run3DFlexTrainJob.get()
        input_group_connect = {"flex_model": "%s.flex_model" % protocolTraining}
        params = {}

        for paramName in self._paramsName:
            if self.getAttributeValue(paramName) is not None:
                params[str(paramName)] = str(self.getAttributeValue(paramName))

        run3DGeneratorJob = enqueueJob(self._className,
                                   self.projectName.get(),
                                   self.workSpaceName.get(),
                                   str(params).replace('\'', '"'),
                                   str(input_group_connect).replace('\'','"'),
                                   self.lane, gpusToUse)

        self.run3DGeneratorJob = String(run3DGeneratorJob.get())
        self.currenJob.set(self.run3DGeneratorJob.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(),
                         self.run3DGeneratorJob.get(),
                         "An error occurred in the 3D Flex Generator process. "
                         "Please, go to cryoSPARC software for more "
                         "details.", self)
        clearIntermediateResults(self.projectName.get(),
                                 self.run3DGeneratorJob.get())
