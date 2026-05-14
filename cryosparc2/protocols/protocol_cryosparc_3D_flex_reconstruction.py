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
from pwem.objects import Volume
from pyworkflow import BETA
from pyworkflow.protocol.params import (PointerParam,  BooleanParam)
from . import ProtCryosparcBase
from ..convert import *
from ..utils import *
from ..constants import *


class ProtCryoSparc3DFlexReconstruction(ProtCryosparcBase):
    """
    Takes in a checkpoint from training as well as prepared high-resolution
    particles and performs high-resolution refinement using L-BFGS under
    the 3DFlex model. This is the stage at which improvements to density
    in high-res regions are computed. Outputs two half-maps that can be used
    for FSC validation, sharpening, and other downstream tasks.
    """
    """
        Performs high-resolution flexible refinement reconstruction using a previously
        trained 3DFlex model in CryoSPARC. The protocol reconstructs flexible and,
        optionally, rigid density maps from prepared particle datasets by applying
        L-BFGS optimization under the learned deformation model. The resulting maps
        can be used for structural interpretation, FSC validation, sharpening, and
        downstream cryo-EM analysis workflows.

        AI Generated:

        3D Flex Reconstruction (ProtCryoSparc3DFlexReconstruction) — User Manual
            Overview

            The 3D Flex Reconstruction protocol performs high-resolution refinement
            using a trained CryoSPARC 3DFlex model together with previously prepared
            particle datasets. Its main purpose is to recover high-resolution density
            information while incorporating conformational variability learned during
            3DFlex training. Unlike conventional rigid refinements, this protocol
            models continuous structural motion and applies flexible deformation
            during reconstruction.

            In practical cryo-EM workflows, this protocol is commonly used after
            3DFlex Training to improve density quality in regions affected by
            conformational heterogeneity. The reconstruction stage generates both
            flexible and optional rigid reconstructions, enabling direct comparison
            between motion-aware refinement and standard rigid refinement approaches.

            Inputs and General Workflow

            The protocol requires as input a completed 3DFlex Training protocol,
            which provides both the trained deformation model and the associated
            prepared particles. The reconstruction process uses these inputs to
            generate high-resolution maps using L-BFGS optimization under the
            learned latent deformation space.

            During execution, CryoSPARC reconstructs two independent half-maps
            following gold-standard refinement principles. These half-maps are later
            combined into a final reconstruction and can be used for FSC validation,
            local resolution estimation, sharpening, and further post-processing.

            The protocol can also generate a rigid reconstruction using the same
            optimization framework. This rigid reconstruction serves as a biological
            and methodological baseline, allowing users to evaluate whether flexible
            refinement genuinely improves structural interpretation.

            Flexible and Rigid Reconstruction

            One of the central aspects of this protocol is the distinction between
            flexible and rigid reconstruction modes. Flexible reconstruction applies
            the learned deformation model to account for structural motion during
            refinement, which often improves density quality in dynamic regions of
            macromolecular complexes.

            The optional rigid reconstruction disables deformation modeling and
            instead reconstructs the structure using a traditional rigid framework.
            Comparing rigid and flexible maps is particularly useful when evaluating
            whether observed improvements correspond to biologically meaningful
            flexibility or potential overfitting.

            From a biological perspective, flexible refinement is especially valuable
            for complexes containing mobile domains, flexible linkers, or continuous
            conformational transitions that are difficult to resolve with standard
            refinement methods.

            Optimization Strategy

            Reconstruction is performed using the L-BFGS optimization algorithm.
            The parameter controlling the maximum number of iterations defines how
            extensively the optimizer refines each half-map during reconstruction.

            In most biological applications, the default number of iterations is
            sufficient to obtain stable and reliable reconstructions. Increasing the
            number of iterations may improve convergence in very high-resolution
            datasets or particularly large molecular assemblies, although it also
            increases computational cost and runtime.

            Excessively large optimization settings are not always beneficial and
            should be interpreted carefully, especially in datasets with limited
            signal-to-noise ratio or strong structural heterogeneity.

            Gold-Standard Half-Map Handling

            The protocol preserves or regenerates the gold-standard particle split
            used during refinement. Maintaining independent half datasets is critical
            for proper FSC validation and for avoiding overfitting during flexible
            refinement.

            Users may optionally force a new gold-standard split. This option is
            useful when input alignments contain inconsistent or imbalanced particle
            distributions between half-sets. Re-splitting ensures statistical
            independence between reconstructions and can improve validation
            reliability.

            Biologically, careful half-map management is essential because flexible
            refinement methods can otherwise introduce misleading high-resolution
            features if overfitting is not properly controlled.

            GPU and Computational Considerations

            The protocol supports GPU acceleration and is designed to integrate
            directly into CryoSPARC scheduling environments. GPU execution is highly
            recommended because flexible high-resolution reconstruction is
            computationally demanding.

            Depending on the execution environment, GPU resources may be selected
            automatically or provided through an external queue system. Runtime is
            influenced by particle count, box size, map resolution, and the number
            of optimization iterations.

            Outputs and Their Interpretation

            After completion, the protocol generates two main reconstruction outputs:
            a flexible reconstruction and, optionally, a rigid reconstruction. Each
            reconstruction includes its corresponding pair of gold-standard half-maps.

            The flexible reconstruction represents the primary scientific output and
            reflects density refinement under the learned deformation model. This map
            is generally expected to improve density continuity and structural detail
            in mobile regions compared to rigid refinement approaches.

            The rigid reconstruction provides a direct comparison reference. When the
            flexible map shows meaningful improvements over the rigid map, it
            suggests that the learned motion model successfully captures biologically
            relevant conformational variability.

            The generated half-maps can be used for FSC analysis, local resolution
            estimation, sharpening, map validation, and further structural modeling
            workflows.

            Practical Recommendations

            In routine cryo-EM practice, it is generally advisable to begin with the
            default reconstruction parameters and visually compare the flexible and
            rigid outputs. Improvements should be assessed carefully, particularly in
            flexible peripheral regions where overfitting risks are higher.

            Increasing the number of L-BFGS iterations may help in difficult
            high-resolution datasets, but excessively aggressive optimization should
            be avoided unless supported by clear validation metrics.

            When working with highly heterogeneous complexes, preserving proper
            gold-standard separation is especially important to ensure reliable FSC
            interpretation and biologically meaningful conclusions.

            Final Perspective

            For many cryo-EM studies, 3D Flex Reconstruction represents the stage
            where conformational variability becomes directly translated into
            high-resolution structural information. Rather than treating flexibility
            as noise, this protocol incorporates molecular motion into the refinement
            process itself.

            Successful interpretation depends not only on reconstruction quality, but
            also on careful biological validation, comparison against rigid
            refinement baselines, and rigorous assessment of map reliability through
            gold-standard FSC procedures.
        """
    _label = '3D flex reconstruction'
    _devStatus = BETA
    _protCompatibility = [V4_1_0, V4_1_1, V4_1_2, V4_2_0, V4_2_1, V4_3_1, V4_4_0, V4_4_1, V4_5_1,
                          V4_5_3, V4_6_0, V4_6_1, V4_6_2, V4_7_0, V4_7_1]

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
        form.addParam('input3DFlexTrainingProt', PointerParam,
                      pointerClass='ProtCryoSparc3DFlexTraining',
                      label="3D flex data prepare protocol",
                      important=True,
                      help='Particle stacks to use.')

        form.addParam('flex_do_noflex_recon', BooleanParam, default=True,
                      label="Do rigid reconstruction",
                      help='If True, the job will also do a rigid '
                           'reconstruction using the same L-BFGS reconstruction '
                           'method as is used for flexible refinement '
                           'reconstruction. This serves as a useful baseline '
                           'for comparisons.')

        form.addParam('flex_bfgs_num_iters', IntParam, default=20,
                      label="Max BFGS iterations",
                      help='The maximum number of L-BFGS iterations that will '
                           'be done during reconstruction of a half-map. '
                           'The default (20) works well in most cases but can '
                           'be increased for very high resolution '
                           'reconstruction or very large volumes potentially.')

        form.addParam('refine_gs_resplit', BooleanParam, default=False,
                      label="Force re-do GS split",
                      help='Force re-splitting the particles into two random '
                           'gold-standard halves. If this is not set, split '
                           'is preserved from input alignments (if connected).'
                           ' If the input alignments do not have equal '
                           'particles in each split, the job will issue a '
                           'warning but will continue.')

        # --------------[Compute settings]---------------------------
        form.addSection(label="Compute settings")
        addComputeSectionParams(form, allowMultipleGPUs=False, needGPU=True)

        """
        # job.param_add('flex_highres', "flex_force_restart", base_value=False,  title="Restart training",  param_type="boolean", desc="Force restart of training even if a model with trained checkpoint is connected.", hidden=True)
        
        job.param_add_section('compute_settings', title='Compute settings', desc='')
        # NB: app checks this param at queue time and tells command to no_check_inputs_ready in enqueue_job:
        job.param_add('compute_settings', "scheduler_no_check_inputs_ready", base_value=False,  title="Override scheduler",  param_type="boolean", desc="Force the scheduler to run this job even if connected inputs are not completed. For example, with this setting on, you can connect a running 3DFlex Training job output to this job and this job will run even though the training is still in progress. This allows visualization of in-progress results.")
        # job.param_add('compute_settings', "compute_use_ssd",        base_value=True,       title="Cache particle images on SSD",    param_type="boolean",   hidden=False,   advanced=False, desc='Use the SSD to cache particles. Speeds up processing significantly.')
        """

    def _insertAllSteps(self):
        self._defineFileNames()
        self._defineParamsName()
        self._initializeCryosparcProject()
        self._insertFunctionStep(self.reconstructionStep)
        self._insertFunctionStep(self.createOutputStep)

    def reconstructionStep(self):
        self.info(pwutils.yellowStr("3D Flex Reconstruction started..."))
        self.doRun3DFlexReconstruction()

    def createOutputStep(self):
        """
         Create the protocol output.  """
        self._initializeUtilsVariables()
        csOutputFolder = os.path.join(self.projectDir.get(),
                                      self.run3DFlexReconstructionJob.get())
        csOutputPattern = "%s%s" % (getOutputPreffix(self.projectName.get()),
                                    self.run3DFlexReconstructionJob.get())

        # Flex volume
        fnFlexVolName = csOutputPattern + "_flex_map.mrc"
        flexHalf1Name = csOutputPattern + "_flex_map_half_A.mrc"
        flexHalf2Name = csOutputPattern + "_flex_map_half_B.mrc"

        # No Flex volume
        fnNoFlexVolName = csOutputPattern + "_noflex_map.mrc"
        flexNoHalf1Name = csOutputPattern + "_noflex_map_half_A.mrc"
        flexNoHalf2Name = csOutputPattern + "_noflex_map_half_B.mrc"

        # Copy the CS output volume and half to extra folder
        copyFiles(csOutputFolder, self._getExtraPath(), files=[fnFlexVolName, flexHalf1Name, flexHalf2Name,
                                                               fnNoFlexVolName, flexNoHalf1Name, flexNoHalf2Name])

        fnVol = os.path.join(self._getExtraPath(), fnFlexVolName)
        half1 = os.path.join(self._getExtraPath(), flexHalf1Name)
        half2 = os.path.join(self._getExtraPath(), flexHalf2Name)

        flexVol = Volume()
        fixVolume([fnVol, half1, half2])
        flexVol.setFileName(fnVol)
        ccp4header = Ccp4Header(fnVol, readHeader=True)
        flexVol.setSamplingRate(ccp4header.getSampling()[0])
        flexVol.setHalfMaps([half1, half2])

        fnVol = os.path.join(self._getExtraPath(), fnNoFlexVolName)
        half1 = os.path.join(self._getExtraPath(), flexNoHalf1Name)
        half2 = os.path.join(self._getExtraPath(), flexNoHalf2Name)

        noFlexVol = Volume()
        fixVolume([fnVol, half1, half2])
        noFlexVol.setFileName(fnVol)
        ccp4header = Ccp4Header(fnVol, readHeader=True)
        noFlexVol.setSamplingRate(ccp4header.getSampling()[0])
        noFlexVol.setHalfMaps([half1, half2])

        self._defineOutputs(flexVolume=flexVol)
        self._defineOutputs(noFlexVolume=noFlexVol)

    def _defineParamsName(self):
        """ Define a list with 3D Flex Reconstruction parameters names"""
        self._paramsName = ['flex_do_noflex_recon', 'flex_bfgs_num_iters',
                            'refine_gs_resplit', 'compute_use_ssd']
        self.lane = str(self.getAttributeValue('compute_lane'))

    def doRun3DFlexReconstruction(self):
        self._className = "flex_highres"
        try:
            if not self.useQueueForSteps() and not self.useQueue():  # not using queue system
                gpusToUse = self.getGpuList()
            else:  # using queue system
                gpusToUse = False
        except Exception:
            gpusToUse = False

        protocolJobTraining = str(self.input3DFlexTrainingProt.get().run3DFlexTrainJob)
        input_group_connect = {"particles": "%s.particles" % protocolJobTraining,
                               "flex_model": "%s.flex_model" % protocolJobTraining}
        params = {}

        for paramName in self._paramsName:
            if self.getAttributeValue(paramName) is not None:
                params[str(paramName)] = str(self.getAttributeValue(paramName))

        run3DReconstructionJob = enqueueJob(self._className,
                                   self.projectName.get(),
                                   self.workSpaceName.get(),
                                   str(params).replace('\'', '"'),
                                   str(input_group_connect).replace('\'','"'),
                                   self.lane, gpusToUse)

        self.run3DFlexReconstructionJob = String(run3DReconstructionJob.get())
        self.currenJob.set(self.run3DFlexReconstructionJob.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(),
                         self.run3DFlexReconstructionJob.get(),
                         "An error occurred in the 3D Flex Reconstruction process. "
                         "Please, go to cryoSPARC software for more "
                         "details.", self)
        clearIntermediateResults(self.projectName.get(),
                                 self.run3DFlexReconstructionJob.get())
