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
    Performs high-resolution flexible reconstruction of cryo-EM datasets using
    a previously trained 3DFlex model. The protocol refines structural detail
    while accounting for continuous conformational variability, producing
    flexible and rigid reference reconstructions together with validated
    half-maps suitable for downstream analysis.

    AI Generated:

    3D Flex Reconstruction (ProtCryoSparc3DFlexReconstruction) — User Manual
        Overview

        The 3D Flex Reconstruction protocol performs the final high-resolution
        reconstruction stage within the cryoSPARC 3DFlex framework. Its purpose
        is to recover detailed cryo-EM density maps while explicitly modeling
        continuous structural flexibility present in the particle population.
        Unlike conventional refinement approaches that assume a single rigid
        structure, this protocol incorporates the learned conformational model
        generated during 3DFlex training and uses it to improve reconstruction
        quality in regions affected by motion or heterogeneity.

        In biological applications, this approach is especially valuable for
        macromolecular assemblies that exhibit continuous domain movements,
        breathing motions, hinge rearrangements, or composational flexibility.
        Flexible refinement can reveal structural features that are blurred or
        weakened in standard consensus reconstructions, allowing more accurate
        interpretation of dynamic molecular processes.

        Inputs and Biological Context

        The protocol requires a previously trained 3DFlex model generated from
        compatible particle datasets. The training stage defines the latent
        representation of conformational variability, and this reconstruction
        stage uses that learned variability to produce high-resolution density
        maps.

        The quality of the reconstruction strongly depends on the biological
        relevance and stability of the training model. If the latent space does
        not capture meaningful structural motions, the resulting reconstruction
        may fail to improve map interpretability. For this reason, users should
        carefully inspect latent distributions and variability analyses before
        launching reconstruction.

        This protocol is particularly appropriate for systems such as ribosomes,
        membrane transporters, molecular motors, spliceosomes, viral assemblies,
        and multi-domain proteins where structural flexibility is expected to
        influence reconstruction quality.

        Flexible Versus Rigid Reconstruction

        One important feature of the protocol is the ability to generate both
        flexible and rigid reconstructions. The flexible reconstruction uses the
        learned conformational model to account for continuous motion during map
        estimation, while the rigid reconstruction serves as a conventional
        baseline generated from the same dataset.

        Comparing these two outputs can provide biologically meaningful insight.
        Improvements observed in the flexible map often indicate regions where
        conformational variability significantly affected the original density.
        Flexible refinement may sharpen mobile domains, improve connectivity,
        or recover secondary-structure detail that was previously obscured by
        motion averaging.

        In contrast, if the rigid and flexible reconstructions appear nearly
        identical, this may indicate that the dataset contains limited
        continuous heterogeneity or that the learned latent model does not
        capture substantial motions.

        Gold-Standard Validation

        The protocol supports independent half-map reconstruction for
        gold-standard Fourier Shell Correlation validation. This is essential
        for reliable resolution estimation and for minimizing overfitting during
        high-resolution refinement.

        In many cryo-EM workflows, preserving the original half-set assignment
        is recommended because it maintains consistency with upstream
        refinements. However, re-splitting the dataset may occasionally be
        useful when previous assignments are imbalanced or biologically
        unsuitable.

        The generated half-maps can later be used for FSC analysis, local
        resolution estimation, map sharpening, density modification, and atomic
        model validation.

        Reconstruction Optimization

        The reconstruction process relies on iterative numerical optimization to
        estimate the most consistent density maps under the flexible model. The
        number of optimization iterations controls the balance between runtime
        and reconstruction convergence.

        For most biological datasets, the default number of iterations provides
        stable and reliable results. Increasing the iteration count may improve
        refinement for particularly high-resolution datasets or very large
        complexes, although this also increases computational cost. Excessive
        optimization may sometimes amplify noise or overfit weak structural
        features, especially when particle quality is limited.

        Practical Interpretation of Results

        The flexible reconstruction output should be interpreted as a structural
        representation that integrates continuous conformational variability
        into the refinement process. Regions showing improved clarity often
        correspond to biologically meaningful motions that were previously
        averaged out in standard reconstruction approaches.

        Users should visually compare the flexible and rigid reconstructions,
        inspect local resolution differences, and evaluate whether improved
        density corresponds to plausible structural behavior. Flexible
        refinement is particularly informative when improvements occur in known
        dynamic regions such as ligand-binding domains, peripheral subunits,
        flexible linkers, or membrane-associated components.

        The protocol does not replace discrete classification methods but rather
        complements them by modeling smooth and continuous variability across
        the particle population.

        Outputs and Downstream Applications

        The protocol produces high-resolution flexible maps together with their
        associated half-maps. Optionally, rigid reference reconstructions are
        also generated for comparison. These outputs can be directly used for
        visualization, sharpening, local refinement assessment, atomic model
        building, flexible fitting, and structural interpretation.

        The half-maps are particularly important for downstream cryo-EM
        validation workflows because they allow independent FSC-based quality
        assessment and facilitate robust map interpretation.

        Practical Recommendations

        For most biological projects, users should first ensure that the
        upstream 3DFlex training stage produced a meaningful latent model with
        interpretable variability. Launching high-resolution reconstruction on
        poorly trained latent spaces is unlikely to improve map quality.

        Flexible refinement is most beneficial when datasets contain clear
        continuous motions rather than strongly discrete conformational states.
        When heterogeneity is dominated by a few well-separated structures,
        traditional 3D classification may still be the preferred approach.

        Users are encouraged to compare flexible and rigid reconstructions side
        by side and to interpret improvements cautiously within the biological
        context of the studied system.

        Final Perspective

        For modern cryo-EM studies focused on structural dynamics, 3D Flex
        Reconstruction provides a powerful strategy for recovering high-quality
        density information from heterogeneous datasets. By integrating learned
        conformational variability directly into high-resolution refinement, the
        protocol enables a more realistic representation of molecular motion and
        offers new opportunities for understanding dynamic biological systems.
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
