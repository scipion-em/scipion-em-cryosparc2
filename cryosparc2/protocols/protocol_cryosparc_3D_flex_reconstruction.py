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
    _label = '3D flex reconstruction'
    _devStatus = BETA
    _protCompatibility = [V4_1_0, V4_1_1, V4_1_2, V4_2_0, V4_2_1, V4_3_1, V4_4_0, V4_4_1, V4_5_1,
                          V4_5_3, V4_6_0, V4_6_1, V4_6_2]

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
            gpusToUse = self.getGpuList()
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
