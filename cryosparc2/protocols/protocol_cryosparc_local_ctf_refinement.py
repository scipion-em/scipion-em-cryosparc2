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

from pwem import ALIGN_PROJ
from pwem.protocols import ProtParticles
import pyworkflow.utils as pwutils
from pyworkflow.object import String
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        LEVEL_ADVANCED, IntParam, Positive)

from .protocol_base import ProtCryosparcBase
from .. import RELIONCOLUMNS
from ..convert import (convertCs2Star, createItemMatrix,
                       setCryosparcAttributes)
from ..utils import (addComputeSectionParams, cryosparcValidate, gpusValidate,
                     enqueueJob, waitForCryosparc, copyFiles)


class ProtCryoSparcLocalCtfRefinement(ProtCryosparcBase, ProtParticles):
    """
    Wrapper protocol for the Cryosparc's per-particle Local CTF refinement.
    Performs per-particle defocus estimation for each particle in a dataset,
    against a given 3D reference structure.
    """

    """
        Local CTF Refinement (ProtCryoSparcLocalCtfRefinement) — User Manual

        Overview

        The Local CTF Refinement protocol performs per-particle Contrast Transfer
        Function (CTF) refinement using cryoSPARC algorithms integrated within the
        Scipion framework. Its primary objective is to improve the estimation of
        local defocus values for individual particles against a reference 3D map,
        allowing more accurate reconstruction and higher-quality structural detail
        recovery in cryo-electron microscopy workflows.

        In practical biological applications, this protocol is especially valuable
        for high-resolution datasets where subtle optical inaccuracies may limit
        the interpretability of reconstructed maps. By refining defocus values at
        the particle level, the protocol helps recover high-frequency information,
        improves density sharpness, and enhances the visibility of fine structural
        features such as side chains, secondary structure elements, or flexible
        regions.

        Inputs and General Workflow

        The protocol requires a set of aligned particles, a reference volume, and
        a soft mask defining the biologically relevant region of the structure.
        The input particles must already contain projection alignment information,
        since the refinement process relies on previously estimated orientations
        to optimize local optical parameters.

        The reference volume acts as the structural model used during refinement,
        while the mask restricts optimization to meaningful molecular regions.
        This prevents solvent noise or flexible peripheral densities from biasing
        the defocus estimation process.

        During execution, the workflow initializes a cryoSPARC refinement job,
        converts input metadata into compatible formats, launches the local CTF
        optimization procedure, and finally imports the refined particle metadata
        back into the Scipion environment.

        If the reference structure contains half maps, the protocol automatically
        connects them to the refinement process in order to estimate resolution
        limits more robustly through FSC-based criteria. This improves refinement
        stability and reduces the risk of overfitting high-frequency noise.

        Masking and Biological Relevance

        Masking is one of the most biologically important aspects of local CTF
        refinement because it determines which structural regions contribute to
        the optimization process. The protocol requires a soft mask either
        provided explicitly by the user or inherited from previous refinement
        protocols.

        In biological datasets containing flexible domains, membrane regions,
        disordered appendages, or solvent noise, appropriate masking becomes
        essential for obtaining reliable defocus estimates. The mask should focus
        on structurally conserved and well-resolved regions while excluding noisy
        or highly mobile areas that may destabilize the refinement.

        Poor masking strategies may lead to inaccurate local CTF estimation,
        reduced map quality, or unstable refinement behavior. For this reason,
        carefully designed masks are particularly important in heterogeneous
        macromolecular assemblies or flexible protein complexes.

        Refinement Parameters and Optimization Strategy

        The protocol provides several advanced parameters controlling the local
        refinement behavior. Users can define the refinement box size, resolution
        limits, plotting behavior, and defocus search range used during the
        optimization process.

        The refinement box size determines the reconstruction region used during
        local optimization. Smaller boxes reduce computational cost, whereas larger
        boxes preserve more structural context. In most standard workflows, the
        default particle box size is sufficient unless sampling differences or
        scaling issues are suspected.

        The minimum and maximum fit resolutions define the frequency range used
        during refinement. Lower resolutions emphasize global structural features,
        while higher resolutions allow the protocol to exploit fine structural
        details when the dataset quality supports it.

        The defocus search range determines how broadly the protocol explores
        possible corrections around the initial defocus estimation. Wider search
        ranges may help when initial CTF parameters are inaccurate, although very
        large ranges may increase computational time and reduce optimization
        stability.

        GPU Integration and Execution

        The protocol is designed for GPU-accelerated execution within cryoSPARC
        environments. During runtime, it automatically determines whether GPU
        resources should be allocated directly or managed through an external queue
        system depending on the execution environment.

        Once the refinement job is launched, the workflow continuously monitors
        the cryoSPARC execution state until completion. If errors occur during the
        refinement process, users are encouraged to inspect the corresponding
        cryoSPARC job logs for detailed diagnostic information.

        Outputs and Their Interpretation

        After execution, the protocol generates an updated particle set containing
        refined local CTF information and updated alignment metadata. The output
        particles preserve the identity of the original dataset while incorporating
        improved optical parameter estimation.

        Biologically, these refined particles often produce sharper reconstructions,
        improved local resolution, and better visualization of structurally relevant
        details. High-resolution refinements particularly benefit from accurate
        local defocus correction because small optical inaccuracies can strongly
        affect map interpretability.

        The protocol also preserves transformation information and cryoSPARC
        metadata, enabling seamless integration with subsequent refinement,
        classification, or reconstruction workflows within Scipion.

        Practical Recommendations

        In most cryo-EM processing pipelines, local CTF refinement is most useful
        after obtaining a stable consensus refinement. Applying local refinement
        too early in the workflow may provide limited improvements because particle
        orientations and structural references may still be unstable.

        For high-resolution biological studies, careful mask design and appropriate
        resolution limits are typically the most important factors affecting the
        quality of the refinement. Flexible complexes usually benefit from masks
        focused on the rigid structural core.

        Datasets collected with tilted acquisition schemes, local defocus
        variability, or subtle optical distortions are particularly good
        candidates for this protocol. In contrast, highly heterogeneous or
        low-resolution datasets may show smaller refinement gains.

        Final Perspective

        Local CTF refinement represents an important optimization stage in modern
        cryo-EM workflows because it improves the accuracy of particle-specific
        optical parameter estimation. Although technically focused on defocus
        correction, its biological impact is substantial, directly influencing
        map sharpness, reconstruction quality, and the interpretability of
        macromolecular structures at high resolution.
        """
    _label = 'local ctf refinement'
    _className = "ctf_refine_local"

    def _initialize(self):
        self._createFilenameTemplates()

    def _createFilenameTemplates(self):
        """ Centralize how files are called. """
        myDict = {
            'input_particles': self._getTmpPath('input_particles.star'),
            'out_particles': self._getExtraPath('output_particle.star')
        }
        self._updateFilenamesDict(myDict)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      label="Input particles", important=True,
                      help='Provide a set of particles for local '
                           'CTF refinement.')
        form.addParam('refVolume', PointerParam, pointerClass='Volume',
                      important=True,
                      label="Input volume",
                      help='Provide a reference volume for local '
                           'CTF refinement.')
        form.addParam('refMask', PointerParam, pointerClass='VolumeMask',
                      label='Mask to be applied to this map',
                      important=True,
                      help="Provide a soft mask. if mask is present, use that, "
                           "otherwise use mask_refine if present, otherwise "
                           "fail")

        # -----------[Global CTF Refinement]------------------------
        form.addSection(label="Local CTF Refinement")
        form.addParam('crl_N', FloatParam, default=None,
                      allowsNull=True,
                      expertLevel=LEVEL_ADVANCED,
                      label='Refinement box size (Voxels)',
                      help='Size of reconstruction/image to use for refinement. '
                           'Blank means to use the particle box size '
                           '(upsampling input maps as needed)')
        form.addParam('crl_num_plots', IntParam, default=50,
                      validators=[Positive],
                      label='Num. groups to plot',
                      help='Number of exposure groups to make plots for. '
                           'After this many, stop plotting to save time.')

        form.addParam('crl_min_res_A', IntParam, default=20,
                      validators=[Positive],
                      label='Minimum Fit Res (A)',
                      help='The minimum resolution to use during refinement of '
                           'image aberrations.')

        form.addParam('crl_max_res_A', FloatParam, default=None,
                      label='Maximum Fit Res (A)',
                      expertLevel=LEVEL_ADVANCED,
                      allowsNull=True,
                      help='The maximum resolution to use during refinement of '
                           'image aberrations. If None, use input half-maps '
                           'to compute FSC and set max to FSC=0.5')

        form.addParam('crl_df_range', IntParam, default=2000,
                      label='Defocus Search Range (A +/-)',
                      help='Defocus search range in Angstroms, searching both '
                           'above and below the input defocus by this amount')

        # --------------[Compute settings]---------------------------
        form.addSection(label="Compute settings")
        addComputeSectionParams(form, allowMultipleGPUs=False)

    # --------------------------- INSERT steps functions -----------------------

    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._defineParamsName()
        self._initializeCryosparcProject()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.processStep)
        self._insertFunctionStep(self.createOutputStep)

    # -------------------------- UTILS functions ------------------------------

    def _defineParamsName(self):
        """ Define a list with all protocol parameters names"""
        self._paramsName = ['crl_N',
                            'crl_num_plots',
                            'crl_min_res_A',
                            'crl_max_res_A',
                            'crl_df_range',
                            'compute_use_ssd']
        self.lane = str(self.getAttributeValue('compute_lane'))

    def _getInputMask(self):
        if self.refMask.get() is not None:
            return self.refMask.get()
        else:
            inputProtocolMask = self._getInputPostProcessProtocol().refMask.get()
            if inputProtocolMask is not None:
                return inputProtocolMask

        return None

    # --------------------------- STEPS functions ------------------------------
    def processStep(self):
        self.info(pwutils.yellowStr("Local Ctf Refinement started..."))
        self.doLocalCtfRefinement()

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to Relion file
        """
        self._initializeUtilsVariables()
        outputStarFn = self._getFileName('out_particles')
        csOutputFolder = os.path.join(self.projectDir.get(),
                                      self.runLocalCtfRefinement.get())
        csFileName = "particles.cs"

        # Copy the CS output particles to extra folder
        copyFiles(csOutputFolder, self._getExtraPath(), files=[csFileName])

        csFile = os.path.join(self._getExtraPath(), csFileName)

        argsList = [csFile, outputStarFn]

        convertCs2Star(argsList)

        imgSet = self._getInputParticles()

        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(imgSet)
        self._fillDataFromIter(outImgSet)

        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(imgSet, outImgSet)

    def _fillDataFromIter(self, imgSet):
        outImgsFn = 'particles@' + self._getFileName('out_particles')
        imgSet.setAlignmentProj()
        imgSet.copyItems(self._getInputParticles(),
                         updateItemCallback=self._createItemMatrix,
                         itemDataIterator=emtable.Table.iterRows(outImgsFn))

    def _createItemMatrix(self, particle, row):
        createItemMatrix(particle, row, align=ALIGN_PROJ)
        setCryosparcAttributes(particle, row,
                               RELIONCOLUMNS.rlnRandomSubset.value)

    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        """ Should be overwritten in subclasses to
            return summary message for NORMAL EXECUTION.
        """
        validateMsgs = cryosparcValidate()
        if not validateMsgs:
            validateMsgs = gpusValidate(self.getGpuList(), checkSingleGPU=True)
            if not validateMsgs:
                self._validateDim(self._getInputParticles(),
                                  self._getInputVolume(),
                                  validateMsgs, 'Input particles',
                                  'Input volume')
        return validateMsgs

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputParticles'):
            summary.append("Output Particles not ready yet.")
        else:
            summary.append("Input Particles: %s" %
                           self.getObjectTag('inputParticles'))
            summary.append("Reference Mask: %s" %
                           self.getObjectTag('refMask'))
            summary.append("--------------------------------------------------")
            summary.append("Output particles %s" %
                           self.getObjectTag('outputParticles'))
        return summary

    def doLocalCtfRefinement(self):
        """
         :return:
         """
        input_group_connect = {"particles": self.particles.get(),
                              "volume": self.volume.get(),
                              "mask": self.mask.get()}

        input_result_connect = None
        if self._getInputVolume().hasHalfMaps():
            input_result_connect = {"volume.0.map_half_A": self.importVolumeHalfA.get(),
                                    "volume.0.map_half_B": self.importVolumeHalfB.get()}

        params = {}

        for paramName in self._paramsName:
            if (paramName != 'crl_max_res_A' and paramName != 'crl_N'):
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif self.crl_max_res_A.get() is not None:
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif self.crl_N.get() is not None:
                params[str(paramName)] = str(self.getAttributeValue(paramName))

        # Determinate the GPUs to use (in dependence of
        # the cryosparc version)
        try:
            if not self.useQueueForSteps() and not self.useQueue():  # not using queue system
                gpusToUse = self.getGpuList()
            else:  # using queue system
                gpusToUse = False
        except Exception:
            gpusToUse = False

        runLocalCtfRefinementJob = enqueueJob(self._className, self.projectName.get(),
                                                self.workSpaceName.get(),
                                                str(params).replace('\'', '"'),
                                                str(input_group_connect).replace('\'', '"'),
                                                self.lane, gpusToUse,
                                                result_connect=input_result_connect)

        self.runLocalCtfRefinement = String(runLocalCtfRefinementJob.get())
        self.currenJob.set(runLocalCtfRefinementJob.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(), self.runLocalCtfRefinement.get(),
                         "An error occurred in the particles subtraction process. "
                         "Please, go to cryoSPARC software for more "
                         "details.", self)