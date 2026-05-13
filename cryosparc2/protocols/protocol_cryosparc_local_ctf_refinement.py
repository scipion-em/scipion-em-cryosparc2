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
    Performs local contrast transfer function refinement on cryo-EM particle
    datasets using a reference reconstruction and an associated mask. The
    protocol improves per-particle optical parameter estimation by refining
    local defocus values against a known 3D structure, helping increase the
    accuracy and interpretability of downstream refinements and reconstructions.

    AI Generated:

    Local CTF Refinement (ProtCryoSparcLocalCtfRefinement) — User Manual
        Overview

        The Local CTF Refinement protocol refines optical parameters for
        individual particles in a cryo-EM dataset using a high-quality 3D
        reference map. Its primary purpose is to improve the consistency
        between experimental particle images and the reconstructed density by
        correcting local defocus variations that may differ from particle to
        particle across a micrograph.

        In practical cryo-EM workflows, local refinement of CTF parameters is
        particularly important for achieving high-resolution reconstructions.
        Even when an accurate global CTF estimation has already been performed,
        local differences caused by uneven ice thickness, sample tilt,
        beam-induced movement, or optical imperfections can reduce map quality.
        Refining these parameters at the particle level frequently improves
        high-resolution signal and produces sharper structural features.

        Inputs and Biological Context

        The protocol requires a set of particles with known projection
        alignments together with a reference volume representing the target
        structure. The reference map acts as the structural guide used to
        estimate improved local optical parameters for every particle image.

        A soft mask is also required to define the structural region that
        should drive the refinement. Biologically, the mask should include the
        stable molecular density while excluding solvent and highly flexible
        regions whenever possible. Appropriate masking improves refinement
        stability because the algorithm focuses on interpretable structural
        signal rather than noise or disordered regions.

        For macromolecular complexes with flexible peripheral domains, membrane
        regions, or compositional variability, careful masking becomes
        especially important. Overly broad masks may introduce solvent noise,
        whereas masks that are too restrictive may exclude useful signal.

        Refinement Strategy

        The protocol estimates local defocus corrections independently for
        individual particles while preserving the previously determined
        particle orientations. This strategy is especially valuable after
        homogeneous or non-uniform refinement steps, where orientation accuracy
        is already high and optical inaccuracies become one of the dominant
        limitations for further resolution improvement.

        The refinement can operate across a specified resolution range. Lower
        resolution limits define the broad structural information used during
        fitting, while higher resolution limits determine how aggressively the
        refinement attempts to capture fine detail. In many biological
        datasets, conservative resolution limits provide more stable behavior,
        especially for smaller particles or noisy datasets.

        The defocus search range determines how broadly the protocol explores
        possible local defocus deviations around the original estimates. Wider
        searches may help when substantial local variation is expected, but
        excessively large ranges can increase runtime and reduce robustness.

        Resolution Considerations

        High-resolution refinement is only meaningful when the reference map
        contains reliable structural information. Using an over-refined or
        inaccurate reference may bias the refinement toward non-physical
        solutions. For this reason, users typically perform local CTF
        refinement only after obtaining a stable and biologically reliable
        reconstruction.

        When half-maps are available, they can provide more robust estimation
        of the effective resolution limit and reduce the risk of overfitting.
        This is especially important in near-atomic resolution workflows where
        subtle optical corrections can significantly influence map quality.

        Outputs and Interpretation

        The protocol produces an updated particle dataset containing refined
        local optical parameters while preserving the original particle
        identities and alignments. These refined particles are generally used
        as improved inputs for subsequent homogeneous refinement,
        non-uniform refinement, local refinement, or high-resolution map
        reconstruction.

        Biologically, successful local CTF refinement often manifests as
        sharper secondary structure features, improved side-chain density,
        clearer ligand visibility, and better-resolved flexible regions.
        However, improvements are typically incremental rather than dramatic,
        and their magnitude depends strongly on data quality and microscope
        stability.

        Practical Recommendations

        In routine cryo-EM processing, local CTF refinement is usually most
        beneficial after obtaining a reasonably high-resolution consensus map.
        Early application to poorly aligned datasets may provide limited
        improvement because orientation inaccuracies dominate the error model.

        Users should begin with moderate defocus search ranges and conservative
        resolution limits, especially for small particles or heterogeneous
        datasets. Excessively aggressive settings can increase computational
        cost without producing biologically meaningful gains.

        For flexible complexes, selecting a stable core region with an
        appropriate mask generally produces more reliable refinements than
        attempting to include highly mobile domains. Visual inspection of the
        resulting maps remains essential to confirm that refinement improves
        true structural detail rather than amplifying noise.

        Final Perspective

        Local CTF refinement represents an important optimization stage in
        modern high-resolution cryo-EM workflows. By correcting subtle
        particle-specific optical differences, the protocol helps maximize the
        structural information recoverable from experimental images. Careful
        choice of the reference map, biologically meaningful masking, and
        realistic refinement limits are the main factors that determine whether
        the refinement leads to robust and interpretable structural
        improvements.
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