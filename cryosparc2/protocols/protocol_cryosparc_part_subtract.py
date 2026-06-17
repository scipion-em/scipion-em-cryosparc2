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
from pwem.protocols import ProtOperateParticles

import pyworkflow.utils as pwutils
from pyworkflow.object import String
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        LEVEL_ADVANCED, Positive, BooleanParam)

from .protocol_base import ProtCryosparcBase
from ..convert import (convertCs2Star, cryosparcToLocation, rowToCtfModel, rowToAlignment)
from ..utils import (addComputeSectionParams, addPreprocessLaneParam, calculateNewSamplingRate,
                     cryosparcValidate, gpusValidate, enqueueJob,
                     waitForCryosparc, clearIntermediateResults, copyFiles)
from ..constants import *


class ProtCryoSparcSubtract(ProtCryosparcBase, ProtOperateParticles):
    """ Signal subtraction protocol of cryoSPARC.
        Subtract projections of a masked volume from particles.
        """
    """
    Particle Subtraction (ProtCryoSparcSubtract) — User Manual

    Overview

    The Particle Subtraction protocol removes selected structural density from
    experimental cryo-EM particles by subtracting projections generated from a
    reference volume. This strategy is commonly used in focused cryo-EM analysis
    when the goal is to isolate specific regions of a macromolecular complex,
    improve the visualization of flexible domains, or perform localized refinement
    on regions that would otherwise be dominated by larger structural features.

    From a biological perspective, particle subtraction is especially valuable in
    studies involving conformational heterogeneity, multi-domain assemblies,
    membrane-associated complexes, or protein–ligand interactions. By
    computationally removing stable or dominant regions of the structure, the
    protocol enables downstream analysis to focus on smaller or more dynamic
    components that are difficult to resolve in conventional refinements.

    The protocol is implemented as a cryoSPARC-based workflow integrated into
    Scipion, handling both the execution of the subtraction job and the conversion
    of outputs into Scipion-compatible particle datasets.

    Inputs and General Workflow

    The protocol requires a set of aligned experimental particles, a reference 3D
    volume, and a subtraction mask. The input particles must already contain
    projection alignment information because the subtraction process depends on
    accurately projecting the reference volume into the orientation of each
    particle image.

    The reference volume represents the density to be projected and removed from
    the particles. In practice, this volume should originate from the same dataset
    and ideally from the same refinement process as the particles themselves.
    Maintaining consistent greyscale normalization and orientation accuracy is
    biologically important because subtraction quality strongly depends on the
    correspondence between the experimental particles and the projected reference
    density.

    The subtraction mask defines which regions of the reference volume will be
    removed. Regions inside the mask are retained for subtraction, while regions
    outside are ignored. Biologically, this means the mask should include the
    structural component intended for removal rather than the region of interest
    that should remain in the particles. This distinction is critical because
    incorrect masking can unintentionally eliminate meaningful density or leave
    strong residual signal that interferes with downstream focused analysis.

    Internally, the protocol prepares cryoSPARC-compatible inputs, launches the
    subtraction job, waits for completion, and converts the resulting particle
    stack into a RELION/Scipion-compatible STAR representation. The transformed
    particles preserve updated alignment and CTF metadata after subtraction.

    Masking and Focused Signal Removal

    Masking is the central biological component of this protocol because it
    determines which structural signal will be computationally removed from the
    particles. The protocol requires a soft mask, typically generated around the
    stable or dominant region of the complex.

    In practical cryo-EM workflows, subtraction is frequently applied to remove
    rigid cores from flexible assemblies, subtract detergent micelles from membrane
    proteins, isolate ligand-binding regions, or separate subcomplexes in large
    molecular machines. A well-designed mask generally follows the boundaries of
    the region intended for subtraction while avoiding sharp edges that may
    introduce Fourier artifacts.

    The protocol also supports optional mask thresholding and hole filling.
    Thresholding allows binarization and expansion of the mask before subtraction,
    while hole filling improves continuity in disconnected regions. These
    operations are particularly useful when masks contain fragmented density or
    internal cavities that could produce unstable subtraction behavior.

    From a biological interpretation standpoint, over-aggressive masking may erase
    relevant signal, whereas overly permissive masks may fail to sufficiently
    isolate the region of interest. Careful visual inspection of the subtraction
    mask is therefore essential before large-scale processing.

    Windowing and Signal Scaling

    The subtraction procedure includes inner and outer reference window radii that
    define how the input particles are spatially windowed before subtraction.
    These parameters control the transition region between preserved and attenuated
    density and influence the stability of scaling during subtraction.

    The protocol additionally supports a premultiplier-based scaling approach that
    improves consistency between the projected reference and experimental particle
    intensities. In most standard workflows this option remains enabled because it
    stabilizes subtraction quality and reduces mismatches in signal amplitude.

    An optional low-pass filter can also be applied to the input reference
    structure before subtraction. This is biologically useful when the reference
    map contains high-resolution noise or overfitted features that should not
    propagate into the subtraction projections. Filtering the reference volume
    often improves subtraction robustness, particularly when working with flexible
    regions or intermediate-resolution reconstructions.

    Gold-Standard Subtraction and Half Maps

    The protocol supports the use of half maps during subtraction in order to
    preserve gold-standard refinement conditions. When enabled, each particle
    subset is subtracted using the corresponding independent half map generated
    during refinement.

    This strategy is biologically and computationally important because it
    minimizes the risk of information leakage between independently refined
    particle halves. Preserving gold-standard conditions is especially critical
    when the subtracted particles will later undergo focused refinement or local
    classification intended for high-resolution structural interpretation.

    If half maps are unavailable, the protocol can still perform subtraction using
    a single consensus map, although this may compromise strict gold-standard
    separation.

    Execution and Data Conversion

    After parameter definition, the protocol prepares the cryoSPARC job
    configuration and determines GPU usage depending on whether the workflow is
    executed locally or through a queue system. The subtraction job is then
    submitted directly to cryoSPARC.

    Once processing is complete, the protocol automatically retrieves the generated
    .cs particle file, converts it into STAR format, updates particle metadata,
    recalculates sampling information when needed, and reconstructs the transformed
    particle set inside Scipion.

    The updated particles preserve alignment transformations, CTF information, and
    image locations, allowing seamless continuation into downstream cryo-EM
    workflows such as focused classification, local refinement, variability
    analysis, or masked reconstruction.

    Outputs and Their Interpretation

    The main output of the protocol is a new set of subtracted particles in which
    the masked reference density has been computationally removed. These particles
    remain aligned and fully compatible with subsequent cryo-EM refinement and
    classification protocols.

    Biologically, the resulting particles should contain enhanced relative signal
    for the remaining structural regions after subtraction. This often improves
    the detectability of flexible domains, weakly occupied ligands, peripheral
    subunits, or compositional variability that would otherwise remain obscured by
    dominant density contributions.

    The protocol summary additionally reports the input particles, reference
    volume, subtraction mask, and the selected subtraction window parameters,
    allowing traceability and reproducibility within complex Scipion workflows.

    Practical Recommendations

    In routine cryo-EM practice, successful particle subtraction depends primarily
    on the quality of the subtraction mask and the consistency between the
    reference volume and the experimental particles. The most reliable results are
    generally obtained when the reference volume originates from the same
    refinement pipeline and shares the same normalization and orientation
    conventions as the particle dataset.

    For flexible assemblies, it is often beneficial to subtract only the most
    stable structural core while preserving regions expected to exhibit
    conformational variability. Applying moderate low-pass filtering to the
    reference volume can further reduce subtraction artifacts caused by
    high-frequency noise.

    When preparing particles for focused refinement, maintaining gold-standard
    half-map subtraction is strongly recommended to avoid artificial resolution
    inflation. Visual inspection of several representative subtracted particles is
    also advisable to verify that the intended signal has been removed without
    introducing strong residual artifacts.

    Final Perspective

    Particle subtraction is not simply a preprocessing step but a biologically
    meaningful strategy for simplifying complex cryo-EM datasets. By selectively
    removing dominant structural features, the protocol enables focused analysis
    of regions that are otherwise difficult to interpret due to flexibility,
    compositional heterogeneity, or weak occupancy.

    Careful definition of the subtraction mask, preservation of gold-standard
    refinement conditions, and consistent reference scaling are the key elements
    for obtaining biologically reliable subtracted particles suitable for
    downstream structural analysis.
    """

    _label = 'subtract projection'
    _className = "particle_subtract"

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
                      help='Select the experimental particles.')
        form.addParam('refVolume', PointerParam, pointerClass='Volume',
                      label="Input map to be projected",
                      important=True,
                      help='Provide the input volume that will be used to '
                           'calculate projections, which will be subtracted '
                           'from the experimental particles. Make sure this '
                           'map was calculated by RELION from the same '
                           'particles as above, and preferably with those '
                           'orientations, as it is crucial that the absolute '
                           'greyscale is the same as in the experimental '
                           'particles.')
        form.addParam('refMask', PointerParam, pointerClass='VolumeMask',
                      label='Mask to be applied to this map',
                      important=True,
                      allowsNull=False,
                      help="Provide a soft mask where the protein density "
                           "you wish to subtract from the experimental "
                           "particles is white (1) and the rest of the "
                           "protein and the solvent is black (0). "
                           "That is: *the mask should INCLUDE the part of the "
                           "volume that you wish to SUBTRACT.*")

        # -----------[Particles Subtraction]------------------------
        form.addSection(label="Particle Subtraction")

        form.addParam('inner_radius', FloatParam, default=0.85,
                      validators=[Positive],
                      label='Inner radius of reference window',
                      help='Inner radius of the windowing applied to the '
                           'particles used to generate the input structure.')

        form.addParam('outer_radius', FloatParam, default=0.99,
                      validators=[Positive],
                      label='Outer radius of reference window',
                      help='Outer radius of the windowing applied to the '
                           'particles used to generate the input structure')

        form.addParam('use_premult', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Use premultiplier for scaling",
                      help="Premultiplier to scale initial real space mode")

        form.addParam('use_halfmaps', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Use halfmaps for Gold-Standard Subtraction",
                      help="Subtract each halfmap from particles used to "
                           "generate it. Disabling this parameter breaks the "
                           "assumptions of gold-standard FSC calculation.")

        form.addParam('lpf_volume', FloatParam, default=None,
                      allowsNull=True,
                      expertLevel=LEVEL_ADVANCED,
                      label='Low-pass Filter Input Structure (A)',
                      help='Apply a lowpass filter to the specified reoslution '
                           'in Angstroms to the input volume before subtraction.'
                           'Leave None to ignore a lowpass filter')

        form.addParam('mask_threshold', FloatParam, default=None,
                      allowsNull=True,
                      expertLevel=LEVEL_ADVANCED,
                      label='Mask threshold',
                      help='The threshold of binarization of the mask. Must be '
                           'set to dilate or pad. Leave None to skip mask processing.')

        form.addParam('mask_fill_holes', BooleanParam, default=False,
                      label="Fill holes",
                      help="Fill the holes in the binarized mask")

        # --------------[Compute settings]---------------------------
        form.addSection(label="Compute settings")
        addComputeSectionParams(form, allowMultipleGPUs=False)
        addPreprocessLaneParam(form)

    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._createFilenameTemplates()
        self._defineParamsName()
        self._initializeCryosparcProject()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.processStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions ------------------------------
    def processStep(self):
        self.info(pwutils.yellowStr("Particles Subtraction started..."))
        self.doPartStract()

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to Relion file
        """
        self._initializeUtilsVariables()
        outputStarFn = self._getFileName('out_particles')
        csOutputFolder = os.path.join(self.projectDir.get(),
                                      self.runPartStract.get())
        csFileName = "subtracted_particles.cs"

        # Create the output folder
        copyFiles(csOutputFolder, os.path.join(self._getExtraPath(),
                                               self.runPartStract.get()))

        csFile = os.path.join(self._getExtraPath(), self.runPartStract.get(),
                              csFileName)
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
        imgSet.copyItems(self._getInputParticles(),
                         updateItemCallback=self._updateItem,
                         itemDataIterator=emtable.Table.iterRows(outImgsFn))

    def _updateItem(self, item, row):
        newFn = row.get(RELIONCOLUMNS.rlnImageName.value)
        index, file = cryosparcToLocation(newFn)
        item.setLocation((index, self._getExtraPath(file)))
        item.setCTF(rowToCtfModel(row))
        pixelSize = item.getSamplingRate()
        item.setTransform(rowToAlignment(row, ALIGN_PROJ, pixelSize))
        item.setSamplingRate(calculateNewSamplingRate(item.getDim(),
                                                      self._getInputParticles().getSamplingRate(),
                                                      self._getInputParticles().getDim()))

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
                                  self.refVolume.get(),
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
            summary.append("Reference Volume: %s" %
                           self.getObjectTag('refVolume'))
            summary.append("Reference Mask: %s" %
                           self.getObjectTag('refMask'))

            summary.append("Inner radius of the window: %s" %
                           str(self.inner_radius.get()))

            summary.append("Outer radius of the window: %s" %
                           str(self.outer_radius.get()))

            summary.append("--------------------------------------------------")
            summary.append("Output particles %s" %
                           self.getObjectTag('outputParticles'))
        return summary

    # ---------------Utils Functions-------------------------------------------

    def _defineParamsName(self):
        """ Define a list with all protocol parameters names"""
        self._paramsName = ['inner_radius',
                            'outer_radius',
                            'use_premult',
                            'use_halfmaps',
                            'lpf_volume',
                            'mask_threshold',
                            'mask_fill_holes',
                            'compute_use_ssd']
        self.lane = str(self.getAttributeValue('compute_lane'))
        self.preprocessLane = str(self.getAttributeValue('preprocess_lane'))

    def doPartStract(self):
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
            if paramName != 'lpf_volume' and paramName != 'mask_threshold':
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif paramName == 'lpf_volume' and self.getAttributeValue(paramName) is not None:
                if float(self.getAttributeValue(paramName)) > 0:
                    params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif paramName == 'mask_threshold' and self.getAttributeValue(paramName) is not None:
                if float(self.getAttributeValue(paramName)) > 0:
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

        runPartStractJob = enqueueJob(self._className, self.projectName.get(),
                                        self.workSpaceName.get(),
                                        str(params).replace('\'', '"'),
                                        str(input_group_connect).replace('\'', '"'),
                                        self.lane, gpusToUse,
                                        result_connect=input_result_connect)

        self.runPartStract = String(runPartStractJob.get())
        self.currenJob.set(self.runPartStract.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(), self.runPartStract.get(),
                         "An error occurred in the particles subtraction process. "
                         "Please, go to cryoSPARC software for more "
                         "details.", self)
        clearIntermediateResults(self.projectName.get(), self.runPartStract.get())
