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
                     waitForCryosparc, clearIntermediateResults, copyFiles, parse_version, getCryosparcVersion)
from ..constants import *


class ProtCryoSparcSubtract(ProtCryosparcBase, ProtOperateParticles):
    """
    Performs signal subtraction on cryo-EM particle images by removing the contribution
    of a selected region of a reference structure from experimental particles. The protocol
    is intended to isolate flexible domains, compositional heterogeneity, or weakly resolved
    regions that may otherwise be obscured by dominant structural features.

    AI Generated:

    Particle Subtraction (ProtCryoSparcSubtract) — User Manual
        Overview

        The Particle Subtraction protocol removes projected density corresponding to a
        selected region of a reconstructed volume from aligned experimental particles.
        In cryo-EM workflows, this operation is widely used to study structural flexibility,
        compositional variability, or local conformational changes that are difficult to
        analyze when the complete particle dominates the signal.

        From a biological perspective, subtraction is particularly useful for large
        macromolecular assemblies containing flexible domains, transient binding partners,
        membrane-associated regions, or heterogeneous subcomplexes. By computationally
        removing part of the signal, the remaining density can be analyzed independently
        through focused classification, local refinement, or variability analysis.

        Inputs and Biological Context

        The protocol requires a set of aligned particles, a reference volume, and a mask
        defining the density to subtract. The biological interpretation of the subtraction
        strongly depends on the quality and consistency of these inputs. The reference map
        should ideally originate from the same dataset and reconstruction workflow as the
        particles being processed. Using inconsistent maps or mismatched scaling may lead
        to subtraction artifacts or biologically misleading residual densities.

        The alignment information associated with the particles is essential because the
        subtraction relies on projecting the reference structure into the orientation of
        each particle image. Accurate angular assignments are therefore critical for
        reliable subtraction quality.

        Defining the Subtraction Mask

        The subtraction mask is the most biologically important parameter of the protocol.
        The mask must include the density intended for removal while excluding the regions
        that should remain available for downstream analysis. In practical cryo-EM work,
        masks are often designed around flexible domains, peripheral subunits, ligand
        densities, or regions undergoing conformational rearrangements.

        Soft masks are generally preferred because they reduce sharp boundary artifacts and
        minimize Fourier ringing effects. Abrupt mask edges can introduce subtraction
        errors that may later appear as artificial densities during focused refinement or
        classification.

        Careful biological interpretation is required when defining the subtraction region.
        Removing too much density may eliminate important structural context, while removing
        too little may leave dominant signal contributions that continue masking the region
        of interest.

        Windowing and Scaling Considerations

        The protocol allows control over the inner and outer radii used during particle
        windowing and normalization. These parameters influence how smoothly the particle
        edges are treated and can affect subtraction stability.

        In most biological applications, default values are appropriate, particularly when
        particles are well centered and reconstructed under standard cryo-EM conditions.
        However, for particles with large box sizes, elongated geometries, or substantial
        peripheral flexibility, adjusting the windowing parameters may improve subtraction
        quality.

        Scaling options are also important because subtraction assumes that the projected
        density and the experimental particle images share consistent intensity scaling.
        Improper scaling may result in over-subtraction or under-subtraction, both of which
        can complicate downstream interpretation.

        Gold-Standard and Half-Map Strategies

        The protocol supports subtraction strategies compatible with gold-standard cryo-EM
        refinement practices. When half-maps are available, subtraction can be performed
        independently for each half dataset, preserving statistical independence during
        subsequent resolution estimation and refinement.

        Maintaining gold-standard separation is particularly important in publication-level
        workflows or whenever focused refinement will be followed by FSC-based resolution
        analysis. Disabling this strategy may introduce correlations that artificially
        inflate reported resolutions.

        Optional Low-Pass Filtering and Mask Processing

        The input structure may optionally be low-pass filtered before subtraction. This is
        biologically useful when the reference map contains high-resolution features that
        are not reliably represented in the particle images. Applying moderate filtering
        can stabilize subtraction and reduce high-frequency artifacts.

        Additional mask processing options allow thresholding and hole filling. These tools
        help produce more continuous subtraction regions, particularly for fragmented masks
        or noisy segmentations. In practice, biologically meaningful masks should remain
        smooth, interpretable, and consistent with the expected molecular boundaries.

        Outputs and Interpretation

        The protocol produces a new particle dataset in which the selected density has been
        computationally removed. These particles preserve their original orientations and
        metadata while emphasizing the remaining structural information.

        Biologically, the resulting particles often reveal weak conformational variability,
        transient interactions, or mobile regions that were previously obscured. The output
        particles are commonly used for focused classification, local refinement, masked
        reconstruction, or heterogeneity analysis.

        Interpretation of subtraction results should always be performed carefully. Residual
        density does not necessarily indicate a true biological feature and may reflect
        imperfect masking, alignment inaccuracies, or subtraction artifacts.

        Practical Recommendations

        In routine cryo-EM workflows, the most reliable subtraction results are obtained
        when the reference map, mask, and particles all originate from a consistent
        refinement strategy. Before large-scale analysis, it is advisable to visually
        inspect a subset of subtracted particles to verify that the targeted density has
        been removed cleanly without introducing strong artifacts.

        For flexible assemblies, subtracting only the dominant rigid core often improves
        focused classification of mobile regions. Conversely, when analyzing ligand binding
        or small peripheral domains, highly precise masks and conservative subtraction
        settings are usually preferable.

        Final Perspective

        Particle subtraction is one of the most powerful focused-analysis techniques in
        modern cryo-EM image processing. Rather than treating the particle as a single
        rigid object, the protocol enables researchers to isolate biologically relevant
        variability and investigate local structural behavior in much greater detail.
        Successful application depends on careful masking, accurate alignments, and
        biologically meaningful interpretation of the remaining signal.
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
        cryosparcVersion = parse_version(getCryosparcVersion())
        csFileName = "subtracted_particles.cs" if cryosparcVersion < parse_version(V5_0_0) else "particles_0000.cs"

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
