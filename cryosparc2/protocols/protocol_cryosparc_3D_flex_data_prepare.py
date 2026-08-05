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
from pyworkflow.protocol.params import PointerParam
from . import ProtCryosparcBase
from ..convert import *
from ..utils import *
from ..constants import *


class ProtCryoSparc3DFlexDataPrepare(ProtCryosparcBase):
    """
    Prepares particle images and consensus density maps for downstream 3DFlex
    training and flexible reconstruction workflows. The protocol standardizes
    particle data, adapts spatial sampling conditions, and generates the
    processed inputs required for continuous heterogeneity analysis.

    AI Generated:

    3D Flex Data Prepare (ProtCryoSparc3DFlexDataPrepare) — User Manual
        Overview

        The 3D Flex Data Prepare protocol is the entry point for workflows
        involving continuous flexibility analysis with 3DFlex. Its purpose is
        to transform an existing cryo-EM particle dataset and an associated
        consensus reconstruction into a format suitable for flexible modeling
        and deformation learning.

        In cryo-EM studies, many biological systems contain continuous
        structural variability that cannot be adequately represented by a small
        number of discrete classes. The preparation stage establishes a stable
        and standardized dataset that can later be used to learn smooth
        conformational transitions and dynamic structural behavior.

        The protocol is especially important because the quality of all later
        stages, including mesh generation, training, reconstruction, and motion
        visualization, depends strongly on the consistency and reliability of
        the prepared data.

        Inputs and Biological Context

        The protocol requires two principal inputs: a set of aligned particles
        and a consensus three-dimensional density map. The consensus map serves
        as the structural reference describing the average state of the
        particle population, while the particle images provide the experimental
        observations used to model structural variability.

        For biological applications, the consensus map should ideally represent
        a well-refined and stable reconstruction with reliable alignment
        parameters. Poor alignments, severe heterogeneity, or inaccurate
        contrast transfer correction may negatively affect the interpretation
        of subsequent flexibility analyses.

        The input particles are expected to contain valid imaging metadata,
        including alignment and CTF information. These elements are essential
        because flexible reconstruction relies on accurately relating each
        particle image to the consensus structural framework.

        Particle Cropping and Spatial Standardization

        One of the most important stages in preparation is the adaptation of
        the particle box size. The protocol allows particles and reference
        volumes to be cropped to a desired spatial extent before flexible
        training begins.

        Biologically, cropping helps focus the analysis on the region of
        interest while reducing unnecessary solvent regions and computational
        overhead. For compact particles, moderate cropping often improves
        efficiency without affecting interpretability. For larger complexes or
        assemblies with flexible peripheral regions, excessive cropping should
        be avoided because biologically relevant motions may become truncated.

        The selected box size also influences the downstream reconstruction
        stage. Since flexible refinement and high-resolution reconstruction
        share the same spatial framework, the preparation stage establishes the
        geometric basis for the entire workflow.

        Downsampling and Training Resolution

        The protocol also defines a training box size that controls the spatial
        sampling used during model training. This operation effectively limits
        the resolution range used for learning conformational variability.

        In biological practice, lower-resolution training is often advantageous
        because large-scale conformational motions are generally dominated by
        low- and medium-resolution structural information. Restricting the
        training resolution also reduces computational demands and helps
        prevent overfitting.

        The training box size should therefore be selected according to the
        biological question and the quality of the consensus reconstruction.
        Systems dominated by large domain motions can often be analyzed
        effectively at moderate resolution, whereas subtle local rearrangements
        may require finer sampling.

        Particle Selection and Dataset Quality

        The protocol provides mechanisms for restricting the subset of particles
        used during flexibility analysis. This allows users to discard
        low-quality particles or limit the dataset size for computational
        efficiency.

        From a biological perspective, careful particle selection can strongly
        influence the stability and interpretability of learned motions.
        Including excessive junk particles, damaged particles, or poorly aligned
        images may introduce artificial variability unrelated to the true
        conformational landscape.

        Conversely, overly aggressive filtering may remove rare but meaningful
        conformational states. A balanced strategy is therefore recommended,
        particularly for systems expected to exhibit broad structural
        heterogeneity.

        Consensus Volume Preparation

        Alongside particle preparation, the protocol also standardizes the
        consensus density map used as the structural reference for later stages
        of the workflow. The prepared map establishes the coordinate system and
        deformation framework that will guide flexible modeling.

        In practical cryo-EM workflows, the consensus map should represent the
        biologically dominant or most reliable state of the system. Large
        reconstruction artifacts, incomplete density, or severe preferred
        orientation effects may complicate later interpretation of flexibility.

        Outputs and Their Interpretation

        After execution, the protocol produces a prepared particle dataset and
        an associated processed consensus volume suitable for downstream 3DFlex
        analysis. These outputs form the standardized foundation used during
        mesh preparation, training, and flexible reconstruction.

        The resulting particles preserve their experimental identity while
        being reorganized into a format compatible with continuous deformation
        workflows. The processed consensus volume defines the spatial reference
        frame used throughout later stages of the analysis.

        Practical Recommendations

        In routine biological workflows, it is generally advisable to begin
        with a carefully refined consensus reconstruction and a well-cleaned
        particle dataset before attempting flexibility analysis. The quality of
        these initial inputs strongly determines the reliability of all later
        results.

        Moderate downsampling is often sufficient for capturing biologically
        relevant large-scale motions while maintaining manageable computational
        requirements. Extremely large box sizes should only be used when the
        biological system genuinely requires fine structural detail during
        flexibility modeling.

        When substantial conformational variability is expected, users should
        avoid excessive particle filtering so that meaningful motions remain
        represented in the dataset. Visual inspection of the consensus map and
        particle quality remains an essential step before beginning training.

        Final Perspective

        The 3D Flex Data Prepare protocol establishes the structural and
        computational foundation for continuous heterogeneity analysis in
        cryo-EM. By standardizing particles and consensus reconstructions into
        a coherent framework, it enables downstream methods to model molecular
        flexibility, conformational landscapes, and biologically meaningful
        structural motion with greater stability and interpretability.
    """
    _label = '3D flex data prepare'
    _devStatus = BETA

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
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      label="Input particles",
                      important=True,
                      help='Particle stacks to use.')
        form.addParam('refVolume', PointerParam, pointerClass='Volume',
                      label='Initial volume',
                      important=True,
                      help='Initial volume raw data')

        form.addSection(label='Data Prepare')
        form.addParam('box_size_pix', IntParam, default=None,
                      allowsNull=True,
                      label="Crop box size (pix)",
                      help="Crop input particles and volume to this box size. "
                           "This is the box size that will be used for high "
                           "resolution reconstruction with 3D Flex. Particles "
                           "are cropped to this box size first, and then "
                           "downsampled to the Training Box Size for training "
                           "time. Default (None) means to keep the original "
                           "box size of the particles.")

        form.addParam('bin_size_pix', IntParam, default=128,
                      allowsNull=True,
                      label="Training box size (pix)",
                      help="Downsample cropped particles (via Fourier cropping)"
                           " to this box size for training the 3D Flex model. "
                           "This should be chosen to limit 3D Flex training to "
                           "a resolution below the gold-standard FSC resolution "
                           "of the consensus reconstruction, in order to ensure "
                           "high resolution reconstructions can be validated. "
                           "Box sizes over 256 may become prohibitively slow")

        form.addParam('alpha_min', IntParam, default=None,
                      allowsNull=True,
                      label="Min. scale to keep",
                      help="Only keep particles with scale factor above this "
                           "value. Useful for discarding particles that might "
                           "be junk.")

        form.addParam('keep_num_particles', IntParam, default=None,
                      allowsNull=True,
                      label="Num. particles to use",
                      help="Only keep the first X particles. The final number"
                           " of particles used during 3D Flex training "
                           "and reconstruction must be divisible by 1000. "
                           "If this is None (default) then the number of input "
                           "particles will be rounded down to the nearest 1000.")

        # --------------[Compute settings]---------------------------
        form.addSection(label="Compute settings")
        addComputeSectionParams(form, allowMultipleGPUs=False, needGPU=False)

    def _insertAllSteps(self):
        self._defineFileNames()
        self._defineParamsPrepareName()
        self._initializeCryosparcProject()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.dataPrepareStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions ------------------------------
    def dataPrepareStep(self):
        self.info(pwutils.yellowStr("3D Flex Data Preparation started..."))
        self.doRun3DFlexDataPrepare()

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to Relion file
        """
        self.info(pwutils.yellowStr("Creating the output..."))
        self._initializeUtilsVariables()
        outputStarFn = self._getFileName('out_particles')
        self.info(pwutils.yellowStr("outputStarFn: %s " % outputStarFn))
        csOutputFolder = os.path.join(self.projectDir.get(), self.run3DFlexDataPrepJob.get())
        self.info("csOutputFolder: %s " % csOutputFolder)
        # csFileName = "subtracted_particles.cs"
        csFileName = "%s_passthrough_particles.cs" % self.run3DFlexDataPrepJob.get()
        self.info("csFileName: %s " % csOutputFolder)
        # Create the output folder
        copyFiles(csOutputFolder,  os.path.join(self._getExtraPath(), self.run3DFlexDataPrepJob.get()))
        self.info("copyFolder: src-> dst %s %s" % (csOutputFolder, os.path.join(self._getExtraPath(), self.run3DFlexDataPrepJob.get())))
        csFile = os.path.join(self._getExtraPath(), self.run3DFlexDataPrepJob.get(), csFileName)
        self.info("csFile (metadata): %s " % csFile)
        argsList = [csFile, outputStarFn]
        self.info("starFile: %s " % outputStarFn)
        convertCs2Star(argsList)
        self.info("Convert to star done.")
        self.info("Creating the particles output")
        imgSet = self._getInputParticles()
        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(imgSet)
        outImgSet.setSamplingRate(imgSet.getSamplingRate())
        imgSetDict = {}
        for img in imgSet.iterItems():
            imgSetDict[str(img.getIndex()) + '@' + os.path.basename(img.getFileName())] = img.clone()
        self._fillDataFromIter(outImgSet, imgSetDict)
        self.info("Creating the consensus volume  output")

        csMapName = "%s_map.mrc" % self.run3DFlexDataPrepJob.get()
        fnVol = os.path.join(self._getExtraPath(), self.run3DFlexDataPrepJob.get(), csMapName)
        vol = Volume()
        fixVolume(fnVol)
        vol.setFileName(fnVol)
        vol.setSamplingRate(calculateNewSamplingRate(vol.getDim(),
                                                     imgSet.getSamplingRate(),
                                                     imgSet.getDim()))

        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(imgSet, outImgSet)
        self._defineOutputs(outputVolume=vol)

    # ------------------------- Utils methods ----------------------------------

    def _fillDataFromIter(self, outImgSet, imgSetDict):
        filename = 'particles@' + self._getFileName('out_particles')
        for imgRow in emtable.Table.iterRows(filename):
            fileName = imgRow.get(RELIONCOLUMNS.rlnImageName.value)
            splitFileName = fileName.split('@')
            key = str(int(splitFileName[0])) + '@' + '_'.join(splitFileName[1].split('_')[1:])
            if key in imgSetDict:
                img = imgSetDict[key]
                outImgSet.append(img)

    def _createItemMatrix(self, particle, row):
        createItemMatrix(particle, row, align=ALIGN_PROJ)
        setCryosparcAttributes(particle, row, RELIONCOLUMNS.rlnRandomSubset.value)

    def _validate(self):
        validateMsgs = cryosparcValidate()
        if not validateMsgs:
            particles = self._getInputParticles()
            if not particles.hasCTF():
                validateMsgs.append(
                    "The Particles has not associated a CTF model")
                if not particles.hasAlignment3D():
                    validateMsgs.append("The Particles has not aligned")
                    if self.keep_num_particles.get() % 1000 != 0:
                        validateMsgs.append("The final number of particles "
                                            "used during 3D Flex training and "
                                            "reconstruction must be divisible by 1000")

        return validateMsgs

    def _defineParamsPrepareName(self):
        """ Define a list with 3D Flex Prepare Data parameters names"""
        self._paramsPrepareName = ['box_size_pix', 'bin_size_pix', 'alpha_min',
                                   'keep_num_particles']
        self.lane = str(self.getAttributeValue('compute_lane'))

    def doRun3DFlexDataPrepare(self):
        self._className = "flex_prep"
        input_group_connect = {"particles": self.particles.get(),
                               "volume": self.volume.get()}
        params = {}

        for paramName in self._paramsPrepareName:
            if self.getAttributeValue(paramName) is not None:
                params[str(paramName)] = str(self.getAttributeValue(paramName))

        run3DFlexDataPrepJob = enqueueJob(self._className, self.projectName.get(),
                                  self.workSpaceName.get(),
                                  str(params).replace('\'', '"'),
                                  str(input_group_connect).replace('\'', '"'),
                                  self.lane, False)

        self.run3DFlexDataPrepJob = String(run3DFlexDataPrepJob.get())
        self.currenJob.set(self.run3DFlexDataPrepJob.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(), self.run3DFlexDataPrepJob.get(),
                         "An error occurred in the 3D Flex Data Preparation process. "
                         "Please, go to cryoSPARC software for more "
                         "details.", self)
        clearIntermediateResults(self.projectName.get(), self.run3DFlexDataPrepJob.get())








