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
    Prepares particles for use in 3DFlex training and reconstruction. At the same
    way, takes in a consensus (rigid) refinement density map, plus optionally
    a segmentation and generates a tetrahedral mesh for 3DFlex.
    """
    """
        3D Flex Data Prepare (ProtCryoSparc3DFlexDataPrepare) — User Manual

        Overview

        The 3D Flex Data Prepare protocol prepares particle datasets and
        consensus maps for subsequent 3D Flex training and reconstruction
        within cryoSPARC. Its main purpose is to transform a rigidly refined
        cryo-EM dataset into an optimized input representation suitable for
        modeling continuous structural heterogeneity. In practice, this
        protocol standardizes particle metadata, crops and downsamples
        particle images when necessary, filters particles according to
        quality criteria, and generates the consensus density map required
        for flexible reconstruction workflows.

        In a typical cryo-EM analysis pipeline, this protocol is used after
        consensus refinement has already been completed and particle poses
        are known. It acts as a preprocessing bridge between traditional
        rigid refinement and neural-network-based flexibility analysis.
        For biological users, this becomes especially important when
        studying molecular motions, conformational continua, or structural
        transitions that cannot be captured by discrete classification alone.

        Inputs and General Workflow

        The protocol requires two essential inputs: a particle set with
        associated CTF information and 3D alignment parameters, and a
        reference consensus volume representing the rigid reconstruction.
        The particle set must already be aligned, since 3D Flex relies on
        known particle orientations as the structural basis for learning
        conformational variability. The reference map serves as the
        structural scaffold from which the flexible model will infer
        continuous deformation fields.

        An optional preprocessing stage allows cropping the particles and
        consensus volume to a smaller box size. This is typically used to
        remove empty solvent regions and reduce computational cost. For
        biological systems where the particle occupies only a fraction of
        the original reconstruction box, cropping can substantially improve
        training efficiency without losing structural information. Care
        should be taken to ensure the biologically relevant density remains
        fully contained within the cropped region.

        The protocol also supports downsampling to a training box size,
        which defines the effective resolution used during model
        optimization. This is one of the most critical parameters for
        successful 3D Flex analysis. Smaller training boxes reduce
        computational burden and often stabilize optimization, while
        excessively large values can dramatically increase runtime and may
        attempt to learn high-resolution detail beyond what the consensus
        reconstruction can reliably support. In most practical workflows,
        the chosen training resolution should remain below the validated FSC
        resolution of the rigid consensus map.

        Particle filtering options allow selecting only subsets of the input
        dataset for training. The minimum scale threshold removes particles
        whose scale factors suggest poor signal quality or contamination,
        which is often useful for excluding junk particles or poorly
        reconstructed views. Similarly, the protocol allows restricting the
        total number of particles used. This is particularly relevant for
        large datasets, since 3D Flex training can become computationally
        expensive. The selected particle count must always be divisible by
        1000, which reflects internal cryoSPARC batching constraints for
        training and reconstruction consistency.

        Internal Processing

        Internally, the protocol launches a cryoSPARC flex_prep job that
        processes the input particles and reference map according to the
        selected parameters. Once execution completes, cryoSPARC outputs are
        automatically copied into the protocol workspace and converted into
        RELION-compatible STAR metadata for integration within the Scipion
        ecosystem. This conversion ensures interoperability with downstream
        Scipion workflows while preserving particle identity and alignment
        information.

        Outputs and Interpretation

        The output consists of two biologically meaningful objects. The
        first is a new prepared particle set, which contains only the
        selected particles reformatted for 3D Flex compatibility while
        preserving their original metadata and sampling information. The
        second is the prepared consensus volume, stored as an MRC map with
        corrected dimensions and sampling rate adjusted to match the
        processed particle box size. Together, these outputs define the
        complete structural input required for subsequent 3D Flex training.

        Validation and Consistency Checks

        Validation checks ensure that the input data satisfy all
        requirements before execution begins. The protocol verifies that
        particles contain associated CTF models and valid 3D alignment
        parameters, since flexible reconstruction cannot proceed without
        these prerequisites. It also confirms that any explicitly requested
        particle count satisfies divisibility requirements. These safeguards
        prevent incompatible datasets from entering computationally
        expensive processing stages.

        Practical Recommendations

        From a biological perspective, careful parameter selection is
        essential for meaningful flexibility analysis. Excessive cropping
        may truncate flexible peripheral regions, while aggressive
        downsampling can suppress subtle motions of interest. Conversely,
        overly large box sizes may increase computational cost without
        improving interpretability. In most workflows, it is advisable to
        begin with conservative downsampling and moderate particle
        filtering, then refine these settings based on the quality and
        resolution of the consensus refinement.

        Final Perspective

        For most cryo-EM practitioners, 3D Flex Data Prepare is more than a
        simple formatting step. It defines the structural representation
        that the flexible neural model will learn from, directly shaping
        the biological motions that can ultimately be resolved. Thoughtful
        preprocessing therefore has a strong influence on the quality,
        interpretability, and reliability of downstream continuous
        heterogeneity analysis.
        """



    _label = '3D flex data prepare'
    _devStatus = BETA
    _protCompatibility = [V4_1_0, V4_1_1, V4_1_2, V4_2_0, V4_2_1, V4_3_1,
                          V4_4_0, V4_4_1, V4_5_1, V4_5_3, V4_6_0, V4_6_1, V4_6_2, V4_7_0, V4_7_1]

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








