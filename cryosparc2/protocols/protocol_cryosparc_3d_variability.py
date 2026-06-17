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
from pyworkflow import BETA
from pyworkflow.protocol.params import (PointerParam, FloatParam, Positive,
                                        BooleanParam, EnumParam)
from pwem.objects import Volume
from pwem.protocols import ProtRefine3D

from . import ProtCryosparcBase
from ..convert import *
from ..utils import *
from ..constants import *


class ProtCryoSparc3DVariability(ProtCryosparcBase, ProtRefine3D):
    """
    Protocol to compute the principle modes of variability with a dataset
    of aligned particles
    """

    """
        Computes the principal modes of structural variability from a set of
        aligned cryo-EM particles using CryoSPARC 3D Variability Analysis.
        The protocol identifies dominant conformational motions by estimating
        the covariance structure of the reconstructed volume, allowing users
        to explore continuous heterogeneity and flexible regions within
        macromolecular complexes.

        AI Generated:

        3D Variability Analysis (ProtCryoSparc3DVariability) — User Manual
            Overview

            The 3D Variability Analysis protocol is designed to identify and
            characterize continuous structural heterogeneity within cryo-EM
            datasets. Using a set of aligned particles together with a spatial
            mask, the protocol computes the principal modes of variability of
            the reconstructed density map.

            Instead of separating particles into discrete conformational classes,
            this approach models structural flexibility as continuous motions.
            Mathematically, the protocol estimates the dominant eigenvectors of
            the 3D covariance matrix, which represent the principal directions
            of structural variability present in the dataset.

            For biological users, this protocol is particularly useful when
            studying flexible proteins, multi-domain assemblies, molecular
            machines, or complexes undergoing conformational transitions.
            It enables visualization and interpretation of continuous motions
            that are often difficult to capture with standard classification
            approaches.

            Inputs and General Workflow

            The protocol requires a set of aligned particles together with an
            input mask defining the region where variability should be analyzed.
            Accurate particle alignment is essential because the protocol assumes
            that observed differences between particles primarily reflect true
            structural variability rather than orientation errors.

            During execution, the aligned particles and mask are connected to
            the CryoSPARC 3D Variability job. The algorithm iteratively estimates
            the dominant modes of structural motion while applying spatial filters
            and regularization constraints to stabilize the optimization process.

            After completion, the protocol generates a variability map and an
            updated particle set containing the refined alignment information and
            associated metadata required for downstream analysis.

            Variability Modes and Biological Interpretation

            One of the most important parameters is the number of modes to solve.
            Each mode represents an independent direction of structural variability
            within the dataset. These modes are mathematically orthogonal and can
            correspond to biologically meaningful motions such as domain rotations,
            hinge movements, breathing motions, ligand-induced rearrangements, or
            partial assembly transitions.

            In many biological applications, a small number of modes is sufficient
            to capture the dominant conformational landscape. Increasing the number
            of modes may reveal additional subtle motions but can also introduce
            noise or non-interpretable variability if the dataset quality is limited.

            Interpretation of variability modes should always be supported by
            visual inspection and biological plausibility. Not every mathematical
            mode necessarily corresponds to a true physical motion.

            Masking and Region Selection

            The input mask plays a central role in determining which structural
            regions contribute to the variability analysis. Proper masking allows
            the protocol to focus on biologically relevant motions while excluding
            solvent noise or unrelated density regions.

            From a biological perspective, the mask should ideally include the
            structurally stable core together with the flexible regions of interest.
            Poor masking may cause variability estimates to become dominated by
            noise, solvent fluctuations, or peripheral artefacts.

            In highly flexible systems, careful masking is often one of the most
            important factors affecting the interpretability of the resulting modes.

            Iterative Optimization and Stability

            The protocol performs an iterative optimization process controlled by
            the number of iterations parameter. Increasing the number of iterations
            may improve convergence and refinement quality, particularly for complex
            datasets with subtle conformational differences.

            To stabilize optimization, the protocol includes a regularization
            parameter referred to as lambda. Larger lambda values suppress unstable
            fluctuations and reduce the likelihood of reconstruction artefacts,
            while lower values allow more flexible variability estimation.

            If optimization diverges or generates unrealistic density changes,
            increasing the regularization strength often improves robustness.

            Frequency Filtering

            Several parameters control the spatial frequency range used during
            variability estimation. Low-pass filtering limits the analysis to
            broader structural motions by suppressing high-resolution noise,
            whereas high-pass filtering removes large-scale low-frequency trends
            that may not correspond to biologically meaningful variability.

            The filter order parameters determine the sharpness of these frequency
            transitions. In most practical cryo-EM workflows, default values provide
            stable behavior, although advanced users may adjust them when working
            with particularly noisy datasets or highly dynamic systems.

            Noise and Scale Modeling

            The protocol allows different noise models to be used during variability
            estimation. By default, a colored noise model is employed, although a
            white noise model can also be selected depending on dataset properties.

            Additionally, the protocol supports several strategies for handling
            per-particle scaling factors. Input scales preserve previously estimated
            scaling values, whereas optimal scaling recomputes particle-specific
            scales dynamically during the variability analysis.

            Proper scale handling can improve stability and consistency, especially
            when datasets originate from multiple refinements or imaging conditions.

            Symmetry Considerations

            Symmetry parameters allow the protocol to enforce structural symmetry
            during variability estimation. Applying the correct symmetry can improve
            signal quality and stability for symmetric particles.

            However, biological users should exercise caution because enforcing
            symmetry may suppress asymmetric conformational motions that are
            biologically relevant. For systems suspected to contain asymmetric
            flexibility, lower symmetry constraints may provide more informative
            variability modes.

            Outputs and Downstream Analysis

            After execution, the protocol produces an output variability volume
            together with an updated particle set containing transformed alignment
            information. Internally, CryoSPARC particle metadata is converted into
            Scipion-compatible STAR files, allowing seamless integration with
            downstream workflows.

            The resulting variability maps can later be visualized as continuous
            trajectories or animated motions, providing an intuitive representation
            of conformational changes across the dataset.

            Biologically, these outputs are especially useful for exploring dynamic
            mechanisms, identifying flexible domains, and understanding structural
            transitions involved in molecular function.

            Practical Recommendations

            In most routine cryo-EM workflows, starting with three variability modes
            and default filtering parameters provides a good balance between
            interpretability and computational stability. If the resulting motions
            appear noisy or unstable, stronger filtering or increased regularization
            may improve robustness.

            Accurate alignment quality is critical for meaningful results. Misaligned
            particles can introduce artificial variability that may be incorrectly
            interpreted as biological motion.

            Careful visual inspection of variability trajectories and reconstructed
            maps is essential before drawing structural or mechanistic conclusions.

            Final Perspective

            The 3D Variability Analysis protocol provides a powerful framework for
            investigating continuous conformational heterogeneity directly from
            cryo-EM particle datasets. By modeling dominant modes of structural
            motion, the protocol enables researchers to explore molecular dynamics
            beyond the limitations of discrete classification methods.

            For biological interpretation, the combination of accurate particle
            alignment, appropriate masking, stable regularization, and careful
            visualization remains essential for obtaining reliable and meaningful
            insights into macromolecular flexibility.
        """
    _label = '3D variability Analysis '
    _devStatus = BETA

    # --------------------------- DEFINE param functions ----------------------
    def _defineFileNames(self):
        """ Centralize how files are called within the protocol. """
        myDict = {
                  'input_particles': self._getTmpPath('input_particles.star'),
                  'out_particles': self._getPath() + '/output_particle.star',
                  'stream_log': self._getPath()+'/stream.log'
                  }
        self._updateFilenamesDict(myDict)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      label="Input particles", important=True,
                      validators=[Positive],
                      help='Select the input images from the project.')
        form.addParam('refMask', PointerParam, pointerClass='VolumeMask',
                      default=None,
                      label='Input Mask',
                      allowsNull=False,
                      help='Mask raw data')

        form.addParallelSection(threads=1, mpi=1)

        # --------------[3D Variability]---------------------------

        form.addSection(label='3D Variability')
        form.addParam('var_K', IntParam, default=3,
                      label="Number of modes to solve",
                      help='The number of orthogonal principle modes (i.e. '
                           'eigenvectors of the 3D covariance) to solve.')

        form.addParam('var_num_iterations', IntParam, default=20,
                      validators=[Positive],
                      label="Number of iterations",
                      help='Number of iterations')

        # form.addParam('var_N', FloatParam, default=-1,
        #               expertLevel=LEVEL_ADVANCED,
        #               label="Refinement box size (Voxels)",
        #               help='The volume size to use for refinement. If this is '
        #                    '-1, use the full image size. Otherwise images '
        #                    'are automatically downsampled')

        addSymmetryParam(form)

        # form.addParam('var_num_particles', IntParam, default=None,
        #               label="Only use this many particles")

        form.addParam('var_filter_res', FloatParam, default=None,
                      validator=[Positive],
                      label="Filter resolution (A)",
                      help='Resolution at which results are filtered')

        form.addParam('var_filter_order', FloatParam, default=1.5,
                      validators=[Positive],
                      label="Filter order",
                      help='Order of filter')

        # form.addParam('var_highpass_res', StringParam, default='',
        #               expertLevel=LEVEL_ADVANCED,
        #               label="Highpass resolution (A)",
        #               help='Resolution below which variability is ignored')

        form.addParam('var_highpass_order', FloatParam, default=8,
                      validators=[Positive],
                      label="Highpass order",
                      help='Order of filter')

        form.addParam('var_use_gramschmidt', BooleanParam, default=True,
                      label="Use Gram-Schmidt",
                      help='Order of filter')

        form.addParam('var_use_white_noise', BooleanParam, default=False,
                      label="Use white noise model",
                      help='Use a white noise model (default) instead of a '
                           'colored noise model. One or the other may work '
                           'better depending on the dataset.')

        form.addParam('var_use_scales', EnumParam,
                      choices=['none', 'input', 'optimal'],
                      default=2,
                      label="Per-particle scale",
                      help='How to treat per-particle scale factors. No scales '
                           'means all particles have scale 1.0 (useful for very '
                           'small particles or strange cases). Input scales '
                           'means to use the input scale factors (which may '
                           'have come from another refinement). Optimal means '
                           'to compute per-particle optimal scales on the fly '
                           'during 3D variability.')

        form.addParam('var_lambda', FloatParam, default=0.01,
                      validators=[Positive],
                      label="Lambda",
                      help='Stabilizing coefficient - try larger values if '
                           'optimization diverges and creates artefacts. '
                           'In v2.12 this was changed to a normalized '
                           'fractional value, where 0.01 should work for '
                           'almost all datasets.')
        # --------------[Compute settings]---------------------------
        form.addSection(label="Compute settings")
        addComputeSectionParams(form, allowMultipleGPUs=False)

    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._defineFileNames()
        self._defineParamsName()
        self._initializeCryosparcProject()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.processStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions ------------------------------
    def processStep(self):
        self.info(pwutils.yellowStr("3D Variability started..."))
        self.doRun3DVariability()

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to Relion file
        """
        self._initializeUtilsVariables()
        csOutputFolder = os.path.join(self.projectDir.get(),
                                      self.run3DVariability.get())
        csParticlesName = (getOutputPreffix(self.projectName.get()) +
                           self.run3DVariability.get() + "_particles.cs")

        csFile = os.path.join(csOutputFolder, csParticlesName)

        # Copy the particles to scipion output folder
        os.system("cp -r " + csFile + " " + self._getExtraPath())
        csFile = os.path.join(self._getExtraPath(), csParticlesName)

        outputStarFn = self._getFileName('out_particles')
        argsList = [csFile, outputStarFn]

        convertCs2Star(argsList)

        fnVolName = (getOutputPreffix(self.projectName.get()) +
                     self.run3DVariability.get() + "_map.mrc")

        fnVol = os.path.join(csOutputFolder, fnVolName)

        # Copy the volumes to extra folder
        os.system("cp -r " + fnVol + " " + self._getExtraPath())
        fnVol = os.path.join(self._getExtraPath(), fnVolName)

        imgSet = self._getInputParticles()
        vol = Volume()
        vol.setFileName(fnVol)
        vol.setSamplingRate(calculateNewSamplingRate(vol.getDim(),
                                                     imgSet.getSamplingRate(),
                                                     imgSet.getDim()))
        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(imgSet)
        self._fillDataFromIter(outImgSet)

        self._defineOutputs(outputVolume=vol)
        self._defineSourceRelation(self.inputParticles.get(), vol)
        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(self.inputParticles.get(), outImgSet)

    #  ----------------------------UTILS functions --------------------------

    def cleanTmp(self):
        """ Delete all files and subdirectories under Tmp folder. For this
        protocol we need to keep the tmp folder content"""
        pass

    def _validate(self):
        validateMsgs = cryosparcValidate()
        if not validateMsgs:
            validateMsgs = gpusValidate(self.getGpuList(), checkSingleGPU=True)
            if not validateMsgs:
                particles = self._getInputParticles()
                if not particles.hasCTF():
                    validateMsgs.append("The Particles has not associated a "
                                        "CTF model")
        return validateMsgs

    def _summary(self):
        summary = []
        if (not hasattr(self, 'outputVolume') or
                not hasattr(self, 'outputParticles') or
                not hasattr(self, 'outputFSC')):
            summary.append("Output objects not ready yet.")
        else:
            summary.append("Input Particles: %s" %
                           self.getObjectTag('inputParticles'))
            summary.append("Input Volume: %s" %
                           self.getObjectTag('refVolume'))
            summary.append("Input Mask: %s" %
                           self.getObjectTag('refMask'))
            summary.append("Symmetry: %s" %
                           getSymmetry(self.symmetryGroup.get(),
                                       self.symmetryOrder.get())
                           )
            summary.append("------------------------------------------")
            summary.append("Output particles %s" %
                           self.getObjectTag('outputParticles'))
            summary.append("Output volume %s" %
                           self.getObjectTag('outputVolume'))
        return summary

    # -------------------------- UTILS functions ------------------------------

    def _fillDataFromIter(self, imgSet):
        outImgsFn = 'particles@' + self._getFileName('out_particles')
        imgSet.setAlignmentProj()
        imgSet.copyItems(self._getInputParticles(),
                         updateItemCallback=self._createItemMatrix,
                         itemDataIterator=emtable.Table.iterRows(
                             fileName=outImgsFn))

    def _createItemMatrix(self, particle, row):
        createItemMatrix(particle, row, align=ALIGN_PROJ)
        setCryosparcAttributes(particle, row,
                               RELIONCOLUMNS.rlnRandomSubset.value)

    def _defineParamsName(self):
        """ Define a list with all protocol parameters names"""
        self._paramsName = ['var_K',
                            'var_filter_res',
                            'var_filter_order',
                            'var_highpass_order',
                            'var_use_gramschmidt',
                            'var_use_white_noise',
                            'var_use_scales',
                            'var_num_iterations',
                            'var_lambda',
                            'compute_use_ssd',
                            'var_symmetry']
        self.lane = str(self.getAttributeValue('compute_lane'))

    def doRun3DVariability(self):
        """
        :return:
        """
        className = "var_3D"
        input_group_conect = {"particles": str(self.particles),
                              "mask": str(self.mask)}
        params = {}

        for paramName in self._paramsName:
            if (paramName != 'var_symmetry' and
                    paramName != 'var_use_scales'):
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif paramName == 'var_symmetry':
                symetryValue = getSymmetry(self.symmetryGroup.get(),
                                           self.symmetryOrder.get())
                params[str(paramName)] = symetryValue
            elif paramName == 'var_use_scales':
                params[str(paramName)] = str(VAR_USE_SCALES[self.var_use_scales.get()])

        # Determinate the GPUs to use (in dependence of
        # the cryosparc version)
        try:
            if not self.useQueueForSteps() and not self.useQueue():  # not using queue system
                gpusToUse = self.getGpuList()
            else:  # using queue system
                gpusToUse = False
        except Exception:
            gpusToUse = False

        self.run3DVariability = enqueueJob(className, self.projectName.get(),
                                           self.workSpaceName.get(),
                                           str(params).replace('\'', '"'),
                                           str(input_group_conect).replace('\'',
                                                                           '"'),
                                           self.lane, gpusToUse)

        self.currenJob.set(self.run3DVariability.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(), self.run3DVariability.get(),
                         "An error occurred in the 3D Variability process. "
                         "Please, go to cryosPARC software for more "
                         "details.", self)



