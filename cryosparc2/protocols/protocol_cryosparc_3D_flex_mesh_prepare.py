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
from pyworkflow.protocol.params import (PointerParam, FloatParam, BooleanParam, FileParam, StringParam)
from . import ProtCryosparcBase
from ..convert import *
from ..utils import *
from ..constants import *


class ProtCryoSparc3DFlexMeshPrepare(ProtCryosparcBase):
    """
    Prepares consensus cryo-EM density maps and associated particle information
    for downstream 3DFlex analysis by generating a tetrahedral mesh representation
    suitable for modeling continuous structural heterogeneity.

    AI Generated:

    3D Flex Mesh Prepare (ProtCryoSparc3DFlexMeshPrepare) — User Manual
        Overview

        The 3D Flex Mesh Prepare protocol generates the structural mesh required
        for CryoSPARC 3DFlex workflows. Its primary purpose is to transform a
        consensus cryo-EM reconstruction into a deformable tetrahedral model
        capable of representing continuous molecular motions. This mesh serves
        as the physical framework that later stages of 3DFlex use to learn and
        reconstruct conformational variability directly from particle images.

        In practical cryo-EM studies, biological macromolecules often exhibit
        flexibility that cannot be adequately described by a small number of
        discrete classes. Instead of separating particles into rigid states,
        3DFlex attempts to model smooth transitions and continuous structural
        landscapes. The mesh preparation stage is therefore essential because
        it defines how the structure can deform during downstream analysis.

        Biological Motivation

        Flexible proteins, ribosomes, membrane complexes, viral assemblies,
        and multi-domain machines frequently contain regions that move relative
        to one another. Traditional refinement approaches may average these
        motions together, causing blurred density and limiting biological
        interpretation. The 3DFlex framework addresses this limitation by
        introducing a physically structured deformation model.

        The mesh generated in this protocol approximates the molecular volume
        using connected tetrahedral elements. During later refinement and
        reconstruction stages, these elements deform smoothly in response to
        latent conformational coordinates. As a result, biologically meaningful
        motions such as domain rotations, hinge bending, breathing motions,
        or gradual conformational transitions can be represented more naturally.

        Inputs and General Workflow

        The protocol requires a previously prepared consensus reconstruction
        originating from a dedicated 3DFlex data preparation workflow. This
        consensus map defines the overall molecular shape from which the mesh
        will be generated. Optionally, a user-defined solvent mask may also
        be provided to constrain the molecular region used during mesh creation.

        When no mask is supplied, the protocol can automatically generate one
        from the consensus density. This automatic masking strategy is often
        sufficient for compact and well-resolved particles. However, for
        challenging datasets containing weak peripheral density, detergent
        micelles, flexible extensions, or surrounding noise, providing a
        carefully designed mask is usually preferable.

        The dimensions of the mask should match the dimensions of the consensus
        reconstruction. Consistent voxel size and box dimensions are important
        because the resulting mesh must accurately correspond to the molecular
        volume used throughout the 3DFlex workflow.

        Solvent Mask Preparation

        Solvent masking is one of the most biologically important aspects of
        mesh generation because it determines which regions of the reconstruction
        participate in the deformation model. Regions outside the mask are
        considered irrelevant background and are excluded from the mesh.

        The protocol allows automatic mask generation by filtering the input
        map, thresholding the density, dilating the selected region, and adding
        soft padding. Together, these operations define a continuous molecular
        envelope suitable for deformation modeling.

        The filtering step suppresses high-frequency noise before thresholding.
        Lower threshold values generally include more peripheral density, while
        higher values produce tighter masks focused on the strongest structural
        regions. Dilation expands the mask slightly to ensure the molecular
        boundaries remain connected, while soft padding smooths transitions
        near the solvent boundary.

        From a biological perspective, the mask should ideally contain all
        regions expected to move while excluding disconnected noise or empty
        solvent areas. Overly restrictive masks may artificially constrain
        biologically meaningful motions, whereas excessively loose masks may
        introduce unstable or unrealistic deformations.

        Mesh Resolution and Tetrahedral Elements

        The central parameter controlling mesh complexity is the number of
        tetrahedral cells spanning the reconstruction volume. This setting
        determines the effective spatial resolution of the deformation model.

        A coarse mesh contains fewer and larger tetrahedra, producing smoother
        and more global motions. Such meshes are computationally efficient and
        often appropriate for large conformational rearrangements or noisy
        datasets. In contrast, finer meshes contain smaller tetrahedral
        elements that can capture more localized flexibility but require
        increased computational resources and may become more sensitive to
        overfitting.

        In biological practice, the optimal mesh density depends on particle
        size, structural complexity, and expected flexibility. Large molecular
        machines with multiple independently moving domains may benefit from
        finer meshes, whereas small or relatively rigid particles are often
        better modeled using coarser representations.

        Segmentation and Subdomain Modeling

        The protocol optionally supports segmentation-guided mesh generation.
        This is particularly useful for complexes composed of multiple domains
        or subunits that move semi-independently.

        Segmentation files can define biologically meaningful structural
        regions, allowing the protocol to construct separate submeshes that
        are later connected into a unified deformation framework. This strategy
        improves the interpretability of flexible motions and can significantly
        stabilize refinement in systems with articulated or modular architecture.

        Segment connection definitions specify how individual regions are linked
        together. Biologically, these connections should reflect realistic
        structural relationships between domains. Incorrect connectivity may
        generate deformation paths that are physically implausible or difficult
        to interpret.

        The protocol also supports rigid segment definitions. Regions marked
        as rigid are constrained to resist deformation more strongly than
        surrounding areas. This is especially useful for highly stable cores,
        membrane-embedded domains, or experimentally validated rigid bodies.

        Rigidity Weighting and Motion Regularization

        Rigidity weighting controls how easily different parts of the structure
        can deform. Dense structural regions generally behave more rigidly,
        while low-density or peripheral regions are allowed greater flexibility.

        The minimum rigidity weight parameter regulates the contrast between
        stable and flexible regions. Lower values permit greater expansion and
        contraction in weak-density areas, whereas higher values produce more
        globally constrained motions.

        An additional option allows peripheral low-density regions to be
        artificially stiffened. This can help prevent unstable or noisy motions
        in poorly resolved regions, particularly in small particles or datasets
        with limited signal-to-noise ratio. However, excessive rigidity may
        oversmooth biologically meaningful transitions and blur the boundaries
        between independently moving domains.

        From a biological standpoint, rigidity regularization should balance
        stability and realism. The goal is to suppress implausible deformations
        while preserving authentic conformational variability.

        Outputs and Their Interpretation

        The protocol produces a tetrahedral mesh representation associated with
        the consensus reconstruction. This mesh becomes the structural foundation
        for downstream 3DFlex training and reconstruction stages.

        A mesh visualization file is also generated for inspection and validation.
        Visual examination of the mesh is strongly recommended before proceeding
        with training. Users should verify that tetrahedral elements adequately
        cover the molecular density, preserve major structural regions, and avoid
        disconnected or excessively distorted areas.

        Biologically meaningful motion modeling depends heavily on the quality
        of this mesh. Poor mesh geometry, incomplete masking, or inappropriate
        segmentation can negatively affect all subsequent stages of flexible
        refinement.

        Practical Recommendations

        For most biological datasets, beginning with automatic mask generation
        and moderate mesh density is a reasonable strategy. If the resulting
        motions appear unstable or physically unrealistic during downstream
        analysis, refining the solvent mask or increasing rigidity constraints
        often improves stability.

        Segmentation-guided meshes are especially valuable for complexes with
        clearly separable domains or hinge-like motions. In such systems,
        incorporating prior biological knowledge into the mesh design can
        substantially improve interpretability.

        Smaller particles or low-resolution datasets may require stronger
        rigidity weighting to avoid overfitting. Conversely, highly flexible
        systems with large conformational transitions may benefit from reduced
        rigidity constraints and finer mesh representations.

        Final Perspective

        The mesh preparation stage is not merely a technical preprocessing step
        but a biologically meaningful definition of how molecular motion will
        be represented throughout the 3DFlex workflow. Careful masking, sensible
        mesh density selection, and biologically informed rigidity constraints
        are essential for obtaining realistic and interpretable continuous
        heterogeneity models in cryo-EM studies.
    """
    _label = '3D flex mesh prepare'
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
        form.addParam('dataPrepare', PointerParam,
                      pointerClass='ProtCryoSparc3DFlexDataPrepare',
                      label="Data prepare protocol",
                      important=True,
                      help='Prepared particles and consensus volume from the 3D flex data prepare protocol.')
        form.addParam('refMask', PointerParam, pointerClass='VolumeMask',
                      default=None,
                      label='Input Mask',
                      allowsNull=True,
                      help='Mask raw data')

        form.addSection(label='Mesh Prepare')

        solventMaskGroup = form.addGroup('Solvent mask preparation',
                                         condition="refMask is None")
        solventMaskGroup.addParam('mask_in_lowpass_A', IntParam, default=10,
                                  condition="refMask is None",
                                  label="Filter input volume res. (A)",
                                  help="Filter the input consensus reconstruction "
                                       "volume to this resolution (A) before thresholding "
                                       "to create the outer solvent mask. Solvent mask is"
                                       " only generated if it is not already input.")

        solventMaskGroup.addParam('mask_in_threshold_level', FloatParam, default=0.5,
                                  condition="refMask is None",
                                  label="Mask threshold",
                                  help="Threshold the input consensus reconstruction "
                                       "volume (after filtering) at this absolute "
                                       "density level to make the solvent mask. "
                                       "Solvent mask is only generated if it is "
                                       "not already input.")

        solventMaskGroup.addParam('mask_dilate_A', IntParam,
                                  default=2,
                                  condition="refMask is None",
                                  label="Mask dilation (A)",
                                  help="After thresholding, dilate by this much"
                                       " (A) to create solvent mask")

        solventMaskGroup.addParam('mask_pad_A', IntParam,
                                  default=5,
                                  condition="refMask is None",
                                  label="Mask soft padding (A)",
                                  help="After thresholding and dilation, soft "
                                       "pad by this much (A) to create solvent mask")

        meshPreparationGroup = form.addGroup('Mesh preparation')
        meshPreparationGroup.addParam('tetra_num_cells', IntParam,
                                  default=20,
                                  label="Base num. tetra cells",
                                  help="Number of tetrahedral cells that fit "
                                       "across the extent of the box. Use this "
                                       "to set the size of mesh elements. If "
                                       "set to e.g. 20, the base tetramesh will"
                                       " have elements that are the right size "
                                       "to create a spacing of 20 tetra elements "
                                       "across the box extent in each x,y,z "
                                       "direction. A higher number makes a finer mesh.")

        meshPreparationGroup.addParam('tetra_segments_path', FileParam,
                                      allowsNull=True,
                                      default=None,
                                      label="Segmentation file path",
                                      help="Absolute path to a segmentation file "
                                           "in either .seg format from UCSF "
                                           "Chimera Segger tool or else .mrc format "
                                           "(see CryoSPARC guide for details), "
                                           "defining subdomain regions that "
                                           "should each have a submesh. "
                                           "Submeshes are fused to make final "
                                           "mesh using the segment connections list.")

        meshPreparationGroup.addParam('tetra_segments_fuse_list', StringParam,
                                      allowsNull=True,
                                      default=None,
                                      label="Segment connections",
                                      help="A comma and '>' separated list of "
                                           "connections between segments to use when "
                                           "fusing sub-meshes to make the final mesh. "
                                           "See CryoSPARC guide for full explanation. "
                                           "For example, '0>3, 0>4, 3>2, 2>1' is a valid connection"
                                           " string. Each pair X>Y denotes that segment Y is "
                                           "joined to segment X. The connections must form "
                                           "a 'tree' structure and cannot have cycles. "
                                           "The first pair X>Y must start with the "
                                           "root of the tree as X. The connections "
                                           "must be in breadth-first order of the tree. "
                                           "When using Chimera Segger segmentation input, "
                                           "the X>Y numbers should be region_ids from Segger (e.g., 948>947)")

        meshPreparationGroup.addParam('tetra_rigid_list', StringParam,
                                      allowsNull=True,
                                      default=None,
                                      label="Rigid segments",
                                      help="A comma separated list of segments to make "
                                           "rigid. This is done by setting the rigidity "
                                           "weight of tetra elements for this region to 20. ")


        rigidityWeighting = form.addGroup('Rigidity weighting')
        rigidityWeighting.addParam('rigidity_penalty_min', FloatParam,
                                      default=0.5,
                                      label="Min. rigidity weight",
                                      help="Rigidity weights of tetra elements"
                                           " are 1.0 in the most dense regions "
                                           "of the input consensus map, and fall "
                                           "off to this value (default 0.5) in the"
                                           " least dense/empty regions. This helps "
                                           "encourage the deformation model to expand/contract"
                                           " empty space without distorting the protein density.")

        rigidityWeighting.addParam('rigidity_penalty_stiffen_low_density', BooleanParam,
                                   default=False,
                                   label="Stiffen low density regions",
                                   help="Turning this on will cause the rigidity"
                                        " weights of tetra elements at the periphery "
                                        "of the input consensus density to be increased to 3.0 . "
                                        "Empty regions will still have low rigidity, but non-empty "
                                        "regions at the periphery of the structure will be rigidified."
                                        " This helps to combat overfitting in smaller particles or"
                                        " poor SNR data where otherwise low density peripheral "
                                        "features start to 'fly around'. However, it can also cause "
                                        "the deformations to be overly smooth and blur motion 'boundaries between domains. ")


        # --------------[Compute settings]---------------------------
        form.addSection(label="Compute settings")
        addComputeSectionParams(form, allowMultipleGPUs=False, needGPU=False)

    def _insertAllSteps(self):
        self._defineFileNames()
        self._defineParamsPrepareName()
        self._defineParamsMeshName()
        self._initializeCryosparcProject()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.mergePrepareStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions ------------------------------
    def mergePrepareStep(self):
        self.info(pwutils.yellowStr("3D Flex Mesh Preparation started..."))
        self.doRun3DFlexMeshPrepare()

    def createOutputStep(self):
        """
        Create the protocol output. Convert cryosparc file to Relion file
        """
        # save pdb file with mesh
        # this will not be scipion output
        # I just want to display it in the viewer
        jobId = self.run3DFlexMeshPrep
        csOutputFolder = os.path.join(self.projectDir.get(),
                                       self.run3DFlexMeshPrep.get())
        pdbMeshName = "%s_mesh_pdb.pdb" % jobId
        csOutputFolder = os.path.join(self.projectDir.get(),
                                      self.run3DFlexMeshPrep.get())
        copyFiles(csOutputFolder, self._getExtraPath(), files=[pdbMeshName])
    # ------------------------- Utils methods ----------------------------------

    def _fillDataFromIter(self, imgSet):
        outImgsFn = 'particles@' + self._getFileName('out_particles')
        imgSet.setAlignmentProj()
        imgSet.copyItems(self._getInputParticles(),
                         updateItemCallback=self._createItemMatrix,
                         itemDataIterator=emtable.Table.iterRows(fileName=outImgsFn))

    def _createItemMatrix(self, particle, row):
        createItemMatrix(particle, row, align=ALIGN_PROJ)
        setCryosparcAttributes(particle, row, RELIONCOLUMNS.rlnRandomSubset.value)

    def _validate(self):
        validateMsgs = cryosparcValidate()
        if not validateMsgs:
            mask = self.refMask.get()
            if mask:
                maskDim = mask.getDim()
                dataPrepareProt = self.dataPrepare.get()
                if dataPrepareProt:
                    if hasattr(dataPrepareProt, 'outputVolume'):
                        outputVolumeDim = dataPrepareProt.outputVolume.getDim()
                        if maskDim != outputVolumeDim:
                            validateMsgs.append('The dimension of the mask must be %s according to the 3D Flex data prepare protocol(Training box size parameter)' % str(outputVolumeDim) )
                else:
                    validateMsgs.append('You need to specify the 3D Flex Data Prepare protocol')

        return validateMsgs

    def _defineParamsPrepareName(self):
        """ Define a list with 3D Flex Prepare Data parameters names"""
        self._paramsPrepareName = ['box_size_pix', 'bin_size_pix', 'alpha_min',
                            'keep_num_particles']
        self.lane = str(self.getAttributeValue('compute_lane'))

    def _defineParamsMeshName(self):
        """ Define a list with 3D Flex Mesh Prepare parameters names"""
        self._maskMeshPrepareName = ['mask_in_lowpass_A', 'mask_in_threshold_level',
                                     'mask_dilate_A', 'mask_pad_A']
        self._paramsMeshName = ['tetra_num_cells', 'tetra_segments_path',
                                'tetra_segments_fuse_list', 'tetra_rigid_list',
                                'rigidity_penalty_min',
                                'rigidity_penalty_stiffen_low_density']
        self.lane = str(self.getAttributeValue('compute_lane'))

    def doRun3DFlexMeshPrepare(self):
        self._className = 'flex_meshprep'
        params = {}
        varDataPrepJob = str(self.dataPrepare.get().run3DFlexDataPrepJob)
        input_group_connect = {"volume": str('%s.volume' % varDataPrepJob)}

        if self.refMask.get() is not None:
            input_group_connect["mask"] = str(self.mask)
        else:
            for paramName in self._maskMeshPrepareName:
                if self.getAttributeValue(paramName) is not None:
                    params[str(paramName)] = str(self.getAttributeValue(paramName))

        for paramName in self._paramsMeshName:
            if self.getAttributeValue(paramName) is not None:
                params[str(paramName)] = str(self.getAttributeValue(paramName))

        run3DFlexMeshPrepJob = enqueueJob(self._className,
                                          self.projectName.get(),
                                          self.workSpaceName.get(),
                                          str(params).replace('\'', '"'),
                                          str(input_group_connect).replace('\'',
                                                                           '"'),
                                          self.lane, False)

        self.run3DFlexMeshPrep = String(run3DFlexMeshPrepJob.get())
        self.currenJob.set(run3DFlexMeshPrepJob.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(),
                         self.run3DFlexMeshPrep.get(),
                         "An error occurred in the 3D Flex Mesh Preparation process. "
                         "Please, go to cryoSPARC software for more "
                         "details.", self)
        clearIntermediateResults(self.projectName.get(),
                                 self.run3DFlexMeshPrep.get())










