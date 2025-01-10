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
    Prepares particles for use in 3DFlex training and reconstruction. At the same
    way, Takes in a consensus (rigid) refinement density map, plus optionally
    a segmentation and generates a tetrahedral mesh for 3DFlex.
    """
    _label = '3D flex mesh prepare'
    _devStatus = BETA
    _protCompatibility = [V4_1_0, V4_1_1, V4_1_2, V4_2_0, V4_2_1, V4_3_1,
                          V4_4_0, V4_4_1, V4_5_1, V4_5_3, V4_6_0, V4_6_1, V4_6_2]

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










