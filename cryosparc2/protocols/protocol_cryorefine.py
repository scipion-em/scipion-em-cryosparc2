# **************************************************************************
# *
# *  Authors:     Szu-Chi Chung (phonchi@stat.sinica.edu.tw)
# *               Yunior C. Fonseca Reyna (cfonseca@cnb.csic.es)
# *
# * SABID Laboratory, Institute of Statistical Science, Academia Sinica
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

import pyworkflow.em as em
import pyworkflow.em.metadata as md
from pyworkflow.protocol.params import (PointerParam, FloatParam, LabelParam,
                                        IntParam, EnumParam, StringParam,
                                        Positive, BooleanParam, PathParam,
                                        LEVEL_ADVANCED)
from pyworkflow.em.data import Volume, FSC
from pyworkflow.em.protocol import ProtRefine3D
from pyworkflow.em import ALIGN_PROJ
from cryosparc2.utils import *
from cryosparc2.constants import *

relionPlugin = pwutils.importFromPlugin("relion.convert", doRaise=True)

import os
import commands
import ast, json


class ProtCryoSparcRefine3D(ProtRefine3D):
    """ Protocol to refine a 3D map using cryosparc.
    """
    _label = '3D refinement'

    # --------------------------- DEFINE param functions ----------------------
    def _defineFileNames(self):
        """ Centralize how files are called within the protocol. """
        myDict = {
                  'input_particles': self._getPath('input_particles.star'),
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
        form.addParam('referenceVolume', PointerParam, pointerClass='Volume',
                       important=True,
                       label="Input volume",
                       help='Initial reference 3D map, it should have the same '
                            'dimensions and the same pixel size as your input '
                            'particles.')

        form.addParallelSection(threads=1, mpi=1)

        # --------------[Homogeneous Refinement]---------------------------

        form.addSection(label='Homogeneous Refinement')
        # form.addParam('refine_N', IntParam, default=None,
        #               expertLevel=LEVEL_ADVANCED,
        #               label="Refinement box size (Voxels)",
        #               help='The volume size to use for refinement. If this is '
        #                    'null, use the full image size. Otherwise images '
        #                    'are automatically downsampled')

        form.addParam('refine_symmetry', StringParam, default='C1',
                      label="Symmetry",
                      help='Symmetry String (C, D, I, O, T). E.g. C1, D7, C4, '
                           'etc')

        form.addParam('refine_symmetry_do_align', BooleanParam, default=True,
                      label="Do symmetry alignment",
                      help='Align the input structure to the symmetry axes')

        form.addParam('refine_do_init_scale_est', BooleanParam, default=True,
                      label="Re-estimate greyscale level of input reference")

        form.addParam('refine_num_final_iterations', IntParam, default=0,
                      expertLevel=LEVEL_ADVANCED,
                      label="Number of extra final passes",
                      help='Number of extra passes through the data to do '
                           'after the GS-FSC resolution has stopped improving')

        # form.addParam('refine_res_align_max', IntParam, default=None,
        #               expertLevel=LEVEL_ADVANCED,
        #               label="Maximum align resolution (A)",
        #               help='Manual override for maximum resolution that is '
        #                    'used for alignment. This value is normally set by '
        #                    'the GS-FSC')

        form.addParam('refine_res_init', IntParam, default=30,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label="Initial lowpass resolution (A)",
                      help='Applied to input structure')

        form.addParam('refine_res_gsfsc_split', IntParam, default=20,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label="GSFSC split resolution (A)",
                      help='Resolution beyond which two GS-FSC halves are '
                           'independent')

        # form.addParam('refine_highpass_res', IntParam, default=None,
        #               expertLevel=LEVEL_ADVANCED,
        #               label="Highpass resolution (A)")

        form.addParam('refine_SPW', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Use SPW")

        # form.addParam('refine_particle_mw_kda', IntParam, default=None,
        #               expertLevel=LEVEL_ADVANCED,
        #               label="Particle MW (KDa)")

        form.addParam('refine_FSC_weight', StringParam, default='fsc_loosemask',
                      expertLevel=LEVEL_ADVANCED,
                      label="FSC Weighting")

        form.addParam('refine_bnb_params', StringParam, default='3D',
                      expertLevel=LEVEL_ADVANCED,
                      label="BnB Params")

        form.addParam('refine_clip', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Enforce non-negativity",
                      help='Clip negative density. Probably should be false')

        form.addParam('refine_window', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Skip interpolant premult",
                      help='Softly window the structure in real space with a '
                           'spherical window. Should be true')

        form.addParam('refine_skip_premult', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Window structure in real space",
                      help='Leave this as true')

        form.addParam('refine_ignore_dc', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Ignore DC component",
                      help='Ignore the DC component of images. Should be true')

        form.addParam('refine_batchsize_init', IntParam, default=0,
                      expertLevel=LEVEL_ADVANCED,
                      label="Initial batchsize",
                      help='Number of images used in the initial iteration. '
                           'Set to zero to autotune')

        form.addParam('refine_batchsize_epsilon', FloatParam, default=0.001,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label="Batchsize epsilon",
                      help='Controls batch size when autotuning batchsizes. '
                           'Set closer to zero for larger batches')

        form.addParam('refine_batchsize_snrfactor', FloatParam, default=40.0,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label="Batchsize snrfactor",
                      help='Specifies the desired improvement in SNR from the '
                           'images when autotuning batchsizes. Directly '
                           'multiplies the number of images in the batch')

        form.addParam('refine_scale_min', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Minimize over per-particle scale")

        form.addParam('refine_scale_align_use_prev', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Use scales from previous iteration during "
                            "alignment")

        form.addParam('refine_scale_ctf_use_current', BooleanParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=False,
                      label="Use scales from current alignment in reconstruction",
                      help='Use scales from current alignment in reconstruction')

        form.addParam('refine_scale_start_iter', IntParam, default=0,
                      expertLevel=LEVEL_ADVANCED,
                      label="Scale min/use start iter",
                      help='Iteration to start minimizing over per-particle scale')

        form.addParam('refine_noise_model', StringParam, default='symmetric',
                      label="Noise model (white, symmetric or coloured)",
                      help='Noise model to be used. Valid options are white, '
                           'coloured or symmetric. Symmetric is the default, '
                           'meaning coloured with radial symmetry')

        form.addParam('refine_noise_priorw', IntParam, default=50,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label="Noise priorw",
                      help='Weight of the prior for estimating noise (units of '
                           '# of images)')

        form.addParam('refine_noise_initw', IntParam, default=200,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label="Noise initw",
                      help='Weight of the initial noise estimate (units of # '
                           'of images)')

        form.addParam('refine_noise_init_sigmascale', IntParam, default=3,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label="Noise initial sigma-scale",
                      help='Scale factor initially applied to the base noise '
                           'estimate')

        form.addParam('refine_minisize', IntParam, default=2000,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label="Computational minibatch size",
                      help='Number of images to use in each minibatch - only '
                           'affects computational performance. 1000 is a good '
                           'number, but try 4000 if you have lots of RAM')

        form.addParam('refine_mask', StringParam, default='dynamic',
                      expertLevel=LEVEL_ADVANCED,
                      label="Mask (dynamic, static, null)",
                      help='Type of masking to use. Either "dynamic", '
                           '"static", or "null"')

        form.addParam('refine_dynamic_mask_thresh_factor', FloatParam, default=0.2,
                      expertLevel=LEVEL_ADVANCED,
                      validators=[Positive],
                      label="Dynamic mask threshold (0-1)",
                      help='Level set threshold for selecting regions that are '
                           'included in the dynamic mask. Probably don\'t need '
                           'to change this')

        form.addParam('refine_dynamic_mask_near_ang', FloatParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=6.0,
                      validators=[Positive],
                      label="Dynamic mask near (A)",
                      help='Controls extent to which mask is expanded. At the '
                           'near distance, the mask value is 1.0 (in A)')

        form.addParam('refine_dynamic_mask_far_ang', FloatParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=14.0,
                      validators=[Positive],
                      label="Dynamic mask far (A)",
                      help='Controls extent to which mask is expanded. At the '
                           'far distance the mask value becomes 0.0 (in A)')

        form.addParam('refine_dynamic_mask_start_res', IntParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=12,
                      validators=[Positive],
                      label="Dynamic mask start resolution (A)",
                      help='Map resolution at which to start dynamic masking '
                           '(in A)')

        form.addParam('refine_dynamic_mask_use_abs', BooleanParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=False,
                      label="Dynamic mask use absolute value",
                      help='Include negative regions if they are more negative '
                           'than the threshold')

        # --------------[Compute settings]---------------------------

        form.addSection(label='Compute settings')

        form.addParam('compute_use_ssd', BooleanParam, default=True,
                      label='Cache particle images on SSD:',
                      help='Use the SSD to cache particles. Speeds up '
                           'processing significantly')


    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._defineFileNames()
        self._defineParamsName()

        objId = self._getInputParticles().getObjId()
        self._insertFunctionStep("convertInputStep", objId)
        self._insertFunctionStep('processStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions ------------------------------
    def convertInputStep(self, particlesId):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file.
        """
        imgSet = self._getInputParticles()
        # Create links to binary files and write the relion .star file
        relionPlugin.writeSetOfParticles(imgSet,
                                         self._getFileName('input_particles'),
                                         outputDir=self._getExtraPath(),
                                         fillMagnification=True)

        self._importParticles()
        self.vol_fn = os.path.join(os.getcwd(),
                                   relionPlugin.convertBinaryVol(self.referenceVolume.get(),
                                                                 self._getExtraPath()))
        self.importVolume = self.doImportVolumes()
        self.importVolume = self.importVolume[-1].split()[-1]

    def processStep(self):
        self.vol = self.importVolume + '.imported_volume.map'

        while getJobStatus(self.projectName, self.importedParticles) != 'completed':
            waitJob(self.projectName, self.importedParticles)

        while getJobStatus(self.projectName, self.importVolume) != 'completed':
            waitJob(self.projectName, self.importVolume)

        print("Refinement started...")
        self.runRefine = self.doRunRefine()[-1].split()[-1]
        while getJobStatus(self.projectName, self.runRefine) != 'completed':
            waitJob(self.projectName, self.importVolume)

    # -------------------------- STEPS functions ------------------------------
    def createOutputStep(self):
        """
        """
        self._program2 = os.path.join(os.environ['PYEM_DIR'], 'csparc2star.py')

        commands.getstatusoutput(self._program + " \'get_job_streamlog(\"" +
                                 self.projectName+"\", \"" + self.runRefine +
                                 "\")\'" + ">" +self._getFileName('stream_log'))

        # Get the metadata information from stream.log
        with open(self._getFileName('stream_log')) as f:
            data = f.readlines()

        x = ast.literal_eval(data[1])

        # Find the ID of last iteration
        for y in x:
            if y.has_key('text'):
                z = str(y['text'])
                if z.startswith('FSC'):
                    idd = y['imgfiles'][2]['fileid']
                    itera = z[-3:]

        self.runJob(self._program2, self.projectPath + '/' + self.projectName + '/' +
                    self.runRefine + "/cryosparc_" +self.projectName + "_" +
                    self.runRefine+"_" + itera + "_particles.cs" + " " +
                    self._getFileName('out_particles') +" -p " + self.projectPath +
                    "/" + self.projectName + '/' + self.runRefine +
                    "/passthrough_particles.cs", numberOfMpi=1)

        # Link the folder on SSD to scipion directory
        os.system("ln -s " + self.projectPath + "/" + self.projectName + '/' +
                  self.runRefine + " " + self._getExtraPath())

        fnVol = (self._getExtraPath() + "/" + self.runRefine + "/cryosparc_" +
                self.projectName + "_" + self.runRefine + "_" + itera +
                 "_volume_map.mrc")
        half1 = (self._getExtraPath() + "/" + self.runRefine + "/cryosparc_" +
                 self.projectName + "_" + self.runRefine + "_" + itera +
                 "_volume_map_half_A.mrc")
        half2 = (self._getExtraPath() + "/" + self.runRefine + "/cryosparc_" +
                 self.projectName + "_" + self.runRefine + "_" + itera +
                 "_volume_map_half_B.mrc")
        imgSet = self._getInputParticles()
        vol = Volume()
        vol.setFileName(fnVol)
        vol.setSamplingRate(imgSet.getSamplingRate())
        vol.setHalfMaps([half1, half2])

        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(imgSet)
        self._fillDataFromIter(outImgSet)

        self._defineOutputs(outputVolume=vol)
        self._defineSourceRelation(self.inputParticles, vol)
        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(self.inputParticles, outImgSet)
        # Need to get the host IP address if it is not stanalone installation
        os.system("wget 127.0.0.1:39000/file/" + idd + " -nd -P" +
                  self._getExtraPath())
        os.system("mv " + self._getExtraPath() + "/" + idd + " " +
                  self._getExtraPath()+"/fsc.txt")
        # Convert into scipion fsc format
        f=open(self._getExtraPath()+"/fsc.txt", "r")
        lines=f.readlines()
        wv=[]
        corr = []
        for x in lines[1:-1]:
            wv.append(str(float(x.split('\t')[0])/(int(self._getInputParticles().getDim()[0])*float(imgSet.getSamplingRate()))))
            corr.append(x.split('\t')[6])
        f.close()

        fsc = FSC(objLabel=self.getRunName())
        fsc.setData(wv, corr)
        wv2, corr2 = fsc.getData()

        self._defineOutputs(outputFSC=fsc)
        self._defineSourceRelation(vol, fsc)

    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        """ Should be overriden in subclasses to
        return summary message for NORMAL EXECUTION.
        """
        validateMsgs = []
        return validateMsgs
    # -------------------------- UTILS functions ------------------------------

    def _getInputParticles(self):
        return self.inputParticles.get()

    def _fillDataFromIter(self, imgSet):
        outImgsFn = self._getFileName('out_particles')
        imgSet.setAlignmentProj()
        imgSet.copyItems(self._getInputParticles(),
                         updateItemCallback=self._createItemMatrix,
                         itemDataIterator=md.iterRows(outImgsFn,
                                                      sortByLabel=md.RLN_IMAGE_ID))

    def _createItemMatrix(self, particle, row):
        relionPlugin.createItemMatrix(particle, row, align=ALIGN_PROJ)
        relionPlugin.setRelionAttributes(particle, row,
                                         md.RLN_PARTICLE_RANDOM_SUBSET)

    def _importParticles(self):
        """
        Initialize all utils cryoSPARC variables
        """
        self._program = getCryosparcProgram()
        self._user = getCryosparcUser()
        self._ssd = getCryosparcSSD()

        # Create a cryoSPARC project dir
        self.projectDirName = suffix + self.getProject().getShortName()
        self.projectPath = pwutils.join(self._ssd, self.projectDirName)
        self.projectDir = createProjectDir(self.projectPath)

        # create empty project or load an exists one
        folderPaths = getProjectPath(self.projectPath)
        if not folderPaths:
            self.a = createEmptyProject(self.projectPath, self.projectDirName)
            self.projectName = self.a[-1].split()[-1]
        else:
            self.projectName = folderPaths[0]

        # create empty workspace
        self.b = createEmptyWorkSpace(self.projectName, self.getRunName(),
                                      self.getObjComment())
        self.workSpaceName = self.b[-1].split()[-1]

        print("Importing Particles")

        # import_particles_star
        self.c = self.doImportParticlesStar()

        self.importedParticles = self.c[-1].split()[-1]
        self.par = self.importedParticles + '.imported_particles'

    def doImportParticlesStar(self):
        """
        do_import_particles_star(puid, wuid, uuid, abs_star_path,
                                 abs_blob_path=None, psize_A=None)
        returns the new uid of the job that was created
        """
        cmd = """ 'do_import_particles_star("%s","%s", "'+%s'", "'%s'", "'%s'", "'%s'")'"""
        import_particles_cmd = (self._program + cmd % (
            self.projectName, self.workSpaceName,
            self._user,
            os.path.join(os.getcwd(),
            self._getFileName('input_particles')),
            os.path.join(os.getcwd(),
            self._getExtraPath()),
            str(self._getInputParticles().getSamplingRate())
        ))
        print(pwutils.greenStr(import_particles_cmd))
        return commands.getstatusoutput(import_particles_cmd)

    def doImportVolumes(self):
        """
        :return:
        """

        className = "import_volumes"
        params = {"volume_blob_path": str(self.vol_fn),
                  "volume_out_name": "map",
                  "volume_psize": str(self._getInputParticles().getSamplingRate())}

        return doJob(className, self.projectName, self.workSpaceName,
                     str(params).replace('\'', '"'), '{}')

    def _defineParamsName(self):
        """ Define a list with all protocol parameters names"""
        self._paramsName = ['refine_symmetry',
                            'refine_symmetry_do_align',
                            'refine_do_init_scale_est',
                            'refine_num_final_iterations',
                            'refine_res_init',
                            'refine_res_gsfsc_split',
                            'refine_FSC_weight', 'refine_bnb_params',
                            'refine_clip',
                            'refine_window', 'refine_skip_premult',
                            'refine_ignore_dc',
                            'refine_batchsize_init',
                            'refine_batchsize_snrfactor',
                            'refine_batchsize_epsilon',
                            'refine_scale_min', 'refine_scale_align_use_prev',
                            'refine_scale_ctf_use_current',
                            'refine_scale_start_iter',
                            'refine_noise_model', 'refine_noise_priorw',
                            'refine_noise_initw',
                            'refine_noise_init_sigmascale',
                            'refine_minisize', 'refine_mask',
                            'refine_dynamic_mask_thresh_factor',
                            'refine_dynamic_mask_near_ang',
                            'refine_dynamic_mask_far_ang',
                            'refine_dynamic_mask_start_res',
                            'refine_dynamic_mask_use_abs',
                            'compute_use_ssd']


    def doRunRefine(self):
        """
        :return:
        """
        className = "homo_refine"
        input_group_conect = {"particles": str(self.par),
                              "volume": str(self.vol)}
        # {'particles' : 'JXX.imported_particles' }
        params = {}

        for paramName in self._paramsName:
            params[str(paramName)] = str(self.getAttributeValue(paramName))

        return doJob(className, self.projectName, self.workSpaceName,
                     str(params).replace('\'', '"'),
                     str(input_group_conect).replace('\'', '"'))





