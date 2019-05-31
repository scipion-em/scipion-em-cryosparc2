# **************************************************************************
# *
# *  Authors:     Szu-Chi Chung (phonchi@stat.sinica.edu.tw) 
# *
# * SABID Laboratory, Institute of Statistical Science, Academia Sinica
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
                                        BooleanParam, PathParam, LEVEL_ADVANCED)
from pyworkflow.em.data import Volume, FSC
from pyworkflow.em.protocol import ProtRefine3D
from pyworkflow.em import ALIGN_PROJ
from pyworkflow.utils import importFromPlugin

relionPlugin = importFromPlugin("relion.convert", doRaise=True)

import os
import commands
import ast


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
                      help='Select the input images from the project.')
        form.addParam('referenceVolume', PointerParam, pointerClass='Volume',
                       important=True,
                       label="Input volume",
                       help='Initial reference 3D map, it should have the same '
                            'dimensions and the same pixel size as your input '
                            'particles.')
        self.addSymmetry(form)
        form.addParallelSection(threads=1, mpi=1)

    def addSymmetry(self, container):
        container.addParam('symmetryGroup', StringParam, default='C1',
                           label="Symmetry",
                           help='If the molecule is asymmetric, set Symmetry ')

    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._defineFileNames()
        objId = self.inputParticles.get().getObjId()
        self._insertFunctionStep("convertInputStep", objId)
        self._insertFunctionStep('processStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions ------------------------------
    def convertInputStep(self, particlesId):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file.
        """
        imgSet = self.inputParticles.get()
        # Create links to binary files and write the relion .star file
        relionPlugin.writeSetOfParticles(imgSet,
                                         self._getFileName('input_particles'),
                                         outputDir=self._getExtraPath(),
                                         fillMagnification=True)

        self._program = os.path.join(os.environ['CRYOSPARC_DIR'],
                                     'cryosparc2_master/bin/cryosparcm cli')
        self._user = os.environ['CRYOSPARC_USER']
        self._ssd = os.environ['CRYOSSD_DIR']
        print("Importing Particles")
        #create_empty_project(owner_user_id, project_container_dir, title=None, desc=None) returns the new uid of the project that was created
        self.a = commands.getstatusoutput(self._program + " \'create_empty_project(\"\'+"+self._user+"\'\", \"\'"+self._ssd+"\'\")\'")
        #create_empty_workspace(project_uid, created_by_user_id, created_by_job_uid=None, title=None, desc=None) returns the new uid of the workspace that was created
        self.b = commands.getstatusoutput(self._program + " \'create_empty_workspace(\""+self.a[-1].split()[-1]+"\", \"\'+"+self._user+"\'\")\'")
        #do_import_particles_star(puid, wuid, uuid, abs_star_path, abs_blob_path=None, psize_A=None) returns the new uid of the job that was created
        self.c = commands.getstatusoutput(self._program + " \'do_import_particles_star(\""+self.a[-1].split()[-1]+"\", \""+self.b[-1].split()[-1]+"\", \"\'+"+self._user+"\'\", \"\'"+ os.path.join(os.getcwd(),self._getFileName('input_particles'))+"\'\", \"\'"+os.path.join(os.getcwd(),self._getExtraPath())+"\'\", \"\'"+str(self.inputParticles.get().getSamplingRate())+"\'\")\'")
        self.par = self.c[-1].split()[-1]+'.imported_particles'

        self.vol_fn = os.path.join(os.getcwd(),relionPlugin.convertBinaryVol(self.referenceVolume.get(), self._getTmpPath()))
        print(self.vol_fn)
        self.e = commands.getstatusoutput(self._program + " \'do_import_volumes(\""+self.a[-1].split()[-1]+"\", \""+self.b[-1].split()[-1]+"\", \"\'+"+self._user+"\'\", \"\'"+ self.vol_fn + "\'\", \"\'map\'\", \"\'"+str(imgSet.getSamplingRate())+"\'\")\'")
        print(self._program + " \'do_import_volumes(\""+self.a[-1].split()[-1]+"\", \""+self.b[-1].split()[-1]+"\", \"\'+"+self._user+"\'\", \"\'"+ self.vol_fn + "\'\", \"\'map\'\", \"\'"+str(imgSet.getSamplingRate())+"\'\")\'")

    def processStep(self):
        self.vol = self.e[-1].split()[-1]+'.imported_volume.map'
        print(self._program + " \'do_run_refine(\""+self.a[-1].split()[-1]+"\", \""+self.b[-1].split()[-1]+"\", \"\'+"+self._user+"\'\", \"" + self.par + "\", \"" + self.vol + "\", None,"+ "\"\'"+self.symmetryGroup.get()+"\'\")\'")
        while commands.getstatusoutput(self._program + " \'get_job_status(\""+self.a[-1].split()[-1]+"\", \""+self.c[-1].split()[-1]+"\")\'")[-1].split()[-1] != 'completed':
            print("waiting...\n")
            commands.getstatusoutput(self._program + " \'wait_job_complete(\""+self.a[-1].split()[-1]+"\", \""+self.c[-1].split()[-1]+"\")\'")

        while commands.getstatusoutput(self._program + " \'get_job_status(\""+self.a[-1].split()[-1]+"\", \""+self.e[-1].split()[-1]+"\")\'")[-1].split()[-1] != 'completed':
            print("waiting...\n")
            commands.getstatusoutput(self._program + " \'wait_job_complete(\""+self.a[-1].split()[-1]+"\", \""+self.e[-1].split()[-1]+"\")\'")

        print("Refinement Started...")
        self.f = commands.getstatusoutput(self._program + " \'do_run_refine(\""+self.a[-1].split()[-1]+"\", \""+self.b[-1].split()[-1]+"\", \"\'+"+self._user+"\'\", \"" + self.par + "\", \"" + self.vol + "\", None,"+ "\"\'"+str(self.symmetryGroup.get())+"\'\")\'")
        while commands.getstatusoutput(self._program + " \'get_job_status(\""+self.a[-1].split()[-1]+"\", \""+self.f[-1].split()[-1]+"\")\'")[-1].split()[-1] != 'completed':
            commands.getstatusoutput(self._program + " \'wait_job_complete(\""+self.a[-1].split()[-1]+"\", \""+self.f[-1].split()[-1]+"\")\'")

    # -------------------------- STEPS functions ------------------------------
    def createOutputStep(self):
        self._program2 = os.path.join(os.environ['PYEM_DIR'], 'csparc2star.py')
        commands.getstatusoutput(self._program + " \'get_job_streamlog(\""+self.a[-1].split()[-1]+"\", \""+self.f[-1].split()[-1]+ "\")\'" + ">" +self._getFileName('stream_log'))
        # Get the metadata information from stream.log
        with open(self._getFileName('stream_log')) as f:
            data = f.readlines()
        x = ast.literal_eval(data[0])

        # Find the ID of last iteration
        for y in x:
            if y.has_key('text'):
                z = str(y['text'])
                if z.startswith('FSC'):
                    idd = y['imgfiles'][2]['fileid']
                    itera = z[-3:]

        self.runJob(self._program2, self._ssd+'/'+self.a[-1].split()[-1]+'/'+self.f[-1].split()[-1]+"/cryosparc_"+self.a[-1].split()[-1]+"_"+self.f[-1].split()[-1]+"_"+itera+"_particles.cs"+" "+self._getFileName('out_particles')+" -p "+self._ssd+"/"+self.a[-1].split()[-1]+'/'+self.f[-1].split()[-1]+"/passthrough_particles.cs", numberOfMpi=1)
        # Link the folder on SSD to scipion directory
        os.system("ln -s "+self._ssd+"/"+self.a[-1].split()[-1]+'/'+self.f[-1].split()[-1]+" "+ self._getExtraPath())

        fnVol = self._getExtraPath() + "/" + self.f[-1].split()[-1] + "/cryosparc_"+self.a[-1].split()[-1]+"_"+self.f[-1].split()[-1]+"_"+itera+"_volume_map.mrc"
        half1 = self._getExtraPath() + "/" + self.f[-1].split()[-1] + "/cryosparc_"+self.a[-1].split()[-1]+"_"+self.f[-1].split()[-1]+"_"+itera+"_volume_map_half_A.mrc"
        half2 = self._getExtraPath() + "/" + self.f[-1].split()[-1] + "/cryosparc_"+self.a[-1].split()[-1]+"_"+self.f[-1].split()[-1]+"_"+itera+"_volume_map_half_B.mrc"
        imgSet = self.inputParticles.get()
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
        os.system("wget 127.0.0.1:39000/file/"+idd+" -nd -P"+ self._getExtraPath())
        os.system("mv "+ self._getExtraPath()+"/"+idd+" "+self._getExtraPath()+"/fsc.txt")
        # Convert into scipion fsc format
        f=open(self._getExtraPath()+"/fsc.txt","r")
        lines=f.readlines()
        wv=[]
        corr = []
        for x in lines[1:-1]:
            wv.append(str(float(x.split('\t')[0])/(int(self.inputParticles.get().getDim()[0])*float(imgSet.getSamplingRate()))))
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

    def _fillDataFromIter(self, imgSet):
        outImgsFn = self._getFileName('out_particles')
        imgSet.setAlignmentProj()
        imgSet.copyItems(self.inputParticles.get(),
                         updateItemCallback=self._createItemMatrix,
                         itemDataIterator=md.iterRows(outImgsFn,
                                                      sortByLabel=md.RLN_IMAGE_ID))

    def _createItemMatrix(self, particle, row):
        relionPlugin.createItemMatrix(particle, row, align=ALIGN_PROJ)
        relionPlugin.setRelionAttributes(particle, row,
                                         md.RLN_PARTICLE_RANDOM_SUBSET)


