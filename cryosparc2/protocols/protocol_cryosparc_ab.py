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

import pyworkflow.em.metadata as md
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        LabelParam, IntParam,
                                        EnumParam, StringParam,
                                        BooleanParam, PathParam,
                                        LEVEL_ADVANCED)
from pyworkflow.em.data import Volume
from pyworkflow.em.protocol import ProtInitialVolume
from pyworkflow.em.packages.relion.convert import (writeSetOfParticles,readSetOfParticles,
                                                   getVersion, relionToLocation)
import os
import commands

class ProtCryoSparcInitialModel(ProtInitialVolume):
    """    
    Generate a 3D initial model _de novo_ from 2D particles using
    CryoSparc Stochastic Gradient Descent (SGD) algorithm.
    """
    _label = 'CryoSparc 3D initial model'


    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineFileNames(self):
        """ Centralize how files are called within the protocol. """
        myDict = {
                  'input_particles': self._getPath('input_particles.star'),
                  'out_particles': self._getPath() + '/output_particle.star'
                  }
        self._updateFilenamesDict(myDict)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      label="Input particles", important=True,
                      help='Select the input images from the project.')
        form.addParallelSection(threads=1, mpi=1)

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
        writeSetOfParticles(imgSet, self._getFileName('input_particles'), outputDir=self._getExtraPath(), fillMagnification=True)

        self._program = os.path.join(os.environ['CRYOSPARC_DIR'], 'cryosparc2_master/bin/cryosparcm cli')    
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
        

    def processStep(self):
        print(self._program + "  \'do_run_abinit(\""+self.a[-1].split()[-1]+"\", \""+self.b[-1].split()[-1]+"\", \"\'"+self._user+"\'\", \"" + self.par + "\",\"\'1\'\")\'")
        while commands.getstatusoutput(self._program +" \'get_job_status(\""+self.a[-1].split()[-1]+"\", \""+self.c[-1].split()[-1]+"\")\'")[-1].split()[-1] != 'completed':
            print("waiting...\n")
            commands.getstatusoutput(self._program + " \'wait_job_complete(\""+self.a[-1].split()[-1]+"\", \""+self.c[-1].split()[-1]+"\")\'")
        print("Ab Initial Model Generation Started...")
        self.d = commands.getstatusoutput(self._program + " \'do_run_abinit(\""+self.a[-1].split()[-1]+"\", \""+self.b[-1].split()[-1]+"\", \"\'"+self._user+"\'\", \"" + self.par + "\",\"\'1\'\")\'")
        while commands.getstatusoutput(self._program + " \'get_job_status(\""+self.a[-1].split()[-1]+"\", \""+self.d[-1].split()[-1]+"\")\'")[-1].split()[-1] != 'completed':
            commands.getstatusoutput(self._program + " \'wait_job_complete(\""+self.a[-1].split()[-1]+"\", \""+self.d[-1].split()[-1]+"\")\'")



    
    def createOutputStep(self):
        self._program2 = os.path.join(os.environ['PYEM_DIR'], 'csparc2star.py')
        self.runJob(self._program2, self._ssd+'/'+self.a[-1].split()[-1]+'/'+self.d[-1].split()[-1]+"/cryosparc_"+self.a[-1].split()[-1]+"_"+self.d[-1].split()[-1]+"_class_00_final_particles.cs"+" "+self._getFileName('out_particles'), numberOfMpi=1)
        # Link the folder on SSD to scipion directory
        os.system("ln -s "+self._ssd+"/"+self.a[-1].split()[-1]+'/'+self.d[-1].split()[-1]+" "+ self._getExtraPath())

       
        imgSet = self.inputParticles.get()
        vol = Volume()
        fnVol = self._getExtraPath() + "/" + self.d[-1].split()[-1] + "/cryosparc_"+self.a[-1].split()[-1]+"_"+self.d[-1].split()[-1]+"_class_00_final_volume.mrc"
        vol.setFileName(fnVol)
        vol.setSamplingRate(imgSet.getSamplingRate())

        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(imgSet)
        self._fillDataFromIter(outImgSet)

        self._defineOutputs(outputVolume=vol)
        self._defineSourceRelation(self.inputParticles, vol)
        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(self.inputParticles, outImgSet)
    
    
    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        validateMsgs = []
        return validateMsgs


    #--------------------------- UTILS functions --------------------------------------------

    def _fillDataFromIter(self, imgSet):
        outImgsFn = self._getFileName('out_particles')
        imgSet.setAlignmentProj()
        imgSet.copyItems(self.inputParticles.get(),
                         updateItemCallback=self._createItemMatrix,
                         itemDataIterator=md.iterRows(outImgsFn, sortByLabel=md.RLN_IMAGE_ID))

    def _createItemMatrix(self, item, row):
        from pyworkflow.em.packages.relion.convert import createItemMatrix
        from pyworkflow.em import ALIGN_PROJ

        createItemMatrix(item, row, align=ALIGN_PROJ)
