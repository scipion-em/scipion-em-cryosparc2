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

import sys
import os
import commands
import pyworkflow.em as em
import pyworkflow.em.metadata as md

from pyworkflow.em.protocol import ProtClassify2D, SetOfClasses2D
from pyworkflow.protocol.params import (PointerParam, BooleanParam,
                                        FloatParam, IntParam, Positive, StringParam)
from pyworkflow.em.packages.relion.convert import (writeSetOfParticles,readSetOfParticles,
                                                   getVersion, relionToLocation, rowToAlignment)

class ProtCryo2D(ProtClassify2D):
    """ Wrapper to CryoSparc 2D clustering program.
    """
    _label = 'perform Cryosparc2D'
    
    def __init__(self, **args):
        ProtClassify2D.__init__(self, **args)
        if self.numberOfMpi.get() < 2:
            self.numberOfMpi.set(2)
    
    def _defineFileNames(self):
        """ Centralize how files are called within the protocol. """
        myDict = {
                  'input_particles': self._getPath('input_particles.star'),
                  'out_particles': self._getPath() + '/output_particle.star',
                  'out_class': self._getPath() + '/output_class.star',
                  'out_class_m2': self._getPath() + '/output_class_m2.star'
                  }
        self._updateFilenamesDict(myDict)
    
    # --------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      label="Input particles", important=True,
                      help='Select the input images from the project.')
        form.addParam('numberOfClasses', IntParam, default=50, validators=[Positive],
                      label='Number of Classes',
                      help='The number of 2D classes into which to sort the dataset. '
                           'Runtime is approximately linear in the number of classes. '
                           'Often, sorting the same dataset into different numbers of '
                           'classes can be helpful in finding junk particles or rare views')
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
        print(self._program + "  \'do_run_class_2D(\""+self.a[-1].split()[-1]+"\", \""+self.b[-1].split()[-1]+"\", \"\'+"+self._user+"\'\", \"" + self.par + "\",\"\'"+str(self.numberOfClasses.get())+"\'\")\'")
        while commands.getstatusoutput(self._program + " \'get_job_status(\""+self.a[-1].split()[-1]+"\", \""+self.c[-1].split()[-1]+"\")\'")[-1].split()[-1] != 'completed':
            print("waiting...\n")
            commands.getstatusoutput(self._program + " \'wait_job_complete(\""+self.a[-1].split()[-1]+"\", \""+self.c[-1].split()[-1]+"\")\'")
        print("2D Classifications Started...")
        #do_run_class_2D(puid, wuid, uuid, particle_group, num_classes=50) returns the new uid of the job that was created
        self.d = commands.getstatusoutput(self._program + " \'do_run_class_2D(\""+self.a[-1].split()[-1]+"\", \""+self.b[-1].split()[-1]+"\", \"\'+"+self._user+"\'\", \"" + self.par + "\",\"\'"+str(self.numberOfClasses.get())+"\'\")\'")
        while commands.getstatusoutput(self._program + " \'get_job_status(\""+self.a[-1].split()[-1]+"\", \""+self.d[-1].split()[-1]+"\")\'")[-1].split()[-1] != 'completed':
            commands.getstatusoutput(self._program + " \'wait_job_complete(\""+self.a[-1].split()[-1]+"\", \""+self.d[-1].split()[-1]+"\")\'")



    
    def createOutputStep(self):
        self._program2 = os.path.join(os.environ['PYEM_DIR'], 'csparc2star.py')
        
        self.runJob(self._program2, self._ssd+'/'+self.a[-1].split()[-1]+'/'+self.d[-1].split()[-1]+"/cryosparc_"+self.a[-1].split()[-1]+"_"+self.d[-1].split()[-1]+"_020_particles.cs"+" "+self._getFileName('out_particles'), numberOfMpi=1)
        self.runJob(self._program2, self._ssd+'/'+self.a[-1].split()[-1]+'/'+self.d[-1].split()[-1]+"/cryosparc_"+self.a[-1].split()[-1]+"_"+self.d[-1].split()[-1]+"_020_class_averages.cs"+" "+self._getFileName('out_class'), numberOfMpi=1)
        # Link the folder on SSD to scipion directory
        os.system("ln -s "+self._ssd+"/"+self.a[-1].split()[-1]+'/'+self.d[-1].split()[-1]+" "+ self._getExtraPath())

        with open(self._getFileName('out_class'), 'r') as input_file, open(self._getFileName('out_class_m2'), 'w') as output_file:
            j = 0 #mutex lock
            i = 0 #start
            k = 1
            l = 0
            for line in input_file:
                if line.startswith("_rln"):
                    output_file.write(line)
                    i = 1
                elif i == 0:
                    output_file.write(line)
                elif j == 0:
                    for n, m in enumerate(line.split()):
                        if '@' in m:
                            break
                    output_file.write(" ".join(line.split()[:n]) + " " + line.split()[n].split('@')[0] + '@'+ self._getExtraPath() +"/"+ line.split()[n].split('@')[1] + " " + " ".join(line.split()[n+1:])+"\n")
                    j = 1
                else:
                    output_file.write(" ".join(line.split()[:n]) + " " + line.split()[n].split('@')[0] + '@'+ self._getExtraPath() +"/"+ line.split()[n].split('@')[1] + " " + " ".join(line.split()[n+1:])+"\n")
        
        self._loadClassesInfo(self._getFileName('out_class_m2'))
        inputParticles = self.inputParticles.get()
        classes2DSet = self._createSetOfClasses2D(inputParticles)
        self._fillClassesFromLevel(classes2DSet)
  
        self._defineOutputs(outputClasses=classes2DSet)
        self._defineSourceRelation(self.inputParticles, classes2DSet)


    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []

        return errors

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputClasses'):
            summary.append("Output classes not ready yet.")
        else:
            summary.append("Input Particles: %s" % self.getObjectTag('inputParticles'))
            summary.append("Classified into *%d* classes." % self.numberOfClasses.get())
            summary.append("Output set: %s" % self.getObjectTag('outputClasses'))

        return summary

    def _methods(self):
        methods = "We classified input particles %s (%d items) " % (
            self.getObjectTag('inputParticles'),
            self._getInputParticles().getSize())
        methods += "into %d classes using CryoSparc " % self.numberOfClasses.get()
        return [methods]
    
    # --------------------------- UTILS functions ------------------------------
    def _getInputParticles(self):
        return self.inputParticles.get()

    def _loadClassesInfo(self, filename):
        """ Read some information about the produced 2D classes
        from the metadata file.
        """
        self._classesInfo = {}  # store classes info, indexed by class id

        mdClasses = md.MetaData(filename)

        for classNumber, row in enumerate(md.iterRows(mdClasses)):
            index, fn = relionToLocation(row.getValue('rlnImageName'))
            # Store info indexed by id, we need to store the row.clone() since
            # the same reference is used for iteration
            self._classesInfo[classNumber + 1] = (index, fn, row.clone())
        self._numClass = index
        print("number of class is "+str(index))


    def _fillClassesFromLevel(self, clsSet):
        """ Create the SetOfClasses2D from a given iteration. """
        xmpMd = self._getFileName("out_particles") #the particle with orientation parameters (all_parameters)
        clsSet.classifyItems(updateItemCallback=self._updateParticle,
                                 updateClassCallback=self._updateClass,
                                 itemDataIterator=md.iterRows(xmpMd,
                                                          sortByLabel=md.MDL_ITEM_ID)) # relion style
    def _updateParticle(self, item, row):
        item.setClassId(row.getValue(md.RLN_PARTICLE_CLASS))
        item.setTransform(rowToAlignment(row, em.ALIGN_2D))
        
        
    def _updateClass(self, item):
        classId = item.getObjId()
        if classId in self._classesInfo:
            index, fn, row = self._classesInfo[classId]
            item.setAlignment2D()
            item.getRepresentative().setLocation(index, fn)
        


