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

import os
import commands
import pyworkflow.em.metadata as md
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        LabelParam, IntParam,
                                        EnumParam, StringParam,
                                        BooleanParam, PathParam,
                                        LEVEL_ADVANCED)
from pyworkflow.em.data import Volume
from pyworkflow.em.protocol import ProtInitialVolume
from pyworkflow.utils import importFromPlugin
from cryosparc2.utils import *

relionPlugin = importFromPlugin("relion.convert", doRaise=True)


class ProtCryoSparcInitialModel(ProtInitialVolume):
    """    
    Generate a 3D initial model _de novo_ from 2D particles using
    CryoSparc Stochastic Gradient Descent (SGD) algorithm.
    """
    _label = 'CryoSparc 3D initial model'

    # --------------------------- DEFINE param functions ----------------------
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
        relionPlugin.writeSetOfParticles(imgSet,
                                         self._getFileName('input_particles'),
                                         outputDir=self._getExtraPath(),
                                         fillMagnification=True)
        self._importParticles()

    def processStep(self):

        while self.getJobStatus(self.importedParticles) != 'completed':
            print("waiting...\n")
            self.waitJob(self.importedParticles)

        print("Ab Initial Model Generation Started...")

        self.runAbinit = self.doRunAbinit()[-1].split()[-1]

        while self.getJobStatus(self.runAbinit) != 'completed':
            self.waitJob(self.runAbinit)

    def createOutputStep(self):

        _program2 = os.path.join(os.environ['PYEM_DIR'], 'csparc2star.py')

        self.runJob(_program2, self._ssd + '/' + self.projectName + '/' +
                    self.runAbinit + "/cryosparc_" +
                    self.projectName + "_" + self.runAbinit +
                    "_class_00_final_particles.cs" + " " +
                    self._getFileName('out_particles'), numberOfMpi=1)

        # Link the folder on SSD to scipion directory
        os.system("ln -s " + self._ssd + "/" + self.projectName + '/' +
                  self.runAbinit + " " + self._getExtraPath())

       
        imgSet = self._getInputParticles()
        vol = Volume()
        fnVol = self._getExtraPath() + "/" + self.runAbinit + "/cryosparc_" +\
                self.projectName+"_"+self.runAbinit+"_class_00_final_volume.mrc"
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

    # --------------------------- UTILS functions ---------------------------

    def _getInputParticles(self):
        return self.inputParticles.get()

    def _fillDataFromIter(self, imgSet):
        outImgsFn = self._getFileName('out_particles')
        imgSet.setAlignmentProj()
        imgSet.copyItems(self._getInputParticles(),
                         updateItemCallback=self._createItemMatrix,
                         itemDataIterator=md.iterRows(outImgsFn,
                                                      sortByLabel=md.RLN_IMAGE_ID))

    def _createItemMatrix(self, item, row):
        from pyworkflow.em import ALIGN_PROJ

        relionPlugin.createItemMatrix(item, row, align=ALIGN_PROJ)

    def _importParticles(self):
        """
        Initialize all utils cryoSPARC variables
        """
        self._program = getCryosparcProgram()
        self._user = getCryosparcUser()
        self._ssd = getCryosparcSSD()
        print("Importing Particles")

        # create empty project
        self.a = createEmptyProject()
        self.projectName = self.a[-1].split()[-1]

        # create empty workspace
        self.b = createEmptyWorkSpace(self.projectName)
        self.workSpaceName = self.b[-1].split()[-1]

        # import_particles_star
        self.c = self.doImportParticlesStar()

        self.importedParticles = self.c[-1].split()[-1]
        self.par = self.importedParticles + '.imported_particles'

    def getJobStatus(self, job):
        """
        Return the job status
        """
        return commands.getstatusoutput(self._program +
                                        " \'get_job_status(\""+
                                        self.projectName+"\", \""+
                                        job+"\")\'")[-1].split()[-1]

    def waitJob(self, job):
        commands.getstatusoutput(self._program +
                                 " \'wait_job_complete(\"" +
                                 self.projectName + "\", \"" +
                                 job + "\")\'")

    def doImportParticlesStar(self):
        """
        do_import_particles_star(puid, wuid, uuid, abs_star_path,
                                 abs_blob_path=None, psize_A=None)
        returns the new uid of the job that was created
        """
        return commands.getstatusoutput(self._program +
                                        " \'do_import_particles_star(\"" +
                                        self.projectName + "\", \"" +
                                        self.workSpaceName + "\", \"\'+" +
                                        self._user + "\'\", \"\'" +
                                        os.path.join(os.getcwd(),
                                                     self._getFileName('input_particles')) +
                                        "\'\", \"\'" + os.path.join(os.getcwd(),
                                                                    self._getExtraPath()) +
                                        "\'\", \"\'" +
                                        str(self._getInputParticles().getSamplingRate()) +
                                        "\'\")\'")

    def doRunAbinit(self):
        """self._program + "  \'do_run_abinit(\"" + self.projectName +
        "\", \"" + self.workSpaceName + "\", \"\'" + self._user + "\'\", \""
        + self.par + "\",\"\'1\'\")\'")
        """
        return commands.getstatusoutput(self._program + " \'do_run_abinit(\"" +
                                        self.projectName + "\", \"" +
                                        self.workSpaceName + "\", \"\'" +
                                        self._user + "\'\", \"" + self.par +
                                        "\",\"\'1\'\")\'")

