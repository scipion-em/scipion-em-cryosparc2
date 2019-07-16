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

from pyworkflow.em.protocol import *
from pyworkflow.tests import *
from pyworkflow.utils import importFromPlugin

from cryosparc2.protocols import *
from cryosparc2.constants import *

relionProtocols = importFromPlugin('relion.protocols', doRaise=True)


class TestCryosparcBase(BaseTest):
    @classmethod
    def setData(cls, dataProject='xmipp_tutorial'):
        cls.dataset = DataSet.getDataSet(dataProject)
        cls.micFn = cls.dataset.getFile('allMics')
        cls.volFn = cls.dataset.getFile('vol2')
        cls.partFn1 = cls.dataset.getFile('particles2')
        cls.partFn2 = cls.dataset.getFile('particles3')
        cls.ctfFn = cls.dataset.getFile('ctf')

    @classmethod
    def runImportMicrograph(cls, pattern, samplingRate, voltage,
                            scannedPixelSize, magnification,
                            sphericalAberration):
        """ Run an Import micrograph protocol. """
        # We have two options: pass the SamplingRate or
        # the ScannedPixelSize + microscope magnification
        kwargs = {
            'filesPath': pattern,
            'magnification': magnification,
            'voltage': voltage,
            'sphericalAberration': sphericalAberration
        }

        if samplingRate is not None:
            kwargs.update({'samplingRateMode': 0,
                           'samplingRate': samplingRate})
        else:
            kwargs.update({'samplingRateMode': 1,
                           'scannedPixelSize': scannedPixelSize})

        cls.protImport = ProtImportMicrographs(**kwargs)
        cls.launchProtocol(cls.protImport)

        # Check that input micrographs have been imported
        if cls.protImport.outputMicrographs is None:
            raise Exception('Import of micrograph: %s, failed. '
                            'outputMicrographs is None.' % pattern)

        return cls.protImport

    @classmethod
    def runImportVolumes(cls, pattern, samplingRate,
                         importFrom=ProtImportParticles.IMPORT_FROM_FILES):
        """ Run an Import particles protocol. """
        cls.protImport = cls.newProtocol(ProtImportVolumes,
                                         filesPath=pattern,
                                         samplingRate=samplingRate
                                         )
        cls.launchProtocol(cls.protImport)
        return cls.protImport

    @classmethod
    def runImportParticles(cls, pattern, samplingRate, checkStack=False,
                           importFrom=ProtImportParticles.IMPORT_FROM_FILES):
        """ Run an Import particles protocol. """
        if importFrom == ProtImportParticles.IMPORT_FROM_SCIPION:
            objLabel = 'from scipion (particles)'
        elif importFrom == ProtImportParticles.IMPORT_FROM_FILES:
            objLabel = 'from file (particles)'

        cls.protImport = cls.newProtocol(ProtImportParticles,
                                         objLabel=objLabel,
                                         filesPath=pattern,
                                         sqliteFile=pattern,
                                         samplingRate=samplingRate,
                                         checkStack=checkStack,
                                         importFrom=importFrom)

        cls.launchProtocol(cls.protImport)
        # Check that input images have been imported (a better way to do this?)
        if cls.protImport.outputParticles is None:
            raise Exception('Import of images: %s, failed. '
                            'outputParticles is None.' % pattern)
        return cls.protImport

    @classmethod
    def runImportMicrographBPV(cls, pattern):
        """ Run an Import micrograph protocol. """
        return cls.runImportMicrograph(pattern,
                                       samplingRate=1.237,
                                       voltage=300,
                                       sphericalAberration=2,
                                       scannedPixelSize=None,
                                       magnification=56000)

    @classmethod
    def runImportMicrographRCT(cls, pattern):
        """ Run an Import micrograph protocol. """
        return cls.runImportMicrograph(pattern,
                                       samplingRate=2.28,
                                       voltage=100,
                                       sphericalAberration=2.9,
                                       scannedPixelSize=None,
                                       magnification=50000)

    @classmethod
    def runImportParticleCryoSPARC(cls, pattern):
        """ Run an Import micrograph protocol. """
        return cls.runImportParticles(pattern,
                                      samplingRate=4.,
                                      checkStack=True,
                                      importFrom=ProtImportParticles.IMPORT_FROM_SCIPION)

    @classmethod
    def runImportVolumesCryoSPARC(cls, pattern):
        """ Run an Import micrograph protocol. """
        return cls.runImportVolumes(pattern,
                                    samplingRate=4.,
                                    importFrom=ProtImportParticles.IMPORT_FROM_FILES)


class TestCryosparcClassify2D(TestCryosparcBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        setupTestProject(cls)
        dataProject = 'grigorieff'
        dataset = DataSet.getDataSet(dataProject)
        TestCryosparcBase.setData()
        particlesPattern = dataset.getFile('particles.sqlite')
        cls.protImport = cls.runImportParticleCryoSPARC(cls.partFn2)

    def testCryosparc2D(self):
        def _runCryosparcClassify2D(label=''):
            prot2D = self.newProtocol(ProtCryo2D,
                                      doCTF=False, maskDiameterA=340,
                                      numberOfMpi=4, numberOfThreads=1)

            # Normalization after the imported particles
            relionProtocol = self.newProtocol(relionProtocols.ProtRelionPreprocessParticles,
                                        doNormalize=True,
                                        doScale=True, scaleSize=140,
                                        doInvert=False)
            relionProtocol.setObjLabel('relion: preprocess particles')
            relionProtocol.inputParticles.set(self.protImport.outputParticles)
            self.launchProtocol(relionProtocol)

            prot2D.inputParticles.set(relionProtocol.outputParticles)
            prot2D.numberOfClasses.set(5)
            prot2D.numberOnlineEMIterator.set(40)
            prot2D.setObjLabel(label)
            prot2D.numberGPU.set(1)
            self.launchProtocol(prot2D)
            return prot2D

        def _checkAsserts(cryosparcProt):

            self.assertIsNotNone(cryosparcProt.outputClasses,
                                 "There was a problem with Cryosparc 2D classify")

            for class2D in cryosparcProt.outputClasses:
                self.assertTrue(class2D.hasAlignment2D())

        cryosparcProtGpu = _runCryosparcClassify2D(label="Cryosparc classify2D GPU")
        _checkAsserts(cryosparcProtGpu)


class TestCryosparc3DInitialModel(TestCryosparcBase):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        setupTestProject(cls)
        dataProject = 'grigorieff'
        dataset = DataSet.getDataSet(dataProject)
        TestCryosparcBase.setData()
        particlesPattern = dataset.getFile('particles.sqlite')
        cls.protImport = cls.runImportParticleCryoSPARC(cls.partFn2)

    def testCryosparcInitialModel(self):
        def _runCryosparcInitialModel(label=''):
            protInitialModel = self.newProtocol(ProtCryoSparcInitialModel,
                                      numberOfMpi=4, numberOfThreads=1)

            # Normalization after the imported particles
            relionProtocol = self.newProtocol(
                relionProtocols.ProtRelionPreprocessParticles,
                doNormalize=True,
                doScale=True, scaleSize=140,
                doInvert=False)
            relionProtocol.setObjLabel('relion: preprocess particles')
            relionProtocol.inputParticles.set(self.protImport.outputParticles)
            self.launchProtocol(relionProtocol)

            protInitialModel.inputParticles.set(relionProtocol.outputParticles)
            protInitialModel.abinit_K.set(1)
            protInitialModel.symmetryGroup.set(SYM_CYCLIC)
            protInitialModel.symmetryOrder.set(1)
            protInitialModel.setObjLabel(label)
            self.launchProtocol(protInitialModel)
            return protInitialModel

        def _checkAsserts(cryosparcProt):
            self.assertIsNotNone(cryosparcProt.outputClasses,
                                 "There was a problem with Cryosparc 3D initial model")

            self.assertIsNotNone(cryosparcProt.outputVolumes,
                                 "There was a problem with Cryosparc 3D initial model")

        cryosparcProtGpu = _runCryosparcInitialModel(label="Cryosparc 3D initial model")
        _checkAsserts(cryosparcProtGpu)


class TestCryosparc3DRefinement(TestCryosparcBase):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        setupTestProject(cls)
        dataProject = 'grigorieff'
        dataset = DataSet.getDataSet(dataProject)
        TestCryosparcBase.setData()
        particlesPattern = dataset.getFile('particles.sqlite')
        cls.protImport = cls.runImportParticleCryoSPARC(cls.partFn2)

    def testCryosparc3DRefinement(self):
        def _runCryosparctest3DRefinement(label=''):
            prot3DRefinement = self.newProtocol(ProtCryoSparcRefine3D,
                                      numberOfMpi=4, numberOfThreads=1)

            # Normalization after the imported particles
            relionProtocol = self.newProtocol(
                relionProtocols.ProtRelionPreprocessParticles,
                doNormalize=True,
                doScale=True, scaleSize=140,
                doInvert=False)
            relionProtocol.setObjLabel('relion: preprocess particles')
            relionProtocol.inputParticles.set(self.protImport.outputParticles)
            self.launchProtocol(relionProtocol)

            importVolumeProt = self.runImportVolumesCryoSPARC(self.volFn)

            prot3DRefinement.inputParticles.set(relionProtocol.outputParticles)
            prot3DRefinement.referenceVolume.set(importVolumeProt.outputVolume)
            prot3DRefinement.symmetryGroup.set(SYM_CYCLIC)
            prot3DRefinement.symmetryOrder.set(1)
            prot3DRefinement.setObjLabel(label)
            self.launchProtocol(prot3DRefinement)
            return prot3DRefinement

        def _checkAsserts(cryosparcProt):
            self.assertIsNotNone(cryosparcProt.outputVolume,
                                 "There was a problem with Cryosparc 3D refinement")

        cryosparcProtGpu = _runCryosparctest3DRefinement(label="Cryosparc 3D refinement")
        _checkAsserts(cryosparcProtGpu)


class TestCryosparcNonUniformRefine3D(TestCryosparcBase):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        setupTestProject(cls)
        dataProject = 'grigorieff'
        dataset = DataSet.getDataSet(dataProject)
        TestCryosparcBase.setData()
        particlesPattern = dataset.getFile('particles.sqlite')
        cls.protImport = cls.runImportParticleCryoSPARC(cls.partFn2)

    def testCryosparcNonUniformRefine3D(self):
        def _runCryosparctestNonUniformRefine3D(label=''):
            protNonUniform3DRefinement = self.newProtocol(ProtCryoSparcNonUniformRefine3D,
                                                          numberOfMpi=4, numberOfThreads=1)

            # Normalization after the imported particles
            relionProtocol = self.newProtocol(
                relionProtocols.ProtRelionPreprocessParticles,
                doNormalize=True,
                doScale=True, scaleSize=140,
                doInvert=False)
            relionProtocol.setObjLabel('relion: preprocess particles')
            relionProtocol.inputParticles.set(self.protImport.outputParticles)
            self.launchProtocol(relionProtocol)

            importVolumeProt = self.runImportVolumesCryoSPARC(self.volFn)

            protNonUniform3DRefinement.inputParticles.set(relionProtocol.outputParticles)
            protNonUniform3DRefinement.referenceVolume.set(importVolumeProt.outputVolume)
            protNonUniform3DRefinement.symmetryGroup.set(SYM_CYCLIC)
            protNonUniform3DRefinement.symmetryOrder.set(1)
            protNonUniform3DRefinement.setObjLabel(label)
            self.launchProtocol(protNonUniform3DRefinement)
            return protNonUniform3DRefinement

        def _checkAsserts(cryosparcProt):
            self.assertIsNotNone(cryosparcProt.outputVolume,
                                 "There was a problem with Cryosparc Non-Uniform 3D refinement")

        cryosparcProtGpu = _runCryosparctestNonUniformRefine3D(label="Cryosparc Non-Uniform 3D refinement")
        _checkAsserts(cryosparcProtGpu)


class TestCryosparcParticlesSubtract(TestCryosparcBase):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        setupTestProject(cls)
        dataProject = 'grigorieff'
        dataset = DataSet.getDataSet(dataProject)
        TestCryosparcBase.setData()
        particlesPattern = dataset.getFile('particles.sqlite')
        cls.protImport = cls.runImportParticleCryoSPARC(cls.partFn2)

    def testCryosparcParticlesSubtract(self):
        def _runCryosparctestParticlesSubtract(label=''):

            protParticlesSubtract = self.newProtocol(ProtCryoSparcSubtract,
                                                numberOfMpi=4,
                                                numberOfThreads=1)

            # Normalization after the imported particles
            relionProtocol = self.newProtocol(
                relionProtocols.ProtRelionPreprocessParticles,
                doNormalize=True,
                doScale=True, scaleSize=140,
                doInvert=False)
            relionProtocol.setObjLabel('relion: preprocess particles')
            relionProtocol.inputParticles.set(self.protImport.outputParticles)
            self.launchProtocol(relionProtocol)

            importVolumeProt = self.runImportVolumesCryoSPARC(self.volFn)

            prot3DRefinement = self.newProtocol(ProtCryoSparcRefine3D,
                                                numberOfMpi=4,
                                                numberOfThreads=1)
            prot3DRefinement.inputParticles.set(relionProtocol.outputParticles)
            prot3DRefinement.referenceVolume.set(importVolumeProt.outputVolume)
            prot3DRefinement.symmetryGroup.set(SYM_CYCLIC)
            prot3DRefinement.symmetryOrder.set(1)
            self.launchProtocol(prot3DRefinement)

            protRelionCreate3DMask = self.newProtocol(
                relionProtocols.ProtRelionCreateMask3D,
                initialLowPassFilterA=20)
            protRelionCreate3DMask.inputVolume.set(prot3DRefinement.outputVolume)
            protRelionCreate3DMask.setObjLabel('relion: create 3d mask')
            self.launchProtocol(protRelionCreate3DMask)

            protParticlesSubtract.inputParticles.set(prot3DRefinement.outputParticles)
            protParticlesSubtract.refVolume.set(prot3DRefinement.outputVolume)
            protParticlesSubtract.refMask.set(protRelionCreate3DMask.outputMask)
            self.launchProtocol(protParticlesSubtract)

            return protParticlesSubtract

        def _checkAsserts(cryosparcProt):
            self.assertIsNotNone(cryosparcProt.outputParticles,
                                 "There was a problem with Cryosparc subtract projection")

        cryosparcProtGpu = _runCryosparctestParticlesSubtract(label="Cryosparc Subtract projection")
        _checkAsserts(cryosparcProtGpu)







