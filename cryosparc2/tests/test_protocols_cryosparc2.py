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

from pwem.protocols import *
from pyworkflow.protocol import PointerList
from pyworkflow.tests import *
from pwem import Domain
xmippProtocols = Domain.importFromPlugin('xmipp3.protocols', doRaise=True)

from ..protocols import *
from ..constants import *


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

        protImportMic = ProtImportMicrographs(**kwargs)
        cls.launchProtocol(protImportMic)

        # Check that input micrographs have been imported
        if protImportMic.outputMicrographs is None:
            raise Exception('Import of micrograph: %s, failed. '
                            'outputMicrographs is None.' % pattern)

        return protImportMic

    @classmethod
    def runImportVolumes(cls, pattern, samplingRate,
                         importFrom=ProtImportParticles.IMPORT_FROM_FILES):
        """ Run an Import volumes protocol. """
        protImportVol = cls.newProtocol(ProtImportVolumes,
                                         filesPath=pattern,
                                         samplingRate=samplingRate
                                         )
        cls.launchProtocol(protImportVol)
        return protImportVol

    @classmethod
    def runImportParticles(cls, pattern, samplingRate, checkStack=False,
                           importFrom=ProtImportParticles.IMPORT_FROM_FILES):
        """ Run an Import particles protocol. """
        if importFrom == ProtImportParticles.IMPORT_FROM_SCIPION:
            objLabel = 'from scipion (particles)'
        elif importFrom == ProtImportParticles.IMPORT_FROM_FILES:
            objLabel = 'from file (particles)'

        protImportPart = cls.newProtocol(ProtImportParticles,
                                         objLabel=objLabel,
                                         filesPath=pattern,
                                         sqliteFile=pattern,
                                         samplingRate=samplingRate,
                                         checkStack=checkStack,
                                         importFrom=importFrom)

        cls.launchProtocol(protImportPart)
        # Check that input images have been imported (a better way to do this?)
        if protImportPart.outputParticles is None:
            raise Exception('Import of images: %s, failed. '
                            'outputParticles is None.' % pattern)
        return protImportPart

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

    @classmethod
    def runCreate3DMask(cls, pattern):
        """ Create a volume mask using xmipp. """
        cls.msk = cls.newProtocol(xmippProtocols.XmippProtCreateMask3D,
                                  inputVolume=pattern,
                                  volumeOperation=0,  # OPERATION_THRESHOLD,
                                  doSmall=False,
                                  smallSize=False,
                                  doBig=False,
                                  doSymmetrize=False,
                                  doMorphological=False,
                                  doInvert=False,
                                  doSmooth=False,
                                  sigmaConvolution=2
                                  )
        cls.launchProtocol(cls.msk)
        return cls.msk

    @classmethod
    def runResizeParticles(cls, particles, doResize, resizeOption, resizeDim):
        cls.protResize = cls.newProtocol(xmippProtocols.XmippProtCropResizeParticles,
                                         doResize=doResize,
                                         resizeOption=resizeOption,
                                         resizeDim=resizeDim)
        cls.protResize.inputParticles.set(particles)
        cls.launchProtocol(cls.protResize)
        return cls.protResize


class TestCryosparcClassify2D(TestCryosparcBase):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        setupTestProject(cls)
        dataProject = 'grigorieff'
        dataset = DataSet.getDataSet(dataProject)
        TestCryosparcBase.setData()
        cls.protImportPart = cls.runImportParticleCryoSPARC(cls.partFn2)

    def testCryosparc2D(self):
        def _runCryosparcClassify2D(label=''):
            prot2D = self.newProtocol(ProtCryo2D,
                                      doCTF=False, maskDiameterA=340)

            prot2D.inputParticles.set(self.protImportPart.outputParticles)
            prot2D.numberOfClasses.set(5)
            prot2D.numberOnlineEMIterator.set(20)
            prot2D.compute_use_ssd.set(False)
            prot2D.setObjLabel(label)
            self.launchProtocol(prot2D)

            # Check the output, the alignment, the sampling rate and the
            # dimensions
            self.assertIsNotNone(prot2D.outputClasses)
            for class2D in prot2D.outputClasses:
                self.assertTrue(class2D.hasAlignment2D())
                self.assertTrue(class2D.getSamplingRate(),
                                self.protImportPart.outputParticles.getSamplingRate())
                self.assertTrue(class2D.getRepresentative().getDimensions(),
                                self.protImportPart.outputParticles.getDimensions())

            # Resize the particles in order to classify them and make a
            # downsampling to original particles
            protResize = self.runResizeParticles(self.protImportPart.outputParticles,
                                                 True,
                                                 xmippProtocols.XmippResizeHelper.RESIZE_DIMENSIONS,
                                                 400)

            prot2D = self.newProtocol(ProtCryo2D, doCTF=False,
                                      maskDiameterA=340)

            prot2D.inputParticles.set(protResize.outputParticles)
            prot2D.numberOfClasses.set(5)
            prot2D.numberOnlineEMIterator.set(20)
            prot2D.compute_use_ssd.set(False)
            prot2D.setObjLabel(label)
            self.launchProtocol(prot2D)

            # Check the output, the alignment, the sampling rate and the
            # dimensions
            self.assertIsNotNone(prot2D.outputClasses)
            for class2D in prot2D.outputClasses:
                self.assertTrue(class2D.hasAlignment2D())
                self.assertEqual(class2D.getSamplingRate(),
                                protResize.outputParticles.getSamplingRate())
                self.assertEqual(class2D.getRepresentative().getDimensions(),
                                protResize.outputParticles.getDimensions())

        _runCryosparcClassify2D(label="Cryosparc classify2D GPU")


class TestCryosparc3DInitialModel(TestCryosparcBase):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        setupTestProject(cls)
        dataProject = 'grigorieff'
        dataset = DataSet.getDataSet(dataProject)
        TestCryosparcBase.setData()
        particlesPattern = dataset.getFile('particles.sqlite')
        cls.protImportPart = cls.runImportParticleCryoSPARC(cls.partFn2)

    def testCryosparcInitialModel(self):
        def _runCryosparcInitialModel(label=''):
            protInitialModel = self.newProtocol(ProtCryoSparcInitialModel)

            protInitialModel.inputParticles.set(self.protImportPart.outputParticles)
            protInitialModel.abinit_K.set(2)
            protInitialModel.symmetryGroup.set(SYM_CYCLIC)
            protInitialModel.symmetryOrder.set(1)
            protInitialModel.compute_use_ssd.set(False)
            protInitialModel.setObjLabel(label)
            self.launchProtocol(protInitialModel)

            # Check the outputs
            self.assertIsNotNone(protInitialModel.outputClasses)
            self.assertIsNotNone(protInitialModel.outputVolumes)

            # Check the number of classes
            self.assertEqual(protInitialModel.outputClasses.getSize(),
                             protInitialModel.abinit_K.get())

            # Check the number of volumes
            self.assertEqual(protInitialModel.outputVolumes.getSize(),
                             protInitialModel.abinit_K.get())

            # Check the volumes alignment and the sampling rate
            for volume in protInitialModel.outputVolumes:
                self.assertEqual(volume.getSamplingRate(),
                                 self.protImportPart.outputParticles.getSamplingRate())

        _runCryosparcInitialModel(label="Cryosparc 3D initial model")


class TestCryosparc3DRefinement(TestCryosparcBase):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        setupTestProject(cls)
        dataProject = 'grigorieff'
        dataset = DataSet.getDataSet(dataProject)
        TestCryosparcBase.setData()
        particlesPattern = dataset.getFile('particles.sqlite')
        cls.protImportPart = cls.runImportParticleCryoSPARC(cls.partFn2)
        cls.protImportVol = cls.runImportVolumesCryoSPARC(cls.volFn)

    def testCryosparc3DRefinement(self):
        def _runCryosparctest3DRefinement(label=''):
            prot3DRefinement = self.newProtocol(ProtCryoSparcRefine3D)

            prot3DRefinement.inputParticles.set(self.protImportPart.outputParticles)
            prot3DRefinement.referenceVolume.set(self.protImportVol.outputVolume)
            prot3DRefinement.symmetryGroup.set(SYM_CYCLIC)
            prot3DRefinement.symmetryOrder.set(1)
            prot3DRefinement.compute_use_ssd.set(False)
            prot3DRefinement.setObjLabel(label)
            self.launchProtocol(prot3DRefinement)

            # Check the outputs, sampling rate,...
            self.assertIsNotNone(prot3DRefinement.outputVolume)
            self.assertEqual(prot3DRefinement.outputVolume.getSamplingRate(),
                             self.protImportPart.outputParticles.getSamplingRate())

            outputParticles = prot3DRefinement.outputParticles
            self.assertIsNotNone(outputParticles)
            self.assertEqual(outputParticles.getSamplingRate(),
                             self.protImportPart.outputParticles.getSamplingRate())
            self.assertTrue(outputParticles.hasAlignmentProj())
            self.assertTrue(outputParticles.hasCTF())
            self.assertEqual(outputParticles.getSize(),
                             self.protImportPart.outputParticles.getSize())

        _runCryosparctest3DRefinement(label="Cryosparc 3D refinement")


class TestProtCryoSparc3DHomogeneousRefine(TestCryosparcBase):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        setupTestProject(cls)
        dataProject = 'grigorieff'
        dataset = DataSet.getDataSet(dataProject)
        TestCryosparcBase.setData()
        particlesPattern = dataset.getFile('particles.sqlite')
        cls.protImportPart = cls.runImportParticleCryoSPARC(cls.partFn2)
        cls.protImportVol = cls.runImportVolumesCryoSPARC(cls.volFn)

    def testCryosparc3DHomogeneousRefinement(self):
        def _runCryosparctest3DHomogeneousRefinement(label=''):
            prot3DHomoRefinement = self.newProtocol(ProtCryoSparc3DHomogeneousRefine)

            prot3DHomoRefinement.inputParticles.set(self.protImportPart.outputParticles)
            prot3DHomoRefinement.referenceVolume.set(self.protImportVol.outputVolume)
            prot3DHomoRefinement.symmetryGroup.set(SYM_CYCLIC)
            prot3DHomoRefinement.symmetryOrder.set(1)
            prot3DHomoRefinement.compute_use_ssd.set(False)
            prot3DHomoRefinement.setObjLabel(label)
            self.launchProtocol(prot3DHomoRefinement)

            # Check the outputs, sampling rate,...
            self.assertIsNotNone(prot3DHomoRefinement.outputVolume)
            self.assertEqual(prot3DHomoRefinement.outputVolume.getSamplingRate(),
                             self.protImportPart.outputParticles.getSamplingRate())

            outputParticles = prot3DHomoRefinement.outputParticles
            self.assertIsNotNone(outputParticles)
            self.assertEqual(outputParticles.getSamplingRate(),
                             self.protImportPart.outputParticles.getSamplingRate())
            self.assertTrue(outputParticles.hasAlignmentProj())
            self.assertTrue(outputParticles.hasCTF())
            self.assertEqual(outputParticles.getSize(),
                             self.protImportPart.outputParticles.getSize())

        _runCryosparctest3DHomogeneousRefinement(label="Cryosparc 3D homogeneous refinement")


class TestCryosparcNonUniformRefine3D(TestCryosparcBase):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        setupTestProject(cls)
        dataProject = 'grigorieff'
        dataset = DataSet.getDataSet(dataProject)
        TestCryosparcBase.setData()
        particlesPattern = dataset.getFile('particles.sqlite')
        cls.protImportPart = cls.runImportParticleCryoSPARC(cls.partFn2)
        cls.protImportVolumeVol = cls.runImportVolumesCryoSPARC(cls.volFn)

    def testCryosparcNonUniformRefine3D(self):
        def _runCryosparctestNonUniformRefine3D(label=''):
            protNonUniform3DRefinement = self.newProtocol(ProtCryoSparcNonUniformRefine3D)
            protNonUniform3DRefinement.inputParticles.set(self.protImportPart.outputParticles)
            protNonUniform3DRefinement.referenceVolume.set(self.protImportVolumeVol.outputVolume)
            protNonUniform3DRefinement.symmetryGroup.set(SYM_CYCLIC)
            protNonUniform3DRefinement.symmetryOrder.set(1)
            protNonUniform3DRefinement.compute_use_ssd.set(False)
            protNonUniform3DRefinement.setObjLabel(label)
            self.launchProtocol(protNonUniform3DRefinement)

            # Check the outputs, dimensions, sampling rate,...
            self.assertIsNotNone(protNonUniform3DRefinement.outputVolume)
            self.assertEqual(protNonUniform3DRefinement.outputVolume.getSamplingRate(),
                             self.protImportPart.outputParticles.getSamplingRate())

            outputParticles = protNonUniform3DRefinement.outputParticles
            self.assertIsNotNone(outputParticles)
            self.assertEqual(outputParticles.getSamplingRate(),
                             self.protImportPart.outputParticles.getSamplingRate())
            self.assertTrue(outputParticles.hasAlignmentProj())
            self.assertTrue(outputParticles.hasCTF())
            self.assertEqual(outputParticles.getSize(),
                             self.protImportPart.outputParticles.getSize())

        _runCryosparctestNonUniformRefine3D(label="Cryosparc Non-Uniform 3D refinement")


class TestCryosparcNewNonUniformRefine3D(TestCryosparcBase):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        setupTestProject(cls)
        dataProject = 'grigorieff'
        dataset = DataSet.getDataSet(dataProject)
        TestCryosparcBase.setData()
        particlesPattern = dataset.getFile('particles.sqlite')
        cls.protImportPart = cls.runImportParticleCryoSPARC(cls.partFn2)
        cls.protImportVolumeVol = cls.runImportVolumesCryoSPARC(cls.volFn)

    def testCryosparcNewNonUniformRefine3D(self):
        def _runCryosparctestNewNonUniformRefine3D(label=''):
            protNewNonUniform3DRefinement = self.newProtocol(ProtCryoSparcNewNonUniformRefine3D)
            protNewNonUniform3DRefinement.inputParticles.set(self.protImportPart.outputParticles)
            protNewNonUniform3DRefinement.referenceVolume.set(self.protImportVolumeVol.outputVolume)
            protNewNonUniform3DRefinement.symmetryGroup.set(SYM_CYCLIC)
            protNewNonUniform3DRefinement.symmetryOrder.set(1)
            protNewNonUniform3DRefinement.compute_use_ssd.set(False)
            protNewNonUniform3DRefinement.setObjLabel(label)
            self.launchProtocol(protNewNonUniform3DRefinement)

            # Check the outputs, dimensions, sampling rate,...
            self.assertIsNotNone(protNewNonUniform3DRefinement.outputVolume)
            self.assertEqual(protNewNonUniform3DRefinement.outputVolume.getSamplingRate(),
                             self.protImportPart.outputParticles.getSamplingRate())

            outputParticles = protNewNonUniform3DRefinement.outputParticles
            self.assertIsNotNone(outputParticles)
            self.assertEqual(outputParticles.getSamplingRate(),
                             self.protImportPart.outputParticles.getSamplingRate())
            self.assertTrue(outputParticles.hasAlignmentProj())
            self.assertTrue(outputParticles.hasCTF())
            self.assertEqual(outputParticles.getSize(),
                             self.protImportPart.outputParticles.getSize())

        _runCryosparctestNewNonUniformRefine3D(label="Cryosparc New Non-Uniform 3D refinement")


class TestCryosparcHelicalRefine3D(TestCryosparcBase):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        setupTestProject(cls)
        dataProject = 'grigorieff'
        dataset = DataSet.getDataSet(dataProject)
        TestCryosparcBase.setData()
        cls.protImportPart = cls.runImportParticleCryoSPARC(cls.partFn2)
        cls.protImportVolumeVol = cls.runImportVolumesCryoSPARC(cls.volFn)

    def testCryosparcHelicalRefine3D(self):
        def _runCryosparctestHelicalRefine3D(label=''):
            protHelical3DRefinement = self.newProtocol(ProtCryoSparcHelicalRefine3D)
            protHelical3DRefinement.inputParticles.set(self.protImportPart.outputParticles)
            protHelical3DRefinement.referenceVolume.set(self.protImportVolumeVol.outputVolume)
            protHelical3DRefinement.symmetryGroup.set(SYM_CYCLIC)
            protHelical3DRefinement.symmetryOrder.set(1)
            protHelical3DRefinement.compute_use_ssd.set(False)
            protHelical3DRefinement.setObjLabel(label)
            self.launchProtocol(protHelical3DRefinement)

            # Check the outputs, dimensions, sampling rate,...
            self.assertIsNotNone(protHelical3DRefinement.outputVolume)
            self.assertEqual(protHelical3DRefinement.outputVolume.getSamplingRate(),
                             self.protImportPart.outputParticles.getSamplingRate())

            outputParticles = protHelical3DRefinement.outputParticles
            self.assertIsNotNone(outputParticles)
            self.assertEqual(outputParticles.getSamplingRate(),
                             self.protImportPart.outputParticles.getSamplingRate())
            self.assertTrue(outputParticles.hasAlignmentProj())
            self.assertTrue(outputParticles.hasCTF())
            self.assertEqual(outputParticles.getSize(),
                             self.protImportPart.outputParticles.getSize())

        _runCryosparctestHelicalRefine3D(label="Cryosparc Helical 3D refinement")


class TestCryosparc3DClassification(TestCryosparcBase):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        setupTestProject(cls)
        dataProject = 'grigorieff'
        dataset = DataSet.getDataSet(dataProject)
        TestCryosparcBase.setData()
        particlesPattern = dataset.getFile('particles.sqlite')
        cls.protImportPart = cls.runImportParticleCryoSPARC(cls.partFn2)
        cls.protImportVolumeVol = cls.runImportVolumesCryoSPARC(cls.volFn)

    def testCryosparc3DClassification(self):
        def _runCryosparctest3DClassification(label=''):
            # Launch a first refinement protocol
            protNonUniform3DRefinement = self.newProtocol(ProtCryoSparcNonUniformRefine3D)


            protNonUniform3DRefinement.inputParticles.set(self.protImportPart.outputParticles)
            protNonUniform3DRefinement.referenceVolume.set(self.protImportVolumeVol.outputVolume)
            protNonUniform3DRefinement.symmetryGroup.set(SYM_CYCLIC)
            protNonUniform3DRefinement.symmetryOrder.set(1)
            protNonUniform3DRefinement.compute_use_ssd.set(False)
            protNonUniform3DRefinement.setObjLabel("protNonUniform3DRefinement_1")
            self.launchProtocol(protNonUniform3DRefinement)
            self.assertIsNotNone(protNonUniform3DRefinement.outputVolume)

            # Launch a second refinement protocol
            protNonUniform3DRefinement1 = self.newProtocol(ProtCryoSparcNonUniformRefine3D)

            protNonUniform3DRefinement1.inputParticles.set(
                self.protImportPart.outputParticles)
            protNonUniform3DRefinement1.referenceVolume.set(
                self.protImportVolumeVol.outputVolume)
            protNonUniform3DRefinement1.symmetryGroup.set(SYM_CYCLIC)
            protNonUniform3DRefinement1.symmetryOrder.set(1)
            protNonUniform3DRefinement1.compute_use_ssd.set(False)
            protNonUniform3DRefinement1.setObjLabel("protNonUniform3DRefinement_2")
            self.launchProtocol(protNonUniform3DRefinement1)
            self.assertIsNotNone(protNonUniform3DRefinement1.outputVolume)

            # Launch a 3D classification protocol
            prot3DClassification = self.newProtocol(ProtCryoSparc3DClassification)
            prot3DClassification.inputParticles.set(protNonUniform3DRefinement.outputParticles)
            volumes = PointerList()
            volumes.append(protNonUniform3DRefinement.outputVolume)
            volumes.append(protNonUniform3DRefinement1.outputVolume)
            prot3DClassification.refVolumes.set(volumes)
            prot3DClassification.compute_use_ssd.set(False)
            self.launchProtocol(prot3DClassification)

            # Check the outputs
            self.assertIsNotNone(prot3DClassification.outputClasses)
            self.assertIsNotNone(prot3DClassification.outputVolumes)

            # Check the number of classes
            self.assertEqual(prot3DClassification.outputClasses.getSize(), 2)

            # Check the number of classes
            self.assertEqual(prot3DClassification.outputVolumes.getSize(), 2)

        _runCryosparctest3DClassification(label="Cryosparc 3D Classification")


class TestCryosparcParticlesSubtract(TestCryosparcBase):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        setupTestProject(cls)
        dataProject = 'grigorieff'
        dataset = DataSet.getDataSet(dataProject)
        TestCryosparcBase.setData()
        particlesPattern = dataset.getFile('particles.sqlite')
        cls.protImportPart = cls.runImportParticleCryoSPARC(cls.partFn2)
        cls.protImportVol = cls.runImportVolumesCryoSPARC(cls.volFn)

    def testCryosparcParticlesSubtract(self):
        def _runCryosparctestParticlesSubtract(label=''):

            protParticlesSubtract = self.newProtocol(ProtCryoSparcSubtract)

            prot3DRefinement = self.newProtocol(ProtCryoSparcRefine3D)
            prot3DRefinement.inputParticles.set(self.protImportPart.outputParticles)
            prot3DRefinement.referenceVolume.set(self.protImportVol.outputVolume)
            prot3DRefinement.symmetryGroup.set(SYM_CYCLIC)
            prot3DRefinement.compute_use_ssd.set(False)
            prot3DRefinement.symmetryOrder.set(1)
            self.launchProtocol(prot3DRefinement)
            self.assertIsNotNone(prot3DRefinement.outputVolume)

            # Create a 3D Mask using xmipp
            protXmippCreate3DMask = self.runCreate3DMask(prot3DRefinement.outputVolume)

            protParticlesSubtract.inputParticles.set(prot3DRefinement.outputParticles)
            protParticlesSubtract.refVolume.set(prot3DRefinement.outputVolume)
            protParticlesSubtract.refMask.set(protXmippCreate3DMask.outputMask)
            protParticlesSubtract.compute_use_ssd.set(False)
            self.launchProtocol(protParticlesSubtract)

            self.assertIsNotNone(protParticlesSubtract.outputParticles)
            self.assertEqual(protParticlesSubtract.outputParticles.getSize(),
                             prot3DRefinement.outputParticles.getSize())

        _runCryosparctestParticlesSubtract(label="Cryosparc Subtract projection")


class TestCryosparcSharppening(TestCryosparcBase):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        setupTestProject(cls)
        dataProject = 'grigorieff'
        dataset = DataSet.getDataSet(dataProject)
        TestCryosparcBase.setData()
        particlesPattern = dataset.getFile('particles.sqlite')
        cls.protImportPart = cls.runImportParticleCryoSPARC(cls.partFn2)
        cls.protImportVol = cls.runImportVolumesCryoSPARC(cls.volFn)

    def testCryosparcParticlesSubtract(self):
        def _runCryosparctestSharppening(label=''):

            protSharppening = self.newProtocol(ProtCryoSparcSharppening)


            prot3DRefinement = self.newProtocol(ProtCryoSparcRefine3D)
            prot3DRefinement.inputParticles.set(self.protImportPart.outputParticles)
            prot3DRefinement.referenceVolume.set(self.protImportVol.outputVolume)
            prot3DRefinement.symmetryGroup.set(SYM_CYCLIC)
            prot3DRefinement.symmetryOrder.set(1)
            prot3DRefinement.compute_use_ssd.set(False)
            self.launchProtocol(prot3DRefinement)
            
            protSharppening.refVolume.set(prot3DRefinement.outputVolume)
            protSharppening.sharp_bfactor.set(-80)
            protSharppening.compute_use_ssd.set(False)
            self.launchProtocol(protSharppening)

            return protSharppening

        def _checkAsserts(cryosparcProt):
            self.assertIsNotNone(cryosparcProt.outputVolume,
                                 "There was a problem with Cryosparc sharppeninh")

        cryosparcProtGpu = _runCryosparctestSharppening(label="Cryosparc Sharppening")
        _checkAsserts(cryosparcProtGpu)


class TestCryosparcGlobalCtfRefinement(TestCryosparcBase):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        setupTestProject(cls)
        dataProject = 'grigorieff'
        dataset = DataSet.getDataSet(dataProject)
        TestCryosparcBase.setData()
        particlesPattern = dataset.getFile('particles.sqlite')
        cls.protImportPart = cls.runImportParticleCryoSPARC(cls.partFn2)
        cls.protImportVol = cls.runImportVolumesCryoSPARC(cls.volFn)

    def testCryosparcGlobalCtfRefinement(self):
        def _runCryosparctestGlobalCtfRefinement(label=''):

            protGlobalCtfRefinement = self.newProtocol(ProtCryoSparcGlobalCtfRefinement)

            prot3DRefinement = self.newProtocol(ProtCryoSparcRefine3D)
            prot3DRefinement.inputParticles.set(self.protImportPart.outputParticles)
            prot3DRefinement.referenceVolume.set(self.protImportVol.outputVolume)
            prot3DRefinement.symmetryGroup.set(SYM_CYCLIC)
            prot3DRefinement.symmetryOrder.set(1)
            prot3DRefinement.compute_use_ssd.set(False)
            self.launchProtocol(prot3DRefinement)

            # Create a 3D Mask using xmipp
            protXmippCreate3DMask = self.runCreate3DMask(prot3DRefinement.outputVolume)

            protGlobalCtfRefinement.inputParticles.set(prot3DRefinement.outputParticles)
            protGlobalCtfRefinement.refVolume.set(prot3DRefinement.outputVolume)
            protGlobalCtfRefinement.refMask.set(protXmippCreate3DMask.outputMask)
            protGlobalCtfRefinement.compute_use_ssd.set(False)
            self.launchProtocol(protGlobalCtfRefinement)

            return protGlobalCtfRefinement

        def _checkAsserts(cryosparcProt):
            self.assertIsNotNone(cryosparcProt.outputParticles,
                                 "There was a problem with Cryosparc Ctf Refinement")

        cryosparcProtGpu = _runCryosparctestGlobalCtfRefinement(label="Cryosparc Global Ctf Refinement")
        _checkAsserts(cryosparcProtGpu)


class TestCryosparcLocalCtfRefinement(TestCryosparcBase):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        setupTestProject(cls)
        dataProject = 'grigorieff'
        dataset = DataSet.getDataSet(dataProject)
        TestCryosparcBase.setData()
        particlesPattern = dataset.getFile('particles.sqlite')
        cls.protImportPart = cls.runImportParticleCryoSPARC(cls.partFn2)
        cls.protImportVol = cls.runImportVolumesCryoSPARC(cls.volFn)

    def testCryosparcLocalCtfRefinement(self):
        def _runCryosparctestLocalCtfRefinement(label=''):

            protLocalCtfRefinement = self.newProtocol(ProtCryoSparcLocalCtfRefinement)

            prot3DRefinement = self.newProtocol(ProtCryoSparcRefine3D)
            prot3DRefinement.inputParticles.set(self.protImportPart.outputParticles)
            prot3DRefinement.referenceVolume.set(self.protImportVol.outputVolume)
            prot3DRefinement.symmetryGroup.set(SYM_CYCLIC)
            prot3DRefinement.symmetryOrder.set(1)
            prot3DRefinement.compute_use_ssd.set(False)
            self.launchProtocol(prot3DRefinement)

            # Create a 3D Mask using xmipp
            protXmippCreate3DMask = self.runCreate3DMask(prot3DRefinement.outputVolume)

            protLocalCtfRefinement.inputParticles.set(prot3DRefinement.outputParticles)
            protLocalCtfRefinement.refVolume.set(prot3DRefinement.outputVolume)
            protLocalCtfRefinement.refMask.set(protXmippCreate3DMask.outputMask)
            protLocalCtfRefinement.compute_use_ssd.set(False)
            self.launchProtocol(protLocalCtfRefinement)

            return protLocalCtfRefinement

        def _checkAsserts(cryosparcProt):
            self.assertIsNotNone(cryosparcProt.outputParticles,
                                 "There was a problem with Cryosparc Ctf Refinement")

        cryosparcProtGpu = _runCryosparctestLocalCtfRefinement(label="Cryosparc Local Ctf Refinement")
        _checkAsserts(cryosparcProtGpu)


class TestCryosparcSymetryExpansion(TestCryosparcBase):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        setupTestProject(cls)
        dataProject = 'grigorieff'
        dataset = DataSet.getDataSet(dataProject)
        TestCryosparcBase.setData()
        particlesPattern = dataset.getFile('particles.sqlite')
        cls.protImportPart = cls.runImportParticleCryoSPARC(cls.partFn2)
        cls.protImportVol = cls.runImportVolumesCryoSPARC(cls.volFn)

    def testCryosparcSymetryExpansion(self):
        def _runCryosparctestSymetryExpansion(label=''):

            protSymExp = self.newProtocol(ProtCryoSparcSymmetryExpansion)

            prot3DRefinement = self.newProtocol(ProtCryoSparcRefine3D)
            prot3DRefinement.inputParticles.set(self.protImportPart.outputParticles)
            prot3DRefinement.referenceVolume.set(self.protImportVol.outputVolume)
            prot3DRefinement.symmetryGroup.set(SYM_CYCLIC)
            prot3DRefinement.symmetryOrder.set(1)
            prot3DRefinement.compute_use_ssd.set(False)
            self.launchProtocol(prot3DRefinement)

            protSymExp.inputParticles.set(prot3DRefinement.outputParticles)
            protSymExp.symmetryGroup.set(SYM_CYCLIC)
            protSymExp.symmetryOrder.set(4)
            protSymExp.compute_use_ssd.set(False)
            self.launchProtocol(protSymExp)

            return protSymExp

        def _checkAsserts(cryosparcProt):
            self.assertIsNotNone(cryosparcProt.outputParticles,
                                 "There was a problem with Cryosparc Symmetry Expansion")
            self.assertEqual(cryosparcProt.outputParticles.getSize(),
                             self.protImportPart.outputParticles.getSize()*4)

        cryosparcProtGpu = _runCryosparctestSymetryExpansion(label="Cryosparc Symmetry Expansion")
        _checkAsserts(cryosparcProtGpu)


class TestCryosparcNaiveLocalRefine(TestCryosparcBase):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        setupTestProject(cls)
        dataProject = 'grigorieff'
        dataset = DataSet.getDataSet(dataProject)
        TestCryosparcBase.setData()
        particlesPattern = dataset.getFile('particles.sqlite')
        cls.protImportPart = cls.runImportParticleCryoSPARC(cls.partFn2)
        cls.protImportVol = cls.runImportVolumesCryoSPARC(cls.volFn)

    def testCryosparcNaiveLocalRefine(self):
        def _runCryosparctestNaiveLocalRefine(label=''):

            protLocalRefine = self.newProtocol(ProtCryoSparcNaiveLocalRefine)

            prot3DRefinement = self.newProtocol(ProtCryoSparcRefine3D)
            prot3DRefinement.inputParticles.set(self.protImportPart.outputParticles)
            prot3DRefinement.referenceVolume.set(self.protImportVol.outputVolume)
            prot3DRefinement.symmetryGroup.set(SYM_CYCLIC)
            prot3DRefinement.symmetryOrder.set(1)
            prot3DRefinement.compute_use_ssd.set(False)
            self.launchProtocol(prot3DRefinement)

            # Create a 3D Mask using xmipp
            protXmippCreate3DMask = self.runCreate3DMask(prot3DRefinement.outputVolume)

            protLocalRefine.inputParticles.set(prot3DRefinement.outputParticles)
            protLocalRefine.refVolume.set(prot3DRefinement.outputVolume)
            protLocalRefine.refMask.set(protXmippCreate3DMask.outputMask)
            protLocalRefine.compute_use_ssd.set(False)
            self.launchProtocol(protLocalRefine)

            return protLocalRefine

        def _checkAsserts(cryosparcProt):
            self.assertIsNotNone(cryosparcProt.outputParticles,
                                 "There was a problem with Cryosparc subtract projection")

        cryosparcProtGpu = _runCryosparctestNaiveLocalRefine(label="Cryosparc Local Refine")
        _checkAsserts(cryosparcProtGpu)


class TestCryosparcLocalRefine(TestCryosparcBase):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        setupTestProject(cls)
        dataProject = 'grigorieff'
        dataset = DataSet.getDataSet(dataProject)
        TestCryosparcBase.setData()
        particlesPattern = dataset.getFile('particles.sqlite')
        cls.protImportPart = cls.runImportParticleCryoSPARC(cls.partFn2)
        cls.protImportVol = cls.runImportVolumesCryoSPARC(cls.volFn)

    def testCryosparcNaiveLocalRefine(self):
        def _runCryosparctestNaiveLocalRefine(label=''):

            protLocalRefine = self.newProtocol(ProtCryoSparcLocalRefine)

            prot3DRefinement = self.newProtocol(ProtCryoSparcRefine3D)
            prot3DRefinement.inputParticles.set(self.protImportPart.outputParticles)
            prot3DRefinement.referenceVolume.set(self.protImportVol.outputVolume)
            prot3DRefinement.symmetryGroup.set(SYM_CYCLIC)
            prot3DRefinement.symmetryOrder.set(1)
            prot3DRefinement.compute_use_ssd.set(False)
            self.launchProtocol(prot3DRefinement)

            # Create a 3D Mask using xmipp
            protXmippCreate3DMask = self.runCreate3DMask(prot3DRefinement.outputVolume)

            protLocalRefine.inputParticles.set(prot3DRefinement.outputParticles)
            protLocalRefine.refVolume.set(prot3DRefinement.outputVolume)
            protLocalRefine.refMask.set(protXmippCreate3DMask.outputMask)
            protLocalRefine.compute_use_ssd.set(False)
            self.launchProtocol(protLocalRefine)

            return protLocalRefine

        def _checkAsserts(cryosparcProt):
            self.assertIsNotNone(cryosparcProt.outputParticles,
                                 "There was a problem with Cryosparc subtract projection")

        cryosparcProtGpu = _runCryosparctestNaiveLocalRefine(label="Cryosparc Local Refine")
        _checkAsserts(cryosparcProtGpu)


class TestCryosparcHomogeneousReconstruction(TestCryosparcBase):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        setupTestProject(cls)
        dataProject = 'grigorieff'
        dataset = DataSet.getDataSet(dataProject)
        TestCryosparcBase.setData()
        cls.protImportPart = cls.runImportParticleCryoSPARC(cls.partFn2)
        cls.protImportVol = cls.runImportVolumesCryoSPARC(cls.volFn)
        cls.protXmippCreate3DMask = cls.runCreate3DMask(cls.protImportVol.outputVolume)

    def testCryosparcHomogeneousReconstruction(self):
        def _runCryosparctestHomogeneousReconstruction(label=''):

            prot3DRefinement = self.newProtocol(ProtCryoSparcRefine3D)
            prot3DRefinement.inputParticles.set(self.protImportPart.outputParticles)
            prot3DRefinement.referenceVolume.set(self.protImportVol.outputVolume)
            prot3DRefinement.symmetryGroup.set(SYM_CYCLIC)
            prot3DRefinement.symmetryOrder.set(1)
            prot3DRefinement.compute_use_ssd.set(False)
            self.launchProtocol(prot3DRefinement)

            protHomogeneousReconst = self.newProtocol(ProtCryoSparcHomogeneousReconstruct)
            protHomogeneousReconst.inputParticles.set(prot3DRefinement.outputParticles)
            protHomogeneousReconst.refMask.set(self.protXmippCreate3DMask.outputMask)
            self.launchProtocol(protHomogeneousReconst)

            return protHomogeneousReconst

        def _checkAsserts(cryosparcProt):
            self.assertIsNotNone(cryosparcProt.outputParticles,
                                 "There was a problem with Cryosparc subtract projection")

        cryosparcProtGpu = _runCryosparctestHomogeneousReconstruction(label="Cryosparc Homogeneous Reconstruction")
        _checkAsserts(cryosparcProtGpu)





