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
import json
import os
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
                                  volumeOperation=2,  # OPERATION_THRESHOLD,
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


# class TestCryosparc3DRefinement(TestCryosparcBase):
#
#     @classmethod
#     def setUpClass(cls):
#         setupTestProject(cls)
#         setupTestProject(cls)
#         dataProject = 'grigorieff'
#         dataset = DataSet.getDataSet(dataProject)
#         TestCryosparcBase.setData()
#         particlesPattern = dataset.getFile('particles.sqlite')
#         cls.protImportPart = cls.runImportParticleCryoSPARC(cls.partFn2)
#         cls.protImportVol = cls.runImportVolumesCryoSPARC(cls.volFn)
#
#     def testCryosparc3DRefinement(self):
#         def _runCryosparctest3DRefinement(label=''):
#             prot3DRefinement = self.newProtocol(ProtCryoSparcRefine3D)
#
#             prot3DRefinement.inputParticles.set(self.protImportPart.outputParticles)
#             prot3DRefinement.referenceVolume.set(self.protImportVol.outputVolume)
#             prot3DRefinement.symmetryGroup.set(SYM_CYCLIC)
#             prot3DRefinement.symmetryOrder.set(1)
#             prot3DRefinement.compute_use_ssd.set(False)
#             prot3DRefinement.setObjLabel(label)
#             self.launchProtocol(prot3DRefinement)
#
#             # Check the outputs, sampling rate,...
#             self.assertIsNotNone(prot3DRefinement.outputVolume)
#             self.assertEqual(prot3DRefinement.outputVolume.getSamplingRate(),
#                              self.protImportPart.outputParticles.getSamplingRate())
#
#             outputParticles = prot3DRefinement.outputParticles
#             self.assertIsNotNone(outputParticles)
#             self.assertEqual(outputParticles.getSamplingRate(),
#                              self.protImportPart.outputParticles.getSamplingRate())
#             self.assertTrue(outputParticles.hasAlignmentProj())
#             self.assertTrue(outputParticles.hasCTF())
#             self.assertEqual(outputParticles.getSize(),
#                              self.protImportPart.outputParticles.getSize())
#
#         _runCryosparctest3DRefinement(label="Cryosparc 3D refinement")


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
            prot3DHomoRefinement.refine_defocus_refine.set(False)
            prot3DHomoRefinement.refine_ctf_global_refine.set(False)
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
            protNewNonUniform3DRefinement.refine_defocus_refine.set(False)
            protNewNonUniform3DRefinement.refine_ctf_global_refine.set(False)
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

        _runCryosparctestHelicalRefine3D(label="Cryosparc Helical 3D refinement")


class TestCryosparc3DHeterogeneousRefinement(TestCryosparcBase):

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

    def testCryosparc3DHeterogeneousRefinement(self):
        def _runCryosparctestHeterogeneousRefinement(label=''):
            # Launch a first refinement protocol
            protNonUniform3DRefinement = self.newProtocol(ProtCryoSparcNewNonUniformRefine3D)
            protNonUniform3DRefinement.inputParticles.set(self.protImportPart.outputParticles)
            protNonUniform3DRefinement.referenceVolume.set(self.protImportVolumeVol.outputVolume)
            protNonUniform3DRefinement.symmetryGroup.set(SYM_CYCLIC)
            protNonUniform3DRefinement.symmetryOrder.set(1)
            protNonUniform3DRefinement.compute_use_ssd.set(False)
            protNonUniform3DRefinement.refine_defocus_refine.set(False)
            protNonUniform3DRefinement.refine_ctf_global_refine.set(False)
            protNonUniform3DRefinement.setObjLabel("protNonUniform3DRefinement_1")
            self.launchProtocol(protNonUniform3DRefinement)
            self.assertIsNotNone(protNonUniform3DRefinement.outputVolume)

            # Launch a second refinement protocol
            protNonUniform3DRefinement1 = self.newProtocol(ProtCryoSparcNewNonUniformRefine3D)

            protNonUniform3DRefinement1.inputParticles.set(self.protImportPart.outputParticles)
            protNonUniform3DRefinement1.referenceVolume.set(self.protImportVolumeVol.outputVolume)
            protNonUniform3DRefinement1.symmetryGroup.set(SYM_CYCLIC)
            protNonUniform3DRefinement1.symmetryOrder.set(1)
            protNonUniform3DRefinement1.compute_use_ssd.set(False)
            protNonUniform3DRefinement1.refine_defocus_refine.set(False)
            protNonUniform3DRefinement1.refine_ctf_global_refine.set(False)
            protNonUniform3DRefinement1.setObjLabel("protNonUniform3DRefinement_2")
            self.launchProtocol(protNonUniform3DRefinement1)
            self.assertIsNotNone(protNonUniform3DRefinement1.outputVolume)

            unionProt = self.newProtocol(ProtUnionSet)
            # Set the input volumes
            unionProt.inputSets.append(protNonUniform3DRefinement.outputParticles)
            unionProt.inputSets.append(protNonUniform3DRefinement1.outputParticles)
            unionProt.setObjLabel("Single particles union")
            self.launchProtocol(unionProt)
            self.assertSetSize(unionProt.outputSet, 746,
                               msg="Union of 2 set of particles does not work.")

            # Launch a 3D Heterogeneous Refinement protocol
            prot3DHeterogeneousRefinement = self.newProtocol(ProtCryoSparc3DClassification)
            prot3DHeterogeneousRefinement.inputParticles.set(unionProt.outputSet)
            volumes = PointerList()
            volumes.append(protNonUniform3DRefinement.outputVolume)
            volumes.append(protNonUniform3DRefinement1.outputVolume)
            prot3DHeterogeneousRefinement.refVolumes.set(volumes)
            prot3DHeterogeneousRefinement.compute_use_ssd.set(False)
            self.launchProtocol(prot3DHeterogeneousRefinement)

            # Check the outputs
            self.assertIsNotNone(prot3DHeterogeneousRefinement.outputClasses)
            self.assertIsNotNone(prot3DHeterogeneousRefinement.outputVolumes)

            # Check the number of classes
            self.assertEqual(prot3DHeterogeneousRefinement.outputClasses.getSize(), 2)

            # Check the number of classes
            self.assertEqual(prot3DHeterogeneousRefinement.outputVolumes.getSize(), 2)

        _runCryosparctestHeterogeneousRefinement(label="Cryosparc 3D Heterogeneous Refinement")


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
            # Create a 3D Mask using xmipp
            protXmippCreate3DMask = self.runCreate3DMask(self.protImportVol.outputVolume)

            # Launch a refinement protocol
            protNonUniform3DRefinement = self.newProtocol(ProtCryoSparcNewNonUniformRefine3D)
            protNonUniform3DRefinement.inputParticles.set(self.protImportPart.outputParticles)
            protNonUniform3DRefinement.referenceVolume.set(self.protImportVol.outputVolume)
            protNonUniform3DRefinement.symmetryGroup.set(SYM_CYCLIC)
            protNonUniform3DRefinement.symmetryOrder.set(1)
            protNonUniform3DRefinement.compute_use_ssd.set(False)
            protNonUniform3DRefinement.refine_defocus_refine.set(False)
            protNonUniform3DRefinement.refine_ctf_global_refine.set(False)
            protNonUniform3DRefinement.setObjLabel("protNonUniform3DRefinement")
            self.launchProtocol(protNonUniform3DRefinement)
            self.assertIsNotNone(protNonUniform3DRefinement.outputVolume)

            protParticlesSubtract = self.newProtocol(ProtCryoSparcSubtract)

            prot3DRefinement = self.newProtocol(ProtCryoSparcLocalRefine)
            prot3DRefinement.inputParticles.set(protNonUniform3DRefinement.outputParticles)
            prot3DRefinement.refVolume.set(protNonUniform3DRefinement.outputVolume)
            prot3DRefinement.refMask.set(protXmippCreate3DMask.outputMask)
            prot3DRefinement.symmetryGroup.set(SYM_CYCLIC)
            prot3DRefinement.compute_use_ssd.set(False)
            prot3DRefinement.symmetryOrder.set(1)
            self.launchProtocol(prot3DRefinement)
            self.assertIsNotNone(prot3DRefinement.outputVolume)

            protParticlesSubtract.inputParticles.set(prot3DRefinement.outputParticles)
            protParticlesSubtract.refVolume.set(prot3DRefinement.outputVolume)
            protParticlesSubtract.refMask.set(protXmippCreate3DMask.outputMask)
            protParticlesSubtract.compute_use_ssd.set(False)
            self.launchProtocol(protParticlesSubtract)

            self.assertIsNotNone(protParticlesSubtract.outputParticles)

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

    def testCryosparcSharpening(self):
        def _runCryosparctestSharppening(label=''):

            protSharppening = self.newProtocol(ProtCryoSparcSharppening)
            prot3DRefinement = self.newProtocol(ProtCryoSparc3DHomogeneousRefine)
            prot3DRefinement.inputParticles.set(self.protImportPart.outputParticles)
            prot3DRefinement.referenceVolume.set(self.protImportVol.outputVolume)
            prot3DRefinement.symmetryGroup.set(SYM_CYCLIC)
            prot3DRefinement.symmetryOrder.set(1)
            prot3DRefinement.compute_use_ssd.set(False)
            prot3DRefinement.refine_defocus_refine.set(False)
            prot3DRefinement.refine_ctf_global_refine.set(False)
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

            prot3DRefinement = self.newProtocol(ProtCryoSparc3DHomogeneousRefine)
            prot3DRefinement.inputParticles.set(self.protImportPart.outputParticles)
            prot3DRefinement.referenceVolume.set(self.protImportVol.outputVolume)
            prot3DRefinement.symmetryGroup.set(SYM_CYCLIC)
            prot3DRefinement.symmetryOrder.set(1)
            prot3DRefinement.compute_use_ssd.set(False)
            prot3DRefinement.refine_defocus_refine.set(False)
            prot3DRefinement.refine_ctf_global_refine.set(False)
            self.launchProtocol(prot3DRefinement)

            # Create a 3D Mask using xmipp
            protXmippCreate3DMask = self.runCreate3DMask(prot3DRefinement.outputVolume)

            protGlobalCtfRefinement.inputParticles.set(prot3DRefinement.outputParticles)
            protGlobalCtfRefinement.refVolume.set(prot3DRefinement.outputVolume)
            protGlobalCtfRefinement.refMask.set(protXmippCreate3DMask.outputMask)
            protGlobalCtfRefinement.crg_min_res_A.set(400)
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

            prot3DRefinement = self.newProtocol(ProtCryoSparc3DHomogeneousRefine)
            prot3DRefinement.inputParticles.set(self.protImportPart.outputParticles)
            prot3DRefinement.referenceVolume.set(self.protImportVol.outputVolume)
            prot3DRefinement.symmetryGroup.set(SYM_CYCLIC)
            prot3DRefinement.symmetryOrder.set(1)
            prot3DRefinement.compute_use_ssd.set(False)
            prot3DRefinement.refine_defocus_refine.set(False)
            prot3DRefinement.refine_ctf_global_refine.set(False)
            self.launchProtocol(prot3DRefinement)

            # Create a 3D Mask using xmipp
            protXmippCreate3DMask = self.runCreate3DMask(prot3DRefinement.outputVolume)

            protLocalCtfRefinement.inputParticles.set(prot3DRefinement.outputParticles)
            protLocalCtfRefinement.refVolume.set(prot3DRefinement.outputVolume)
            protLocalCtfRefinement.refMask.set(protXmippCreate3DMask.outputMask)
            protLocalCtfRefinement.crl_min_res_A.set(400)
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

            prot3DRefinement = self.newProtocol(ProtCryoSparc3DHomogeneousRefine)
            prot3DRefinement.inputParticles.set(self.protImportPart.outputParticles)
            prot3DRefinement.referenceVolume.set(self.protImportVol.outputVolume)
            prot3DRefinement.symmetryGroup.set(SYM_CYCLIC)
            prot3DRefinement.symmetryOrder.set(1)
            prot3DRefinement.compute_use_ssd.set(False)
            prot3DRefinement.refine_defocus_refine.set(False)
            prot3DRefinement.refine_ctf_global_refine.set(False)
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

    def testCryosparcLocalRefine(self):
        def _runCryosparctestLocalRefine(label=''):

            protLocalRefine = self.newProtocol(ProtCryoSparcLocalRefine)

            prot3DRefinement = self.newProtocol(ProtCryoSparc3DHomogeneousRefine)
            prot3DRefinement.inputParticles.set(self.protImportPart.outputParticles)
            prot3DRefinement.referenceVolume.set(self.protImportVol.outputVolume)
            prot3DRefinement.symmetryGroup.set(SYM_CYCLIC)
            prot3DRefinement.symmetryOrder.set(1)
            prot3DRefinement.compute_use_ssd.set(False)
            prot3DRefinement.refine_defocus_refine.set(False)
            prot3DRefinement.refine_ctf_global_refine.set(False)
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
                                 "There was a problem with Cryosparc local refine")
            transform = cryosparcProt.outputParticles.getFirstItem().getTransform()
            shifts = transform.getShifts()
            shiftsCeros = shifts[0] == shifts[1] == shifts[2] == 0
            self.assertFalse(shiftsCeros, "The transform matrix must have shifts")

        cryosparcProtGpu = _runCryosparctestLocalRefine(label="Cryosparc Local Refine")
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

            prot3DRefinement = self.newProtocol(ProtCryoSparc3DHomogeneousRefine)
            prot3DRefinement.inputParticles.set(self.protImportPart.outputParticles)
            prot3DRefinement.referenceVolume.set(self.protImportVol.outputVolume)
            prot3DRefinement.symmetryGroup.set(SYM_CYCLIC)
            prot3DRefinement.symmetryOrder.set(1)
            prot3DRefinement.compute_use_ssd.set(False)
            prot3DRefinement.refine_defocus_refine.set(False)
            prot3DRefinement.refine_ctf_global_refine.set(False)
            self.launchProtocol(prot3DRefinement)

            protHomogeneousReconst = self.newProtocol(ProtCryoSparcHomogeneousReconstruct)
            protHomogeneousReconst.inputParticles.set(prot3DRefinement.outputParticles)
            protHomogeneousReconst.refMask.set(self.protXmippCreate3DMask.outputMask)
            self.launchProtocol(protHomogeneousReconst)

            return protHomogeneousReconst

        def _checkAsserts(cryosparcProt):
            self.assertIsNotNone(cryosparcProt.outputParticles,
                                 "There was a problem with Cryosparc Homogeneous Reconstruction")

        cryosparcProtGpu = _runCryosparctestHomogeneousReconstruction(label="Cryosparc Homogeneous Reconstruction")
        _checkAsserts(cryosparcProtGpu)


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

    def testTestCryosparc3DClassification(self):
        def _runCryosparctestTestCryosparc3DClassification(label=''):
            # Launch a first refinement protocol
            protNonUniform3DRefinement = self.newProtocol(ProtCryoSparcNewNonUniformRefine3D)
            protNonUniform3DRefinement.inputParticles.set(self.protImportPart.outputParticles)
            protNonUniform3DRefinement.referenceVolume.set(self.protImportVolumeVol.outputVolume)
            protNonUniform3DRefinement.symmetryGroup.set(SYM_CYCLIC)
            protNonUniform3DRefinement.symmetryOrder.set(1)
            protNonUniform3DRefinement.compute_use_ssd.set(False)
            protNonUniform3DRefinement.refine_defocus_refine.set(False)
            protNonUniform3DRefinement.refine_ctf_global_refine.set(False)
            protNonUniform3DRefinement.setObjLabel("protNonUniform3DRefinement_1")
            self.launchProtocol(protNonUniform3DRefinement)
            self.assertIsNotNone(protNonUniform3DRefinement.outputVolume)

            # Launch a second refinement protocol
            protNonUniform3DRefinement1 = self.newProtocol(ProtCryoSparcNewNonUniformRefine3D)

            protNonUniform3DRefinement1.inputParticles.set(self.protImportPart.outputParticles)
            protNonUniform3DRefinement1.referenceVolume.set(self.protImportVolumeVol.outputVolume)
            protNonUniform3DRefinement1.symmetryGroup.set(SYM_CYCLIC)
            protNonUniform3DRefinement1.symmetryOrder.set(1)
            protNonUniform3DRefinement1.compute_use_ssd.set(False)
            protNonUniform3DRefinement1.refine_defocus_refine.set(False)
            protNonUniform3DRefinement1.refine_ctf_global_refine.set(False)
            protNonUniform3DRefinement1.setObjLabel("protNonUniform3DRefinement_2")
            self.launchProtocol(protNonUniform3DRefinement1)
            self.assertIsNotNone(protNonUniform3DRefinement1.outputVolume)

            # Launch a second refinement protocol
            protNonUniform3DRefinement2 = self.newProtocol(ProtCryoSparcNewNonUniformRefine3D)

            protNonUniform3DRefinement2.inputParticles.set(self.protImportPart.outputParticles)
            protNonUniform3DRefinement2.referenceVolume.set(self.protImportVolumeVol.outputVolume)
            protNonUniform3DRefinement2.symmetryGroup.set(SYM_CYCLIC)
            protNonUniform3DRefinement2.symmetryOrder.set(1)
            protNonUniform3DRefinement2.compute_use_ssd.set(False)
            protNonUniform3DRefinement2.refine_defocus_refine.set(False)
            protNonUniform3DRefinement2.refine_ctf_global_refine.set(False)
            protNonUniform3DRefinement2.setObjLabel("protNonUniform3DRefinement_3")
            self.launchProtocol(protNonUniform3DRefinement2)
            self.assertIsNotNone(protNonUniform3DRefinement2.outputVolume)

            # Launch a join set
            unionProt = self.newProtocol(ProtUnionSet)
            # Set the input volumes
            unionProt.inputSets.append(protNonUniform3DRefinement.outputParticles)
            unionProt.inputSets.append(protNonUniform3DRefinement1.outputParticles)
            unionProt.inputSets.append(protNonUniform3DRefinement2.outputParticles)
            unionProt.setObjLabel("Single particles union")
            self.launchProtocol(unionProt)
            self.assertSetSize(unionProt.outputSet, 1119, msg="Union of 3 set of particles does not work.")

            # Launch a 3D Classification protocol
            prot3DClassification = self.newProtocol(ProtCryoSparcNew3DClassification)
            prot3DClassification.inputParticles.set(unionProt.outputSet)
            prot3DClassification.class3D_init_mode.set(0)
            prot3DClassification.class3D_N_ik.set(400)
            prot3DClassification.compute_use_ssd.set(False)
            self.launchProtocol(prot3DClassification)

            # Check the outputs
            self.assertIsNotNone(prot3DClassification.outputClasses)
            self.assertIsNotNone(prot3DClassification.outputVolumes)

            # Check the number of classes
            self.assertEqual(prot3DClassification.outputClasses.getSize(), 2)

            # Check the number of classes
            self.assertEqual(prot3DClassification.outputVolumes.getSize(), 2)

        _runCryosparctestTestCryosparc3DClassification(label="Cryosparc 3D Classification")


class TestCryosparcJobsOfAStep(BaseTest):
    """How a protocol keeps track of its cryoSPARC jobs across executions (Continue, Restart, Abort), with the
    cryoSPARC tools simulated: no server is needed."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    def _newProt(self):
        import pyworkflow.object as pwobj
        prot = self.newProtocol(ProtCryo2D)
        prot.projectName = pwobj.String('P1')
        prot.workSpaceName = pwobj.String('W1')
        prot.particles = pwobj.String('J1.imported_particles')
        prot._defineParamsName()
        self.proj.saveProtocol(prot)
        return prot

    def _initialize(self, prot):
        """What every launch does before inserting the steps (also on Continue and Restart)."""
        from unittest import mock
        from ..protocols import protocol_base as base

        def utilsVariables():
            prot.projectDirName, prot.projectPath = 'scipion', '/cs/scipion'
            prot.projectContainerDir = self.getOutputPath()

        with mock.patch.object(prot, '_initializeUtilsVariables', side_effect=utilsVariables), \
                mock.patch.object(base, 'getCryosparcProjectsList',
                                  return_value=[{'title': 'scipion', 'project_dir': '/cs/P1', 'uid': 'P1'}]), \
                mock.patch.object(base, 'getProjectPath', return_value=['/cs/P1']), \
                mock.patch.object(base, 'getCryosparcVersion', return_value='4.4.0'), \
                mock.patch.object(base, 'getCryosparcWorkSpaces', return_value=[{'uid': 'W1'}]):
            prot._initializeCryosparcProject()

    def _run2D(self, prot, status, newJob='J10'):
        """Run the 2D classification step with cryoSPARC reporting status for the jobs it already knows."""
        from unittest import mock
        import pyworkflow.object as pwobj
        from ..protocols import protocol_base as base
        from ..protocols import protocol_cryosparc2d as module2D
        enqueue = mock.Mock(return_value=pwobj.String(newJob))
        with mock.patch.object(module2D, 'enqueueJob', enqueue), \
                mock.patch.object(base, 'enqueueJob', enqueue), \
                mock.patch.object(base, 'getJobStatus', return_value=status), \
                mock.patch.object(module2D, 'waitForCryosparc') as wait, \
                mock.patch.object(module2D, 'clearIntermediateResults'), \
                mock.patch.object(module2D, 'isCryosparcStandalone', return_value=True):
            prot.doRunClass2D()
        return enqueue, wait

    def testContinueKeepsTheJobItWasRunning(self):
        prot = self._newProt()
        self._initialize(prot)
        prot.currenJob.set('J5')

        self._initialize(prot)  # Continue

        self.assertEqual(prot.currenJob.get(), 'J5', 'Abort could not kill the job anymore.')

    def testContinueKeepsItsCryosparcProjectAndWorkspace(self):
        from unittest import mock
        from ..protocols import protocol_base as base
        prot = self._newProt()
        self._initialize(prot)
        prot.projectDir.set('/cs/P1')

        # Meanwhile, another cryoSPARC project with the same title and another workspace were created
        with mock.patch.object(prot, '_initializeUtilsVariables'), \
                mock.patch.object(base, 'getCryosparcProjectsList',
                                  return_value=[{'title': 'scipion', 'project_dir': '/cs/P1', 'uid': 'P1'},
                                                {'title': 'scipion', 'project_dir': '/cs/P7', 'uid': 'P7'}]), \
                mock.patch.object(base, 'getProjectPath', return_value=['/cs/P7']), \
                mock.patch.object(base, 'getCryosparcVersion', return_value='4.4.0'), \
                mock.patch.object(base, 'getCryosparcWorkSpaces', return_value=[{'uid': 'W1'}, {'uid': 'W9'}]):
            prot.projectDirName, prot.projectPath = 'scipion', '/cs/scipion'
            prot.projectContainerDir = self.getOutputPath()
            prot._initializeCryosparcProject()  # Continue

        self.assertEqual((prot.projectName.get(), prot.workSpaceName.get(), prot.projectDir.get()),
                         ('P1', 'W1', '/cs/P1'), 'The jobs stored by the run belong to P1.')

    def testAbortingWhileCryosparcIsDownStillAborts(self):
        from unittest import mock
        from ..protocols import protocol_base as base
        prot = self._newProt()
        self._initialize(prot)
        prot.currenJob.set('J5')

        with mock.patch.object(base, 'getJobStatus', side_effect=Exception('Connection refused')):
            prot.setAborted()

        self.assertTrue(prot.isAborted())

    def testAnImportJobCanBeKilledWhileItRuns(self):
        from unittest import mock
        import pyworkflow.object as pwobj
        from .. import utils
        from ..protocols import protocol_base as base
        prot = self._newProt()
        self._initialize(prot)
        stored = []
        with mock.patch.object(utils, 'enqueueJob', return_value=pwobj.String('J7')), \
                mock.patch.object(base, 'enqueueJob', return_value=pwobj.String('J7')), \
                mock.patch.object(utils, 'waitForCryosparc',
                                  side_effect=lambda *args, **kwargs: stored.append(prot.currenJob.get())), \
                mock.patch.object(prot, '_getFileName', return_value='input_particles.star'), \
                mock.patch.object(prot, '_getInputParticles', return_value=mock.Mock(getSamplingRate=lambda: 1.0)):
            utils.doImportParticlesStar(prot)

        self.assertEqual(stored, ['J7'], 'While the import runs, Abort does not know which job to kill.')

    def testContinueTakesBackTheJobOfThePreviousExecution(self):
        prot = self._newProt()
        self._initialize(prot)
        self._run2D(prot, status=None, newJob='J9')  # The execution that crashed while waiting for J9

        self._initialize(prot)  # Continue
        enqueue, wait = self._run2D(prot, status='running')

        enqueue.assert_not_called()
        self.assertEqual(str(wait.call_args[0][1]), 'J9')

    def testACompletedJobOfThePreviousExecutionIsReplaced(self):
        # The step runs again because something changed (e.g. its input particles):
        # the job completed before may hold the results of other inputs.
        prot = self._newProt()
        self._initialize(prot)
        self._run2D(prot, status=None, newJob='J9')

        self._initialize(prot)
        enqueue, wait = self._run2D(prot, status='completed', newJob='J13')

        enqueue.assert_called_once()
        self.assertEqual(str(wait.call_args[0][1]), 'J13')

    def testAFailedJobOfThePreviousExecutionIsReplaced(self):
        prot = self._newProt()
        self._initialize(prot)
        self._run2D(prot, status=None, newJob='J9')

        self._initialize(prot)
        enqueue, wait = self._run2D(prot, status='failed', newJob='J11')

        enqueue.assert_called_once()
        self.assertEqual(str(wait.call_args[0][1]), 'J11')

    def testRestartKillsTheRunningJobsAndDoesNotTakeThemBack(self):
        from unittest import mock
        from pyworkflow.protocol.constants import MODE_RESTART
        from ..protocols import protocol_base as base
        prot = self._newProt()
        self._initialize(prot)
        self._run2D(prot, status=None, newJob='J9')

        prot.runMode.set(MODE_RESTART)
        with mock.patch.object(base, 'getJobStatus', return_value='running'), \
                mock.patch.object(base, 'killJob') as kill, mock.patch.object(base, 'clearJob'):
            self._initialize(prot)
        enqueue, wait = self._run2D(prot, status='completed', newJob='J12')

        kill.assert_called_once_with('P1', 'J9')
        self.assertEqual(str(wait.call_args[0][1]), 'J12')


class TestCryosparcStepDetails(BaseTest):
    """Details of single steps that broke a run, sent wrong values to cryoSPARC or lost data, with cryoSPARC
    simulated: no server is needed."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    def _newProt(self, protClass, **kwargs):
        prot = self.newProtocol(protClass, **kwargs)
        self.proj.saveProtocol(prot)
        prot.makePathsAndClean()
        return prot

    def test2DReadsTheLastIteration(self):
        from ..protocols import ProtCryo2D
        for online, final, suffix in [(5, 3, '_007'), (2, 1, '_002'), (20, 2, '_021'), (40, 60, '_099'),
                                      (90, 20, '_109')]:
            prot = self.newProtocol(ProtCryo2D, numberOnlineEMIterator=online, numberFinalIterator=final)
            self.assertEqual(prot._getNumberOfIterSuffix(), suffix, 'online=%s, final=%s' % (online, final))

    def testHeterogeneousRefinementContinuesFromItsJob(self):
        from unittest import mock
        import pyworkflow.object as pwobj
        from ..protocols import protocol_base as base
        prot = self._newProt(ProtCryoSparc3DClassification)
        prot.importVolumes = pwobj.CsvList()
        prot.importVolumes.append('J3')

        # Continue: the conversion step ran in a previous execution
        with mock.patch.object(base, 'getCryosparcVersion', return_value='4.4.0'), \
                mock.patch.object(prot, 'do3DClasification'):
            prot.processStep()

        self.assertEqual(prot.volumes, ['J3.imported_volume_1.map'])

    def testLocalCtfRefinementSendsOnlyTheParametersItHas(self):
        from unittest import mock
        import pyworkflow.object as pwobj
        from ..protocols import protocol_cryosparc_local_ctf_refinement as module
        prot = self._newProt(ProtCryoSparcLocalCtfRefinement, crl_max_res_A=3.5)
        prot.crl_N.set(None)
        prot._defineParamsName()
        prot.particles, prot.volume, prot.mask = (pwobj.String('J1.particles'), pwobj.String('J2.volume'),
                                                  pwobj.String('J3.mask'))
        prot.projectName, prot.workSpaceName = pwobj.String('P1'), pwobj.String('W1')
        prot.currenJob = pwobj.String()
        enqueue = mock.Mock(return_value=pwobj.String('J4'))
        with mock.patch.object(prot, '_enqueueOrResumeJob', enqueue), \
                mock.patch.object(prot, '_getInputVolume', return_value=mock.Mock(hasHalfMaps=lambda: False)), \
                mock.patch.object(module, 'waitForCryosparc'), mock.patch.object(prot, '_store'):
            prot.doLocalCtfRefinement()

        params = json.loads(enqueue.call_args[0][4])
        self.assertNotIn('None', params.values(), 'cryoSPARC receives the text None: %s' % params)
        self.assertEqual(params.get('crl_max_res_A'), '3.5')

    def testVariabilityDisplayAcceptsAnEmptyIntermediateWidth(self):
        from unittest import mock
        import pyworkflow.object as pwobj
        from ..protocols import protocol_cryosparc_3d_variability_display as module
        prot = self._newProt(ProtCryoSparc3DVariabilityDisplay)
        prot.var_intermediate_width.set(None)
        prot._defineParamsName()
        prot.projectName, prot.workSpaceName = pwobj.String('P1'), pwobj.String('W1')
        prot.currenJob = pwobj.String()
        enqueue = mock.Mock(return_value=pwobj.String('J4'))
        with mock.patch.object(prot, '_enqueueOrResumeJob', enqueue), \
                mock.patch.object(prot, 'input3DVariablityAnalisysProt',
                                  mock.Mock(get=lambda: mock.Mock(run3DVariability='J3'))), \
                mock.patch.object(module, 'waitForCryosparc'), mock.patch.object(prot, '_store'):
            prot.doRun3DVariabilityDisplay()

        self.assertEqual(json.loads(enqueue.call_args[0][4])['var_intermediate_width'], '0')

    def testTwoVariabilityDisplaysDoNotShareTheirParticleClasses(self):
        first = self.newProtocol(ProtCryoSparc3DVariabilityDisplay)
        second = self.newProtocol(ProtCryoSparc3DVariabilityDisplay)
        first._initialize()
        second._initialize()

        self.assertIsNot(first.partClassDict, second.partClassDict)

    def testPatchCtfGivesEachMicrographItsOwnCtf(self):
        from pwem.objects import SetOfMicrographs, Micrograph
        micsFn = self.getOutputPath('mics.sqlite')
        if os.path.exists(micsFn):
            os.remove(micsFn)
        mics = SetOfMicrographs(filename=micsFn)
        mics.setSamplingRate(1.0)
        for name in ('mic_1', 'mic_2'):
            mic = Micrograph(location=self.getOutputPath(name + '.mrc'))
            mic.setMicName(name)
            mics.append(mic)
        mics.write()
        prot = self._newProt(ProtCryoSparcPatchCTFEstimate)
        prot.inputMicrographs.set(mics)
        starFn = prot._getExtraPath('ctf.star')
        with open(starFn, 'w') as f:
            f.write('data_micrograph\n\nloop_\n_rlnDefocusU #1\n_rlnDefocusV #2\n_rlnPhaseShift #3\n'
                    '_rlnCtfMaxResolution #4\n_rlnDefocusAngle #5\n'
                    '10000 10100 0 4 10\n20000 20100 0 4 20\n')
        micList = {mic.getMicName(): mic.clone() for mic in mics}

        ctfs = prot._fillSetOfCTF(starFn, micList)

        self.assertEqual(sorted(ctf.getDefocusU() for ctf in ctfs), [10000, 20000])

    def testFlexReconstructionWithoutRigidReconstructionHasItsFlexVolume(self):
        from unittest import mock
        import mrcfile
        import numpy as np
        import pyworkflow.object as pwobj
        from ..protocols import protocol_cryosparc_3D_flex_reconstruction as module
        prot = self._newProt(ProtCryoSparc3DFlexReconstruction, flex_do_noflex_recon=False)
        jobDir = self.getOutputPath('P1', 'J8')
        os.makedirs(jobDir, exist_ok=True)
        for name in ('J8_flex_map.mrc', 'J8_flex_map_half_A.mrc', 'J8_flex_map_half_B.mrc'):
            with mrcfile.new(os.path.join(jobDir, name), overwrite=True) as mrc:
                mrc.set_data(np.zeros((8, 8, 8), dtype=np.float32))
                mrc.voxel_size = 1.5
        prot.projectDir, prot.projectName = pwobj.String(os.path.dirname(jobDir)), pwobj.String('P1')
        prot.run3DFlexReconstructionJob = pwobj.String('J8')
        with mock.patch.object(prot, '_initializeUtilsVariables'), \
                mock.patch.object(module, 'getOutputPreffix', return_value=''):
            prot.createOutputStep()

        self.assertTrue(hasattr(prot, 'flexVolume'))
        self.assertFalse(hasattr(prot, 'noFlexVolume'), 'A rigid reconstruction that was not requested.')

    def testPositiveParametersAreValidated(self):
        from pyworkflow.protocol.params import Positive
        for protClass, params in [(ProtCryoSparc3DVariability, ['var_filter_res']),
                                  (ProtCryoSparcNewNonUniformRefine3D, ['refine_nu_order', 'refine_nu_awf'])]:
            prot = self.newProtocol(protClass)
            for name in params:
                self.assertIn(Positive, prot.getParam(name).validators, '%s.%s' % (protClass.__name__, name))

    def _reconstructionFscId(self, logEvents):
        from unittest import mock
        import pyworkflow.object as pwobj
        from ..protocols import protocol_cryosparc_homogeneous_reconstruction as module
        prot = self._newProt(ProtCryoSparcHomogeneousReconstruct)
        prot._defineFileNames()
        prot.projectName = pwobj.String('P1')

        def streamLog(project, job, outputFile):
            with open(outputFile, 'w') as f:
                f.write(repr(logEvents))

        with mock.patch.object(module, 'get_job_streamlog', side_effect=streamLog):
            return prot.findLastIteration('J5')

    def testReconstructionFscLookupToleratesALogWithoutIt(self):
        self.assertIsNone(self._reconstructionFscId([{'text': 'Done.'}]))

    def testReconstructionFscLookupTakesTheTxtFile(self):
        event = {'text': 'FSC, after mask auto-tightening',
                 'imgfiles': [{'fileid': 'png1', 'filetype': 'png'}, {'fileid': 'txt1', 'filetype': 'txt'}]}

        self.assertEqual(self._reconstructionFscId([event]), 'txt1')

    def _validated(self, protClass, particles, **kwargs):
        from unittest import mock
        from ..protocols import protocol_cryosparc_homogeneous_reconstruction as reconstruction
        from ..protocols import protocol_cryosparc_3D_flex_data_prepare as flexPrepare
        prot = self.newProtocol(protClass, **kwargs)
        with mock.patch.object(reconstruction, 'cryosparcValidate', return_value=[]), \
                mock.patch.object(reconstruction, 'gpusValidate', return_value=[]), \
                mock.patch.object(flexPrepare, 'cryosparcValidate', return_value=[]), \
                mock.patch.object(prot, '_getInputParticles', return_value=particles):
            return prot._validate()

    def testReconstructionNeedsAlignedParticles(self):
        from unittest import mock
        unaligned = mock.Mock(hasCTF=lambda: True, hasAlignmentProj=lambda: False)
        aligned = mock.Mock(hasCTF=lambda: True, hasAlignmentProj=lambda: True)

        self.assertTrue(any('align' in msg for msg in self._validated(ProtCryoSparcHomogeneousReconstruct, unaligned)))
        self.assertEqual(self._validated(ProtCryoSparcHomogeneousReconstruct, aligned), [])

    def testFlexDataPreparationChecksItsInputWithoutFailing(self):
        from unittest import mock
        unaligned = mock.Mock(hasCTF=lambda: True, hasAlignmentProj=lambda: False)
        aligned = mock.Mock(hasCTF=lambda: True, hasAlignmentProj=lambda: True)

        self.assertTrue(any('align' in msg for msg in
                            self._validated(ProtCryoSparc3DFlexDataPrepare, unaligned, keep_num_particles=None)))
        self.assertEqual(self._validated(ProtCryoSparc3DFlexDataPrepare, aligned, keep_num_particles=None), [])
        self.assertTrue(any('1000' in msg for msg in
                            self._validated(ProtCryoSparc3DFlexDataPrepare, aligned, keep_num_particles=1500)))

    def testFlexGeneratorWithATrainingThatNeverRanSaysSo(self):
        from unittest import mock
        from ..protocols import protocol_cryosparc_3D_flex_generator as module
        training = self.newProtocol(ProtCryoSparc3DFlexTraining)
        prot = self.newProtocol(ProtCryoSparc3DFlexGenerator)
        prot.flexTraining.set(training)

        with mock.patch.object(module, 'cryosparcValidate', return_value=[]):
            messages = prot._validate()

        self.assertTrue(messages, 'It would run on a training without a model.')

    def testSharpeningSummaryShowsItsInput(self):
        import pyworkflow.object as pwobj
        prot = self.newProtocol(ProtCryoSparcSharppening)
        prot.outputVolume = pwobj.String('volume.mrc')

        self.assertTrue(any('Input' in line for line in prot._summary()))

    def testDownloadedFscFilesStayInTheRunFolder(self):
        from unittest import mock
        from ..protocols import protocol_base as base
        prot = self._newProt(ProtCryoSparc3DHomogeneousRefine)
        response = mock.Mock(content=b'0.0 1.0\n')
        with mock.patch.object(base, 'getSystemInfo', return_value=(0, repr({'master_hostname': 'h', 'port_webapp': 1}))), \
                mock.patch.object(base, 'getCryosparcVersion', return_value='4.0.0'), \
                mock.patch.object(base.requests, 'get', return_value=response), \
                mock.patch.object(prot, '_looksLikeFscTxt', return_value=True):
            fscFile = prot._downloadLegacyFscFileById('f1')

        self.assertTrue(os.path.abspath(fscFile).startswith(os.path.abspath(prot._getTmpPath())), fscFile)
