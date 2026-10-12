"""Regression: reference volumes picked from one SetOfVolumes stay distinct.

A Pointer to a Set item resolves through set[id], which native Scipion
Sets serve from one reused object. Volumes kept beyond their own lookup
must be cloned, or every reference ends up as the last one resolved.
"""
import os
import shutil
import tempfile
import unittest

from pwem import objects as emobj
from pyworkflow.object import Pointer, PointerList

from cryosparc2.protocols.protocol_cryosparc_3D_classification import ProtCryoSparc3DClassification
from cryosparc2.protocols.protocol_cryosparc_new_3D_classification import ProtCryoSparcNew3DClassification


class _Protocol:
    def __init__(self, refVolumes):
        self.refVolumes = refVolumes

    def hasAttribute(self, name):
        return hasattr(self, name)


class TestReferenceVolumesFromOneSet(unittest.TestCase):
    def setUp(self):
        self.workDir = tempfile.mkdtemp()
        fileName = os.path.join(self.workDir, 'volumes.sqlite')
        volumes = emobj.SetOfVolumes(filename=fileName)
        volumes.setSamplingRate(1.0)
        for itemId in (1, 2, 3):
            volume = emobj.Volume()
            volume.setObjId(itemId)
            volume.setFileName('class%03d.mrc' % itemId)
            volume.setSamplingRate(1.0)
            volumes.append(volume)
        volumes.write()
        volumes.close()
        self.volumes = emobj.SetOfVolumes(filename=fileName)

    def tearDown(self):
        self.volumes.close()
        shutil.rmtree(self.workDir)

    def _references(self):
        references = PointerList()
        for itemId in (1, 2, 3):
            references.append(Pointer(self.volumes, extended=str(itemId)))
        return references

    def test_EveryClassificationKeepsEachSelectedReference(self):
        for protocolClass in (ProtCryoSparc3DClassification, ProtCryoSparcNew3DClassification):
            protocol = _Protocol(self._references())
            volumes = protocolClass._getInputVolume(protocol)
            self.assertEqual(['class001.mrc', 'class002.mrc', 'class003.mrc'],
                             [volume.getFileName() for volume in volumes], protocolClass.__name__)
            self.assertEqual([1, 2, 3], [volume.getObjId() for volume in volumes])


if __name__ == '__main__':
    unittest.main()
