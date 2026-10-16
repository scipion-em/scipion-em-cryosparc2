import getpass
import os
import tempfile
import unittest
from unittest.mock import patch, Mock

from cryosparc2 import V_UNKNOWN, V3_0_0
from cryosparc2.utils import (cryosparcValidate, cryosparcExists,
                              isCryosparcRunning, calculateNewSamplingRate,
                              getProjectName, getCryosparcVersion,
                              getCryosparcPreprocessLane)

import cryosparc2.utils as csutils

class TestUtils(unittest.TestCase):

    def testProjectName(self):

        scipionProjectName = "TestProject"

        projectName = getProjectName(scipionProjectName)
        user = getpass.getuser()
        expected = scipionProjectName + "-" + user
        self.assertEqual(projectName, expected, "Project name does not equals %s." % expected)

    def testExists(self):

        with patch('cryosparc2.utils.getCryosparcDir') as csDir:
            csDir.return_value = None
            self.assertFalse(cryosparcExists(), "cryosparcExists fail to detect a wrong installation")

            csDir.return_value = "."
            self.assertTrue(cryosparcExists(), "cryosparcExists fail to detect installation")

    def testIsRunning(self):

        with patch("cryosparc2.utils.getCryosparcProgram") as getProg:

            getProg.return_value = None
            self.assertFalse(isCryosparcRunning(), "Cryosparc status RUNNING when program is None")

            getProg.return_value = "something"

            with patch("subprocess.getstatusoutput") as cmdoutput:
                cmdoutput.return_value = (0, "Ok")
                self.assertTrue(isCryosparcRunning(), "isCryosparcRunning running but not detected")

                cmdoutput.return_value = (1, "Ok")
                self.assertFalse(isCryosparcRunning(), "isCryosparcRunning not running but not detected")

    def testValidate(self):

        with patch('cryosparc2.utils.cryosparcExists') as exists:
            # Case, CS not found
            exists.return_value = False
            result = cryosparcValidate()
            self.assertTrue("cryoSPARC software not found" in result[0], "Validation did not detect CS not found")

            exists.return_value = True

            # Running?
            with patch('cryosparc2.utils.isCryosparcRunning') as running:
                running.return_value = False
                result = cryosparcValidate()
                self.assertTrue("connect" in result[0], "Validation did not detect CS not running")

                running.return_value = True

                # Installed version:
                with patch('cryosparc2.utils.getCryosparcVersion') as getVersion:

                    # Low version
                    getVersion.return_value = "1.0.0"
                    result = cryosparcValidate()
                    self.assertTrue('not compatible' in result[0], "Validation did not detect CS low version")

                    # Supported version
                    getVersion.return_value = V3_0_0
                    result = cryosparcValidate()
                    self.assertEqual(0, len(result), "Validation did not detectCS correct version.")

                    # Patch printing Warning function--> redStr
                    with patch('pyworkflow.utils.yellowStr') as yellowStr:
                        # Higher version
                        highV = "100.1.1"
                        getVersion.return_value = highV
                        result = cryosparcValidate()
                        self.assertEqual(0, len(result), "Validation did not allows higher versions than supported.")
                        self.assertTrue(yellowStr.called, "Warning not printed when working with higher versions")

    def testSamplingRateConvertion(self):

        sr = calculateNewSamplingRate((2, 2, 2), 4, (4, 4, 4))
        self.assertEqual(sr, 8, "Wrong sampling rate conversion 1")

        sr = calculateNewSamplingRate((2, 2, 2), 1.5, (4, 4, 4))
        self.assertEqual(sr, 3, "Wrong sampling rate conversion 2")

        sr = calculateNewSamplingRate((3, 3, 3), 1.5, (4, 4, 4))
        self.assertEqual(sr, 2, "Wrong sampling rate conversion 3")

    def testGetVersion(self):

        VERSION_FROM_FILE = "1.2.3"
        VERSION_FROM_ENV = "3.2.1"

        # Installed version:
        with patch('cryosparc2.utils.getCryosparcEnvInformation') as getEnvInfo:
            with patch('cryosparc2.utils._getCryosparcVersionFromFile') as getFromFile:

                # Set up the mocks
                getEnvInfo.return_value = VERSION_FROM_ENV
                getFromFile.return_value = VERSION_FROM_FILE

                # Call getting version from file
                version = getCryosparcVersion()
                self.assertEqual(version, VERSION_FROM_FILE)
                version = getCryosparcVersion()
                self.assertEqual(version, VERSION_FROM_FILE)
                getFromFile.assert_called_once()

                # Reset
                csutils._csVersion = None
                getFromFile.reset_mock()
                getFromFile.side_effect = Exception("version file missing")

                # Call getting version from env
                getFromFile.assert_not_called()
                version = getCryosparcVersion()
                self.assertEqual(version, VERSION_FROM_ENV)
                getFromFile.assert_called_once()
                getEnvInfo.assert_called_once()

                version = getCryosparcVersion()
                self.assertEqual(version, VERSION_FROM_ENV)
                getFromFile.assert_called_once()
                getEnvInfo.assert_called_once()

                # Reset
                csutils._csVersion = None
                getFromFile.reset_mock()
                getFromFile.side_effect = Exception('version file missing')
                getEnvInfo.reset_mock()
                getEnvInfo.side_effect = Exception('Env fails')

                # Call getting version from env
                version = getCryosparcVersion()
                self.assertEqual(version, V_UNKNOWN)
                getFromFile.assert_called_once()
                getEnvInfo.assert_called_once()

    def testGetPreprocessLane(self):

        with patch.dict('os.environ', {}, clear=True):
            self.assertIsNone(getCryosparcPreprocessLane())

        with patch.dict('os.environ', {'CRYOSPARC_DEFAULT_LANE': 'heavy-lane'}, clear=True):
            self.assertEqual(getCryosparcPreprocessLane(), 'heavy-lane')

        with patch.dict('os.environ', {'CRYOSPARC_DEFAULT_LANE': 'heavy-lane',
                                       'CRYOSPARC_PREPROCESS_LANE': 'light-lane'}, clear=True):
            self.assertEqual(getCryosparcPreprocessLane(), 'light-lane')



class SleepGuard:
    """Stands for time.sleep: fails the test when a wait loops for ever."""

    def __init__(self, maxCalls=200):
        self.calls = 0
        self.maxCalls = maxCalls

    def __call__(self, seconds):
        self.calls += 1
        if self.calls > self.maxCalls:
            raise AssertionError('The wait never gives up')


class TestJobLifecycle(unittest.TestCase):
    """How the plugin waits for cryoSPARC jobs and reads what the cryoSPARC tools return."""

    def _wait(self, status):
        with patch.object(csutils, 'getJobStatus', side_effect=status), \
                patch.object(csutils, '_logCryosparcJobEvents'), \
                patch.object(csutils, 'waitJob'), \
                patch.object(csutils.time, 'sleep', SleepGuard()):
            return csutils.waitForCryosparc('P1', 'J5', 'The job failed.')

    def testTheWaitGivesUpWhenCryosparcKeepsFailing(self):
        with self.assertRaisesRegex(Exception, 'J5'):
            self._wait(Exception('Connection refused'))

    def testTheWaitFailsOnAStatusItDoesNotKnow(self):
        with self.assertRaisesRegex(Exception, 'mystery'):
            self._wait(lambda *args: 'mystery')

    def testTheWaitUnderstandsEnumLikeStatuses(self):
        self.assertEqual(self._wait(lambda *args: csutils._normalizeJobStatus('JobStatus.completed')), 'completed')

    def testOnlyTheLastLineOfTheCliOutputIsTheValue(self):
        with patch.object(csutils, '_runCliExpression', return_value=(0, "WARNING: a banner line\n'J5'")):
            self.assertEqual(csutils._runCliValue("api.jobs.create(...).uid"), 'J5')

    def testAnUnknownVersionGivesAClearError(self):
        with patch.object(csutils, '_isCryosparcV5OrNewer', return_value=False), \
                patch.object(csutils, 'getCryosparcVersion', return_value=V_UNKNOWN), \
                patch.object(csutils, 'isCryosparcStandalone', return_value=True), \
                patch.object(csutils, 'runCmd') as runCmd:
            with self.assertRaisesRegex(Exception, 'version'):
                csutils.enqueueJob('class_2D', 'P1', 'W1', '{}', '{}', 'default')
        runCmd.assert_not_called()


class TestParticleKeys(unittest.TestCase):
    """The key a particle is found by, from Scipion and from what cryoSPARC lists."""

    def testAParticleIsFoundWhateverExtensionItsStackWasConvertedTo(self):
        item = Mock(getIndex=lambda: 3, getFileName=lambda: '/data/Runs/000123/extra/parts.stk')

        self.assertEqual(csutils.particleKey(item), csutils.csParticleKey('000003@J1/imported/012345678901234567890_parts.mrcs'))

    def testAStackNameWithDotsKeepsAllOfItsName(self):
        item = Mock(getIndex=lambda: 1, getFileName=lambda: '/data/mic_001.aligned.mrcs')

        self.assertEqual(csutils.particleKey(item), csutils.csParticleKey('1@J1/imported/012345678901234567890_mic_001.aligned.mrcs'))
        self.assertNotEqual(csutils.particleKey(item), csutils.csParticleKey('1@J1/imported/012345678901234567890_mic_001.mrcs'))


class TestClearingIntermediateResults(unittest.TestCase):

    def testClearingAJobNeverClearsTheWholeProject(self):
        tried = []

        def cli(expression, printCmd=True):
            tried.append(expression)
            raise Exception('No such API')

        with patch.object(csutils, '_isCryosparcV5OrNewer', return_value=True), \
                patch.object(csutils, '_runCliValue', side_effect=cli):
            csutils.clearIntermediateResults('P1', 'J5')

        self.assertTrue(tried)
        self.assertEqual([e for e in tried if 'J5' not in e], [],
                         'Other runs of the project would lose their intermediate results.')

    def testAFailedClearingDoesNotFailTheJob(self):
        with patch.object(csutils, '_isCryosparcV5OrNewer', return_value=False), \
                patch.object(csutils, '_runCliExpression', side_effect=Exception('cli error')):
            csutils.clearIntermediateResults('P1', 'J5')


class TestFilesFromCryosparc(unittest.TestCase):
    """Copies and conversions of the files a cryoSPARC job leaves."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.src = os.path.join(self.tmp, 'J5')
        self.dst = os.path.join(self.tmp, 'extra')
        os.makedirs(self.src)
        os.makedirs(self.dst)

    def _write(self, folder, name, text='x'):
        with open(os.path.join(folder, name), 'w') as f:
            f.write(text)

    def testCopyingAFileThatDoesNotExistFails(self):
        with self.assertRaises(Exception):
            csutils.copyFiles(self.src, self.dst, files=['J5_particles.cs'])

    def testCopyingAJobFolderOverAPreviousCopyUpdatesIt(self):
        self._write(self.dst, 'J5_volume.mrc', 'old')  # Left by a previous execution of the step
        self._write(self.src, 'J5_volume.mrc', 'new')

        csutils.copyFiles(self.src, self.dst)

        with open(os.path.join(self.dst, 'J5_volume.mrc')) as f:
            self.assertEqual(f.read(), 'new')

    def testAFailedConversionFailsAndLeavesNoStaleStar(self):
        from cryosparc2.convert import convert
        csFile = os.path.join(self.tmp, 'J5_particles.cs')
        starFile = os.path.join(self.tmp, 'particles.star')
        self._write(self.tmp, 'J5_particles.cs')
        self._write(self.tmp, 'particles.star', 'from a previous execution')
        process = Mock(returncode=1)
        process.communicate.return_value = (b'', b'Traceback: mapping error')

        with patch.object(convert.subprocess, 'Popen', return_value=process), \
                patch.object(convert.Plugin, 'getCondaActivationCmd', return_value=''), \
                patch.object(convert.Plugin, 'getPyemEnvActivation', return_value='true'):
            with self.assertRaises(Exception):
                convert.convertCs2Star([csFile, starFile])
        self.assertFalse(os.path.exists(starFile), 'A star from a previous execution would be read as the result.')


class TestParticlesConversion(unittest.TestCase):

    def testParticlesWithoutMagnificationCanBeWritten(self):
        from cryosparc2.convert import convert
        imgSet = Mock()
        imgSet.getSamplingRate.return_value = 1.0
        imgSet.getAcquisition.return_value.getMagnification.return_value = None
        with patch.object(convert, 'convertBinaryFiles', return_value={}), \
                patch.object(convert, 'setOfImagesToMd'), patch.object(convert.md, 'MetaData'):
            convert.cryosPARCwriteSetOfParticles(imgSet, 'particles.star', 'out', fillMagnification=True)

    def testAnEmptyResultSaysSo(self):
        from cryosparc2.convert import convert
        starFile = os.path.join(tempfile.mkdtemp(), 'particles.star')
        with open(starFile, 'w') as f:
            f.write('data_particles\n\nloop_\n_rlnImageName #1\n')

        with self.assertRaisesRegex(Exception, 'No particles'):
            convert.readSetOfParticles(starFile, Mock(), alignType=None)


class TestProtocolsMenu(unittest.TestCase):

    def testTheMenuListsEveryProtocolAndOnlyThem(self):
        import json
        import re
        import cryosparc2
        import cryosparc2.protocols as protocols
        from cryosparc2.protocols import ProtCryosparcBase
        with open(os.path.join(os.path.dirname(cryosparc2.__file__), 'protocols.conf')) as f:
            listed = set(re.findall(r'"value":\s*"(\w+)"', f.read()))
        available = {name for name, obj in vars(protocols).items()
                     if isinstance(obj, type) and issubclass(obj, ProtCryosparcBase) and obj is not ProtCryosparcBase}

        self.assertEqual(sorted(listed - available), [], 'The menu offers protocols that do not exist.')
        self.assertEqual(sorted(available - listed), [], 'Protocols only reachable by searching them.')


if __name__ == '__main__':
    unittest.main()
