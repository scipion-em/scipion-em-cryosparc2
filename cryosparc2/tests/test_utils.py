import getpass
import unittest
from unittest.mock import patch

from cryosparc2 import V2_14_2, V_UNKNOWN
from cryosparc2.utils import (cryosparcValidate, cryosparcExists,
                              isCryosparcRunning, calculateNewSamplingRate,
                              getProjectName, getCryosparcVersion)

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
                    getVersion.return_value = V2_14_2
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

                version = getCryosparcVersion()
                self.assertEqual(version, V_UNKNOWN)
                getFromFile.assert_called_once()
                getEnvInfo.assert_called_once()


if __name__ == '__main__':
    unittest.main()
