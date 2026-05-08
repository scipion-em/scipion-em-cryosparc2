import unittest
from unittest.mock import patch

from cryosparc2.wizards import ProtCryosparcLanesWizard


class _ProtocolWithPreprocessLane:
    preprocess_lane = True


class _DummyForm:
    def __init__(self, protocol):
        self.protocol = protocol
        self.root = None
        self.values = {}

    def setVar(self, key, value):
        self.values[key] = value


class _DummyDialog:
    values = ['lane-a']

    def resultYes(self):
        return True


class _DummyDialogView:
    def __init__(self, *args, **kwargs):
        pass

    def show(self):
        return _DummyDialog()


class TestProtCryosparcLanesWizard(unittest.TestCase):

    @patch('cryosparc2.wizards.LanesDialogView', _DummyDialogView)
    @patch('cryosparc2.wizards.cryosparcValidate', return_value=[])
    def test_show_sets_only_requested_parameter(self, _validate):
        form = _DummyForm(_ProtocolWithPreprocessLane())

        wizard = ProtCryosparcLanesWizard()
        wizard.show(form, 'preprocess_lane')

        self.assertEqual(form.values.get('preprocess_lane'), 'lane-a')
        self.assertNotIn('compute_lane', form.values)


if __name__ == '__main__':
    unittest.main()
