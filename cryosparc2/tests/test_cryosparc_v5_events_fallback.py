"""CryoSPARC 5 plain-text job-events fallback regression contracts.

Exercises actual functions extracted from utils.py; never calls CryoSPARC or
loads the Scipion runtime. The CLI fixture uses real 5.0.7 output structure.
"""
import ast
from pathlib import Path
import unittest
from unittest.mock import Mock, call


FUNCTIONS = {"_getCryosparcJobLogEvents", "_getCryosparcLogText", "_logCryosparcJobEvents"}
CLI_TEXT = (
    "[2026-10-16 16:16:01] [761 MB] Compiling job outputs...\n"
    "[2026-10-16 16:16:01] [768 MB] Job complete. Total time 343.77s\n"
)
CLI_LINES = CLI_TEXT.splitlines()


def actual_functions():
    path = Path(__file__).resolve().parents[1] / "utils.py"
    syntax = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    methods = [node for node in syntax.body
               if isinstance(node, ast.FunctionDef) and node.name in FUNCTIONS]
    missing = FUNCTIONS.difference(node.name for node in methods)
    if missing:
        raise AssertionError("Missing production functions: %r" % sorted(missing))
    space = {}
    exec(compile(ast.Module(body=methods, type_ignores=[]), str(path), "exec"), space)
    return space


class Protocol:
    def __init__(self):
        self.line = 0

    def getLogLine(self):
        return self.line

    def setLogLine(self, line):
        self.line = line


class CryosparcV5EventsFallbackTest(unittest.TestCase):
    def setUp(self):
        self.ns = actual_functions()
        self.api = Mock(return_value=None)
        self.cli = Mock(return_value=(0, CLI_TEXT))
        self.logger = Mock()
        self.ns.update({
            "_isCryosparcV5OrNewer": Mock(return_value=True),
            "_tryV5StructuredEventLog": self.api,
            "_runCryosparcmSubcommand": self.cli,
            "_pythonLiteral": lambda value, default=None: ast.literal_eval(value) if isinstance(value, str) else value,
            "logger": self.logger,
        })

    def get_texts(self):
        events = self.ns["_getCryosparcJobLogEvents"]("P1", "J32")
        return [self.ns["_getCryosparcLogText"](e) for e in events]

    def assert_cli_used(self):
        self.cli.assert_called_once_with("job", "events", "P1", "J32", printCmd=False)

    def test_falls_back_when_structured_api_returns_an_empty_list(self):
        self.api.return_value = []
        self.assertEqual(self.get_texts(), CLI_LINES)
        self.assert_cli_used()

    def test_falls_back_when_structured_api_wraps_an_empty_events_list(self):
        self.api.return_value = {"events": []}
        self.assertEqual(self.get_texts(), CLI_LINES)
        self.assert_cli_used()

    def test_falls_back_when_structured_api_events_have_no_printable_text(self):
        self.api.return_value = [{"unrelated_metadata": 4}]
        self.assertEqual(self.get_texts(), CLI_LINES)
        self.assert_cli_used()

    def test_warns_when_fallback_cli_fails_after_empty_api_response(self):
        self.api.return_value = []
        self.cli.side_effect = RuntimeError("CryoSPARC CLI unavailable")
        self.assertEqual(self.get_texts(), [])
        self.assert_cli_used()
        self.logger.warning.assert_called()

    def test_prefers_valid_structured_api_events_without_running_cli(self):
        self.api.return_value = [{"text": "Queued"}, {"message": "Running"}]
        self.assertEqual(self.get_texts(), ["Queued", "Running"])
        self.cli.assert_not_called()

    def test_delivers_real_cli_event_lines_once_when_api_is_unavailable(self):
        self.api.return_value = None
        protocol = Protocol()
        self.ns["_logCryosparcJobEvents"]("P1", "J32", protocol=protocol)
        self.ns["_logCryosparcJobEvents"]("P1", "J32", protocol=protocol)
        self.assertEqual(self.logger.info.call_args_list, [call(line) for line in CLI_LINES])


if __name__ == "__main__":
    unittest.main()
