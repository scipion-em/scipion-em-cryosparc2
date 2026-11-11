"""Isolated contract tests for cryoSPARC -> Scipion log streaming.

Only the real functions are extracted from utils.py via Python AST, so no
cryoSPARC server, Scipion runtime, GPU, or plugin imports are required.
"""

import ast
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, call


FUNCTIONS = {
    "waitForCryosparc",
    "_getCryosparcJobLogEvents",
    "_getCryosparcLogText",
    "_logCryosparcJobEvents",
}


def load_real_functions():
    source_path = Path(__file__).resolve().parents[1] / "utils.py"
    source = ast.parse(source_path.read_text(encoding="utf-8"), filename=str(source_path))
    selected = [node for node in source.body if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name in FUNCTIONS]
    missing = FUNCTIONS.difference(node.name for node in selected)
    if missing:
        raise AssertionError("Missing cryoSPARC log functions: %s" % sorted(missing))
    namespace = {}
    exec(compile(ast.Module(body=selected, type_ignores=[]), str(source_path), "exec"), namespace)
    return namespace


class ProtocolCursor:
    def __init__(self):
        self._logLastLine = 0

    def getLogLine(self):
        return self._logLastLine

    def setLogLine(self, count):
        self._logLastLine = count


class CryosparcLogStreamingTest(unittest.TestCase):
    def setUp(self):
        self.functions = load_real_functions()
        self.logger = Mock()
        self.functions.update({
            "logger": self.logger,
            "time": SimpleNamespace(sleep=Mock()),
            "STOP_STATUSES": ("completed", "failed", "aborted", "killed"),
            "ACTIVE_STATUSES": ("queued", "running", "started", "launched", "building", "waiting"),
            "STATUS_COMPLETED": "completed",
            "_pythonLiteral": lambda value, default=None: ast.literal_eval(value) if isinstance(value, str) else value,
        })

    def test_polls_events_while_v5_job_is_still_running(self):
        # A blocking waitJob call would suppress all intermediate log refreshes.
        self.functions["getJobStatus"] = Mock(side_effect=["running", "running", "completed"])
        self.functions["waitJob"] = Mock()
        events = Mock()
        self.functions["_logCryosparcJobEvents"] = events
        protocol = ProtocolCursor()

        result = self.functions["waitForCryosparc"]("P1", "J1", "Job failed", protocol=protocol, sleepTime=0)

        self.assertEqual(result, "completed")
        self.functions["waitJob"].assert_not_called()
        self.assertEqual(events.call_args_list, [call("P1", "J1", protocol=protocol)] * 3)
        self.assertEqual(self.functions["time"].sleep.call_count, 2)

    def test_does_not_repeat_the_last_event_when_no_new_event_arrives(self):
        self.functions["_getCryosparcJobLogEvents"] = Mock(return_value=[
            {"text": "Starting"}, {"message": "Processing"},
        ])
        protocol = ProtocolCursor()

        self.functions["_logCryosparcJobEvents"]("P1", "J1", protocol=protocol)
        self.functions["_logCryosparcJobEvents"]("P1", "J1", protocol=protocol)

        self.assertEqual(self.logger.info.call_args_list, [call("Starting"), call("Processing")])

    def test_tracks_event_offsets_separately_for_consecutive_jobs(self):
        events_by_job = {
            "J1": [{"text": "job1-first"}, {"text": "job1-second"}],
            "J2": [{"text": "job2-first"}],
        }
        self.functions["_getCryosparcJobLogEvents"] = Mock(side_effect=lambda project, job: list(events_by_job[job]))
        protocol = ProtocolCursor()

        self.functions["_logCryosparcJobEvents"]("P1", "J1", protocol=protocol)
        self.functions["_logCryosparcJobEvents"]("P1", "J2", protocol=protocol)
        events_by_job["J2"].extend([{"text": "job2-second"}, {"text": "job2-third"}])
        self.functions["_logCryosparcJobEvents"]("P1", "J2", protocol=protocol)

        self.assertEqual(self.logger.info.call_args_list, [
            call("job1-first"), call("job1-second"), call("job2-first"),
            call("job2-second"), call("job2-third"),
        ])

    def test_uses_v5_cli_events_when_structured_api_is_unavailable(self):
        self.functions["_isCryosparcV5OrNewer"] = Mock(return_value=True)
        self.functions["_tryV5StructuredEventLog"] = Mock(return_value=None)
        self.functions["_runCryosparcmSubcommand"] = Mock(return_value=(0, "Event one\nEvent two\n"))

        events = self.functions["_getCryosparcJobLogEvents"]("P1", "J1")

        self.assertEqual([self.functions["_getCryosparcLogText"](event) for event in events], ["Event one", "Event two"])
        self.functions["_runCryosparcmSubcommand"].assert_called_once_with("job", "events", "P1", "J1", printCmd=False)

    def test_keeps_structured_v5_events_compatible(self):
        expected = [{"text": "Queued"}, {"message": "Running"}]
        self.functions["_isCryosparcV5OrNewer"] = Mock(return_value=True)
        self.functions["_tryV5StructuredEventLog"] = Mock(return_value=expected)

        self.assertEqual(self.functions["_getCryosparcJobLogEvents"]("P1", "J1"), expected)

    def test_preserves_legacy_event_retrieval(self):
        expected = [{"text": "Legacy event"}]
        self.functions["_isCryosparcV5OrNewer"] = Mock(return_value=False)
        self.functions["getJobStreamlog"] = Mock(return_value=(0, repr(expected)))

        self.assertEqual(self.functions["_getCryosparcJobLogEvents"]("P1", "J1"), expected)


if __name__ == "__main__":
    unittest.main()
