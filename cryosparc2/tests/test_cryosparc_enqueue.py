"""Queuing a job on cryoSPARC 5: chosen GPUs go with the lane's machine (5.1), else with the lane (5.0), else
the job goes to the lane alone."""
import unittest
from unittest.mock import patch

from cryosparc2 import utils


class _Cli:
    """cryoSPARC's cli: answers the targets, and refuses the enqueue calls listed in refuse."""

    def __init__(self, targets, refuse=()):
        self.targets = targets
        self.refuse = refuse
        self.enqueued = []

    def __call__(self, expression, printCmd=True):
        if "find_targets" in expression:
            return self.targets
        if "api.jobs.enqueue" in expression:
            self.enqueued.append(expression)
            if any(word in expression for word in self.refuse):
                raise Exception("api.jobs.enqueue failed --> Exit code 1, message code 422")
            return "J9"
        raise AssertionError(expression)


ONE_NODE = [("gpu1", "gpu1.cnb.csic.es", "node")]


class TestEnqueueV5(unittest.TestCase):

    def _enqueue(self, cli, gpus, standalone=True):
        with patch.object(utils, "_runCliValue", cli):
            return utils._enqueueV5("P1", "J9", "default", gpus, standalone)

    def test_ChosenGpusGoWithTheLanesMachine(self):
        cli = _Cli(ONE_NODE)
        self._enqueue(cli, [0, 1])
        self.assertEqual(len(cli.enqueued), 1)
        self.assertIn("hostname='gpu1', gpus=[0, 1]", cli.enqueued[0])
        self.assertNotIn("lane=", cli.enqueued[0])

    def test_WithoutChosenGpusTheJobGoesToTheLane(self):
        cli = _Cli(ONE_NODE)
        self._enqueue(cli, [])
        self.assertEqual(len(cli.enqueued), 1)
        self.assertIn("lane='default'", cli.enqueued[0])
        self.assertNotIn("gpus", cli.enqueued[0])

    def test_AnEarlierReleaseTakesTheGpusWithTheLane(self):
        cli = _Cli(ONE_NODE, refuse=("hostname=",))
        self._enqueue(cli, [0])
        self.assertEqual(len(cli.enqueued), 2)
        self.assertIn("lane='default', gpus=[0]", cli.enqueued[1])

    def test_AtWorstTheJobGoesToTheLaneAndCryosparcChoosesTheGpus(self):
        cli = _Cli(ONE_NODE, refuse=("gpus=",))
        with self.assertLogs(utils.logger, level="WARNING"):
            self._enqueue(cli, [0])
        self.assertEqual(cli.enqueued[-1], "api.jobs.enqueue('P1', 'J9', lane='default', no_check_inputs_ready=False)")

    def test_ALaneOfSeveralMachinesUsesThisOne(self):
        cli = _Cli([("a", "a.example.org", "node"), ("b", "b.example.org", "node")])
        with patch("socket.gethostname", return_value="b"), patch("socket.getfqdn", return_value="b.example.org"):
            self._enqueue(cli, [0], standalone=False)
        self.assertIn("hostname='b'", cli.enqueued[0])

    def test_AClusterLaneHasNoMachineForChosenGpus(self):
        cli = _Cli([("slurm", "login", "cluster")], refuse=("gpus=",))
        with self.assertLogs(utils.logger, level="WARNING"):
            self._enqueue(cli, [0], standalone=False)
        self.assertNotIn("hostname=", " ".join(cli.enqueued))
        self.assertIn("lane='default', no_check_inputs_ready=False", cli.enqueued[-1])


if __name__ == "__main__":
    unittest.main()
