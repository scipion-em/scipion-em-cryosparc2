"""Regression tests for resolving the cryoSPARC master CLI.

The subject function is compiled from the repository source, so this test suite
requires neither an installed cryoSPARC server nor Scipion's runtime imports.
"""

import ast
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


UTILS = Path(__file__).resolve().parents[1] / "utils.py"


def load_command_resolver(home):
    tree = ast.parse(UTILS.read_text(encoding="utf-8"), filename=str(UTILS))
    functions = [
        node for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == "getCryosparcProgram"
    ]
    if len(functions) != 1:
        raise AssertionError("Expected exactly one getCryosparcProgram in cryosparc2/utils.py")

    # Inject the same constants and dependencies that getCryosparcProgram uses,
    # without importing Scipion or connecting to cryoSPARC.
    namespace = {
        "os": os,
        "CRYOSPARCM_CMD": "CRYOSPARCM_CMD",
        "CRYOSPARC_MASTER": "cryosparc_master",
        "getCryosparcDir": lambda *parts: os.path.join(home, *parts) if home is not None else None,
    }
    isolated = ast.Module(body=functions, type_ignores=[])
    exec(compile(ast.fix_missing_locations(isolated), str(UTILS), "exec"), namespace)
    return namespace["getCryosparcProgram"]


class CryosparcCommandResolutionTest(unittest.TestCase):
    def setUp(self):
        # An empty override is the ordinary installation path. Prevent an
        # existing shell variable from making these tests machine-specific.
        self.env = patch.dict(os.environ, {"CRYOSPARCM_CMD": ""})
        self.env.start()
        self.addCleanup(self.env.stop)
        self.folder = tempfile.TemporaryDirectory(prefix="cryosparc-command-test-")
        self.addCleanup(self.folder.cleanup)
        self.root = Path(self.folder.name)

    def make_master_command(self, executable=True):
        binary = self.root / "cryosparc_master" / "bin" / "cryosparcm"
        binary.parent.mkdir(parents=True, exist_ok=True)
        binary.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
        binary.chmod(0o755 if executable else 0o644)
        return binary

    def assert_configuration_error(self, resolve, must_mention):
        with self.assertRaises((RuntimeError, FileNotFoundError, ValueError)) as caught:
            resolve()
        message = str(caught.exception)
        self.assertTrue(
            all(term in message for term in must_mention),
            "Configuration error should mention %s; got %r" % (must_mention, message),
        )
        self.assertNotIn("cryosparc2_master", message)

    def test_resolves_real_master_executable_for_cli_and_direct_commands(self):
        binary = self.make_master_command()
        resolve = load_command_resolver(str(self.root))
        self.assertEqual(resolve(), str(binary) + " cli")
        self.assertEqual(resolve(""), str(binary))
        self.assertEqual(resolve("status"), str(binary) + " status")

    def test_preserves_explicit_wrapper_override_without_home_directory(self):
        resolve = load_command_resolver("")
        with patch.dict(os.environ, {"CRYOSPARCM_CMD": "/usr/local/bin/cryosparc-wrapper"}):
            self.assertEqual(resolve(), "/usr/local/bin/cryosparc-wrapper cli")
            self.assertEqual(resolve(""), "/usr/local/bin/cryosparc-wrapper")
            self.assertEqual(resolve("status"), "/usr/local/bin/cryosparc-wrapper status")

    def test_rejects_empty_home_instead_of_inventing_legacy_command(self):
        resolve = load_command_resolver("")
        self.assert_configuration_error(resolve, ("CRYOSPARC_HOME",))

    def test_rejects_missing_master_even_if_old_legacy_directory_exists(self):
        legacy = self.root / "cryosparc2_master" / "bin" / "cryosparcm"
        legacy.parent.mkdir(parents=True, exist_ok=True)
        legacy.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
        legacy.chmod(0o755)
        resolve = load_command_resolver(str(self.root))
        self.assert_configuration_error(resolve, ("cryosparc_master",))

    def test_rejects_missing_binary_in_existing_master_bin_directory(self):
        (self.root / "cryosparc_master" / "bin").mkdir(parents=True)
        resolve = load_command_resolver(str(self.root))
        self.assert_configuration_error(resolve, ("cryosparcm",))

    def test_rejects_nonexecutable_master_binary(self):
        self.make_master_command(executable=False)
        resolve = load_command_resolver(str(self.root))
        self.assert_configuration_error(resolve, ("cryosparcm",))


if __name__ == "__main__":
    unittest.main()
