"""Test configuration and global fixtures.

Provides automatic cleanup of temporary test artifact directories that some
tests create via environment-variable based path resolution.

Rationale:
  Certain tests previously used relative paths like ./test_datasets which
  caused persistent directories (test_datasets, test_batches, test_outputs,
  test_data, relative) to remain in the repository root after the test run.
  This fixture removes those directories at the end of the test session to
  keep the workspace clean.

If you need to inspect these directories after a run, set the environment
variable KEEP_TEST_ARTIFACTS=1 before running pytest.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path

import pytest

_ARTIFACT_DIR_NAMES = [
    "test_batches",
    "test_data",
    "test_datasets",
    "test_outputs",
    "relative",  # created by some path normalization tests
]


@pytest.fixture(scope="session", autouse=True)
def _cleanup_test_artifacts():  # pragma: no cover - housekeeping
    """Session fixture that cleans up known temporary test directories.

    Runs once after the full test session. Safe-guards:
      * Only deletes directories whose names are in the allowlist above.
      * Only deletes if located directly under the repository root.
    """

    yield  # Run tests first

    if os.getenv("KEEP_TEST_ARTIFACTS"):
        return

    repo_root = Path.cwd()
    for name in _ARTIFACT_DIR_NAMES:
        candidate = repo_root / name
        try:
            if candidate.exists() and candidate.is_dir():
                # Extra safety: ensure we are not pointing outside repo
                if repo_root == candidate.parent:
                    shutil.rmtree(candidate, ignore_errors=True)
        except Exception:  # noqa: BLE001
            # Best effort cleanup; ignore failures
            pass


"""Config pytest: garante que diretório raiz esteja no sys.path para importar 'algorithms'."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# Força auto-discovery dos algoritmos
try:
    import algorithms  # noqa: F401
except Exception:  # pragma: no cover
    pass
