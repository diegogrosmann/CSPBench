import time
from typing import Any, Dict

from src.infrastructure.orchestration.algorithm_runner import run_algorithm
from src.domain.algorithms import CSPAlgorithm, register_algorithm, global_registry
from src.domain.distance import create_distance_calculator
from src.domain.status import BaseStatus


# --- Stubs -----------------------------------------------------------------
class StubMonitor:
    def __init__(self):
        self.progress = []
        self.warnings = []
        self.errors = []
        self._cancelled = False

    def on_progress(self, p: float, msg: str, *_, **__):
        self.progress.append((p, msg))

    def on_warning(self, msg: str, *_, **__):
        self.warnings.append(msg)

    def on_error(self, msg: str, exc: Exception | None = None, *_, **__):
        self.errors.append((msg, exc.__class__.__name__ if exc else None))

    def cancel(self):
        self._cancelled = True

    def is_cancelled(self) -> bool:
        return self._cancelled


class StubExecutionController:
    def __init__(self, status_sequence=None, internal_jobs: int = 1):
        self._status_sequence = status_sequence or [BaseStatus.RUNNING]
        self._idx = 0
        self.internal_jobs = internal_jobs
        # Add missing timeout_per_item attribute for compatibility
        self.timeout_per_item = 30.0

    def check_status(self):
        if self._idx < len(self._status_sequence) - 1:
            current = self._status_sequence[self._idx]
            self._idx += 1
            return current
        return self._status_sequence[-1]
    
    # Add missing method for compatibility
    def item_timeout(self):
        # Return a no-op context manager
        from contextlib import nullcontext
        return nullcontext()


# --- Fake Algorithms --------------------------------------------------------
@register_algorithm
class _AlgSuccess(CSPAlgorithm):  # type: ignore
    name = "_alg_success"

    def run(self):  # type: ignore[override]
        return {
            "success": True,
            "center_string": self.strings[0],
            "max_distance": 0,
            "parameters": self.params,
            "error": None,
            "metadata": {"runs": 1},
        }


@register_algorithm
class _AlgSleep(CSPAlgorithm):  # type: ignore
    name = "_alg_sleep"

    def run(self):  # type: ignore[override]
        # Loop que respeita cancelamento
        start = time.time()
        while time.time() - start < 2:  # 2s planned
            if self._monitor and self._monitor.is_cancelled():  # type: ignore[attr-defined]
                break
            time.sleep(0.05)
        return {
            "success": True,
            "center_string": self.strings[0],
            "max_distance": 0,
            "parameters": self.params,
            "error": None,
            "metadata": {"slept": True},
        }


@register_algorithm
class _AlgFail(CSPAlgorithm):  # type: ignore
    name = "_alg_fail"

    def run(self):  # type: ignore[override]
        raise RuntimeError("boom")


# --- Helpers ----------------------------------------------------------------
def _basic_inputs():
    strings = ["AAAA", "AAAA"]
    alphabet = "A"
    dist = create_distance_calculator("hamming", strings, use_cache=False)
    return strings, alphabet, dist


# --- Tests ------------------------------------------------------------------
def test_run_algorithm_success():
    strings, alphabet, dist = _basic_inputs()
    monitor = StubMonitor()
    controller = StubExecutionController()
    res = run_algorithm(
        algorithm_name="_alg_success",
        strings=strings,
        alphabet=alphabet,
        distance_calculator=dist,
        execution_controller=controller,
        monitor=monitor,
        params={},
    )
    assert res["status"] == BaseStatus.COMPLETED
    assert res["algorithm_result"]["success"] is True


def test_run_algorithm_timeout():
    strings, alphabet, dist = _basic_inputs()
    monitor = StubMonitor()
    controller = StubExecutionController()
    res = run_algorithm(
        algorithm_name="_alg_sleep",
        strings=strings,
        alphabet=alphabet,
        distance_calculator=dist,
        execution_controller=controller,
        monitor=monitor,
        params={"max_time": 0.2},  # força timeout rápido
    )
    assert res["timeout_triggered"] is True
    assert res["status"] in (BaseStatus.CANCELED, BaseStatus.ERROR)


def test_run_algorithm_canceled_externally():
    strings, alphabet, dist = _basic_inputs()
    monitor = StubMonitor()
    # Sequência: primeiro running, depois canceled
    controller = StubExecutionController([BaseStatus.RUNNING, BaseStatus.CANCELED])
    res = run_algorithm(
        algorithm_name="_alg_sleep",
        strings=strings,
        alphabet=alphabet,
        distance_calculator=dist,
        execution_controller=controller,
        monitor=monitor,
        params={"max_time": 5},
    )
    assert res["status"] == BaseStatus.CANCELED


def test_run_algorithm_failure():
    strings, alphabet, dist = _basic_inputs()
    monitor = StubMonitor()
    controller = StubExecutionController()
    res = run_algorithm(
        algorithm_name="_alg_fail",
        strings=strings,
        alphabet=alphabet,
        distance_calculator=dist,
        execution_controller=controller,
        monitor=monitor,
        params={},
    )
    assert res["status"] in (BaseStatus.ERROR, BaseStatus.FAILED)
    assert res["error"]
