from src.domain.distance import create_distance_calculator
from src.domain.algorithms import global_registry


class MemoryStore:
    def __init__(self):
        self.events = []

    def report_algorithm_progress(self, progress, message, **extra):
        self.events.append(
            {
                "progress": progress,
                "message": message,
                **extra,
            }
        )

    def warning(self, message):
        self.events.append({"warning": message})


def build_algorithm(strings, alphabet, store=None, **params):
    calc = create_distance_calculator("hamming", strings)
    AlgCls = global_registry["CSC"]
    return AlgCls(strings, alphabet, distance_calculator=calc, store=store, **params)


def test_progress_phases_cover_expected():
    strings = ["ACGT", "ACGA", "ACGG", "ACGC"]
    store = MemoryStore()
    alg = build_algorithm(strings, "ACGT", store=store)
    alg.run()
    phases = {e.get("phase") for e in store.events if "phase" in e}
    # Pelo menos as fases principais
    expected = {
        "start",
        "parameters",
        "clustering",
        "consensus",
        "candidates",
        "evaluation",
        "local_search",
        "finish",
    } - {
        "local_search"
    }  # local_search pode estar embutida
    assert expected.issubset(phases)
    # progresso final 100
    assert any(e.get("progress") == 100.0 for e in store.events)
