"""AlgorithmRunner: executa uma unidade (repetition/trial/sample).
Assume retorno (center, distance, metadata) para CSPAlgorithm.
"""

from __future__ import annotations
from typing import Any, Dict, Tuple
import time
import traceback
from src.domain.algorithms import global_registry, CSPAlgorithm


class AlgorithmExecutionError(Exception):
    pass


def run_algorithm(
    algorithm_name: str,
    strings: list[str],
    alphabet: str,
    params: Dict[str, Any],
) -> Dict[str, Any]:
    t0 = time.time()
    try:
        cls = global_registry.get(algorithm_name)
        if cls is None:
            raise AlgorithmExecutionError(f"Algorithm not registered: {algorithm_name}")
        alg: CSPAlgorithm = cls(strings, alphabet, **params)
        center, dist, meta = alg.run()
        duration = time.time() - t0
        return {
            "status": "ok",
            "objective": float(dist),  # padrão: distância int -> float
            "center": center,
            "distance": dist,
            "duration_s": duration,
            "metadata": meta,
        }
    except Exception as e:  # noqa: BLE001
        duration = time.time() - t0
        return {
            "status": "error",
            "error_type": e.__class__.__name__,
            "error_message": str(e),
            "traceback": traceback.format_exc(limit=5),
            "duration_s": duration,
            "objective": float("inf"),
        }
