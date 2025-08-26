"""DP-CSP exact Closest String via iterative decision procedure.

Features:
- Unified AlgorithmResult contract
- Search strategies: linear (0..max_d) or binary (first feasible d)
- Single complexity guard: (d+1)^n <= complexity_abort_threshold
- Structured progress + metadata (optimal d, iterations, strategy)
"""

from __future__ import annotations

import time
from typing import Any

from src.domain.algorithms import AlgorithmResult, CSPAlgorithm, register_algorithm

from .config import DP_CSP_DEFAULTS
from .implementation import _dp_decision


@register_algorithm
class DPCSPAlgorithm(CSPAlgorithm):
    name = "DP-CSP"
    default_params = DP_CSP_DEFAULTS
    supports_internal_parallel = False
    deterministic = True

    def __init__(
        self,
        strings: list[str],
        alphabet: str,
        distance_calculator,
        monitor=None,
        seed: int | None = None,
        internal_jobs: int = 1,
        **params: Any,
    ) -> None:
        if strings:
            lens = {len(s) for s in strings}
            if len(lens) > 1:
                raise ValueError(
                    f"DP-CSP requires strings of equal length. Found lengths: {sorted(lens)}"
                )
        super().__init__(
            strings,
            alphabet,
            distance_calculator=distance_calculator,
            monitor=monitor,
            seed=seed,
            internal_jobs=internal_jobs,
            **params,
        )

    def run(self) -> AlgorithmResult:  # type: ignore[override]
        start_time = time.time()
        try:
            if not self.strings:
                raise ValueError("String list cannot be empty")
            if not self.alphabet:
                raise ValueError("Alphabet cannot be empty")

            L = len(self.strings[0])
            n = len(self.strings)
            if self._monitor:
                self._monitor.on_progress(
                    0.0,
                    "Starting DP-CSP",
                    phase="start",
                    n=n,
                    L=L,
                    alphabet_size=len(self.alphabet),
                )

            # Verificação inicial de cancelamento
            if self._monitor and self._monitor.is_cancelled():
                return self._build_cancelled_result(start_time)

            max_d_param = self.params.get("max_d")
            baseline_upper = self.max_distance(self.strings[0])
            if max_d_param is None:
                self.params["max_d"] = baseline_upper
            max_d = self.params["max_d"]
            complexity_abort_threshold = self.params.get(
                "complexity_abort_threshold",
                DP_CSP_DEFAULTS["complexity_abort_threshold"],
            )
            search_strategy = self.params.get(
                "search_strategy", DP_CSP_DEFAULTS.get("search_strategy", "linear")
            )
            if self._monitor:
                self._monitor.on_progress(
                    0.05,
                    f"Parameters set max_d={max_d}",
                    phase="parameters",
                    max_d=max_d,
                    baseline_upper_bound=baseline_upper,
                    complexity_abort_threshold=complexity_abort_threshold,
                    search_strategy=search_strategy,
                )

            solution: str | None = None
            optimal_d: int | None = None
            iterations_tested = 0

            def abort_complexity(d: int, complexity_est: int) -> AlgorithmResult:
                self._report_warning(
                    f"Complexity exceeded at d={d}: (d+1)^n={complexity_est:,} > {complexity_abort_threshold:,}"
                )
                execution_time = time.time() - start_time
                return AlgorithmResult(
                    success=False,
                    center_string="",
                    max_distance=-1,
                    parameters=self.get_actual_params(),
                    error=f"Aborted due to complexity: (d+1)^n={complexity_est:,} exceeds {complexity_abort_threshold:,}",
                    metadata={
                        "algorithm_name": self.name,
                        "execution_time": execution_time,
                        "num_strings": n,
                        "string_length": L,
                        "alphabet_size": len(self.alphabet),
                        "baseline_upper_bound": baseline_upper,
                        "iterations_tested": iterations_tested,
                        "complexity_estimate_last": complexity_est,
                        "termination_reason": "complexity_abort",
                        "deterministic": True,
                        "seed": self.seed,
                        "internal_jobs": self.internal_jobs,
                        "search_strategy": search_strategy,
                    },
                )

            if search_strategy == "linear":
                for d in range(max_d + 1):
                    # Verificar cancelamento a cada iteração
                    if self._monitor and self._monitor.is_cancelled():
                        return self._build_cancelled_result(start_time)

                    iterations_tested += 1
                    complexity_est = (d + 1) ** n
                    if complexity_est > complexity_abort_threshold:
                        return abort_complexity(d, complexity_est)
                    progress = 0.10 + (d / max_d) * 0.80 if max_d > 0 else 0.90
                    if self._monitor:
                        self._monitor.on_progress(
                            progress,
                            f"Testing d={d}",
                            phase="search",
                            current_d=d,
                            max_d=max_d,
                            complexity_estimate=complexity_est,
                            iterations_tested=iterations_tested,
                            search_strategy=search_strategy,
                        )
                    center = _dp_decision(self.strings, self.alphabet, d)
                    if center is not None:
                        solution = center
                        optimal_d = d
                        break
            elif search_strategy == "binary":
                low, high = 0, max_d
                feasible_center: str | None = None
                feasible_d: int | None = None
                total_initial = max_d + 1
                while low <= high:
                    # Verificar cancelamento a cada iteração do binary search
                    if self._monitor and self._monitor.is_cancelled():
                        return self._build_cancelled_result(start_time)

                    mid = (low + high) // 2
                    iterations_tested += 1
                    complexity_est = (mid + 1) ** n
                    if complexity_est > complexity_abort_threshold:
                        return abort_complexity(mid, complexity_est)
                    interval_size = high - low + 1
                    searched_fraction = 1 - (interval_size / total_initial)
                    progress = 0.10 + min(1.0, searched_fraction) * 0.80
                    if self._monitor:
                        self._monitor.on_progress(
                            progress,
                            f"Testing d={mid} (binary)",
                            phase="search",
                            current_d=mid,
                            max_d=max_d,
                            complexity_estimate=complexity_est,
                            iterations_tested=iterations_tested,
                            search_strategy=search_strategy,
                            interval_low=low,
                            interval_high=high,
                        )
                    center = _dp_decision(self.strings, self.alphabet, mid)
                    if center is not None:
                        feasible_center = center
                        feasible_d = mid
                        high = mid - 1
                    else:
                        low = mid + 1
                if feasible_center is not None:
                    solution = feasible_center
                    optimal_d = feasible_d
            else:  # pragma: no cover
                raise ValueError(f"invalid search_strategy: {search_strategy}")

            if solution is None:
                execution_time = time.time() - start_time
                if self._monitor:
                    self._monitor.on_warning("No solution found within max_d")
                return AlgorithmResult(
                    success=False,
                    center_string="",
                    max_distance=-1,
                    parameters=self.get_actual_params(),
                    error="No solution found up to max_d",
                    metadata={
                        "algorithm_name": self.name,
                        "execution_time": execution_time,
                        "num_strings": n,
                        "string_length": L,
                        "alphabet_size": len(self.alphabet),
                        "baseline_upper_bound": baseline_upper,
                        "iterations_tested": iterations_tested,
                        "complexity_estimate_last": (max_d + 1) ** n,
                        "termination_reason": "no_solution_within_max_d",
                        "deterministic": True,
                        "seed": self.seed,
                        "internal_jobs": self.internal_jobs,
                        "search_strategy": search_strategy,
                    },
                )

            if self._monitor:
                self._monitor.on_progress(
                    0.95, "Validating solution", phase="validation", optimal_d=optimal_d
                )
            assert solution is not None and optimal_d is not None
            actual_max_distance = self.max_distance(solution)
            if actual_max_distance != optimal_d:
                self._report_warning(
                    f"Inconsistency: optimal_d={optimal_d} but actual distance={actual_max_distance}"
                )
            avg_dist = self.average_distance(solution)
            total_dist = self.total_distance(solution)
            execution_time = time.time() - start_time
            if self._monitor:
                self._monitor.on_progress(1.0, "DP-CSP finished", phase="finish")
            return AlgorithmResult(
                success=True,
                center_string=solution,
                max_distance=actual_max_distance,
                parameters=self.get_actual_params(),
                error=None,
                metadata={
                    "algorithm_name": self.name,
                    "execution_time": execution_time,
                    "num_strings": n,
                    "string_length": L,
                    "alphabet_size": len(self.alphabet),
                    "baseline_upper_bound": baseline_upper,
                    "optimal_d": optimal_d,
                    "iterations_tested": iterations_tested,
                    "complexity_estimate_last": (
                        (optimal_d + 1) ** n if optimal_d is not None else None
                    ),
                    "validation_max_distance": actual_max_distance,
                    "avg_distance": avg_dist,
                    "total_distance": total_dist,
                    "termination_reason": "optimal_found",
                    "deterministic": True,
                    "seed": self.seed,
                    "internal_jobs": self.internal_jobs,
                    "search_strategy": search_strategy,
                },
            )
        except Exception as e:  # general error
            execution_time = time.time() - start_time
            msg = f"Error during DP-CSP execution: {e}"
            if self._monitor:
                self._monitor.on_warning(msg)
                self._monitor.on_progress(
                    0.0, msg, phase="error", error_type=type(e).__name__
                )
            return AlgorithmResult(
                success=False,
                center_string="",
                max_distance=-1,
                parameters=self.get_actual_params(),
                error=msg,
                metadata={
                    "algorithm_name": self.name,
                    "execution_time": execution_time,
                    "num_strings": len(self.strings),
                    "string_length": len(self.strings[0]) if self.strings else 0,
                    "alphabet_size": len(self.alphabet) if self.alphabet else 0,
                    "termination_reason": "error",
                    "deterministic": True,
                    "seed": self.seed,
                    "internal_jobs": self.internal_jobs,
                    "error_type": type(e).__name__,
                    "search_strategy": self.params.get("search_strategy"),
                },
            )

    def _build_cancelled_result(self, start_time: float) -> AlgorithmResult:
        """Constrói resultado para execução cancelada."""
        execution_time = time.time() - start_time
        return AlgorithmResult(
            success=False,
            center_string="",
            max_distance=-1,
            parameters=self.get_actual_params(),
            error="Algorithm execution was cancelled",
            metadata={
                "algorithm_name": self.name,
                "execution_time": execution_time,
                "num_strings": len(self.strings),
                "string_length": len(self.strings[0]) if self.strings else 0,
                "alphabet_size": len(self.alphabet) if self.alphabet else 0,
                "deterministic": True,
                "seed": self.seed,
                "internal_jobs": self.internal_jobs,
                "termination_reason": "cancelled",
            },
        )
