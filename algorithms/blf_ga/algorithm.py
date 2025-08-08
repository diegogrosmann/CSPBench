"""
BLF-GA: Blockwise Learning Fusion + Genetic Algorithm for CSP.

Classes:
    BLFGAAlgorithm: Wrapper for BLF-GA integration into CSP framework.
"""

from collections.abc import Callable

from src.domain.algorithms import CSPAlgorithm, register_algorithm

from .config import BLF_GA_DEFAULTS
from .implementation import BLFGA


@register_algorithm
class BLFGAAlgorithm(CSPAlgorithm):
    """
    BLF-GA: Blockwise Learning Fusion + Genetic Algorithm for the Closest String Problem.

    A hybrid metaheuristic that combines:
    - Blockwise learning
    - Global genetic algorithm
    - Advanced adaptive mechanisms

    Features:
    - Internal parallelism support
    - Multiple genetic operators
    - Adaptive parameter control
    - Local refinement of elites
    - Flexible stopping criteria

    Args:
        strings (list[str]): List of input strings.
        alphabet (str): Alphabet used.
        **params: Algorithm parameters.

    Methods:
        run(): Executes BLF-GA and returns (center, maximum distance, metadata).
        run_with_history(): Executes and returns complete history.
        get_metadata(): Returns detailed algorithm metadata.
    """

    name = "BLF-GA"
    default_params = BLF_GA_DEFAULTS
    supports_internal_parallel = True  # BLF-GA can use internal parallelism
    is_deterministic = False  # Is stochastic (can have seed for reproducibility)

    def __init__(self, strings: list[str], alphabet: str, **params):
        super().__init__(strings, alphabet, **params)

        # Map TEMPLATE.yaml parameter names to BLF-GA implementation names
        param_mapping = {
            "crossover_method": "crossover_type",
            "mutation_method": "mutation_type",
            "refinement_method": "refinement_type",
            # selection_method is not directly supported, ignore it
        }

        # Filter and map parameters to pass only those that BLFGA implementation knows
        blfga_params = {}
        for k, v in self.params.items():
            # Skip framework-specific history parameters
            if k in ["save_history", "history_frequency"]:
                continue
            # Skip unsupported selection_method parameter
            if k == "selection_method":
                continue
            # Map parameter name if needed
            param_name = param_mapping.get(k, k)
            blfga_params[param_name] = v

        self.blf_ga_instance = BLFGA(self.strings, self.alphabet, **blfga_params)

        # Configure history callback if enabled
        if self.save_history:
            self.blf_ga_instance.set_history_callback(self._save_dynamic_history_entry)

    def set_params(self, **params) -> None:
        """
        Set new parameters for the algorithm.

        Args:
            **params: Parameters to be updated
        """
        super().set_params(**params)
        # Update BLF-GA instance directly
        self.blf_ga_instance.update_params(**params)

    def set_progress_callback(self, callback: Callable[[str], None]) -> None:
        """
        Pass callback to BLFGA instance.

        Args:
            callback (Callable[[str], None]): Progress callback function.
        """
        super().set_progress_callback(callback)
        self.blf_ga_instance.set_progress_callback(callback)

    def set_warning_callback(self, callback: Callable[[str], None]) -> None:
        """
        Set callback for warnings and pass to BLFGA instance.

        Args:
            callback (Callable[[str], None]): Warning callback function.
        """
        super().set_warning_callback(callback)
        # BLF-GA can generate warnings about convergence, parameters, etc.
        # For now, just store the callback
        self._blfga_warning_callback = callback

    def run(self) -> tuple[str, int, dict]:
        """
        Execute BLF-GA and return center string, maximum distance and detailed metadata.

        Returns:
            tuple: (best_string, best_fitness, metadata)
                - best_string: Best string found
                - best_fitness: Maximum distance of best string
                - metadata: Dictionary with detailed execution information
        """
        import time

        start_time = time.time()

        try:
            # Clear previous history
            self.clear_history()

            self._report_progress("Starting BLF-GA...")

            # Save initial state
            self._save_history_entry(
                0,
                phase="initialization",
                population_size=self.blf_ga_instance.pop_size,
                max_generations=self.blf_ga_instance.max_gens,
                initial_fitness=None,
            )

            best, best_val, ga_history = self.blf_ga_instance.run()
            execution_time = time.time() - start_time

            # Process GA history to our standard format
            if ga_history and self.save_history:
                for i, fitness in enumerate(ga_history):
                    self._save_history_entry(
                        i + 1,
                        phase="evolution",
                        generation=i,
                        best_fitness=fitness,
                        improvement=(
                            ga_history[0] - fitness
                            if i == 0
                            else ga_history[i - 1] - fitness
                        ),
                    )

            # Save final state
            self._save_history_entry(
                len(ga_history) if ga_history else 1,
                phase="completion",
                final_solution=best,
                final_fitness=best_val,
                total_generations=len(ga_history) if ga_history else 0,
                execution_time=execution_time,
            )

            # Detailed execution metadata
            metadata = {
                "algorithm": "BLF-GA",
                "status": "completed",
                "execution_time_seconds": execution_time,
                "generations_executed": len(ga_history) if ga_history else 0,
                "best_fitness": best_val,
                "convergence_info": {
                    "initial_fitness": ga_history[0] if ga_history else None,
                    "final_fitness": best_val,
                    "improvement": (ga_history[0] - best_val) if ga_history else 0,
                    "fitness_history_length": len(ga_history) if ga_history else 0,
                },
                "algorithm_config": {
                    "population_size": self.blf_ga_instance.pop_size,
                    "max_generations": self.blf_ga_instance.max_gens,
                    "crossover_prob": self.blf_ga_instance.cross_prob,
                    "mutation_prob": self.blf_ga_instance.mut_prob,
                    "elite_rate": self.blf_ga_instance.elite_rate,
                },
                "adaptive_mechanisms_used": {
                    "immigrants": self.blf_ga_instance.immigrant_freq > 0,
                    "adaptive_mutation": True,
                    "block_redivision": self.blf_ga_instance.rediv_freq > 0,
                    "elite_refinement": self.blf_ga_instance.refine_elites != "none",
                },
                # Detailed history if enabled
                "history": self.get_history() if self.save_history else [],
                "ga_fitness_history": ga_history if ga_history else [],
                # Compatibility with previous version
                "iterations": len(ga_history) if ga_history else 0,
                "best_distance": best_val,
                "complete_history": len(ga_history) if ga_history else 0,
            }

            return best, best_val, metadata

        except Exception as e:
            execution_time = time.time() - start_time
            error_metadata = {
                "algorithm": "BLF-GA",
                "status": "error",
                "error_message": str(e),
                "execution_time_seconds": execution_time,
                "generations_executed": 0,
                "best_fitness": float("inf"),
                "history": self.get_history() if self.save_history else [],
            }
            # Re-raise the exception but with metadata available if needed
            raise RuntimeError(f"BLF-GA execution failed: {e}") from e

    def _save_blfga_history(
        self, generation: int, best_fitness: float, **kwargs
    ) -> None:
        """Callback to save history during BLF-GA execution."""
        self._save_history_entry(
            generation,
            phase="evolution",
            generation=generation,
            best_fitness=best_fitness,
            **kwargs,
        )

    def run_with_history(self) -> tuple[str, int, list]:
        """
        Execute BLF-GA and return the center string, maximum distance and distance history.
        """
        return self.blf_ga_instance.run()

    def get_metadata(self) -> dict:
        """
        Return detailed metadata of the BLF-GA algorithm.

        Returns:
            dict: Metadata including configuration and current state
        """
        base_metadata = super().get_metadata()

        # Add BLF-GA specific metadata
        blfga_metadata = {
            "algorithm_type": "hybrid_metaheuristic",
            "components": [
                "blockwise_learning",
                "genetic_algorithm",
                "adaptive_mechanisms",
            ],
            "population_size": self.blf_ga_instance.pop_size,
            "min_population_size": self.blf_ga_instance.min_pop_size,
            "initial_blocks": self.blf_ga_instance.initial_blocks,
            "crossover_type": self.blf_ga_instance.crossover_type,
            "mutation_type": self.blf_ga_instance.mutation_type,
            "refinement_type": self.blf_ga_instance.refinement_type,
            "adaptive_features": {
                "immigrant_injection": self.blf_ga_instance.immigrant_freq > 0,
                "adaptive_mutation": True,
                "adaptive_blocking": self.blf_ga_instance.rediv_freq > 0,
                "elite_refinement": self.blf_ga_instance.refine_elites != "none",
                "niching_enabled": self.blf_ga_instance.niching,
            },
            "stopping_criteria": {
                "max_generations": self.blf_ga_instance.max_gens,
                "max_time_seconds": self.blf_ga_instance.max_time,
                "early_stopping": self.blf_ga_instance.no_improve_patience > 0,
                "restart_enabled": self.blf_ga_instance.restart_patience > 0,
            },
        }

        # Merge metadata
        base_metadata.update(blfga_metadata)
        return base_metadata

    def _save_dynamic_history_entry(self, generation: int, event_data: dict) -> None:
        """
        Callback to register dynamic events during BLF-GA execution.

        Args:
            generation: Current generation
            event_data: Dictionary with event data including 'event' and other information
        """
        if self.save_history:
            # Extract event type and other data
            event_type = event_data.pop("event", "unknown_event")
            self._save_history_entry(
                generation + 1,  # +1 because iteration 0 is initialization
                phase="dynamic_event",
                generation=generation,
                event_type=event_type,
                **event_data,
            )
