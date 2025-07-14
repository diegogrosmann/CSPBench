"""
BLF-GA: Blockwise Learning Fusion + Genetic Algorithm para CSP.

Classes:
    BLFGAAlgorithm: Wrapper para integração do BLF-GA ao framework CSP.
"""

from collections.abc import Callable

from src.domain.algorithms import CSPAlgorithm, register_algorithm

from .config import BLF_GA_DEFAULTS
from .implementation import BLFGA


@register_algorithm
class BLFGAAlgorithm(CSPAlgorithm):
    """
    BLF-GA: Blockwise Learning Fusion + Genetic Algorithm para o Closest String Problem.

    Uma metaheurística híbrida que combina:
    - Aprendizado por blocos (blockwise learning)
    - Algoritmo genético global
    - Mecanismos adaptativos avançados

    Características:
    - Suporte a paralelismo interno
    - Múltiplos operadores genéticos
    - Controle adaptativos de parâmetros
    - Refinamento local dos elites
    - Critérios de parada flexíveis

    Args:
        strings (list[str]): Lista de strings de entrada.
        alphabet (str): Alfabeto utilizado.
        **params: Parâmetros do algoritmo.

    Métodos:
        run(): Executa o BLF-GA e retorna (centro, distância máxima, metadata).
        run_with_history(): Executa e retorna histórico completo.
        get_metadata(): Retorna metadados detalhados do algoritmo.
    """

    name = "BLF-GA"
    default_params = BLF_GA_DEFAULTS
    supports_internal_parallel = True  # BLF-GA pode usar paralelismo interno
    is_deterministic = False  # É estocástico (pode ter seed para reprodutibilidade)

    def __init__(self, strings: list[str], alphabet: str, **params):
        super().__init__(strings, alphabet, **params)

        # Filtrar parâmetros para passar apenas os que a implementação BLFGA conhece
        # Remover parâmetros específicos do framework de histórico
        blfga_params = {
            k: v
            for k, v in self.params.items()
            if k not in ["save_history", "history_frequency"]
        }

        self.blf_ga_instance = BLFGA(self.strings, self.alphabet, **blfga_params)

        # Configurar callback de histórico se habilitado
        if self.save_history:
            self.blf_ga_instance.set_history_callback(self._save_dynamic_history_entry)

    def set_params(self, **params) -> None:
        """
        Define novos parâmetros para o algoritmo.

        Args:
            **params: Parâmetros a serem atualizados
        """
        super().set_params(**params)
        # Atualizar instância do BLF-GA diretamente
        self.blf_ga_instance.update_params(**params)

    def set_progress_callback(self, callback: Callable[[str], None]) -> None:
        """
        Passa o callback para a instância do BLFGA.

        Args:
            callback (Callable[[str], None]): Função de callback de progresso.
        """
        super().set_progress_callback(callback)
        self.blf_ga_instance.set_progress_callback(callback)

    def set_warning_callback(self, callback: Callable[[str], None]) -> None:
        """
        Define callback para warnings e passa para a instância do BLFGA.

        Args:
            callback (Callable[[str], None]): Função de callback de warning.
        """
        super().set_warning_callback(callback)
        # BLF-GA pode gerar warnings sobre convergência, parâmetros, etc.
        # Por enquanto, apenas armazena o callback
        self._blfga_warning_callback = callback

    def run(self) -> tuple[str, int, dict]:
        """
        Executa o BLF-GA e retorna a string central, a distância máxima e metadata detalhada.

        Returns:
            tuple: (best_string, best_fitness, metadata)
                - best_string: Melhor string encontrada
                - best_fitness: Distância máxima da melhor string
                - metadata: Dicionário com informações detalhadas da execução
        """
        import time

        start_time = time.time()

        try:
            # Limpar histórico anterior
            self.clear_history()

            self._report_progress("Iniciando BLF-GA...")

            # Salvar estado inicial
            self._save_history_entry(
                0,
                phase="initialization",
                population_size=self.blf_ga_instance.pop_size,
                max_generations=self.blf_ga_instance.max_gens,
                initial_fitness=None,
            )

            best, best_val, ga_history = self.blf_ga_instance.run()
            execution_time = time.time() - start_time

            # Processar histórico do GA para nosso formato padrão
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

            # Salvar estado final
            self._save_history_entry(
                len(ga_history) if ga_history else 1,
                phase="completion",
                final_solution=best,
                final_fitness=best_val,
                total_generations=len(ga_history) if ga_history else 0,
                execution_time=execution_time,
            )

            # Metadata detalhada da execução
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
                # Histórico detalhado se habilitado
                "history": self.get_history() if self.save_history else [],
                "ga_fitness_history": ga_history if ga_history else [],
                # Compatibilidade com versão anterior
                "iteracoes": len(ga_history) if ga_history else 0,
                "melhor_distancia": best_val,
                "historico_completo": len(ga_history) if ga_history else 0,
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
        """Callback para salvar histórico durante execução do BLF-GA."""
        self._save_history_entry(
            generation,
            phase="evolution",
            generation=generation,
            best_fitness=best_fitness,
            **kwargs,
        )

    def run_with_history(self) -> tuple[str, int, list]:
        """
        Executa o BLF-GA e retorna a string central, a distância máxima e o histórico de distâncias.
        """
        return self.blf_ga_instance.run()

    def get_metadata(self) -> dict:
        """
        Retorna metadados detalhados do algoritmo BLF-GA.

        Returns:
            dict: Metadados incluindo configuração e estado atual
        """
        base_metadata = super().get_metadata()

        # Adiciona metadados específicos do BLF-GA
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

        # Merge dos metadados
        base_metadata.update(blfga_metadata)
        return base_metadata

    def _save_dynamic_history_entry(self, generation: int, event_data: dict) -> None:
        """
        Callback para registrar eventos dinâmicos durante execução do BLF-GA.

        Args:
            generation: Geração atual
            event_data: Dicionário com dados do evento incluindo 'event' e demais informações
        """
        if self.save_history:
            # Extrair o tipo de evento e demais dados
            event_type = event_data.pop("event", "unknown_event")
            self._save_history_entry(
                generation + 1,  # +1 porque iteration 0 é initialization
                phase="dynamic_event",
                generation=generation,
                event_type=event_type,
                **event_data,
            )
