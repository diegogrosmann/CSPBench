"""
Exportador para formato JSON com relatórios avançados.

Especialização do FileExporter que adiciona geração automática de relatórios
baseado no tipo de dados detectado.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional

from .file_exporter import FileExporter


class JsonExporter(FileExporter):
    """Exportador para formato JSON com relatórios avançados."""

    def __init__(self, output_path: str, config: Optional[Dict[str, Any]] = None):
        """Inicializa exportador JSON com configuração opcional."""
        super().__init__(output_path)
        self.config = config or {}

    def export_batch_results(
        self, batch_results: List[Dict[str, Any]], format_type: str, destination: str
    ) -> str:
        """Exporta resultados de batch com relatórios avançados."""
        # Exportar dados JSON básicos
        json_path = super().export_batch_results(
            batch_results, format_type, destination
        )

        # Gerar relatório avançado se configurado
        if self.config and self._should_generate_report():
            try:
                # Verificar tipos de dados específicos
                is_sensitivity_analysis = self._is_sensitivity_data(batch_results)
                is_optimization_data = self._is_optimization_data(batch_results)

                if is_sensitivity_analysis:
                    self._generate_sensitivity_report(batch_results)
                elif is_optimization_data:
                    self._generate_optimization_report(batch_results)
                else:
                    self._generate_execution_report(batch_results)

            except Exception as e:
                print(f"⚠️ Erro ao gerar relatório avançado: {e}")

        return json_path

    def _generate_sensitivity_report(self, batch_results: List[Dict[str, Any]]) -> None:
        """Gera relatório específico para análise de sensibilidade."""
        try:
            from ..report_generators.sensitivity_report_generator import (
                SensitivityReportGenerator,
            )

            # Construir dados de sensibilidade
            sensitivity_data = {"batch_results": batch_results}

            session_path = Path(self.output_path)
            sensitivity_report_generator = SensitivityReportGenerator(session_path)
            sensitivity_report_generator.generate_sensitivity_report(sensitivity_data)

        except ImportError as e:
            print(
                f"⚠️ Não foi possível importar gerador de relatórios de sensibilidade: {e}"
            )
        except Exception as e:
            print(f"⚠️ Erro ao gerar relatório de sensibilidade: {e}")

    def _generate_optimization_report(
        self, batch_results: List[Dict[str, Any]]
    ) -> None:
        """Gera relatório específico para otimização."""
        print(
            "📊 Dados de otimização detectados - usando relatório específico de otimização"
        )
        # TODO: Implementar gerador específico de otimização
        # try:
        #     from ..report_generators.optimization_report_generator import (
        #         OptimizationReportGenerator,
        #     )
        #
        #     optimization_data = {"batch_results": batch_results}
        #     session_path = Path(self.output_path)
        #     optimization_report_generator = OptimizationReportGenerator(session_path)
        #     optimization_report_generator.generate_optimization_report(optimization_data)
        # except ImportError as e:
        #     print(f"⚠️ Não foi possível importar gerador de relatórios de otimização: {e}")
        # except Exception as e:
        #     print(f"⚠️ Erro ao gerar relatório de otimização: {e}")

    def _generate_execution_report(self, batch_results: List[Dict[str, Any]]) -> None:
        """Gera relatório para execução normal."""
        try:
            from ..report_generators.execution_report_generator import (
                ExecutionReportGenerator,
            )

            # Construir dados do batch para o relatório
            batch_data = {
                "batch_results": batch_results,
                "summary": {"total_experiments": len(batch_results)},
            }

            session_path = Path(self.output_path)
            execution_report_generator = ExecutionReportGenerator(
                self.config, session_path
            )
            execution_report_generator.generate_report(batch_data)

        except ImportError as e:
            print(f"⚠️ Não foi possível importar gerador de relatórios de execução: {e}")
        except Exception as e:
            print(f"⚠️ Erro ao gerar relatório de execução: {e}")

    def _is_sensitivity_data(self, batch_results: List[Dict[str, Any]]) -> bool:
        """Verifica se os dados são de análise de sensibilidade."""
        if not batch_results:
            return False

        # Verificar estruturas aninhadas típicas de dados de sensibilidade
        for result in batch_results:
            # Verificar se tem estrutura de batch com dados de sensibilidade
            if "batch_summary" in result:
                summary_results = result.get("batch_summary", {}).get("results", [])
                for summary_result in summary_results:
                    if any(
                        key in summary_result
                        for key in [
                            "sensitivity_indices",
                            "analysis_method",
                            "parameters_analyzed",
                        ]
                    ):
                        return True

            # Verificar se tem dados de sensibilidade diretos
            if "detailed_results" in result:
                detailed_results = result.get("detailed_results", [])
                for detailed_result in detailed_results:
                    if any(
                        key in detailed_result
                        for key in [
                            "sensitivity_indices",
                            "analysis_method",
                            "parameters_analyzed",
                        ]
                    ):
                        return True

            # Verificar diretamente no resultado
            if any(
                key in result
                for key in [
                    "sensitivity_indices",
                    "analysis_method",
                    "parameters_analyzed",
                    "n_samples",
                ]
            ):
                return True

        return False

    def _is_optimization_data(self, batch_results: List[Dict[str, Any]]) -> bool:
        """Verifica se os dados são de otimização (Optuna)."""
        if not batch_results:
            return False

        # Verificar estruturas aninhadas típicas de dados de otimização
        for result in batch_results:
            # Verificar se tem estrutura de batch com dados de otimização
            if "batch_summary" in result:
                summary_results = result.get("batch_summary", {}).get("results", [])
                for summary_result in summary_results:
                    if any(
                        key in summary_result
                        for key in [
                            "study_name",
                            "n_trials",
                            "best_value",
                            "trials",
                            "direction",
                        ]
                    ):
                        return True

            # Verificar se tem dados de otimização diretos
            if "detailed_results" in result:
                detailed_results = result.get("detailed_results", [])
                for detailed_result in detailed_results:
                    if any(
                        key in detailed_result
                        for key in [
                            "study_name",
                            "n_trials",
                            "best_value",
                            "trials",
                            "direction",
                        ]
                    ):
                        return True

            # Verificar diretamente no resultado
            if any(
                key in result
                for key in [
                    "study_name",
                    "n_trials",
                    "best_value",
                    "trials",
                    "direction",
                ]
            ):
                return True

        return False

    def _should_generate_report(self) -> bool:
        """Verifica se deve gerar relatório baseado na configuração."""
        try:
            report_config = (
                self.config.get("infrastructure", {})
                .get("result", {})
                .get("report", {})
            )
            return report_config.get("save_detailed_results", False)
        except Exception:
            return False

    def get_supported_formats(self) -> list[str]:
        """Lista formatos suportados."""
        return ["json"]
