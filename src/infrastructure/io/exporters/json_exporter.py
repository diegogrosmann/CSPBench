"""
JSON Format Exporter with Advanced Reports.

FileExporter specialization that adds automatic report generation
based on detected data type.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional

from .file_exporter import FileExporter


class JsonExporter(FileExporter):
    """JSON format exporter with advanced reports."""

    def __init__(self, output_path: str, config: Optional[Dict[str, Any]] = None):
        """Initialize JSON exporter with optional configuration."""
        super().__init__(output_path)
        self.config = config or {}

    def export_batch_results(
        self, batch_results: List[Dict[str, Any]], format_type: str, destination: str
    ) -> str:
        """Export batch results with advanced reports."""
        # Export basic JSON data
        json_path = super().export_batch_results(
            batch_results, format_type, destination
        )

        # Generate advanced report if configured
        if self.config and self._should_generate_report():
            try:
                if self._is_sensitivity_data(batch_results):
                    print(
                        "📊 Sensitivity data detected - generating sensitivity report"
                    )
                    self._generate_sensitivity_report(batch_results)
                elif self._is_optimization_data(batch_results):
                    print(
                        "📊 Optimization data detected - generating optimization report"
                    )
                    self._generate_optimization_report(batch_results)
                else:
                    print("📊 Standard execution data - generating execution report")
                    self._generate_execution_report(batch_results)

            except Exception as e:
                print(f"⚠️ Error generating advanced report: {e}")

        return json_path

    def _generate_sensitivity_report(self, batch_results: List[Dict[str, Any]]) -> None:
        """Generate specific report for sensitivity analysis."""
        try:
            from ..report_generators.sensitivity_report_generator import (
                SensitivityReportGenerator,
            )

            # Build sensitivity data
            sensitivity_data = {"batch_results": batch_results}

            session_path = Path(self.output_path)
            sensitivity_report_generator = SensitivityReportGenerator(session_path)
            sensitivity_report_generator.generate_sensitivity_report(sensitivity_data)

        except ImportError as e:
            print(f"⚠️ Could not import sensitivity report generator: {e}")
        except Exception as e:
            print(f"⚠️ Error generating sensitivity report: {e}")

    def _generate_optimization_report(
        self, batch_results: List[Dict[str, Any]]
    ) -> None:
        """Generate specific report for optimization."""
        print("📊 Optimization data detected - using specific optimization report")

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
