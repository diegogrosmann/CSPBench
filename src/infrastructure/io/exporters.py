"""
Exporters de dados

Implementa interfaces para exportar resultados em diferentes formatos.
"""

import csv
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from src.application.ports import ExportPort
from src.infrastructure.io.report_generator import AdvancedReportGenerator


class FileExporter(ExportPort):
    """Exportador base para arquivos."""

    def __init__(self, output_path: str):
        """Inicializa exportador com caminho de saída."""
        self.output_path = Path(output_path)
        self.output_path.mkdir(parents=True, exist_ok=True)

    def export_results(
        self, results: Dict[str, Any], format_type: str, destination: str
    ) -> str:
        """Exporta resultados em formato específico."""
        dest_path = self.output_path / destination

        # Se destination é um diretório, gerar nome de arquivo
        if dest_path.is_dir() or not dest_path.suffix:
            from datetime import datetime

            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"results_{timestamp}.{format_type.lower()}"
            dest_path = dest_path / filename

        dest_path.parent.mkdir(parents=True, exist_ok=True)

        if format_type.lower() == "json":
            self._write_json(results, dest_path)
        elif format_type.lower() == "csv":
            self._write_csv(results, dest_path)
        elif format_type.lower() == "txt":
            self._write_txt(results, dest_path)
        else:
            # Default para JSON
            self._write_json(results, dest_path)

        return str(dest_path)

    def export_batch_results(
        self, batch_results: List[Dict[str, Any]], format_type: str, destination: str
    ) -> str:
        """Exporta resultados de batch."""
        dest_path = self.output_path / destination

        # Se destination é um diretório, usar nome de arquivo fixo
        if dest_path.is_dir() or not dest_path.suffix:
            # Usar nome fixo para resultado principal
            filename = f"results.{format_type.lower()}"
            dest_path = dest_path / filename

        dest_path.parent.mkdir(parents=True, exist_ok=True)

        batch_data = {
            "batch_results": batch_results,
            "summary": {
                "total_experiments": len(batch_results),
                "successful": sum(
                    1 for r in batch_results if r.get("status") == "success"
                ),
                "failed": sum(1 for r in batch_results if r.get("status") == "error"),
            },
        }

        if format_type.lower() == "json":
            self._write_json(batch_data, dest_path)
        elif format_type.lower() == "csv":
            self._write_csv(batch_results, dest_path)  # Para CSV, só os resultados
        elif format_type.lower() == "txt":
            self._write_txt(batch_data, dest_path)
        else:
            self._write_json(batch_data, dest_path)

        return str(dest_path)

    def get_supported_formats(self) -> List[str]:
        """Lista formatos suportados."""
        return ["json", "csv", "txt"]

    def export_optimization_results(
        self, optimization_data: Dict[str, Any], destination: str
    ) -> str:
        """Exporta resultados de otimização."""
        return self.export_results(optimization_data, "json", destination)

    def export(self, data: Any, filename: Optional[str] = None) -> None:
        """Exporta dados para arquivo (método legacy)."""
        if filename:
            path = self.output_path / filename
        else:
            path = self.output_path

        # Detecta formato pelo arquivo
        if path.suffix.lower() == ".json":
            self._write_json(data, path)
        elif path.suffix.lower() == ".csv":
            self._write_csv(data, path)
        else:
            self._write_txt(data, path)

    def _write_json(self, data: Any, path: Path) -> None:
        """Escreve dados em formato JSON."""
        with open(path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, ensure_ascii=False)

    def _write_csv(self, data: Any, path: Path) -> None:
        """Escreve dados em formato CSV."""
        if isinstance(data, dict):
            # Converte dict para lista de dicts
            if all(isinstance(v, dict) for v in data.values()):
                # Dict de dicts -> lista de dicts com chave extra
                rows = []
                for key, value in data.items():
                    row = {"name": key}
                    row.update(value)
                    rows.append(row)
                data = rows
            else:
                # Dict simples -> uma linha
                data = [data]

        if isinstance(data, list) and data:
            with open(path, "w", newline="", encoding="utf-8") as f:
                if isinstance(data[0], dict):
                    fieldnames = data[0].keys()
                    writer = csv.DictWriter(f, fieldnames=fieldnames)
                    writer.writeheader()
                    writer.writerows(data)
                else:
                    writer = csv.writer(f)
                    writer.writerows(data)

    def _write_txt(self, data: Any, path: Path) -> None:
        """Escreve dados em formato texto."""
        with open(path, "w", encoding="utf-8") as f:
            if isinstance(data, (list, tuple)):
                for item in data:
                    f.write(f"{item}\n")
            elif isinstance(data, dict):
                for key, value in data.items():
                    f.write(f"{key}: {value}\n")
            else:
                f.write(str(data))


class CsvExporter(FileExporter):
    """Exportador para formato CSV."""

    pass


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
                    # Gerar relatório específico para análise de sensibilidade
                    try:
                        from .sensitivity_report_generator import (
                            SensitivityReportGenerator,
                        )

                        # Construir dados de sensibilidade
                        sensitivity_data = {"batch_results": batch_results}

                        session_path = Path(self.output_path)
                        sensitivity_report_generator = SensitivityReportGenerator(
                            session_path
                        )
                        sensitivity_report_generator.generate_sensitivity_report(
                            sensitivity_data
                        )

                    except ImportError as e:
                        print(
                            f"⚠️ Não foi possível importar gerador de relatórios de sensibilidade: {e}"
                        )
                    except Exception as e:
                        print(f"⚠️ Erro ao gerar relatório de sensibilidade: {e}")

                    return json_path

                elif is_optimization_data:
                    print(
                        "📊 Dados de otimização detectados - usando relatório específico de otimização"
                    )
                    print(
                        "📄 Os resultados detalhados estão disponíveis em:", json_path
                    )
                    # TODO: Implementar relatório específico de otimização
                    return json_path

                # Construir dados do batch para o relatório (apenas para execução normal)
                # Os batch_results já contêm a estrutura completa
                batch_data = {
                    "batch_results": batch_results,
                    "summary": {"total_experiments": len(batch_results)},
                }

                # Usar o output_path como session_path, que já é o diretório da sessão
                session_path = Path(self.output_path)
                report_generator = AdvancedReportGenerator(self.config, session_path)
                report_generator.generate_report(batch_data)
            except Exception as e:
                print(f"⚠️ Erro ao gerar relatório avançado: {e}")

        return json_path

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


class TxtExporter(FileExporter):
    """Exportador para formato texto simples."""

    pass
