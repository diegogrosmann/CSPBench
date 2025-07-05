"""
Sistema de relatórios para execução em lote.

Classes:
    BatchReporter: Gerador de relatórios detalhados
    ReportFormat: Enum com formatos de relatório
    ReportSection: Enum com seções de relatório

Funcionalidades:
    - Geração de relatórios em múltiplos formatos
    - Análise estatística detalhada
    - Comparação de algoritmos
    - Visualizações (quando disponível)
    - Exportação para CSV, JSON, HTML
"""

import json
import logging
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any

import pandas as pd

from src.core.io.results_formatter import ResultsFormatter

logger = logging.getLogger(__name__)


class ReportFormat(Enum):
    """Enum para formatos de relatório."""

    CSV = "csv"
    JSON = "json"
    HTML = "html"
    TXT = "txt"


class ReportSection(Enum):
    """Enum para seções de relatório."""

    SUMMARY = "summary"
    DETAILED = "detailed"
    COMPARISON = "comparison"
    STATISTICS = "statistics"
    CONFIGURATION = "configuration"


class BatchReporter:
    """
    Gerador de relatórios detalhados para execução em lote.

    Analisa resultados de execução em lote e gera relatórios
    em diferentes formatos com estatísticas e comparações.
    """

    def __init__(self, output_dir: str | Path | None = None):
        """
        Inicializa o gerador de relatórios.

        Args:
            output_dir: Diretório de saída para relatórios
        """
        self.output_dir = Path(output_dir) if output_dir else Path("reports")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.results_formatter = ResultsFormatter()

    def generate_batch_report(
        self,
        batch_results: dict[str, Any],
        config: dict[str, Any],
        formats: list[ReportFormat] | None = None,
        sections: list[ReportSection] | None = None,
    ) -> dict[str, Path]:
        """
        Gera relatório completo de execução em lote.

        Args:
            batch_results: Resultados da execução em lote
            config: Configuração utilizada
            formats: Formatos de relatório desejados
            sections: Seções a incluir no relatório

        Returns:
            Dicionário com paths dos arquivos gerados
        """
        if formats is None:
            formats = [ReportFormat.CSV, ReportFormat.JSON, ReportFormat.TXT]

        if sections is None:
            sections = list(ReportSection)

        logger.info(f"Gerando relatório em lote para: {config.get('nome', 'sem nome')}")

        # Análise dos resultados
        analysis = self._analyze_batch_results(batch_results, config)

        # Geração dos relatórios
        generated_files = {}

        for fmt in formats:
            try:
                file_path = self._generate_report_format(analysis, fmt, sections)
                generated_files[fmt.value] = file_path
                logger.info(f"Relatório {fmt.value.upper()} gerado: {file_path}")
            except Exception as e:
                logger.error(f"Erro ao gerar relatório {fmt.value}: {e}")

        return generated_files

    def _analyze_batch_results(self, batch_results: dict[str, Any], config: dict[str, Any]) -> dict[str, Any]:
        """
        Analisa resultados de execução em lote.

        Args:
            batch_results: Resultados da execução
            config: Configuração utilizada

        Returns:
            Análise estruturada dos resultados
        """
        analysis = {
            "config": config,
            "summary": {},
            "algorithm_stats": {},
            "base_stats": {},
            "execution_details": [],
            "comparison": {},
            "metadata": {
                "generated_at": datetime.now().isoformat(),
                "total_executions": 0,
                "successful_executions": 0,
                "failed_executions": 0,
            },
        }

        # Análise por algoritmo
        for alg_name, alg_results in batch_results.items():
            if not isinstance(alg_results, dict):
                continue

            alg_stats = self._analyze_algorithm_results(alg_results)
            analysis["algorithm_stats"][alg_name] = alg_stats

            # Atualizar contadores globais
            analysis["metadata"]["total_executions"] += alg_stats["total_executions"]
            analysis["metadata"]["successful_executions"] += alg_stats["successful_executions"]
            analysis["metadata"]["failed_executions"] += alg_stats["failed_executions"]

        # Análise comparativa
        analysis["comparison"] = self._generate_comparison_analysis(analysis["algorithm_stats"])

        # Resumo geral
        analysis["summary"] = self._generate_summary(analysis)

        return analysis

    def _analyze_algorithm_results(self, alg_results: dict[str, Any]) -> dict[str, Any]:
        """
        Analisa resultados de um algoritmo específico.

        Args:
            alg_results: Resultados do algoritmo

        Returns:
            Estatísticas do algoritmo
        """
        stats = {
            "total_executions": 0,
            "successful_executions": 0,
            "failed_executions": 0,
            "distances": [],
            "times": [],
            "iterations": [],
            "success_rate": 0.0,
            "avg_distance": 0.0,
            "avg_time": 0.0,
            "avg_iterations": 0.0,
            "min_distance": float("inf"),
            "max_distance": 0.0,
            "bases": {},
        }

        # Processar execuções por base
        for base_key, base_results in alg_results.items():
            if not isinstance(base_results, dict) or "execucoes" not in base_results:
                continue

            base_stats = {
                "executions": [],
                "success_count": 0,
                "avg_distance": 0.0,
                "avg_time": 0.0,
            }

            executions = base_results["execucoes"]
            for exec_result in executions:
                stats["total_executions"] += 1
                base_stats["executions"].append(exec_result)

                if exec_result.get("sucesso", False):
                    stats["successful_executions"] += 1
                    base_stats["success_count"] += 1

                    distance = exec_result.get("distancia", float("inf"))
                    time_taken = exec_result.get("tempo", 0.0)
                    iterations = exec_result.get("iteracoes", 0)

                    stats["distances"].append(distance)
                    stats["times"].append(time_taken)
                    stats["iterations"].append(iterations)

                    stats["min_distance"] = min(stats["min_distance"], distance)
                    stats["max_distance"] = max(stats["max_distance"], distance)
                else:
                    stats["failed_executions"] += 1

            # Calcular estatísticas da base
            if base_stats["success_count"] > 0:
                successful_executions = [e for e in base_stats["executions"] if e.get("sucesso", False)]
                base_stats["avg_distance"] = sum(e.get("distancia", 0) for e in successful_executions) / len(
                    successful_executions
                )
                base_stats["avg_time"] = sum(e.get("tempo", 0) for e in successful_executions) / len(
                    successful_executions
                )

            stats["bases"][base_key] = base_stats

        # Calcular estatísticas finais
        if stats["successful_executions"] > 0:
            stats["success_rate"] = stats["successful_executions"] / stats["total_executions"]
            stats["avg_distance"] = sum(stats["distances"]) / len(stats["distances"])
            stats["avg_time"] = sum(stats["times"]) / len(stats["times"])
            stats["avg_iterations"] = sum(stats["iterations"]) / len(stats["iterations"])

        return stats

    def _generate_comparison_analysis(self, algorithm_stats: dict[str, Any]) -> dict[str, Any]:
        """
        Gera análise comparativa entre algoritmos.

        Args:
            algorithm_stats: Estatísticas por algoritmo

        Returns:
            Análise comparativa
        """
        if not algorithm_stats:
            return {}

        comparison = {
            "best_algorithm": None,
            "best_distance": float("inf"),
            "fastest_algorithm": None,
            "fastest_time": float("inf"),
            "most_reliable": None,
            "highest_success_rate": 0.0,
            "ranking": [],
        }

        # Encontrar melhores em cada categoria
        for alg_name, stats in algorithm_stats.items():
            if stats["avg_distance"] < comparison["best_distance"]:
                comparison["best_distance"] = stats["avg_distance"]
                comparison["best_algorithm"] = alg_name

            if stats["avg_time"] < comparison["fastest_time"]:
                comparison["fastest_time"] = stats["avg_time"]
                comparison["fastest_algorithm"] = alg_name

            if stats["success_rate"] > comparison["highest_success_rate"]:
                comparison["highest_success_rate"] = stats["success_rate"]
                comparison["most_reliable"] = alg_name

        # Criar ranking ponderado
        ranking = []
        for alg_name, stats in algorithm_stats.items():
            score = self._calculate_algorithm_score(stats)
            ranking.append(
                {
                    "algorithm": alg_name,
                    "score": score,
                    "distance": stats["avg_distance"],
                    "time": stats["avg_time"],
                    "success_rate": stats["success_rate"],
                }
            )

        ranking.sort(key=lambda x: x["score"], reverse=True)
        comparison["ranking"] = ranking

        return comparison

    def _calculate_algorithm_score(self, stats: dict[str, Any]) -> float:
        """
        Calcula score ponderado para um algoritmo.

        Args:
            stats: Estatísticas do algoritmo

        Returns:
            Score ponderado
        """
        if stats["total_executions"] == 0:
            return 0.0

        # Fatores de ponderação
        success_weight = 0.4
        distance_weight = 0.4
        time_weight = 0.2

        # Normalizar valores (menores são melhores)
        success_score = stats["success_rate"]
        distance_score = 1.0 / (1.0 + stats["avg_distance"]) if stats["avg_distance"] > 0 else 0.0
        time_score = 1.0 / (1.0 + stats["avg_time"]) if stats["avg_time"] > 0 else 0.0

        return success_weight * success_score + distance_weight * distance_score + time_weight * time_score

    def _generate_summary(self, analysis: dict[str, Any]) -> dict[str, Any]:
        """
        Gera resumo geral da execução.

        Args:
            analysis: Análise completa

        Returns:
            Resumo estruturado
        """
        metadata = analysis["metadata"]
        comparison = analysis["comparison"]

        return {
            "total_algorithms": len(analysis["algorithm_stats"]),
            "total_executions": metadata["total_executions"],
            "success_rate": (
                metadata["successful_executions"] / metadata["total_executions"]
                if metadata["total_executions"] > 0
                else 0.0
            ),
            "best_algorithm": comparison.get("best_algorithm"),
            "best_distance": comparison.get("best_distance"),
            "fastest_algorithm": comparison.get("fastest_algorithm"),
            "execution_time": metadata["generated_at"],
            "configuration": analysis["config"].get("nome", "sem nome"),
        }

    def _generate_report_format(
        self, analysis: dict[str, Any], fmt: ReportFormat, sections: list[ReportSection]
    ) -> Path:
        """
        Gera relatório em formato específico.

        Args:
            analysis: Análise dos resultados
            fmt: Formato do relatório
            sections: Seções a incluir

        Returns:
            Path do arquivo gerado
        """
        config_name = analysis["config"].get("nome", "batch_report")
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"{config_name}_{timestamp}.{fmt.value}"
        output_path = self.output_dir / filename

        if fmt == ReportFormat.CSV:
            self._generate_csv_report(analysis, output_path, sections)
        elif fmt == ReportFormat.JSON:
            self._generate_json_report(analysis, output_path, sections)
        elif fmt == ReportFormat.HTML:
            self._generate_html_report(analysis, output_path, sections)
        elif fmt == ReportFormat.TXT:
            self._generate_txt_report(analysis, output_path, sections)

        return output_path

    def _generate_csv_report(self, analysis: dict[str, Any], output_path: Path, sections: list[ReportSection]):
        """Gera relatório em CSV."""
        # Criar DataFrame com resultados detalhados
        rows = []
        for alg_name, stats in analysis["algorithm_stats"].items():
            for base_name, base_stats in stats["bases"].items():
                for i, execution in enumerate(base_stats["executions"]):
                    row = {
                        "algorithm": alg_name,
                        "base": base_name,
                        "execution": i + 1,
                        "success": execution.get("sucesso", False),
                        "distance": execution.get("distancia", float("inf")),
                        "time": execution.get("tempo", 0.0),
                        "iterations": execution.get("iteracoes", 0),
                        "error": execution.get("erro", ""),
                    }
                    rows.append(row)

        df = pd.DataFrame(rows)
        df.to_csv(output_path, index=False)

    def _generate_json_report(self, analysis: dict[str, Any], output_path: Path, sections: list[ReportSection]):
        """Gera relatório em JSON."""
        # Filtrar seções
        filtered_analysis = {}
        for section in sections:
            if section.value in analysis:
                filtered_analysis[section.value] = analysis[section.value]

        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(filtered_analysis, f, indent=2, ensure_ascii=False, default=str)

    def _generate_html_report(self, analysis: dict[str, Any], output_path: Path, sections: list[ReportSection]):
        """Gera relatório em HTML."""
        html_content = self._create_html_template(analysis, sections)
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(html_content)

    def _generate_txt_report(self, analysis: dict[str, Any], output_path: Path, sections: list[ReportSection]):
        """Gera relatório em texto."""
        with open(output_path, "w", encoding="utf-8") as f:
            f.write("RELATÓRIO DE EXECUÇÃO EM LOTE\n")
            f.write("=" * 50 + "\n\n")

            # Resumo
            if ReportSection.SUMMARY in sections:
                f.write("RESUMO GERAL\n")
                f.write("-" * 20 + "\n")
                summary = analysis["summary"]
                f.write(f"Configuração: {summary['configuration']}\n")
                f.write(f"Algoritmos testados: {summary['total_algorithms']}\n")
                f.write(f"Execuções totais: {summary['total_executions']}\n")
                f.write(f"Taxa de sucesso: {summary['success_rate']:.2%}\n")
                f.write(f"Melhor algoritmo: {summary['best_algorithm']}\n")
                f.write(f"Melhor distância: {summary['best_distance']}\n\n")

            # Estatísticas por algoritmo
            if ReportSection.STATISTICS in sections:
                f.write("ESTATÍSTICAS POR ALGORITMO\n")
                f.write("-" * 30 + "\n")
                for alg_name, stats in analysis["algorithm_stats"].items():
                    f.write(f"\n{alg_name}:\n")
                    f.write(f"  Execuções: {stats['total_executions']}\n")
                    f.write(f"  Sucessos: {stats['successful_executions']}\n")
                    f.write(f"  Taxa de sucesso: {stats['success_rate']:.2%}\n")
                    f.write(f"  Distância média: {stats['avg_distance']:.3f}\n")
                    f.write(f"  Tempo médio: {stats['avg_time']:.3f}s\n")

            # Comparação
            if ReportSection.COMPARISON in sections:
                f.write("\nCOMPARAÇÃO DE ALGORITMOS\n")
                f.write("-" * 25 + "\n")
                comparison = analysis["comparison"]
                f.write(f"Melhor distância: {comparison['best_algorithm']} ({comparison['best_distance']:.3f})\n")
                f.write(f"Mais rápido: {comparison['fastest_algorithm']} ({comparison['fastest_time']:.3f}s)\n")
                f.write(f"Mais confiável: {comparison['most_reliable']} ({comparison['highest_success_rate']:.2%})\n")

                f.write("\nRanking:\n")
                for i, entry in enumerate(comparison["ranking"], 1):
                    f.write(f"  {i}. {entry['algorithm']} (score: {entry['score']:.3f})\n")

    def _create_html_template(self, analysis: dict[str, Any], sections: list[ReportSection]) -> str:
        """Cria template HTML para relatório."""
        html = """
        <!DOCTYPE html>
        <html>
        <head>
            <title>Relatório CSP-BLFGA</title>
            <style>
                body { font-family: Arial, sans-serif; margin: 20px; }
                h1, h2 { color: #333; }
                table { border-collapse: collapse; width: 100%; margin: 20px 0; }
                th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
                th { background-color: #f2f2f2; }
                .summary { background-color: #f9f9f9; padding: 15px; border-radius: 5px; }
            </style>
        </head>
        <body>
            <h1>Relatório de Execução em Lote</h1>
        """

        if ReportSection.SUMMARY in sections:
            summary = analysis["summary"]
            html += f"""
            <div class="summary">
                <h2>Resumo</h2>
                <p><strong>Configuração:</strong> {summary['configuration']}</p>
                <p><strong>Algoritmos testados:</strong> {summary['total_algorithms']}</p>
                <p><strong>Execuções totais:</strong> {summary['total_executions']}</p>
                <p><strong>Taxa de sucesso:</strong> {summary['success_rate']:.2%}</p>
                <p><strong>Melhor algoritmo:</strong> {summary['best_algorithm']}</p>
            </div>
            """

        if ReportSection.STATISTICS in sections:
            html += "<h2>Estatísticas por Algoritmo</h2><table>"
            html += "<tr><th>Algoritmo</th><th>Execuções</th><th>Sucessos</th><th>Taxa</th><th>Dist. Média</th><th>Tempo Médio</th></tr>"
            for alg_name, stats in analysis["algorithm_stats"].items():
                html += f"""
                <tr>
                    <td>{alg_name}</td>
                    <td>{stats['total_executions']}</td>
                    <td>{stats['successful_executions']}</td>
                    <td>{stats['success_rate']:.2%}</td>
                    <td>{stats['avg_distance']:.3f}</td>
                    <td>{stats['avg_time']:.3f}s</td>
                </tr>
                """
            html += "</table>"

        html += "</body></html>"
        return html
