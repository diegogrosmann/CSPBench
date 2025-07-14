"""
Gerador de Relatórios de Análise de Sensibilidade

Responsável por gerar relatórios específicos para análises de sensibilidade,
incluindo visualizações dos índices de Morris, Sobol, etc.
"""

import json
from pathlib import Path
from typing import Any, Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


class SensitivityReportGenerator:
    """
    Gerador de relatórios específicos para análises de sensibilidade.
    """

    def __init__(self, session_path: Path):
        """
        Inicializa o gerador de relatórios de sensibilidade.

        Args:
            session_path: Caminho da sessão onde salvar o relatório
        """
        self.session_path = session_path

        # Configurar matplotlib
        plt.style.use("seaborn-v0_8")
        sns.set_palette("husl")

    def generate_sensitivity_report(self, sensitivity_data: Dict[str, Any]) -> None:
        """
        Gera relatório específico para análise de sensibilidade.

        Args:
            sensitivity_data: Dados da análise de sensibilidade
        """
        # Criar diretório de relatório
        report_dir = self.session_path / "sensitivity_report"
        report_dir.mkdir(exist_ok=True)

        # Processar dados de sensibilidade
        analysis_results = self._extract_sensitivity_results(sensitivity_data)

        if not analysis_results:
            print("⚠️ Nenhum dado de sensibilidade encontrado para relatório")
            return

        # Gerar gráficos de sensibilidade
        self._generate_sensitivity_plots(analysis_results, report_dir)

        # Gerar relatório HTML específico
        self._generate_sensitivity_html_report(analysis_results, report_dir)

        # Exportar dados estruturados
        self._export_sensitivity_csv(analysis_results, report_dir)

        print(f"📊 Relatório de sensibilidade gerado em: {report_dir}")

    def _extract_sensitivity_results(
        self, data: Dict[str, Any]
    ) -> List[Dict[str, Any]]:
        """Extrai resultados de sensibilidade da estrutura de dados."""
        results = []

        # Navegar pela estrutura aninhada
        batch_results = data.get("batch_results", [])

        for batch_result in batch_results:
            # Verificar batch_summary
            batch_summary = batch_result.get("batch_summary", {})
            summary_results = batch_summary.get("results", [])

            for result in summary_results:
                if "sensitivity_indices" in result:
                    results.append(result)

            # Verificar detailed_results
            detailed_results = batch_result.get("detailed_results", [])
            for result in detailed_results:
                if "sensitivity_indices" in result:
                    results.append(result)

        return results

    def _generate_sensitivity_plots(
        self, results: List[Dict[str, Any]], report_dir: Path
    ):
        """Gera gráficos específicos para análise de sensibilidade."""
        plots_dir = report_dir / "plots"
        plots_dir.mkdir(exist_ok=True)

        for i, result in enumerate(results):
            algorithm = result.get("algorithm", f"Algorithm_{i}")
            method = result.get("analysis_method", "unknown")
            sensitivity_indices = result.get("sensitivity_indices", {})

            # Gerar gráfico para cada métrica
            for metric, indices in sensitivity_indices.items():
                if not indices or "parameter_names" not in indices:
                    continue

                self._plot_morris_indices(indices, algorithm, metric, method, plots_dir)

    def _plot_morris_indices(
        self,
        indices: Dict[str, Any],
        algorithm: str,
        metric: str,
        method: str,
        plots_dir: Path,
    ):
        """Gera gráfico dos índices de Morris."""
        if method.lower() != "morris":
            return

        parameter_names = indices.get("parameter_names", [])
        mu_star = indices.get("mu_star", [])
        sigma = indices.get("sigma", [])

        if not parameter_names or not mu_star or not sigma:
            return

        # Criar figura
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

        # Gráfico de barras μ*
        ax1.bar(parameter_names, mu_star, color="skyblue", alpha=0.7)
        ax1.set_title(f"{algorithm} - Morris μ* ({metric})")
        ax1.set_ylabel("μ* (Influência)")
        ax1.tick_params(axis="x", rotation=45)
        ax1.grid(True, alpha=0.3)

        # Gráfico de dispersão μ* vs σ
        ax2.scatter(mu_star, sigma, s=100, alpha=0.7, color="coral")
        for i, param in enumerate(parameter_names):
            ax2.annotate(
                param,
                (mu_star[i], sigma[i]),
                xytext=(5, 5),
                textcoords="offset points",
                fontsize=9,
            )
        ax2.set_xlabel("μ* (Efeito Principal)")
        ax2.set_ylabel("σ (Interações/Não-linearidade)")
        ax2.set_title(f"{algorithm} - Morris μ* vs σ ({metric})")
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()

        # Salvar gráfico
        filename = f"morris_{algorithm}_{metric}.png".replace(" ", "_").lower()
        plt.savefig(plots_dir / filename, dpi=300, bbox_inches="tight")
        plt.close()

    def _generate_sensitivity_html_report(
        self, results: List[Dict[str, Any]], report_dir: Path
    ):
        """Gera relatório HTML específico para sensibilidade."""
        html_content = self._build_sensitivity_html(results)

        html_path = report_dir / "sensitivity_report.html"
        with open(html_path, "w", encoding="utf-8") as f:
            f.write(html_content)

        print(f"📋 Relatório HTML de sensibilidade: {html_path}")

    def _build_sensitivity_html(self, results: List[Dict[str, Any]]) -> str:
        """Constrói conteúdo HTML para relatório de sensibilidade."""
        html = """
<!DOCTYPE html>
<html>
<head>
    <title>CSPBench - Relatório de Análise de Sensibilidade</title>
    <meta charset="utf-8">
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; background: #f9f9f9; }
        .header { background: #e3f2fd; padding: 20px; border-radius: 8px; border-left: 5px solid #2196f3; }
        .section { margin: 25px 0; background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
        .analysis-table { border-collapse: collapse; width: 100%; margin: 15px 0; }
        .analysis-table th, .analysis-table td { border: 1px solid #ddd; padding: 12px; text-align: left; }
        .analysis-table th { background-color: #f5f5f5; font-weight: bold; }
        .plot { margin: 20px 0; text-align: center; }
        .plot img { max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 4px; }
        .metric-section { margin: 15px 0; padding: 15px; border-left: 3px solid #4caf50; background: #f1f8e9; }
        .parameter-importance { display: inline-block; margin: 5px; padding: 8px 12px; background: #e8eaf6; border-radius: 4px; }
        .high-importance { background: #ffcdd2; }
        .medium-importance { background: #fff3e0; }
        .low-importance { background: #e8f5e8; }
    </style>
</head>
<body>
    <div class="header">
        <h1>🔬 CSPBench - Análise de Sensibilidade</h1>
        <p>Relatório detalhado da análise de sensibilidade de parâmetros dos algoritmos</p>
    </div>
    
    <div class="section">
        <h2>📊 Resumo Geral</h2>
        <table class="analysis-table">
            <tr><th>Métrica</th><th>Valor</th></tr>
"""

        # Adicionar estatísticas gerais
        total_analyses = len(results)
        algorithms = list(set(r.get("algorithm", "Unknown") for r in results))
        methods = list(set(r.get("analysis_method", "Unknown") for r in results))

        html += f"""
            <tr><td>Total de Análises</td><td>{total_analyses}</td></tr>
            <tr><td>Algoritmos Analisados</td><td>{', '.join(algorithms)}</td></tr>
            <tr><td>Métodos Utilizados</td><td>{', '.join(methods)}</td></tr>
        </table>
    </div>
"""

        # Adicionar seção para cada análise
        for i, result in enumerate(results):
            algorithm = result.get("algorithm", f"Algorithm_{i}")
            method = result.get("analysis_method", "unknown")
            n_samples = result.get("n_samples", 0)
            parameters = result.get("parameters_analyzed", [])
            sensitivity_indices = result.get("sensitivity_indices", {})

            html += f"""
    <div class="section">
        <h2>🧠 {algorithm} - Método {method.upper()}</h2>
        <table class="analysis-table">
            <tr><td>Amostras Geradas</td><td>{n_samples}</td></tr>
            <tr><td>Parâmetros Analisados</td><td>{len(parameters)}</td></tr>
            <tr><td>Métricas Avaliadas</td><td>{len(sensitivity_indices)}</td></tr>
        </table>
        
        <h3>📈 Parâmetros Analisados</h3>
        <p>
"""
            for param in parameters:
                html += f'<span class="parameter-importance">{param}</span>'

            html += "</p>"

            # Adicionar resultados por métrica
            for metric, indices in sensitivity_indices.items():
                if not indices or "parameter_names" not in indices:
                    continue

                html += f"""
        <div class="metric-section">
            <h4>📊 Métrica: {metric}</h4>
"""

                if method.lower() == "morris" and "mu_star" in indices:
                    mu_star = indices.get("mu_star", [])
                    parameter_names = indices.get("parameter_names", [])

                    html += """
            <table class="analysis-table">
                <tr><th>Parâmetro</th><th>μ* (Influência)</th><th>Classificação</th></tr>
"""
                    for j, (param, influence) in enumerate(
                        zip(parameter_names, mu_star)
                    ):
                        if influence > 1.0:
                            classification = "Alta"
                            css_class = "high-importance"
                        elif influence > 0.5:
                            classification = "Média"
                            css_class = "medium-importance"
                        else:
                            classification = "Baixa"
                            css_class = "low-importance"

                        html += f"""
                <tr class="{css_class}">
                    <td>{param}</td>
                    <td>{influence:.4f}</td>
                    <td>{classification}</td>
                </tr>
"""
                    html += "</table>"

                # Adicionar gráfico se existe
                plot_filename = f"morris_{algorithm}_{metric}.png".replace(
                    " ", "_"
                ).lower()
                html += f"""
            <div class="plot">
                <img src="plots/{plot_filename}" alt="Gráfico {metric}" onerror="this.style.display='none'">
            </div>
"""

                html += "</div>"

            html += "</div>"

        html += """
</body>
</html>
"""
        return html

    def _export_sensitivity_csv(self, results: List[Dict[str, Any]], report_dir: Path):
        """Exporta dados de sensibilidade em formato CSV."""
        csv_data = []

        for result in results:
            algorithm = result.get("algorithm", "Unknown")
            method = result.get("analysis_method", "unknown")
            n_samples = result.get("n_samples", 0)
            sensitivity_indices = result.get("sensitivity_indices", {})

            for metric, indices in sensitivity_indices.items():
                if not indices or "parameter_names" not in indices:
                    continue

                parameter_names = indices.get("parameter_names", [])

                if method.lower() == "morris" and "mu_star" in indices:
                    mu_star = indices.get("mu_star", [])
                    sigma = indices.get("sigma", [])

                    for i, param in enumerate(parameter_names):
                        csv_data.append(
                            {
                                "algorithm": algorithm,
                                "method": method,
                                "metric": metric,
                                "parameter": param,
                                "n_samples": n_samples,
                                "mu_star": mu_star[i] if i < len(mu_star) else 0,
                                "sigma": sigma[i] if i < len(sigma) else 0,
                            }
                        )

        if csv_data:
            df = pd.DataFrame(csv_data)
            csv_path = report_dir / "sensitivity_indices.csv"
            df.to_csv(csv_path, index=False)
            print(f"📄 Dados CSV exportados: {csv_path}")
