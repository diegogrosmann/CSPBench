"""
Gerador de Relatórios Avançados

Responsável por gerar relatórios detalhados com gráficos, estatísticas e análises
dos resultados dos experimentos do CSPBench.
"""

import json
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from jinja2 import Environment, FileSystemLoader
from scipy import stats

from .history_plotter import HistoryPlotGenerator


class AdvancedReportGenerator:
    """
    Gerador de relatórios avançados com gráficos e estatísticas.
    """

    def __init__(self, config: Dict[str, Any], session_path: Path):
        """
        Inicializa o gerador de relatórios.

        Args:
            config: Configuração do sistema
            session_path: Caminho da sessão onde salvar o relatório
        """
        self.config = config
        self.session_path = session_path
        self.report_config = config["infrastructure"]["result"]["report"]
        self.export_config = config["infrastructure"]["result"]["export"]

        # Inicializar gerador de gráficos de histórico
        self.history_plotter = HistoryPlotGenerator(config, session_path)

        # Configurar matplotlib
        plt.style.use("seaborn-v0_8")
        sns.set_palette("husl")

    def generate_report(self, results_data: Dict[str, Any]) -> None:
        """
        Gera relatório completo com gráficos e estatísticas.

        Args:
            results_data: Dados dos resultados dos experimentos
        """
        # Criar diretório de relatório
        report_dir = self.session_path / "report"
        report_dir.mkdir(exist_ok=True)

        # Processar dados
        df_results = self._process_results_data(results_data)

        if df_results.empty:
            print("⚠️ Nenhum dado disponível para gerar relatório")
            return

        # Gerar gráficos de histórico primeiro
        self.history_plotter.generate_history_plots(results_data.get("results", []))

        # Gerar gráficos se habilitado
        if self.report_config.get("save_plots", True):
            self._generate_plots(df_results, report_dir)

        # Calcular estatísticas se habilitado
        statistics = {}
        if self.report_config.get("calculate_statistics", True):
            statistics = self._calculate_statistics(df_results)

        # Exportar dados se habilitado
        if self.export_config.get("export_csv", True):
            self._export_csv(df_results, report_dir)

        if self.export_config.get("export_json", True):
            self._export_json(results_data, statistics, report_dir)

        # Gerar relatório HTML/PDF
        if self.report_config.get("include_summary", True):
            self._generate_html_report(df_results, statistics, report_dir)

        print(f"📊 Relatório completo gerado em: {report_dir}")

    def _process_results_data(self, results_data: Dict[str, Any]) -> pd.DataFrame:
        """Converte dados de resultados em DataFrame do pandas."""
        try:
            if "results" not in results_data:
                return pd.DataFrame()

            # Extrair dados dos experimentos
            experiments = []
            for result in results_data["results"]:
                exp_data = {
                    "algorithm": result.get("algorithm", "Unknown"),
                    "dataset": result.get("dataset_id", "Unknown"),
                    "max_distance": result.get("max_distance", 0),
                    "runtime": result.get("execution_time", 0),
                    "success": result.get("status") == "success",
                    "best_string": result.get("result_string", ""),
                }

                # Consolidar parâmetros em uma única coluna estruturada
                if "params" in result and result["params"]:
                    exp_data["params"] = result["params"]

                # Consolidar metadados em uma única coluna estruturada
                if "metadata" in result and result["metadata"]:
                    exp_data["metadata"] = result["metadata"]

                # Adicionar informações do dataset se disponíveis
                if "dataset_info" in result and result["dataset_info"]:
                    exp_data["dataset_info"] = result["dataset_info"]

                experiments.append(exp_data)

            return pd.DataFrame(experiments)

        except Exception as e:
            print(f"❌ Erro ao processar dados: {e}")
            return pd.DataFrame()

    def _generate_plots(self, df: pd.DataFrame, report_dir: Path) -> None:
        """Gera todos os gráficos configurados."""
        plot_format = self.report_config.get("plot_format", "png")
        plots_dir = report_dir / "plots"
        plots_dir.mkdir(exist_ok=True)

        try:
            # Gráfico de comparação entre algoritmos
            if self.report_config.get("plot_comparison", True):
                self._plot_algorithm_comparison(df, plots_dir, plot_format)

            # Box plots
            if self.report_config.get("plot_boxplots", True):
                self._plot_boxplots(df, plots_dir, plot_format)

            # Scatter plots
            if self.report_config.get("plot_scatter", True):
                self._plot_scatter(df, plots_dir, plot_format)

            # Heatmap
            if self.report_config.get("plot_heatmap", True):
                self._plot_heatmap(df, plots_dir, plot_format)

            # Gráfico de tempo de execução
            if self.report_config.get("plot_runtime", True):
                self._plot_runtime(df, plots_dir, plot_format)

            # Taxa de sucesso
            if self.report_config.get("plot_success_rate", True):
                self._plot_success_rate(df, plots_dir, plot_format)

            # Qualidade vs Tempo
            if self.report_config.get("plot_quality_vs_time", True):
                self._plot_quality_vs_time(df, plots_dir, plot_format)

        except Exception as e:
            print(f"❌ Erro ao gerar gráficos: {e}")

    def _plot_algorithm_comparison(
        self, df: pd.DataFrame, plots_dir: Path, fmt: str
    ) -> None:
        """Gráfico de comparação entre algoritmos."""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

        # Comparação de distância máxima
        df.groupby("algorithm")["max_distance"].mean().plot(kind="bar", ax=ax1)
        ax1.set_title("Distância Máxima Média por Algoritmo")
        ax1.set_ylabel("Distância Máxima")
        ax1.tick_params(axis="x", rotation=45)

        # Comparação de tempo de execução
        df.groupby("algorithm")["runtime"].mean().plot(kind="bar", ax=ax2)
        ax2.set_title("Tempo de Execução Médio por Algoritmo")
        ax2.set_ylabel("Tempo (segundos)")
        ax2.tick_params(axis="x", rotation=45)

        plt.tight_layout()
        plt.savefig(
            plots_dir / f"algorithm_comparison.{fmt}", dpi=300, bbox_inches="tight"
        )
        plt.close()

    def _plot_boxplots(self, df: pd.DataFrame, plots_dir: Path, fmt: str) -> None:
        """Box plots das métricas por algoritmo."""
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))

        # Box plot da distância máxima
        sns.boxplot(data=df, x="algorithm", y="max_distance", ax=axes[0, 0])
        axes[0, 0].set_title("Distribuição da Distância Máxima")
        axes[0, 0].tick_params(axis="x", rotation=45)

        # Box plot do tempo de execução
        sns.boxplot(data=df, x="algorithm", y="runtime", ax=axes[0, 1])
        axes[0, 1].set_title("Distribuição do Tempo de Execução")
        axes[0, 1].tick_params(axis="x", rotation=45)

        # Violin plot da distância
        sns.violinplot(data=df, x="algorithm", y="max_distance", ax=axes[1, 0])
        axes[1, 0].set_title("Densidade da Distância Máxima")
        axes[1, 0].tick_params(axis="x", rotation=45)

        # Violin plot do tempo
        sns.violinplot(data=df, x="algorithm", y="runtime", ax=axes[1, 1])
        axes[1, 1].set_title("Densidade do Tempo de Execução")
        axes[1, 1].tick_params(axis="x", rotation=45)

        plt.tight_layout()
        plt.savefig(plots_dir / f"boxplots.{fmt}", dpi=300, bbox_inches="tight")
        plt.close()

    def _plot_scatter(self, df: pd.DataFrame, plots_dir: Path, fmt: str) -> None:
        """Scatter plot de qualidade vs tempo."""
        plt.figure(figsize=(10, 8))

        algorithms = df["algorithm"].unique()
        colors = plt.cm.tab10(range(len(algorithms)))

        for i, algo in enumerate(algorithms):
            algo_data = df[df["algorithm"] == algo]
            plt.scatter(
                algo_data["runtime"],
                algo_data["max_distance"],
                label=algo,
                alpha=0.7,
                c=[colors[i]],
                s=60,
            )

        plt.xlabel("Tempo de Execução (segundos)")
        plt.ylabel("Distância Máxima")
        plt.title("Qualidade vs Tempo de Execução")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig(
            plots_dir / f"scatter_quality_time.{fmt}", dpi=300, bbox_inches="tight"
        )
        plt.close()

    def _plot_heatmap(self, df: pd.DataFrame, plots_dir: Path, fmt: str) -> None:
        """Heatmap de correlação entre métricas."""
        # Selecionar apenas colunas numéricas
        numeric_cols = df.select_dtypes(include=["number"]).columns
        if len(numeric_cols) < 2:
            return

        plt.figure(figsize=(10, 8))
        correlation_matrix = df[numeric_cols].corr()
        sns.heatmap(
            correlation_matrix,
            annot=True,
            cmap="coolwarm",
            center=0,
            square=True,
            fmt=".2f",
        )
        plt.title("Matriz de Correlação entre Métricas")
        plt.tight_layout()
        plt.savefig(
            plots_dir / f"correlation_heatmap.{fmt}", dpi=300, bbox_inches="tight"
        )
        plt.close()

    def _plot_runtime(self, df: pd.DataFrame, plots_dir: Path, fmt: str) -> None:
        """Gráfico detalhado de tempo de execução."""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

        # Histograma de tempos
        df["runtime"].hist(bins=20, ax=ax1, alpha=0.7)
        ax1.set_title("Distribuição dos Tempos de Execução")
        ax1.set_xlabel("Tempo (segundos)")
        ax1.set_ylabel("Frequência")

        # Linha do tempo por experimento
        df.reset_index().plot(x="index", y="runtime", ax=ax2, marker="o")
        ax2.set_title("Tempo de Execução por Experimento")
        ax2.set_xlabel("Número do Experimento")
        ax2.set_ylabel("Tempo (segundos)")

        plt.tight_layout()
        plt.savefig(plots_dir / f"runtime_analysis.{fmt}", dpi=300, bbox_inches="tight")
        plt.close()

    def _plot_success_rate(self, df: pd.DataFrame, plots_dir: Path, fmt: str) -> None:
        """Gráfico de taxa de sucesso por algoritmo."""
        success_rate = df.groupby("algorithm")["success"].mean()

        plt.figure(figsize=(10, 6))
        success_rate.plot(kind="bar", color="skyblue", alpha=0.8)
        plt.title("Taxa de Sucesso por Algoritmo")
        plt.ylabel("Taxa de Sucesso")
        plt.xlabel("Algoritmo")
        plt.xticks(rotation=45)
        plt.ylim(0, 1)

        # Adicionar valores nas barras
        for i, v in enumerate(success_rate.values):
            plt.text(i, v + 0.02, f"{v:.1%}", ha="center", va="bottom")

        plt.tight_layout()
        plt.savefig(plots_dir / f"success_rate.{fmt}", dpi=300, bbox_inches="tight")
        plt.close()

    def _plot_quality_vs_time(
        self, df: pd.DataFrame, plots_dir: Path, fmt: str
    ) -> None:
        """Análise de qualidade vs tempo."""
        plt.figure(figsize=(12, 8))

        # Criar bins de tempo para análise
        df["time_bin"] = pd.cut(
            df["runtime"],
            bins=5,
            labels=["Muito Rápido", "Rápido", "Médio", "Lento", "Muito Lento"],
        )

        # Box plot de qualidade por bin de tempo
        sns.boxplot(data=df, x="time_bin", y="max_distance")
        plt.title("Qualidade da Solução vs Velocidade de Execução")
        plt.xlabel("Categoria de Tempo")
        plt.ylabel("Distância Máxima")
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(
            plots_dir / f"quality_vs_time_bins.{fmt}", dpi=300, bbox_inches="tight"
        )
        plt.close()

    def _calculate_statistics(self, df: pd.DataFrame) -> Dict[str, Any]:
        """Calcula estatísticas descritivas e testes estatísticos."""
        statistics = {}

        try:
            # Estatísticas gerais
            statistics["general"] = {
                "total_experiments": len(df),
                "unique_algorithms": df["algorithm"].nunique(),
                "unique_datasets": df["dataset"].nunique(),
                "overall_success_rate": df["success"].mean(),
                "mean_runtime": df["runtime"].mean(),
                "mean_max_distance": df["max_distance"].mean(),
            }

            # Estatísticas por algoritmo
            statistics["by_algorithm"] = {}
            for algo in df["algorithm"].unique():
                algo_data = df[df["algorithm"] == algo]
                statistics["by_algorithm"][algo] = {
                    "count": len(algo_data),
                    "success_rate": algo_data["success"].mean(),
                    "runtime_mean": algo_data["runtime"].mean(),
                    "runtime_std": algo_data["runtime"].std(),
                    "distance_mean": algo_data["max_distance"].mean(),
                    "distance_std": algo_data["max_distance"].std(),
                    "distance_median": algo_data["max_distance"].median(),
                }

            # Testes estatísticos se habilitado
            if self.report_config.get("statistical_tests", True):
                statistics["statistical_tests"] = self._perform_statistical_tests(df)

            # Análise de outliers se habilitado
            if self.report_config.get("include_outliers", True):
                statistics["outliers"] = self._detect_outliers(df)

        except Exception as e:
            print(f"❌ Erro ao calcular estatísticas: {e}")
            statistics["error"] = str(e)

        return statistics

    def _perform_statistical_tests(self, df: pd.DataFrame) -> Dict[str, Any]:
        """Realiza testes estatísticos entre algoritmos."""
        tests = {}

        try:
            algorithms = df["algorithm"].unique()
            if len(algorithms) < 2:
                return {"message": "Menos de 2 algoritmos, testes não aplicáveis"}

            # Teste ANOVA para distância máxima
            groups_distance = [
                df[df["algorithm"] == algo]["max_distance"] for algo in algorithms
            ]
            f_stat, p_value = stats.f_oneway(*groups_distance)
            tests["anova_distance"] = {"f_statistic": f_stat, "p_value": p_value}

            # Teste ANOVA para tempo de execução
            groups_runtime = [
                df[df["algorithm"] == algo]["runtime"] for algo in algorithms
            ]
            f_stat, p_value = stats.f_oneway(*groups_runtime)
            tests["anova_runtime"] = {"f_statistic": f_stat, "p_value": p_value}

            # Teste de normalidade (Shapiro-Wilk) para distâncias
            _, p_normal = stats.shapiro(df["max_distance"])
            tests["normality_distance"] = {"p_value": p_normal}

        except Exception as e:
            tests["error"] = str(e)

        return tests

    def _detect_outliers(self, df: pd.DataFrame) -> Dict[str, Any]:
        """Detecta outliers usando IQR."""
        outliers = {}

        try:
            for metric in ["max_distance", "runtime"]:
                Q1 = df[metric].quantile(0.25)
                Q3 = df[metric].quantile(0.75)
                IQR = Q3 - Q1
                lower_bound = Q1 - 1.5 * IQR
                upper_bound = Q3 + 1.5 * IQR

                outlier_indices = df[
                    (df[metric] < lower_bound) | (df[metric] > upper_bound)
                ].index.tolist()

                outliers[metric] = {
                    "count": len(outlier_indices),
                    "percentage": len(outlier_indices) / len(df) * 100,
                    "indices": outlier_indices,
                }

        except Exception as e:
            outliers["error"] = str(e)

        return outliers

    def _export_csv(self, df: pd.DataFrame, report_dir: Path) -> None:
        """Exporta dados para CSV."""
        try:
            csv_path = report_dir / "results_data.csv"
            df.to_csv(csv_path, index=False)
            print(f"📄 Dados exportados para CSV: {csv_path}")
        except Exception as e:
            print(f"❌ Erro ao exportar CSV: {e}")

    def _export_json(
        self, results_data: Dict[str, Any], statistics: Dict[str, Any], report_dir: Path
    ) -> None:
        """Exporta dados e estatísticas para JSON."""
        try:
            # Dados completos
            if self.export_config.get("export_raw_data", True):
                with open(report_dir / "raw_data.json", "w") as f:
                    json.dump(results_data, f, indent=2, default=str)

            # Estatísticas
            with open(report_dir / "statistics.json", "w") as f:
                json.dump(statistics, f, indent=2, default=str)

            print(f"📄 Dados exportados para JSON: {report_dir}")
        except Exception as e:
            print(f"❌ Erro ao exportar JSON: {e}")

    def _generate_html_report(
        self, df: pd.DataFrame, statistics: Dict[str, Any], report_dir: Path
    ) -> None:
        """Gera relatório HTML completo."""
        try:
            # Template HTML básico
            html_template = """
<!DOCTYPE html>
<html>
<head>
    <title>CSPBench - Relatório de Resultados</title>
    <meta charset="utf-8">
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        .header { background: #f4f4f4; padding: 20px; border-radius: 5px; }
        .section { margin: 20px 0; }
        .stats-table { border-collapse: collapse; width: 100%; }
        .stats-table th, .stats-table td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        .stats-table th { background-color: #f2f2f2; }
        .plot { margin: 20px 0; text-align: center; }
        .plot img { max-width: 100%; height: auto; }
    </style>
</head>
<body>
    <div class="header">
        <h1>CSPBench - Relatório de Resultados</h1>
        <p>Relatório gerado automaticamente com análise completa dos experimentos</p>
    </div>
    
    <div class="section">
        <h2>Resumo Geral</h2>
        <table class="stats-table">
            <tr><th>Métrica</th><th>Valor</th></tr>
            <tr><td>Total de Experimentos</td><td>{{ stats.general.total_experiments }}</td></tr>
            <tr><td>Algoritmos Únicos</td><td>{{ stats.general.unique_algorithms }}</td></tr>
            <tr><td>Datasets Únicos</td><td>{{ stats.general.unique_datasets }}</td></tr>
            <tr><td>Taxa de Sucesso Geral</td><td>{{ "%.1f%%" | format(stats.general.overall_success_rate * 100) }}</td></tr>
            <tr><td>Tempo Médio de Execução</td><td>{{ "%.2f segundos" | format(stats.general.mean_runtime) }}</td></tr>
            <tr><td>Distância Máxima Média</td><td>{{ "%.2f" | format(stats.general.mean_max_distance) }}</td></tr>
        </table>
    </div>
    
    <div class="section">
        <h2>Gráficos de Análise</h2>
        {% if plots_exist %}
        <div class="plot">
            <h3>Comparação entre Algoritmos</h3>
            <img src="plots/algorithm_comparison.png" alt="Comparação de Algoritmos">
        </div>
        <div class="plot">
            <h3>Distribuições (Box Plots)</h3>
            <img src="plots/boxplots.png" alt="Box Plots">
        </div>
        <div class="plot">
            <h3>Qualidade vs Tempo</h3>
            <img src="plots/scatter_quality_time.png" alt="Scatter Plot">
        </div>
        {% endif %}
    </div>
    
    <div class="section">
        <h2>Estatísticas por Algoritmo</h2>
        {% for algo, data in stats.by_algorithm.items() %}
        <h3>{{ algo }}</h3>
        <table class="stats-table">
            <tr><td>Número de Experimentos</td><td>{{ data.count }}</td></tr>
            <tr><td>Taxa de Sucesso</td><td>{{ "%.1f%%" | format(data.success_rate * 100) }}</td></tr>
            <tr><td>Tempo Médio</td><td>{{ "%.2f ± %.2f segundos" | format(data.runtime_mean, data.runtime_std) }}</td></tr>
            <tr><td>Distância Média</td><td>{{ "%.2f ± %.2f" | format(data.distance_mean, data.distance_std) }}</td></tr>
            <tr><td>Distância Mediana</td><td>{{ "%.2f" | format(data.distance_median) }}</td></tr>
        </table>
        {% endfor %}
    </div>
    
</body>
</html>
            """

            # Renderizar template
            from jinja2 import Template

            template = Template(html_template)
            html_content = template.render(
                stats=statistics, plots_exist=(report_dir / "plots").exists()
            )

            # Salvar HTML
            html_path = report_dir / "report.html"
            with open(html_path, "w", encoding="utf-8") as f:
                f.write(html_content)

            print(f"📋 Relatório HTML gerado: {html_path}")

        except Exception as e:
            print(f"❌ Erro ao gerar relatório HTML: {e}")
