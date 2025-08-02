"""
Execution Report Generator

Responsible for generating specific reports for normal algorithm executions,
including comparison charts, statistics and analyses.
"""

import json
import warnings
from pathlib import Path
from typing import Any, Dict

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from .history_plotter import HistoryPlotter

# Configurar warnings para suprimir avisos do seaborn sobre NaN
warnings.filterwarnings("ignore", category=RuntimeWarning, module="seaborn")

# Importar logger
from src.infrastructure.logging_config import get_logger


class ExecutionReportGenerator:
    """
    Specific report generator for algorithm executions.
    """

    def __init__(self, config: Dict[str, Any], session_path: Path):
        """
        Initialize the execution report generator.

        Args:
            config: System configuration
            session_path: Session path where to save the report
        """
        self.config = config
        self.session_path = session_path

        # Use new unified output configuration
        output_config = config.get("output", {})
        results_config = output_config.get("results", {})
        plots_config = output_config.get("plots", {})

        self.report_config = results_config.get("content", {})
        self.export_config = results_config.get("formats", {})

        # Get plot format from new config
        plot_formats = plots_config.get("formats", ["png"])

        self.logger = get_logger(__name__)
        self.plot_format = plot_formats[0] if plot_formats else "png"

        # Initialize history plot generator
        self.history_plotter = HistoryPlotter(config, session_path)

        # Configure matplotlib
        plt.style.use("seaborn-v0_8")
        sns.set_palette("husl")

    def generate_report(self, results_data: Dict[str, Any]) -> None:
        """
        Generate complete report for algorithm execution.

        Args:
            results_data: Execution results data
        """
        try:
            self.logger.info("Iniciando geração de relatório de execução")

            # Criar diretório de relatório
            report_dir = self.session_path / "report"
            report_dir.mkdir(exist_ok=True)

            # Processar dados
            df = self._process_results_data(results_data)
            if df.empty:
                self.logger.warning(
                    "Nenhum dado válido encontrado para gerar relatório"
                )
                return

            # Gerar estatísticas
            statistics = self._calculate_statistics(df)

            # Salvar dados processados
            self._save_processed_data(df, report_dir)

            # Gerar gráficos
            plots_dir = report_dir / "plots"
            plots_dir.mkdir(exist_ok=True)
            self._generate_plots(df, plots_dir)

            # Gerar gráficos de histórico
            if "batch_results" in results_data:
                self.history_plotter.generate_history_plots(
                    results_data["batch_results"]
                )

            # Gerar relatório HTML
            self._generate_html_report(df, statistics, report_dir)

            self.logger.info(f"Relatório de execução gerado em: {report_dir}")

        except Exception as e:
            self.logger.error(f"Erro ao gerar relatório de execução: {e}")
            raise

    def _process_results_data(self, results_data: Dict[str, Any]) -> pd.DataFrame:
        """
        Processa os dados de resultados para DataFrame.

        Args:
            results_data: Dados brutos dos resultados

        Returns:
            DataFrame processado
        """
        try:
            # Extrair batch_results se existir
            if "batch_results" in results_data:
                batch_results = results_data["batch_results"]
            else:
                batch_results = [results_data]

            processed_data = []

            for batch_result in batch_results:
                if "batch_summary" in batch_result:
                    # Estrutura de batch
                    for result in batch_result["batch_summary"]["results"]:
                        processed_result = self._flatten_result(result)
                        processed_data.append(processed_result)
                else:
                    # Resultado direto
                    processed_result = self._flatten_result(batch_result)
                    processed_data.append(processed_result)

            if not processed_data:
                return pd.DataFrame()

            df = pd.DataFrame(processed_data)

            # Adicionar colunas derivadas
            self._add_derived_columns(df)

            return df

        except Exception as e:
            self.logger.error(f"Erro ao processar dados: {e}")
            return pd.DataFrame()

    def _flatten_result(self, result: Dict[str, Any]) -> Dict[str, Any]:
        """
        Achata resultado para formato tabular.

        Args:
            result: Resultado individual

        Returns:
            Resultado achatado
        """
        flattened = {}

        # Campos básicos
        basic_fields = [
            "algorithm",
            "best_string",
            "max_distance",
            "execution_time",
            "status",
            "execution_name",
            "dataset_id",
            "algorithm_id",
            "algorithm_name",
            "repetition",
            "total_repetitions",
        ]

        for field in basic_fields:
            flattened[field] = result.get(field)

        # Informações do dataset
        if "dataset" in result:
            dataset_info = result["dataset"]
            flattened["dataset_size"] = dataset_info.get("size")
            flattened["dataset_length"] = dataset_info.get("length")
            flattened["dataset_alphabet"] = dataset_info.get("alphabet")

        # Parâmetros (convertidos para string para facilitar visualização)
        if "params" in result:
            params = result["params"]
            for key, value in params.items():
                flattened[f"param_{key}"] = value

        # Metadados específicos
        if "metadata" in result:
            metadata = result["metadata"]
            flattened["iterations"] = metadata.get(
                "iteracoes", metadata.get("iterations")
            )
            flattened["execution_time_internal"] = metadata.get("execution_time")

            # Para algoritmos evolutivos
            if "generations_executed" in metadata:
                flattened["generations"] = metadata["generations_executed"]
            if "best_fitness" in metadata:
                flattened["best_fitness"] = metadata["best_fitness"]

        return flattened

    def _add_derived_columns(self, df: pd.DataFrame) -> None:
        """
        Adiciona colunas derivadas ao DataFrame.

        Args:
            df: DataFrame para modificar
        """
        # Categoria de tempo de execução
        if "execution_time" in df.columns:
            df["time_bin"] = pd.cut(
                df["execution_time"],
                bins=[0, 0.001, 0.01, 0.1, 1.0, float("inf")],
                labels=["Muito Rápido", "Rápido", "Médio", "Lento", "Muito Lento"],
                include_lowest=True,
            )

        # Sucesso (status == 'success')
        if "status" in df.columns:
            df["success"] = df["status"] == "success"

        # Dataset como categoria
        if "dataset_id" in df.columns:
            df["dataset"] = df["dataset_id"]

    def _calculate_statistics(self, df: pd.DataFrame) -> Dict[str, Any]:
        """
        Calcula estatísticas descritivas dos dados.

        Args:
            df: DataFrame com os dados

        Returns:
            Dicionário com estatísticas
        """
        stats_dict = {}

        # Estatísticas gerais
        stats_dict["general"] = {
            "total_experiments": len(df),
            "unique_algorithms": (
                df["algorithm"].nunique() if "algorithm" in df.columns else 0
            ),
            "unique_datasets": (
                df["dataset"].nunique() if "dataset" in df.columns else 0
            ),
            "overall_success_rate": (
                df["success"].mean() if "success" in df.columns else 0
            ),
        }

        # Estatísticas por algoritmo
        if "algorithm" in df.columns:
            algo_stats = {}
            for algo in df["algorithm"].unique():
                algo_df = df[df["algorithm"] == algo]
                algo_stats[algo] = {
                    "count": len(algo_df),
                    "success_rate": (
                        algo_df["success"].mean() if "success" in algo_df.columns else 0
                    ),
                }

                # Estatísticas de runtime
                if "execution_time" in algo_df.columns:
                    runtime_data = algo_df["execution_time"].dropna()
                    if len(runtime_data) > 0:
                        algo_stats[algo]["runtime_mean"] = runtime_data.mean()
                        algo_stats[algo]["runtime_std"] = runtime_data.std()
                        algo_stats[algo]["runtime_median"] = runtime_data.median()

                # Estatísticas de distância
                if "max_distance" in algo_df.columns:
                    distance_data = algo_df["max_distance"].dropna()
                    if len(distance_data) > 0:
                        algo_stats[algo]["distance_mean"] = distance_data.mean()
                        algo_stats[algo]["distance_std"] = distance_data.std()
                        algo_stats[algo]["distance_median"] = distance_data.median()

            stats_dict["by_algorithm"] = algo_stats

        return stats_dict

    def _save_processed_data(self, df: pd.DataFrame, report_dir: Path) -> None:
        """
        Salva dados processados em diferentes formatos.

        Args:
            df: DataFrame para salvar
            report_dir: Diretório de destino
        """
        try:
            # CSV para análise manual
            csv_path = report_dir / "results_data.csv"
            df.to_csv(csv_path, index=False, encoding="utf-8")

            # JSON para integração
            json_path = report_dir / "raw_data.json"
            df.to_json(json_path, orient="records", indent=2, force_ascii=False)

            self.logger.info(f"Dados salvos em: {csv_path} e {json_path}")

        except Exception as e:
            self.logger.error(f"Erro ao salvar dados processados: {e}")

    def _generate_plots(self, df: pd.DataFrame, plots_dir: Path) -> None:
        """
        Gera gráficos de análise.

        Args:
            df: DataFrame com os dados
            plots_dir: Diretório para salvar gráficos
        """
        try:
            # Gráfico de comparação entre algoritmos
            if "algorithm" in df.columns and "max_distance" in df.columns:
                self._plot_algorithm_comparison(df, plots_dir)

            # Boxplots de métricas
            if "algorithm" in df.columns:
                self._plot_boxplots(df, plots_dir)

            # Análise de tempo de execução
            if "execution_time" in df.columns:
                self._plot_runtime_analysis(df, plots_dir)

        except Exception as e:
            self.logger.error(f"Erro ao gerar gráficos: {e}")

    def _plot_algorithm_comparison(self, df: pd.DataFrame, plots_dir: Path) -> None:
        """Gera gráfico de comparação entre algoritmos."""
        try:
            plt.figure(figsize=(12, 8))

            # Gráfico de barras com média e desvio padrão
            algo_stats = (
                df.groupby("algorithm")["max_distance"]
                .agg(["mean", "std"])
                .reset_index()
            )

            ax = sns.barplot(data=algo_stats, x="algorithm", y="mean")
            plt.errorbar(
                range(len(algo_stats)),
                algo_stats["mean"],
                yerr=algo_stats["std"],
                fmt="none",
                color="black",
                capsize=5,
            )

            plt.title("Comparação de Algoritmos - Distância Máxima")
            plt.xlabel("Algoritmo")
            plt.ylabel("Distância Máxima (média ± desvio padrão)")
            plt.xticks(rotation=45)
            plt.tight_layout()

            plt.savefig(
                plots_dir / f"algorithm_comparison.{self.plot_format}",
                dpi=300,
                bbox_inches="tight",
            )
            plt.close()

        except Exception as e:
            self.logger.error(f"Erro ao gerar gráfico de comparação: {e}")

    def _plot_boxplots(self, df: pd.DataFrame, plots_dir: Path) -> None:
        """Gera boxplots para métricas principais."""
        try:
            if "max_distance" in df.columns:
                plt.figure(figsize=(10, 6))
                sns.boxplot(data=df, x="algorithm", y="max_distance")
                plt.title("Distribuição de Distância Máxima por Algoritmo")
                plt.xticks(rotation=45)
                plt.tight_layout()
                plt.savefig(
                    plots_dir / f"distance_boxplots.{self.plot_format}",
                    dpi=300,
                    bbox_inches="tight",
                )
                plt.close()

            if "execution_time" in df.columns:
                plt.figure(figsize=(10, 6))
                # Usar escala logarítmica para tempo de execução
                plt.yscale("log")
                sns.boxplot(data=df, x="algorithm", y="execution_time")
                plt.title(
                    "Distribuição de Tempo de Execução por Algoritmo (escala log)"
                )
                plt.xticks(rotation=45)
                plt.tight_layout()
                plt.savefig(
                    plots_dir / f"runtime_boxplots.{self.plot_format}",
                    dpi=300,
                    bbox_inches="tight",
                )
                plt.close()

        except Exception as e:
            self.logger.error(f"Erro ao gerar boxplots: {e}")

    def _plot_runtime_analysis(self, df: pd.DataFrame, plots_dir: Path) -> None:
        """Gera análise de tempo de execução."""
        try:
            if "algorithm" in df.columns and "execution_time" in df.columns:
                plt.figure(figsize=(10, 6))

                # Histograma de tempos por algoritmo
                for algo in df["algorithm"].unique():
                    algo_data = df[df["algorithm"] == algo]["execution_time"]
                    plt.hist(algo_data, alpha=0.7, label=algo, bins=20)

                plt.xlabel("Tempo de Execução (s)")
                plt.ylabel("Frequência")
                plt.title("Distribuição de Tempos de Execução")
                plt.legend()
                plt.yscale("log")
                plt.tight_layout()

                plt.savefig(
                    plots_dir / f"runtime_distribution.{self.plot_format}",
                    dpi=300,
                    bbox_inches="tight",
                )
                plt.close()

        except Exception as e:
            self.logger.error(f"Erro ao gerar análise de runtime: {e}")

    def _generate_html_report(
        self, df: pd.DataFrame, statistics: Dict[str, Any], report_dir: Path
    ) -> None:
        """
        Gera relatório HTML final.

        Args:
            df: DataFrame com os dados
            statistics: Estatísticas calculadas
            report_dir: Diretório de destino
        """
        try:
            # Salvar estatísticas em JSON
            stats_path = report_dir / "statistics.json"
            with open(stats_path, "w", encoding="utf-8") as f:
                json.dump(statistics, f, indent=2, ensure_ascii=False, default=str)

            # Template HTML simples
            html_content = self._create_html_template(statistics, df)

            html_path = report_dir / "report.html"
            with open(html_path, "w", encoding="utf-8") as f:
                f.write(html_content)

            self.logger.info(f"Relatório HTML gerado em: {html_path}")

        except Exception as e:
            self.logger.error(f"Erro ao gerar relatório HTML: {e}")

    def _create_html_template(
        self, statistics: Dict[str, Any], df: pd.DataFrame
    ) -> str:
        """Cria template HTML simples para o relatório."""
        html = f"""
<!DOCTYPE html>
<html lang="pt-BR">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Relatório de Execução - CSPBench</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; }}
        .header {{ background-color: #f8f9fa; padding: 20px; border-radius: 8px; }}
        .section {{ margin: 20px 0; }}
        .stats-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 20px; }}
        .stat-card {{ background-color: #e9ecef; padding: 15px; border-radius: 8px; }}
        .plot-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(400px, 1fr)); gap: 20px; }}
        table {{ border-collapse: collapse; width: 100%; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #f2f2f2; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Relatório de Execução - CSPBench</h1>
        <p>Relatório gerado automaticamente com análise dos resultados dos algoritmos.</p>
    </div>

    <div class="section">
        <h2>Estatísticas Gerais</h2>
        <div class="stats-grid">
            <div class="stat-card">
                <h4>Total de Experimentos</h4>
                <p><strong>{statistics.get('general', {}).get('total_experiments', 0)}</strong></p>
            </div>
            <div class="stat-card">
                <h4>Algoritmos Únicos</h4>
                <p><strong>{statistics.get('general', {}).get('unique_algorithms', 0)}</strong></p>
            </div>
            <div class="stat-card">
                <h4>Datasets Únicos</h4>
                <p><strong>{statistics.get('general', {}).get('unique_datasets', 0)}</strong></p>
            </div>
            <div class="stat-card">
                <h4>Taxa de Sucesso</h4>
                <p><strong>{statistics.get('general', {}).get('overall_success_rate', 0):.2%}</strong></p>
            </div>
        </div>
    </div>

    <div class="section">
        <h2>Gráficos</h2>
        <div class="plot-grid">
            <div>
                <h4>Comparação de Algoritmos</h4>
                <img src="plots/algorithm_comparison.png" alt="Comparação de Algoritmos" style="max-width: 100%;">
            </div>
            <div>
                <h4>Distribuição de Distâncias</h4>
                <img src="plots/distance_boxplots.png" alt="Boxplots de Distância" style="max-width: 100%;">
            </div>
            <div>
                <h4>Distribuição de Tempos</h4>
                <img src="plots/runtime_boxplots.png" alt="Boxplots de Runtime" style="max-width: 100%;">
            </div>
            <div>
                <h4>Análise de Runtime</h4>
                <img src="plots/runtime_distribution.png" alt="Distribuição de Runtime" style="max-width: 100%;">
            </div>
        </div>
    </div>

    <div class="section">
        <h2>Arquivos de Dados</h2>
        <ul>
            <li><a href="results_data.csv">Dados em CSV</a></li>
            <li><a href="raw_data.json">Dados em JSON</a></li>
            <li><a href="statistics.json">Estatísticas Detalhadas</a></li>
        </ul>
    </div>
</body>
</html>
"""
        return html
