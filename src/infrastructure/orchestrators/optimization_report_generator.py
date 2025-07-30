"""
Optimization Report Generator - Gerador de Relatórios Avançados de Otimização

Gera 7 tipos de gráficos, estatísticas avançadas e relatórios em múltiplos formatos.
"""

import csv
import json
import os
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import optuna
import pandas as pd
import seaborn as sns
from scipy import stats
from sklearn.preprocessing import StandardScaler

from src.infrastructure.logging_config import get_logger

# Configurar matplotlib para não mostrar warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning, module="seaborn")
plt.style.use("seaborn-v0_8")


class OptimizationReportGenerator:
    """Gerador de relatórios avançados para otimização."""

    def __init__(self, destination: str, config: Dict[str, Any]):
        """
        Inicializa o gerador de relatórios.

        Args:
            destination: Diretório de destino
            config: Configuração completa
        """
        self.destination = Path(destination)
        self.config = config
        self.logger = get_logger(__name__)

        # Configuração de exportação
        self.export_config = config.get("export", {})
        self.formats = self.export_config.get("formats", {})
        self.include = self.export_config.get("include", [])

        # Configuração de gráficos
        self.plot_config = config.get("plots", {})
        self.plot_format = self.plot_config.get("formats", ["png"])[0] if self.plot_config.get("formats") else "png"

        # Criar diretórios
        self.plots_dir = self.destination / "plots"
        self.reports_dir = self.destination / "reports"
        self.data_dir = self.destination / "data"

        for dir_path in [self.plots_dir, self.reports_dir, self.data_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)

        # Configuração de cores e estilo
        self.colors = {
            "primary": "#2E86AB",
            "secondary": "#A23B72",
            "accent": "#F18F01",
            "success": "#C73E1D",
            "warning": "#F4A261",
            "info": "#264653",
            "light": "#F1FAEE",
            "dark": "#1D3557",
        }

    def generate_all_reports(self, study: optuna.Study, results: Dict[str, Any]):
        """
        Gera todos os relatórios configurados.

        Args:
            study: Estudo do Optuna
            results: Resultados da otimização
        """
        try:
            self.logger.info("Iniciando geração de relatórios")

            # Preparar dados
            df_trials = self._prepare_trials_dataframe(study)

            # Gerar gráficos
            self._generate_plots(study, df_trials, results)

            # Gerar relatórios de dados
            self._generate_data_reports(study, df_trials, results)

            # Gerar relatório HTML
            self._generate_html_report(study, df_trials, results)

            self.logger.info("Relatórios gerados com sucesso")

        except Exception as e:
            self.logger.error(f"Erro ao gerar relatórios: {e}")
            raise

    def _prepare_trials_dataframe(self, study: optuna.Study) -> pd.DataFrame:
        """Prepara DataFrame com dados dos trials."""
        trials_data = []

        for trial in study.trials:
            trial_data = {
                "trial_number": trial.number,
                "value": trial.value,
                "state": trial.state.name,
                "datetime_start": trial.datetime_start,
                "datetime_complete": trial.datetime_complete,
                "duration": trial.duration.total_seconds() if trial.duration else None,
            }

            # Adicionar parâmetros
            for param_name, param_value in trial.params.items():
                trial_data[f"param_{param_name}"] = param_value

            trials_data.append(trial_data)

        return pd.DataFrame(trials_data)

    def _generate_plots(
        self, study: optuna.Study, df_trials: pd.DataFrame, results: Dict[str, Any]
    ):
        """Gera todos os tipos de gráficos."""
        plot_types = [
            ("convergence", self._plot_convergence),
            ("comparison", self._plot_comparison),
            ("boxplots", self._plot_boxplots),
            ("scatter", self._plot_scatter),
            ("heatmap", self._plot_heatmap),
            ("runtime", self._plot_runtime),
            ("success_rate", self._plot_success_rate),
        ]

        for plot_name, plot_func in plot_types:
            if self.plot_config.get(plot_name, True):
                try:
                    self.logger.debug(f"Generating plot: {plot_name}")
                    plot_func(study, df_trials, results)
                except Exception as e:
                    self.logger.error(f"Error generating plot {plot_name}: {e}")

    def _plot_convergence(
        self, study: optuna.Study, df_trials: pd.DataFrame, results: Dict[str, Any]
    ):
        """Gera gráfico de convergência."""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

        # Gráfico de valores por trial
        completed_trials = df_trials[df_trials["state"] == "COMPLETE"]
        if not completed_trials.empty:
            x = completed_trials["trial_number"]
            y = completed_trials["value"]

            ax1.plot(
                x,
                y,
                "o-",
                color=self.colors["primary"],
                alpha=0.7,
                linewidth=2,
                markersize=4,
            )
            ax1.set_xlabel("Trial Number")
            ax1.set_ylabel("Objective Value")
            ax1.set_title("Convergência por Trial")
            ax1.grid(True, alpha=0.3)

            # Linha do melhor valor
            if study.best_value is not None:
                ax1.axhline(
                    y=study.best_value,
                    color=self.colors["success"],
                    linestyle="--",
                    linewidth=2,
                    label=f"Melhor: {study.best_value:.4f}",
                )
                ax1.legend()

        # Gráfico de melhores valores acumulados
        if not completed_trials.empty:
            if study.direction == optuna.study.StudyDirection.MINIMIZE:
                cumulative_best = completed_trials["value"].cummin()
            else:
                cumulative_best = completed_trials["value"].cummax()

            ax2.plot(
                completed_trials["trial_number"],
                cumulative_best,
                color=self.colors["secondary"],
                linewidth=3,
            )
            ax2.set_xlabel("Trial Number")
            ax2.set_ylabel("Best Value So Far")
            ax2.set_title("Convergência do Melhor Valor")
            ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(
            self.plots_dir / f"convergence.{self.plot_format}",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close()

    def _plot_comparison(
        self, study: optuna.Study, df_trials: pd.DataFrame, results: Dict[str, Any]
    ):
        """Gera gráfico de comparação de parâmetros."""
        completed_trials = df_trials[df_trials["state"] == "COMPLETE"]
        if completed_trials.empty:
            return

        # Identificar colunas de parâmetros
        param_cols = [
            col for col in completed_trials.columns if col.startswith("param_")
        ]

        if not param_cols:
            return

        # Determinar número de subplots
        n_params = len(param_cols)
        n_cols = min(3, n_params)
        n_rows = (n_params + n_cols - 1) // n_cols

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
        if n_params == 1:
            axes = [axes]
        elif n_rows == 1:
            axes = axes.reshape(1, -1)

        for i, param_col in enumerate(param_cols):
            row = i // n_cols
            col = i % n_cols
            ax = axes[row, col] if n_rows > 1 else axes[col]

            param_name = param_col.replace("param_", "")

            # Determinar tipo de parâmetro
            param_values = completed_trials[param_col].dropna()
            if param_values.dtype == "object" or len(param_values.unique()) <= 10:
                # Parâmetro categórico
                value_counts = param_values.value_counts()
                ax.bar(
                    range(len(value_counts)),
                    value_counts.values,
                    color=self.colors["accent"],
                    alpha=0.8,
                )
                ax.set_xticks(range(len(value_counts)))
                ax.set_xticklabels(value_counts.index, rotation=45)
                ax.set_ylabel("Frequency")
            else:
                # Parâmetro numérico
                ax.hist(
                    param_values,
                    bins=20,
                    color=self.colors["primary"],
                    alpha=0.7,
                    edgecolor="black",
                )
                ax.set_ylabel("Frequency")

            ax.set_title(f"Distribuição: {param_name}")
            ax.grid(True, alpha=0.3)

        # Remover subplots vazios
        for i in range(n_params, n_rows * n_cols):
            row = i // n_cols
            col = i % n_cols
            if n_rows > 1:
                fig.delaxes(axes[row, col])
            else:
                fig.delaxes(axes[col])

        plt.tight_layout()
        plt.savefig(
            self.plots_dir / f"comparison.{self.plot_format}",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close()

    def _plot_boxplots(
        self, study: optuna.Study, df_trials: pd.DataFrame, results: Dict[str, Any]
    ):
        """Gera boxplots dos parâmetros."""
        completed_trials = df_trials[df_trials["state"] == "COMPLETE"]
        if completed_trials.empty:
            return

        # Identificar parâmetros categóricos
        param_cols = [
            col for col in completed_trials.columns if col.startswith("param_")
        ]
        categorical_params = []

        for param_col in param_cols:
            param_values = completed_trials[param_col].dropna()
            if param_values.dtype == "object" or len(param_values.unique()) <= 10:
                categorical_params.append(param_col)

        if not categorical_params:
            return

        n_params = len(categorical_params)
        n_cols = min(2, n_params)
        n_rows = (n_params + n_cols - 1) // n_cols

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(8 * n_cols, 6 * n_rows))
        if n_params == 1:
            axes = [axes]
        elif n_rows == 1:
            axes = axes.reshape(1, -1)

        for i, param_col in enumerate(categorical_params):
            row = i // n_cols
            col = i % n_cols
            ax = axes[row, col] if n_rows > 1 else axes[col]

            param_name = param_col.replace("param_", "")

            # Criar dados para boxplot
            groups = completed_trials.groupby(param_col)["value"].apply(list).to_dict()

            if len(groups) > 1:
                data = list(groups.values())
                labels = list(groups.keys())

                bp = ax.boxplot(data, labels=labels, patch_artist=True)

                # Colorir boxes
                for patch, color in zip(
                    bp["boxes"], plt.cm.Set3(np.linspace(0, 1, len(labels)))
                ):
                    patch.set_facecolor(color)
                    patch.set_alpha(0.7)

                ax.set_ylabel("Objective Value")
                ax.set_title(f"Boxplot: {param_name}")
                ax.grid(True, alpha=0.3)

        # Remover subplots vazios
        for i in range(n_params, n_rows * n_cols):
            row = i // n_cols
            col = i % n_cols
            if n_rows > 1:
                fig.delaxes(axes[row, col])
            else:
                fig.delaxes(axes[col])

        plt.tight_layout()
        plt.savefig(
            self.plots_dir / f"boxplots.{self.plot_format}",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close()

    def _plot_scatter(
        self, study: optuna.Study, df_trials: pd.DataFrame, results: Dict[str, Any]
    ):
        """Gera scatter plots dos parâmetros vs objetivo."""
        completed_trials = df_trials[df_trials["state"] == "COMPLETE"]
        if completed_trials.empty:
            return

        # Identificar parâmetros numéricos
        param_cols = [
            col for col in completed_trials.columns if col.startswith("param_")
        ]
        numeric_params = []

        for param_col in param_cols:
            param_values = completed_trials[param_col].dropna()
            if param_values.dtype != "object" and len(param_values.unique()) > 10:
                numeric_params.append(param_col)

        if not numeric_params:
            return

        n_params = len(numeric_params)
        n_cols = min(3, n_params)
        n_rows = (n_params + n_cols - 1) // n_cols

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
        if n_params == 1:
            axes = [axes]
        elif n_rows == 1:
            axes = axes.reshape(1, -1)

        for i, param_col in enumerate(numeric_params):
            row = i // n_cols
            col = i % n_cols
            ax = axes[row, col] if n_rows > 1 else axes[col]

            param_name = param_col.replace("param_", "")

            x = completed_trials[param_col]
            y = completed_trials["value"]

            # Scatter plot
            ax.scatter(
                x,
                y,
                alpha=0.6,
                color=self.colors["primary"],
                s=50,
                edgecolors="black",
                linewidth=0.5,
            )

            # Linha de tendência
            try:
                z = np.polyfit(x, y, 1)
                p = np.poly1d(z)
                ax.plot(x, p(x), color=self.colors["secondary"], linewidth=2, alpha=0.8)

                # Correlação
                correlation = np.corrcoef(x, y)[0, 1]
                ax.text(
                    0.05,
                    0.95,
                    f"r = {correlation:.3f}",
                    transform=ax.transAxes,
                    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
                )
            except:
                pass

            ax.set_xlabel(param_name)
            ax.set_ylabel("Objective Value")
            ax.set_title(f"Scatter: {param_name}")
            ax.grid(True, alpha=0.3)

        # Remover subplots vazios
        for i in range(n_params, n_rows * n_cols):
            row = i // n_cols
            col = i % n_cols
            if n_rows > 1:
                fig.delaxes(axes[row, col])
            else:
                fig.delaxes(axes[col])

        plt.tight_layout()
        plt.savefig(
            self.plots_dir / f"scatter.{self.plot_format}", dpi=300, bbox_inches="tight"
        )
        plt.close()

    def _plot_heatmap(
        self, study: optuna.Study, df_trials: pd.DataFrame, results: Dict[str, Any]
    ):
        """Gera heatmap de correlações entre parâmetros."""
        completed_trials = df_trials[df_trials["state"] == "COMPLETE"]
        if completed_trials.empty:
            return

        # Identificar parâmetros numéricos
        param_cols = [
            col for col in completed_trials.columns if col.startswith("param_")
        ]
        numeric_params = []

        for param_col in param_cols:
            param_values = completed_trials[param_col].dropna()
            # Verificar se é numérico e tem variância
            if (
                param_values.dtype != "object"
                and len(param_values) > 1
                and param_values.var() > 0
            ):
                numeric_params.append(param_col)

        if len(numeric_params) < 2:
            self.logger.warning(
                "Heatmap: Parâmetros numéricos insuficientes para correlação"
            )
            return

        # Criar matriz de correlação
        correlation_data = completed_trials[numeric_params + ["value"]].corr()

        # Verificar se a matriz contém dados válidos
        valid_data = correlation_data.dropna(axis=0, how="all").dropna(
            axis=1, how="all"
        )
        if valid_data.empty or valid_data.isna().all().all():
            self.logger.warning("Heatmap: Matriz de correlação contém apenas NaN")
            return

        # Renomear colunas
        correlation_data.columns = [
            col.replace("param_", "") for col in correlation_data.columns
        ]
        correlation_data.index = [
            idx.replace("param_", "") for idx in correlation_data.index
        ]

        # Criar heatmap
        plt.figure(figsize=(10, 8))
        mask = np.triu(np.ones_like(correlation_data, dtype=bool))

        sns.heatmap(
            correlation_data,
            mask=mask,
            annot=True,
            cmap="RdBu_r",
            center=0,
            square=True,
            fmt=".3f",
            cbar_kws={"shrink": 0.5},
        )

        plt.title("Correlação entre Parâmetros e Objetivo")
        plt.tight_layout()
        plt.savefig(
            self.plots_dir / f"heatmap.{self.plot_format}", dpi=300, bbox_inches="tight"
        )
        plt.close()

    def _plot_runtime(
        self, study: optuna.Study, df_trials: pd.DataFrame, results: Dict[str, Any]
    ):
        """Gera gráfico de tempo de execução."""
        completed_trials = df_trials[df_trials["state"] == "COMPLETE"]
        if completed_trials.empty or "duration" not in completed_trials.columns:
            return

        duration_data = completed_trials["duration"].dropna()
        if duration_data.empty:
            return

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

        # Histograma de durações
        ax1.hist(
            duration_data,
            bins=20,
            color=self.colors["info"],
            alpha=0.7,
            edgecolor="black",
        )
        ax1.set_xlabel("Duration (seconds)")
        ax1.set_ylabel("Frequency")
        ax1.set_title("Distribuição de Tempo de Execução")
        ax1.grid(True, alpha=0.3)

        # Tempo vs trial number
        ax2.scatter(
            completed_trials["trial_number"],
            duration_data,
            alpha=0.6,
            color=self.colors["warning"],
            s=50,
            edgecolors="black",
            linewidth=0.5,
        )
        ax2.set_xlabel("Trial Number")
        ax2.set_ylabel("Duration (seconds)")
        ax2.set_title("Tempo de Execução por Trial")
        ax2.grid(True, alpha=0.3)

        # Linha de tendência
        try:
            z = np.polyfit(completed_trials["trial_number"], duration_data, 1)
            p = np.poly1d(z)
            ax2.plot(
                completed_trials["trial_number"],
                p(completed_trials["trial_number"]),
                color=self.colors["secondary"],
                linewidth=2,
                alpha=0.8,
            )
        except:
            pass

        plt.tight_layout()
        plt.savefig(
            self.plots_dir / f"runtime.{self.plot_format}", dpi=300, bbox_inches="tight"
        )
        plt.close()

    def _plot_success_rate(
        self, study: optuna.Study, df_trials: pd.DataFrame, results: Dict[str, Any]
    ):
        """Gera gráfico de taxa de sucesso."""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

        # Pizza de estados
        state_counts = df_trials["state"].value_counts()
        colors = [
            self.colors["primary"],
            self.colors["secondary"],
            self.colors["warning"],
        ][: len(state_counts)]

        ax1.pie(
            state_counts.values,
            labels=state_counts.index,
            autopct="%1.1f%%",
            colors=colors,
            startangle=90,
        )
        ax1.set_title("Distribuição de Estados dos Trials")

        # Taxa de sucesso ao longo do tempo
        window_size = max(10, len(df_trials) // 20)
        success_rates = []
        trial_numbers = []

        for i in range(window_size, len(df_trials) + 1):
            window_trials = df_trials.iloc[i - window_size : i]
            success_rate = (window_trials["state"] == "COMPLETE").sum() / len(
                window_trials
            )
            success_rates.append(success_rate * 100)
            trial_numbers.append(i)

        if success_rates:
            ax2.plot(
                trial_numbers, success_rates, color=self.colors["success"], linewidth=2
            )
            ax2.set_xlabel("Trial Number")
            ax2.set_ylabel("Success Rate (%)")
            ax2.set_title(f"Taxa de Sucesso (janela de {window_size} trials)")
            ax2.grid(True, alpha=0.3)
            ax2.set_ylim(0, 100)

        plt.tight_layout()
        plt.savefig(
            self.plots_dir / f"success_rate.{self.plot_format}",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close()

    def _generate_data_reports(
        self, study: optuna.Study, df_trials: pd.DataFrame, results: Dict[str, Any]
    ):
        """Gera relatórios de dados em diferentes formatos."""
        # JSON completo
        if "study_results" in self.include:
            json_format = self.formats.get("study_results", "json")
            if json_format == "json":
                json_file = self.data_dir / "study_results.json"
                with open(json_file, "w", encoding="utf-8") as f:
                    json.dump(results, f, indent=2, ensure_ascii=False)

        # CSV de trials
        if "optimization_history" in self.include:
            csv_format = self.formats.get("optimization_history", "csv")
            if csv_format == "csv":
                csv_file = self.data_dir / "optimization_history.csv"
                df_trials.to_csv(csv_file, index=False)

        # Melhores parâmetros
        if "best_parameters" in self.include:
            params_format = self.formats.get("best_params", "yaml")
            if params_format == "yaml":
                import yaml

                yaml_file = self.data_dir / "best_parameters.yaml"
                with open(yaml_file, "w", encoding="utf-8") as f:
                    yaml.dump(results["best_params"], f, default_flow_style=False)
            elif params_format == "json":
                json_file = self.data_dir / "best_parameters.json"
                with open(json_file, "w", encoding="utf-8") as f:
                    json.dump(results["best_params"], f, indent=2, ensure_ascii=False)

        # Estatísticas avançadas
        if "trial_details" in self.include:
            self._generate_advanced_statistics(df_trials, results)

    def _generate_advanced_statistics(
        self, df_trials: pd.DataFrame, results: Dict[str, Any]
    ):
        """Gera estatísticas avançadas."""
        completed_trials = df_trials[df_trials["state"] == "COMPLETE"]
        if completed_trials.empty:
            return

        statistics = {
            "summary": {
                "total_trials": len(df_trials),
                "completed_trials": len(completed_trials),
                "success_rate": len(completed_trials) / len(df_trials) * 100,
                "best_value": results["best_value"],
                "mean_value": completed_trials["value"].mean(),
                "std_value": completed_trials["value"].std(),
                "median_value": completed_trials["value"].median(),
                "min_value": completed_trials["value"].min(),
                "max_value": completed_trials["value"].max(),
            },
            "parameter_analysis": {},
            "statistical_tests": {},
            "outliers": {},
        }

        # Análise de parâmetros
        param_cols = [
            col for col in completed_trials.columns if col.startswith("param_")
        ]
        for param_col in param_cols:
            param_name = param_col.replace("param_", "")
            param_values = completed_trials[param_col].dropna()

            if param_values.dtype != "object":
                # Parâmetro numérico
                statistics["parameter_analysis"][param_name] = {
                    "mean": param_values.mean(),
                    "std": param_values.std(),
                    "min": param_values.min(),
                    "max": param_values.max(),
                    "median": param_values.median(),
                    "correlation_with_objective": param_values.corr(
                        completed_trials["value"]
                    ),
                }

                # Teste de normalidade
                if len(param_values) > 8:
                    _, p_value = stats.shapiro(param_values)
                    statistics["statistical_tests"][f"{param_name}_normality"] = {
                        "test": "Shapiro-Wilk",
                        "p_value": p_value,
                        "is_normal": p_value > 0.05,
                    }

                # Detecção de outliers
                Q1 = param_values.quantile(0.25)
                Q3 = param_values.quantile(0.75)
                IQR = Q3 - Q1
                lower_bound = Q1 - 1.5 * IQR
                upper_bound = Q3 + 1.5 * IQR

                outliers = param_values[
                    (param_values < lower_bound) | (param_values > upper_bound)
                ]
                statistics["outliers"][param_name] = {
                    "count": len(outliers),
                    "percentage": len(outliers) / len(param_values) * 100,
                    "values": outliers.tolist(),
                }

        # ANOVA para parâmetros categóricos
        categorical_params = [
            col
            for col in param_cols
            if completed_trials[col].dtype == "object"
            or len(completed_trials[col].unique()) <= 10
        ]

        for param_col in categorical_params:
            param_name = param_col.replace("param_", "")
            groups = [
                group["value"].values
                for name, group in completed_trials.groupby(param_col)
            ]

            if len(groups) > 1 and all(len(group) > 1 for group in groups):
                try:
                    f_stat, p_value = stats.f_oneway(*groups)
                    statistics["statistical_tests"][f"{param_name}_anova"] = {
                        "test": "One-way ANOVA",
                        "f_statistic": f_stat,
                        "p_value": p_value,
                        "significant": p_value < 0.05,
                    }
                except:
                    pass

        # Salvar estatísticas
        stats_file = self.data_dir / "advanced_statistics.json"
        with open(stats_file, "w", encoding="utf-8") as f:
            json.dump(statistics, f, indent=2, ensure_ascii=False)

    def _generate_html_report(
        self, study: optuna.Study, df_trials: pd.DataFrame, results: Dict[str, Any]
    ):
        """Gera relatório HTML completo."""
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>Relatório de Otimização - {results['study_name']}</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                .header {{ background-color: #f0f0f0; padding: 20px; border-radius: 5px; }}
                .section {{ margin: 20px 0; }}
                .metric {{ background-color: #e6f3ff; padding: 10px; margin: 5px; border-radius: 3px; }}
                .plot {{ text-align: center; margin: 20px 0; }}
                table {{ border-collapse: collapse; width: 100%; }}
                th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
                th {{ background-color: #f2f2f2; }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>Relatório de Otimização</h1>
                <p><strong>Estudo:</strong> {results['study_name']}</p>
                <p><strong>Algoritmo:</strong> {results['algorithm']}</p>
                <p><strong>Dataset:</strong> {results['dataset']}</p>
                <p><strong>Direção:</strong> {results['direction']}</p>
                <p><strong>Data:</strong> {results.get('start_time', 'N/A')}</p>
            </div>
            
            <div class="section">
                <h2>Resumo Executivo</h2>
                <div class="metric">
                    <strong>Melhor Valor:</strong> {results['best_value']}
                </div>
                <div class="metric">
                    <strong>Número de Trials:</strong> {results['n_trials']}
                </div>
                <div class="metric">
                    <strong>Tempo Total:</strong> {results['total_time']:.2f} segundos
                </div>
                <div class="metric">
                    <strong>Melhor Trial:</strong> {results['best_trial']}
                </div>
            </div>
            
            <div class="section">
                <h2>Melhores Parâmetros</h2>
                <table>
                    <tr><th>Parâmetro</th><th>Valor</th></tr>
        """

        for param, value in results["best_params"].items():
            html_content += f"<tr><td>{param}</td><td>{value}</td></tr>"

        html_content += """
                </table>
            </div>
            
            <div class="section">
                <h2>Gráficos</h2>
        """

        # Adicionar gráficos
        plot_files = list(self.plots_dir.glob(f"*.{self.plot_format}"))
        for plot_file in plot_files:
            plot_name = plot_file.stem.replace("_", " ").title()
            html_content += f"""
                <div class="plot">
                    <h3>{plot_name}</h3>
                    <img src="plots/{plot_file.name}" style="max-width: 100%;">
                </div>
            """

        html_content += """
            </div>
        </body>
        </html>
        """

        # Salvar relatório HTML
        html_file = self.reports_dir / "optimization_report.html"
        with open(html_file, "w", encoding="utf-8") as f:
            f.write(html_content)

        self.logger.info(f"Relatório HTML salvo: {html_file}")
