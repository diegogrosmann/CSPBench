"""
Módulo de visualização para otimização e análise de sensibilidade.

Classes:
    OptimizationVisualizer: Visualizador para resultados de otimização
    SensitivityVisualizer: Visualizador para análise de sensibilidade

Funções:
    plot_optimization_history: Plota histórico de otimização
    plot_parameter_importance: Plota importância dos parâmetros
    save_visualization: Salva visualização em arquivo
"""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from src.optimization.optuna_optimizer import OptimizationResult
from src.optimization.sensitivity_analyzer import SensitivityResult

logger = logging.getLogger(__name__)


class OptimizationVisualizer:
    """Visualizador para resultados de otimização."""

    def __init__(self, result: OptimizationResult):
        self.result = result

    def plot_optimization_history(
        self, save_path: str | None = None, interactive: bool = True
    ):
        """Plota histórico de otimização."""

        if interactive:
            return self._plot_history_interactive(save_path)
        else:
            return self._plot_history_static(save_path)

    def _plot_history_interactive(self, save_path: str | None = None):
        """Plota histórico interativo com Plotly."""

        # Extrair dados
        trials = [t for t in self.result.all_trials if t["value"] is not None]
        if not trials:
            logger.warning("Nenhum trial válido para plotar")
            return None

        trial_numbers = [t["number"] for t in trials]
        values = [t["value"] for t in trials]

        # Calcular melhor valor até agora
        best_values = []
        current_best = float("inf")
        for value in values:
            if value < current_best:
                current_best = value
            best_values.append(current_best)

        # Criar subplot
        fig = make_subplots(
            rows=2,
            cols=1,
            subplot_titles=("Valores dos Trials", "Melhor Valor até o Momento"),
            vertical_spacing=0.1,
        )

        # Plot 1: Valores dos trials
        fig.add_trace(
            go.Scatter(
                x=trial_numbers,
                y=values,
                mode="markers",
                name="Trial Values",
                marker={"color": "blue", "size": 6},
                hovertemplate="Trial: %{x}<br>Valor: %{y:.6f}<extra></extra>",
            ),
            row=1,
            col=1,
        )

        # Plot 2: Melhor valor
        fig.add_trace(
            go.Scatter(
                x=trial_numbers,
                y=best_values,
                mode="lines",
                name="Melhor Valor",
                line={"color": "red", "width": 2},
                hovertemplate="Trial: %{x}<br>Melhor: %{y:.6f}<extra></extra>",
            ),
            row=2,
            col=1,
        )

        # Layout
        fig.update_layout(
            title=f"Histórico de Otimização - {self.result.study_name}",
            height=600,
            showlegend=False,
        )

        fig.update_xaxes(title_text="Trial", row=2, col=1)
        fig.update_yaxes(title_text="Valor Objetivo", row=1, col=1)
        fig.update_yaxes(title_text="Melhor Valor", row=2, col=1)

        if save_path:
            fig.write_html(save_path)
            logger.info("Gráfico salvo em %s", save_path)

        return fig

    def _plot_history_static(self, save_path: str | None = None):
        """Plota histórico estático com Matplotlib."""

        # Extrair dados
        trials = [t for t in self.result.all_trials if t["value"] is not None]
        if not trials:
            logger.warning("Nenhum trial válido para plotar")
            return None

        trial_numbers = [t["number"] for t in trials]
        values = [t["value"] for t in trials]

        # Calcular melhor valor até agora
        best_values = []
        current_best = float("inf")
        for value in values:
            if value < current_best:
                current_best = value
            best_values.append(current_best)

        # Criar figura
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

        # Plot 1: Valores dos trials
        ax1.scatter(trial_numbers, values, alpha=0.6, s=30)
        ax1.set_ylabel("Valor Objetivo")
        ax1.set_title(f"Histórico de Otimização - {self.result.study_name}")
        ax1.grid(True, alpha=0.3)

        # Plot 2: Melhor valor
        ax2.plot(trial_numbers, best_values, "r-", linewidth=2)
        ax2.set_xlabel("Trial")
        ax2.set_ylabel("Melhor Valor")
        ax2.set_title("Melhor Valor até o Momento")
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches="tight")
            logger.info("Gráfico salvo em %s", save_path)

        return fig

    def plot_parameter_importance(
        self, save_path: str | None = None, interactive: bool = True
    ):
        """Plota importância dos parâmetros baseado na variância."""

        # Calcular variância de cada parâmetro
        param_variance = {}
        trials = [t for t in self.result.all_trials if t["value"] is not None]

        if not trials:
            logger.warning("Nenhum trial válido para calcular importância")
            return None

        # Coletar todos os parâmetros
        all_params = set()
        for trial in trials:
            all_params.update(trial["params"].keys())

        for param in all_params:
            values = []
            for trial in trials:
                if param in trial["params"]:
                    values.append(trial["params"][param])

            if values:
                param_variance[param] = np.var(values)

        if not param_variance:
            logger.warning("Nenhum parâmetro encontrado")
            return None

        # Ordenar por importância
        sorted_params = sorted(param_variance.items(), key=lambda x: x[1], reverse=True)
        param_names = [p[0] for p in sorted_params]
        variances = [p[1] for p in sorted_params]

        if interactive:
            # Plot interativo
            fig = go.Figure(
                data=[
                    go.Bar(
                        x=variances,
                        y=param_names,
                        orientation="h",
                        marker_color="steelblue",
                    )
                ]
            )

            fig.update_layout(
                title="Importância dos Parâmetros (Variância)",
                xaxis_title="Variância",
                yaxis_title="Parâmetros",
                height=max(400, len(param_names) * 30),
            )

            if save_path:
                fig.write_html(save_path)
                logger.info("Gráfico salvo em %s", save_path)

            return fig

        else:
            # Plot estático
            fig, ax = plt.subplots(figsize=(10, max(6, len(param_names) * 0.4)))

            bars = ax.barh(param_names, variances)
            ax.set_xlabel("Variância")
            ax.set_title("Importância dos Parâmetros (Variância)")
            ax.grid(True, alpha=0.3, axis="x")

            plt.tight_layout()

            if save_path:
                plt.savefig(save_path, dpi=300, bbox_inches="tight")
                logger.info("Gráfico salvo em %s", save_path)

            return fig


class SensitivityVisualizer:
    """Visualizador para análise de sensibilidade."""

    def __init__(self, result: SensitivityResult):
        self.result = result

    def plot_sensitivity_indices(
        self, save_path: str | None = None, interactive: bool = True
    ):
        """Plota índices de sensibilidade."""

        if interactive:
            return self._plot_sensitivity_interactive(save_path)
        else:
            return self._plot_sensitivity_static(save_path)

    def _plot_sensitivity_interactive(self, save_path: str | None = None):
        """Plota índices de sensibilidade interativo."""

        param_names = self.result.parameter_names
        s1_values = [self.result.first_order[p] for p in param_names]
        st_values = [self.result.total_order[p] for p in param_names]

        # Criar subplot
        fig = make_subplots(
            rows=1,
            cols=2,
            subplot_titles=("Índices de Primeira Ordem (S1)", "Índices Totais (ST)"),
            horizontal_spacing=0.1,
        )

        # S1
        fig.add_trace(
            go.Bar(
                x=s1_values,
                y=param_names,
                orientation="h",
                name="S1",
                marker_color="lightblue",
                hovertemplate="%{y}: %{x:.4f}<extra></extra>",
            ),
            row=1,
            col=1,
        )

        # ST
        fig.add_trace(
            go.Bar(
                x=st_values,
                y=param_names,
                orientation="h",
                name="ST",
                marker_color="orange",
                hovertemplate="%{y}: %{x:.4f}<extra></extra>",
            ),
            row=1,
            col=2,
        )

        fig.update_layout(
            title=f"Análise de Sensibilidade - {self.result.method.upper()}",
            height=max(400, len(param_names) * 30),
            showlegend=False,
        )

        fig.update_xaxes(title_text="Índice S1", row=1, col=1)
        fig.update_xaxes(title_text="Índice ST", row=1, col=2)

        if save_path:
            fig.write_html(save_path)
            logger.info("Gráfico salvo em %s", save_path)

        return fig

    def _plot_sensitivity_static(self, save_path: str | None = None):
        """Plota índices de sensibilidade estático."""

        param_names = self.result.parameter_names
        s1_values = [self.result.first_order[p] for p in param_names]
        st_values = [self.result.total_order[p] for p in param_names]

        # Criar figura
        fig, (ax1, ax2) = plt.subplots(
            1, 2, figsize=(15, max(6, len(param_names) * 0.4))
        )

        # S1
        ax1.barh(param_names, s1_values, color="lightblue")
        ax1.set_xlabel("Índice S1")
        ax1.set_title("Índices de Primeira Ordem")
        ax1.grid(True, alpha=0.3, axis="x")

        # ST
        ax2.barh(param_names, st_values, color="orange")
        ax2.set_xlabel("Índice ST")
        ax2.set_title("Índices Totais")
        ax2.grid(True, alpha=0.3, axis="x")

        plt.suptitle(f"Análise de Sensibilidade - {self.result.method.upper()}")
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches="tight")
            logger.info("Gráfico salvo em %s", save_path)

        return fig

    def plot_parameter_ranking(
        self, save_path: str | None = None, interactive: bool = True
    ):
        """Plota ranking dos parâmetros por importância."""

        # Ordenar por índices totais
        sorted_params = sorted(
            self.result.total_order.items(), key=lambda x: x[1], reverse=True
        )

        param_names = [p[0] for p in sorted_params]
        st_values = [p[1] for p in sorted_params]
        s1_values = [self.result.first_order[p[0]] for p in sorted_params]

        if interactive:
            # Plot interativo
            fig = go.Figure()

            fig.add_trace(
                go.Bar(
                    x=param_names,
                    y=s1_values,
                    name="S1 (Primeira Ordem)",
                    marker_color="lightblue",
                )
            )

            fig.add_trace(
                go.Bar(
                    x=param_names,
                    y=st_values,
                    name="ST (Total)",
                    marker_color="orange",
                    opacity=0.7,
                )
            )

            fig.update_layout(
                title="Ranking de Importância dos Parâmetros",
                xaxis_title="Parâmetros",
                yaxis_title="Índice de Sensibilidade",
                barmode="overlay",
            )

            if save_path:
                fig.write_html(save_path)
                logger.info("Gráfico salvo em %s", save_path)

            return fig

        else:
            # Plot estático
            fig, ax = plt.subplots(figsize=(12, 6))

            x = np.arange(len(param_names), dtype=float)
            width = 0.35

            ax.bar(
                x - width / 2,
                s1_values,
                width,
                label="S1 (Primeira Ordem)",
                color="lightblue",
            )
            ax.bar(
                x + width / 2,
                st_values,
                width,
                label="ST (Total)",
                color="orange",
                alpha=0.7,
            )

            ax.set_xlabel("Parâmetros")
            ax.set_ylabel("Índice de Sensibilidade")
            ax.set_title("Ranking de Importância dos Parâmetros")
            ax.set_xticks(x)
            ax.set_xticklabels(param_names, rotation=45, ha="right")
            ax.legend()
            ax.grid(True, alpha=0.3, axis="y")

            plt.tight_layout()

            if save_path:
                plt.savefig(save_path, dpi=300, bbox_inches="tight")
                logger.info("Gráfico salvo em %s", save_path)

            return fig


def plot_optimization_history(
    result: OptimizationResult, save_path: str | None = None, interactive: bool = True
):
    """Função conveniente para plotar histórico de otimização."""
    visualizer = OptimizationVisualizer(result)
    return visualizer.plot_optimization_history(save_path, interactive)


def plot_parameter_importance(
    result: SensitivityResult, save_path: str | None = None, interactive: bool = True
):
    """Função conveniente para plotar importância dos parâmetros."""
    visualizer = SensitivityVisualizer(result)
    return visualizer.plot_sensitivity_indices(save_path, interactive)


def save_visualization(fig, save_path: str, format: str = "html"):
    """
    Salva visualização em arquivo.

    Args:
        fig: Figura (Plotly ou Matplotlib)
        save_path: Caminho para salvar
        format: Formato do arquivo (html, png, pdf)
    """

    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    if hasattr(fig, "write_html"):  # Plotly
        if format == "html":
            fig.write_html(str(path))
        elif format == "png":
            fig.write_image(str(path))
        elif format == "pdf":
            fig.write_image(str(path))
    else:  # Matplotlib
        fig.savefig(str(path), dpi=300, bbox_inches="tight")

    logger.info("Visualização salva em %s", path)
