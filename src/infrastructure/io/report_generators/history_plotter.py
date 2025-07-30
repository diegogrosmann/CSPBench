"""
History Algorithm Graphics Generator and Export

Responsible for generating specific graphics from algorithm execution history
and exporting raw data, working with generic data independent of algorithm type.
"""

import csv
import json
from pathlib import Path
from typing import Any, Dict, List

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


class HistoryPlotter:
    """
    Graphics generator and data export for algorithm history.
    Works with generic data independent of algorithm type.
    """

    def __init__(self, config: Dict[str, Any], output_dir: Path):
        """
        Initialize the history graphics generator.

        Args:
            config: System configuration
            output_dir: Directory where to save graphics
        """
        self.config = config
        self.output_dir = output_dir
        self.history_config = config.get("infrastructure", {}).get("history", {})
        self.plot_config = config.get("plots", {})
        self.plot_format = self.plot_config.get("formats", ["png"])[0] if self.plot_config.get("formats") else "png"

        # Configure matplotlib
        plt.style.use("seaborn-v0_8")
        sns.set_palette("husl")

    def generate_history_plots(self, experiment_results: List[Dict[str, Any]]) -> None:
        """
        Generate all configured history plots and export raw data.

        Args:
            experiment_results: List of experiment results
        """
        if not self.history_config.get("generate_plots", False):
            return

        print("📈 Generating history plots and export...")

        # Create history directory
        history_dir = self.output_dir / "history"
        history_dir.mkdir(exist_ok=True)

        # Process history data
        history_data = self._extract_history_data(experiment_results)

        if not history_data:
            print("⚠️ No history data found")
            return

        # Export raw data first
        self._export_raw_history_data(history_data, history_dir)

        print(f"📊 History processed and saved in: {history_dir}")

    def _extract_history_data(
        self, experiment_results: List[Dict[str, Any]]
    ) -> Dict[str, Dict[str, Any]]:
        """Extract history data from experiment results."""
        history_data = {}

        for i, result in enumerate(experiment_results):
            metadata = result.get("metadata", {})
            history = metadata.get("history", [])
            algorithm = result.get("algorithm", f"Unknown_{i}")

            if history:
                history_data[f"{algorithm}_{i}"] = {
                    "algorithm": algorithm,
                    "experiment_id": i,
                    "history": history,
                    "final_fitness": result.get("max_distance", 0),
                    "runtime": result.get("execution_time", 0),
                }

        return history_data

    def _export_raw_history_data(
        self, history_data: Dict[str, Dict[str, Any]], history_dir: Path
    ) -> None:
        """Export raw history data in CSV and JSON formats."""
        print("💾 Exporting raw history data...")

        # Export complete JSON
        json_file = history_dir / "history_raw_data.json"
        with open(json_file, "w", encoding="utf-8") as f:
            json.dump(history_data, f, indent=2, ensure_ascii=False, default=str)

        # Prepare data for CSV (tabular format)
        csv_rows = []
        for exp_key, exp_data in history_data.items():
            base_info = {
                "experiment_key": exp_key,
                "algorithm": exp_data.get("algorithm", "Unknown"),
                "experiment_id": exp_data.get("experiment_id", 0),
                "final_fitness": exp_data.get("final_fitness", 0),
                "runtime": exp_data.get("runtime", 0),
            }

            # Process each history entry
            for entry in exp_data.get("history", []):
                row = base_info.copy()

                # Add generic entry data
                for key, value in entry.items():
                    if isinstance(value, (dict, list)):
                        row[f"history_{key}"] = json.dumps(value)
                    else:
                        row[f"history_{key}"] = value

                csv_rows.append(row)

        # Export detailed CSV
        if csv_rows:
            csv_file = history_dir / "history_data.csv"
            fieldnames = set()
            for row in csv_rows:
                fieldnames.update(row.keys())

            with open(csv_file, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=sorted(fieldnames))
                writer.writeheader()
                writer.writerows(csv_rows)

        # Exportar summary CSV (apenas métricas principais)
        summary_rows = []
        for exp_key, exp_data in history_data.items():
            summary_row = {
                "experiment_key": exp_key,
                "algorithm": exp_data.get("algorithm", "Unknown"),
                "experiment_id": exp_data.get("experiment_id", 0),
                "final_fitness": exp_data.get("final_fitness", 0),
                "runtime": exp_data.get("runtime", 0),
                "history_entries": len(exp_data.get("history", [])),
            }

            # Calcular estatísticas básicas se houver dados de fitness
            history = exp_data.get("history", [])
            fitness_values = [
                entry.get("best_fitness")
                for entry in history
                if "best_fitness" in entry
            ]

            if fitness_values:
                summary_row.update(
                    {
                        "initial_fitness": (
                            fitness_values[0] if fitness_values else None
                        ),
                        "best_fitness_achieved": max(fitness_values),
                        "fitness_improvement": (
                            fitness_values[0] - fitness_values[-1]
                            if len(fitness_values) > 1
                            else 0
                        ),
                        "convergence_iterations": len(fitness_values),
                    }
                )

            summary_rows.append(summary_row)

        if summary_rows:
            summary_file = history_dir / "history_summary.csv"
            with open(summary_file, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=summary_rows[0].keys())
                writer.writeheader()
                writer.writerows(summary_rows)

        print(f"✅ Dados exportados: JSON, CSV detalhado e resumo")
