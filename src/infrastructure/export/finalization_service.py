"""
Finalization Service for Comprehensive Export Operations.

Provides finalization services for complete work export including raw database
dumps, manifest generation, full results compilation, and summary generation.
This service handles the final phase of work execution by creating comprehensive
export artifacts.

Features:
    - Raw database table exports to CSV format
    - Full results JSON compilation with metadata
    - Manifest generation with artifact inventory
    - Summary report generation in Markdown format
    - Optimization results export and visualization
    - Sensitivity analysis results export and plotting

Export Structure:
    - raw_db/: CSV exports of database tables
    - optuna/: Optimization trials and plots
    - sensitivity/native/: Sensitivity analysis results and visualizations
    - full_results.json: Complete results compilation
    - manifest.json: Export inventory and metadata
    - summary.md: Human-readable summary report

This implementation intentionally drops legacy/export-config driven outputs
in favor of a standardized export format.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

import matplotlib
import tomllib

matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt

try:  # Optuna is optional
    import optuna  # type: ignore
except Exception:  # pragma: no cover - absence should not break
    optuna = None

from src.infrastructure.logging_config import get_logger

logger = get_logger("CSPBench.Finalization")


@dataclass
class FinalizationConfig:
    """
    Configuration for finalization service operations.

    Attributes:
        work_id (str): Unique work identifier for export.
        output_dir (Path): Base directory for export artifacts.
        tool_version (Optional[str]): Version of the tool generating exports.
    """

    work_id: str
    output_dir: Path
    tool_version: str | None = None


class FinalizationService:
    """
    Comprehensive finalization and export service.

    Handles the complete finalization process for work items including
    database exports, result compilation, visualization generation, and
    artifact manifest creation.

    Features:
        - Database table exports in CSV format
        - Optimization trial exports with plots
        - Sensitivity analysis exports with visualizations
        - Full results JSON compilation
        - Export manifest with metadata
        - Summary report generation

    Attributes:
        config (FinalizationConfig): Service configuration.
        work_store: Work-scoped persistence store.
        output_dir (Path): Base output directory.
        raw_dir (Path): Raw database exports directory.
        optuna_dir (Path): Optimization artifacts directory.
        sensitivity_dir (Path): Sensitivity analysis artifacts directory.
    """

    def __init__(self, config: FinalizationConfig, work_store: Any = None):
        """
        Initialize finalization service with configuration and storage.

        Args:
            config (FinalizationConfig): Service configuration.
            work_store: WorkScopedPersistence instance for data access.
        """
        self.config = config
        self.work_store = work_store  # WorkScopedPersistence instance
        self.output_dir = config.output_dir
        self.raw_dir = self.output_dir / "raw_db"
        self.optuna_dir = self.output_dir / "optuna"
        self.sensitivity_dir = self.output_dir / "sensitivity" / "native"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.raw_dir.mkdir(exist_ok=True)

    def _get_db_connection(self):
        """
        Get database connection (DEPRECATED).

        Direct database access should not be used. Use work_store methods instead.

        Raises:
            ValueError: If work_store is not available.
        """
        logger.warning(
            "Direct database access is deprecated. Use work_store methods instead."
        )
        if self.work_store is None:
            raise ValueError("work_store is required for database operations")
        # Get the raw database connection from SQLAlchemy engine
        engine = self.work_store.store._engine
        return engine.raw_connection()

    # -------------- PUBLIC --------------
    def run(self) -> None:
        """
        Execute the complete finalization process.

        Performs all finalization operations including database exports,
        specialized exports (Optuna/Sensitivity), result compilation,
        manifest generation, and summary creation.
        """
        logger.info(
            "Starting finalization export (raw db + manifest + full_results + summary)"
        )
        tables_info = self._dump_sqlite_tables()

        # Exportações especiais (Optuna / Sensitivity)
        optuna_exports = self._export_optuna_trials_and_plots()
        sensitivity_exports = self._export_sensitivity_results_and_plots()

        full_results_path = self._generate_full_results(
            tables_info, optuna_exports, sensitivity_exports
        )
        manifest_path = self._generate_manifest(
            tables_info, full_results_path, optuna_exports, sensitivity_exports
        )
        summary_path = self._generate_summary(
            tables_info, full_results_path, manifest_path
        )
        logger.info(
            f"Finalization complete: manifest={manifest_path.name}, full_results={full_results_path.name}, summary={summary_path.name}"
        )

    # -------------- INTERNAL --------------
    def _dump_sqlite_tables(self) -> Dict[str, Dict[str, Any]]:
        """
        Export database tables using work_store methods instead of direct SQL.

        Exports all relevant database tables to CSV format using the work_store
        interface rather than direct database access for better abstraction.

        Returns:
            Dict[str, Dict[str, Any]]: Dictionary containing export metadata
                for each table including row counts and file paths.
        """
        info: Dict[str, Dict[str, Any]] = {}

        if self.work_store is None:
            raise ValueError("work_store is required for export operations")

        try:
            # Export work data
            work_data = self.work_store.store.get_work_export_data(self.config.work_id)
            if work_data:
                self._write_csv_data("work", [work_data], info)

            # Export combinations
            combinations_data = self.work_store.store.get_combinations_for_export(
                self.config.work_id
            )
            if combinations_data:
                self._write_csv_data("combinations", combinations_data, info)

            # Export executions
            executions_data = self.work_store.store.get_executions_for_export(
                self.config.work_id
            )
            if executions_data:
                self._write_csv_data("executions", executions_data, info)

            # Export execution progress
            progress_data = self.work_store.store.get_execution_progress_for_export(
                self.config.work_id
            )
            if progress_data:
                self._write_csv_data("execution_progress", progress_data, info)

            # Export events
            events_data = self.work_store.store.get_events_for_export(
                self.config.work_id
            )
            if events_data:
                self._write_csv_data("events", events_data, info)

            # Export datasets
            datasets_data = self.work_store.store.get_datasets_for_export(
                self.config.work_id
            )
            if datasets_data:
                self._write_csv_data("datasets", datasets_data, info)

            # Export dataset sequences
            sequences_data = self.work_store.store.get_dataset_sequences_for_export(
                self.config.work_id
            )
            if sequences_data:
                self._write_csv_data("dataset_sequences", sequences_data, info)

        except Exception as e:
            logger.error(f"Error dumping tables: {e}")

        return info

    def _write_csv_data(
        self,
        table_name: str,
        data: List[Dict[str, Any]],
        info: Dict[str, Dict[str, Any]],
    ) -> None:
        """
        Write data to CSV file with proper formatting.

        Handles CSV writing with proper escaping for special characters
        and updates the export metadata information.

        Args:
            table_name (str): Name of the table being exported.
            data (List[Dict[str, Any]]): Data rows to export.
            info (Dict[str, Dict[str, Any]]): Export metadata to update.
        """
        if not data:
            return

        csv_path = self.raw_dir / f"{table_name}.csv"

        # Get column names from first row
        col_names = list(data[0].keys())

        # Write CSV manually (avoid pandas dependency here)
        with csv_path.open("w", encoding="utf-8") as fh:
            fh.write(",".join(col_names) + "\n")
            for row in data:
                # Basic escaping: replace newlines/commas inside str
                formatted = []
                for col in col_names:
                    cell = row.get(col)
                    if cell is None:
                        formatted.append("")
                    else:
                        text = str(cell).replace("\n", " ")
                        if "," in text:
                            text = '"' + text.replace('"', "''") + '"'
                        formatted.append(text)
                fh.write(",".join(formatted) + "\n")

        info[table_name] = {
            "rows": len(data),
            "file": str(csv_path.relative_to(self.output_dir)),
        }

    def _generate_full_results(
        self,
        tables_info: Dict[str, Dict[str, Any]],
        optuna_exports: Dict[str, Any],
        sensitivity_exports: Dict[str, Any],
    ) -> Path:
        """
        Generate comprehensive full results JSON file.

        Compiles all execution results into a single JSON file with
        metadata, table information, and specialized analysis results.

        Args:
            tables_info (Dict[str, Dict[str, Any]]): Database table export metadata.
            optuna_exports (Dict[str, Any]): Optimization export results.
            sensitivity_exports (Dict[str, Any]): Sensitivity analysis export results.

        Returns:
            Path: Path to the generated full results JSON file.
        """
        # Build executions expansion (lightweight) using work_store
        executions: List[Dict[str, Any]] = []
        try:
            # Use the export method instead of direct DB access
            executions_data = self.work_store.store.get_executions_for_export(
                self.config.work_id
            )

            for exec_row in executions_data:
                executions.append(
                    {
                        "execution_id": exec_row.get("id"),
                        "combination_id": exec_row.get("combination_id"),
                        "status": exec_row.get("status"),
                        "objective": exec_row.get("objective"),
                        "params": self._safe_parse_json(exec_row.get("params_json")),
                        "result": self._safe_parse_json(exec_row.get("result_json")),
                    }
                )
        except Exception as e:  # pragma: no cover
            logger.error(f"Error building executions section: {e}")

        data: Dict[str, Any] = {
            "work_id": self.config.work_id,
            "export_timestamp": self._now_iso(),
            "tool_version": self.config.tool_version,
            "tables": tables_info,
            "results": {
                "executions": executions,
            },
            "optimization": optuna_exports,
            "sensitivity": sensitivity_exports,
        }

        path = self.output_dir / "full_results.json"
        with path.open("w", encoding="utf-8") as fh:
            json.dump(data, fh, ensure_ascii=False, indent=2)
        return path

    def _generate_manifest(
        self,
        tables_info: Dict[str, Dict[str, Any]],
        full_results_path: Path,
        optuna_exports: Dict[str, Any],
        sensitivity_exports: Dict[str, Any],
    ) -> Path:
        """
        Generate export manifest with artifact inventory.

        Creates a comprehensive manifest file listing all export artifacts
        with metadata and version information.

        Args:
            tables_info (Dict[str, Dict[str, Any]]): Database table export metadata.
            full_results_path (Path): Path to full results file.
            optuna_exports (Dict[str, Any]): Optimization export results.
            sensitivity_exports (Dict[str, Any]): Sensitivity analysis export results.

        Returns:
            Path: Path to the generated manifest JSON file.
        """
        manifest = {
            "version": "1.0",
            "work_id": self.config.work_id,
            "tool_version": self.config.tool_version,
            "export_timestamp": self._now_iso(),
            "tables": tables_info,
            "artefacts": {
                "full_results": full_results_path.name,
                "raw_db_dir": str(self.raw_dir.relative_to(self.output_dir)),
                "summary": "summary.md",
                "optimization": optuna_exports.get("studies", []),
                "sensitivity": sensitivity_exports.get("analyses", []),
            },
        }
        path = self.output_dir / "manifest.json"
        with path.open("w", encoding="utf-8") as fh:
            json.dump(manifest, fh, ensure_ascii=False, indent=2)
        return path

    def _generate_summary(
        self,
        tables_info: Dict[str, Dict[str, Any]],
        full_results_path: Path,
        manifest_path: Path,
    ) -> Path:
        """
        Generate human-readable summary report in Markdown format.

        Creates a summary report with key statistics and artifact listings
        for easy human consumption of export results.

        Args:
            tables_info (Dict[str, Dict[str, Any]]): Database table export metadata.
            full_results_path (Path): Path to full results file.
            manifest_path (Path): Path to manifest file.

        Returns:
            Path: Path to the generated summary Markdown file.
        """
        executions_rows = tables_info.get("executions", {}).get("rows", 0)
        work_rows = tables_info.get("work", {}).get("rows", 0)
        summary_lines = [
            "# CSPBench Export Summary",
            f"Work ID: {self.config.work_id}",
            f"Export Timestamp: {self._now_iso()}",
            f"Tool Version: {self.config.tool_version or 'unknown'}",
            "",
            "## Tables",
        ]
        for table, meta in sorted(tables_info.items()):
            summary_lines.append(f"- {table}: {meta['rows']} rows -> {meta['file']}")
        summary_lines += [
            "",
            "## Artefacts",
            f"- Manifest: {manifest_path.name}",
            f"- Full Results: {full_results_path.name}",
            f"- Raw DB Directory: {self.raw_dir.name}",
            "",
            "(Optimization and sensitivity artefacts will appear when those phases run.)",
        ]
        path = self.output_dir / "summary.md"
        path.write_text("\n".join(summary_lines), encoding="utf-8")
        return path

    # -------------- OPTUNA EXPORT --------------
    def _export_optuna_trials_and_plots(self) -> Dict[str, Any]:
        """
        Export optimization trials and generate plots.

        Exports optimization trial data and generates visualizations
        without necessarily depending on Optuna storage format.

        Strategy:
            1. Read executions whose unit_id starts with 'optimization:'
            2. Group by (task_id, dataset_id, preset_id, algorithm_id)
            3. Generate consolidated trials.csv file
            4. Generate plots for each study (objective_vs_trial, best_objective_vs_trial)

        Returns:
            Dict[str, Any]: Dictionary containing export metadata including
                trials CSV path, study information, and plot paths.
        """
        try:
            # Use work_store method instead of direct DB access
            executions_data = (
                self.work_store.store.get_optimization_executions_for_export(
                    self.config.work_id
                )
            )
        except Exception as e:  # pragma: no cover
            logger.warning(f"Optuna export skipped (query error): {e}")
            return {"studies": []}

        if not executions_data:
            return {"studies": []}

        self.optuna_dir.mkdir(exist_ok=True)
        trials_csv = self.optuna_dir / "trials.csv"
        studies: Dict[str, Dict[str, Any]] = {}

        with trials_csv.open("w", encoding="utf-8") as fh:
            fh.write(
                "study_id,task_id,dataset_id,preset_id,algorithm_id,trial_number,status,objective,params_json,result_json\n"
            )
            for exec_data in executions_data:
                unit_id = exec_data.get("unit_id")
                seq = exec_data.get("sequencia")
                status = exec_data.get("status")
                params_json = exec_data.get("params_json")
                result_json = exec_data.get("result_json")
                objective = exec_data.get("objective")

                # unit pattern: optimization:task:dataset:preset:algorithm:trial
                parts = unit_id.split(":")
                if len(parts) < 6:
                    continue
                _, task_id, dataset_id, preset_id, algorithm_id, trial_no = parts[:6]
                study_key = f"{task_id}|{dataset_id}|{preset_id}|{algorithm_id}"
                study = studies.setdefault(
                    study_key,
                    {
                        "study_name": study_key,
                        "task_id": task_id,
                        "dataset_id": dataset_id,
                        "preset_id": preset_id,
                        "algorithm_id": algorithm_id,
                        "trials": [],
                        "plots": {},
                    },
                )
                trial_record = {
                    "trial_number": int(trial_no),
                    "status": status,
                    "objective": objective,
                    "params": self._safe_parse_json(params_json),
                    "result": self._safe_parse_json(result_json),
                }
                study["trials"].append(trial_record)
                fh.write(
                    f"{study_key},{task_id},{dataset_id},{preset_id},{algorithm_id},{trial_no},{status},{objective},{json.dumps(trial_record['params'])},{json.dumps(trial_record['result'])}\n"
                )

        # Gerar plots simples
        for study_key, meta in studies.items():
            trials = sorted(meta["trials"], key=lambda t: t["trial_number"])
            if not trials:
                continue
            numbers = [t["trial_number"] for t in trials]
            objectives = [t["objective"] for t in trials]
            # Filtrar None
            numbers_obj, objectives_clean = (
                zip(*[(n, o) for n, o in zip(numbers, objectives) if o is not None])
                if any(o is not None for o in objectives)
                else ([], [])
            )
            if objectives_clean:
                best_running = []
                current_best = math.inf
                for o in objectives_clean:
                    if o < current_best:
                        current_best = o
                    best_running.append(current_best)

                safe_dir = self.optuna_dir / study_key.replace("|", "__")
                safe_dir.mkdir(exist_ok=True)

                # objective_vs_trial
                plt.figure()
                plt.plot(numbers_obj, objectives_clean, marker="o")
                plt.xlabel("Trial")
                plt.ylabel("Objective")
                plt.title(f"Objective vs Trial ({study_key})")
                p1 = safe_dir / "objective_vs_trial.png"
                plt.tight_layout()
                plt.savefig(p1)
                plt.close()

                # best_objective_vs_trial
                plt.figure()
                plt.plot(numbers_obj, best_running, marker="o")
                plt.xlabel("Trial")
                plt.ylabel("Best Objective")
                plt.title(f"Best Objective vs Trial ({study_key})")
                p2 = safe_dir / "best_objective_vs_trial.png"
                plt.tight_layout()
                plt.savefig(p2)
                plt.close()

                meta["plots"] = {
                    "objective_vs_trial": str(p1.relative_to(self.output_dir)),
                    "best_objective_vs_trial": str(p2.relative_to(self.output_dir)),
                }
            meta["trial_count"] = len(trials)
        studies_list = list(studies.values())
        return {
            "trials_csv": str(trials_csv.relative_to(self.output_dir)),
            "studies": studies_list,
        }

    # -------------- SENSITIVITY EXPORT --------------
    def _export_sensitivity_results_and_plots(self) -> Dict[str, Any]:
        """
        Export sensitivity analysis results and generate plots.

        Exports sensitivity analysis results based on events and generates
        visualizations. Future enhancement will integrate native SALib support.

        Returns:
            Dict[str, Any]: Dictionary containing export metadata including
                analyses file path, analysis data, and plot paths.
        """
        try:
            # Use work_store method instead of direct DB access
            events_data = self.work_store.store.get_sensitivity_events_for_export(
                self.config.work_id
            )
        except Exception as e:  # pragma: no cover
            logger.warning(f"Sensitivity export skipped (query error): {e}")
            return {"analyses": []}

        if not events_data:
            return {"analyses": []}

        self.sensitivity_dir.mkdir(parents=True, exist_ok=True)
        analyses_csv = self.sensitivity_dir / "analyses.csv"
        analyses: List[Dict[str, Any]] = []
        with analyses_csv.open("w", encoding="utf-8") as fh:
            fh.write("unit_id,method,param,score,num_samples,analysis_timestamp\n")
            for event_data in events_data:
                entity_data_json = event_data.get("entity_data_json")
                ts = event_data.get("timestamp")
                try:
                    data = json.loads(entity_data_json)
                except Exception:
                    continue
                ctx = data.get("context", {})
                method = ctx.get("method")
                num_samples = ctx.get("num_samples")
                param_scores = ctx.get("param_scores", {})
                unit_id = data.get("unit_id")
                for param, score in param_scores.items():
                    rec = {
                        "unit_id": unit_id,
                        "method": method,
                        "param": param,
                        "score": score,
                        "num_samples": num_samples,
                        "timestamp": ts,
                    }
                    analyses.append(rec)
                    fh.write(f"{unit_id},{method},{param},{score},{num_samples},{ts}\n")

        # Gerar um plot por unit_id agregando scores
        plots: List[str] = []
        from collections import defaultdict

        grouped: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
        for a in analyses:
            grouped[a["unit_id"]].append(a)
        for unit_id, items in grouped.items():
            # bar plot dos scores
            params = [i["param"] for i in items]
            scores = [i["score"] for i in items]
            plt.figure(figsize=(max(4, len(params) * 0.6), 3))
            plt.bar(params, scores)
            plt.ylabel("Score")
            plt.xlabel("Parameter")
            plt.title(f"Sensitivity Scores ({unit_id})")
            plt.xticks(rotation=45, ha="right")
            plt.tight_layout()
            filename = f"{unit_id.replace(':', '__')}_scores.png"
            out_path = self.sensitivity_dir / filename
            plt.savefig(out_path)
            plt.close()
            plots.append(str(out_path.relative_to(self.output_dir)))

        return {
            "analyses_file": str(analyses_csv.relative_to(self.output_dir)),
            "analyses": analyses,
            "plots": plots,
        }

    @staticmethod
    def _safe_parse_json(text: str | None) -> Any:
        """
        Safely parse JSON text with fallback.

        Args:
            text (Optional[str]): JSON text to parse.

        Returns:
            Any: Parsed JSON object, original text if parsing fails, or None.
        """
        if not text:
            return None
        try:
            return json.loads(text)
        except Exception:
            return text

    @staticmethod
    def _now_iso() -> str:
        """
        Get current timestamp in ISO format.

        Returns:
            str: Current timestamp in ISO 8601 format with timezone.
        """
        return datetime.now(timezone.utc).isoformat()

    @staticmethod
    def detect_tool_version(pyproject_path: Path | None = None) -> str | None:
        """
        Attempt to read version from pyproject.toml (PEP 621).

        Tries to extract the tool version from the project configuration
        file for inclusion in export metadata.

        Args:
            pyproject_path (Optional[Path]): Path to pyproject.toml file.
                If None, attempts to resolve relative to repository root.

        Returns:
            Optional[str]: Tool version string if found, None on any failure.
        """
        try:
            if pyproject_path is None:
                # Attempt to resolve relative to repository root (three levels up from this file)
                root = Path(__file__).resolve().parents[3]
                pyproject_path = root / "pyproject.toml"
            if not pyproject_path.exists():
                return None
            with pyproject_path.open("rb") as fh:  # tomllib needs binary
                data = tomllib.load(fh)
            return data.get("project", {}).get("version")
        except Exception:  # pragma: no cover - best effort
            return None
