"""
Base file exporter.

Implements the ExportPort interface and provides common functionality for all exporters.
"""

import csv
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from src.application.ports import ExportPort


class FileExporter(ExportPort):
    """Base file exporter."""

    def __init__(self, output_path: str):
        """Initialize exporter with output path."""
        self.output_path = Path(output_path)
        self.output_path.mkdir(parents=True, exist_ok=True)

    def export_results(
        self, results: Dict[str, Any], format_type: str, destination: str
    ) -> str:
        """Export results in specific format."""
        dest_path = self.output_path / destination

        # If destination is a directory, generate filename
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
            # Default to JSON
            self._write_json(results, dest_path)

        return str(dest_path)

    def export_batch_results(
        self, batch_results: List[Dict[str, Any]], format_type: str, destination: str
    ) -> str:
        """Export batch results."""
        dest_path = self.output_path / destination

        # If destination is a directory, use fixed filename
        if dest_path.is_dir() or not dest_path.suffix:
            # Use fixed name for main result
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
            self._write_csv(batch_results, dest_path)  # For CSV, only results
        elif format_type.lower() == "txt":
            self._write_txt(batch_data, dest_path)
        else:
            self._write_json(batch_data, dest_path)

        return str(dest_path)

    def get_supported_formats(self) -> List[str]:
        """List supported formats."""
        return ["json", "csv", "txt"]

    def export_optimization_results(
        self, optimization_data: Dict[str, Any], destination: str
    ) -> str:
        """Export optimization results."""
        return self.export_results(optimization_data, "json", destination)

    def export(self, data: Any, filename: Optional[str] = None) -> None:
        """Export data to file (legacy method)."""
        if filename:
            file_path = self.output_path / filename
        else:
            from datetime import datetime

            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            file_path = self.output_path / f"export_{timestamp}.json"

        file_path.parent.mkdir(parents=True, exist_ok=True)
        self._write_json(data, file_path)

    def _write_json(self, data: Any, dest_path: Path) -> None:
        """Write data in JSON format."""
        with open(dest_path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, ensure_ascii=False, default=str)

    def _write_csv(self, data: Any, dest_path: Path) -> None:
        """Write data in CSV format."""
        if isinstance(data, list) and data:
            # Assume list of dictionaries
            fieldnames = set()
            for item in data:
                if isinstance(item, dict):
                    fieldnames.update(item.keys())

            fieldnames = sorted(fieldnames)

            with open(dest_path, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                for item in data:
                    if isinstance(item, dict):
                        writer.writerow(item)
        elif isinstance(data, dict):
            # Try flatten for CSV
            with open(dest_path, "w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                for key, value in data.items():
                    writer.writerow([key, value])
        else:
            # Fallback: convert to string
            with open(dest_path, "w", encoding="utf-8") as f:
                f.write(str(data))

    def _write_txt(self, data: Any, dest_path: Path) -> None:
        """Write data in text format."""
        with open(dest_path, "w", encoding="utf-8") as f:
            if isinstance(data, (list, tuple)):
                for item in data:
                    f.write(f"{item}\n")
            elif isinstance(data, dict):
                for key, value in data.items():
                    f.write(f"{key}: {value}\n")
            else:
                f.write(str(data))
