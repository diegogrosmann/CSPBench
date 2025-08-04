"""
Base file exporter.

Implements the ExportPort interface and provides common functionality for all exporters.
"""

import csv
import json
import pickle
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
        self,
        results: Any,
        format_type: str,
        destination: Optional[str] = None,
        options: Optional[Dict[str, Any]] = None,
        session_id: Optional[str] = None,
    ) -> str:
        """Export optimization results to the specified format."""
        if session_id:
            filename = f"results_{session_id}.{format_type.lower()}"
        else:
            from datetime import datetime
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"results_{timestamp}.{format_type.lower()}"
            
        # Determine output path
        if destination:
            dest_path = self.output_path / destination
            if dest_path.is_dir() or not dest_path.suffix:
                dest_path = dest_path / filename
        else:
            dest_path = self.output_path / filename
            
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Export based on format
        if format_type.lower() == "json":
            self._write_json(results, dest_path, options)
        elif format_type.lower() == "csv":
            # For CSV, handle both single results and batch results
            if isinstance(results, list):
                self._write_csv(results, dest_path, options)
            else:
                # Convert single result to list for CSV export
                self._write_csv([results], dest_path, options)
        elif format_type.lower() == "parquet":
            self._write_parquet(results, dest_path)
        elif format_type.lower() == "pickle":
            self._write_pickle(results, dest_path)
        elif format_type.lower() == "txt":
            self._write_txt(results, dest_path)
        else:
            self._write_json(results, dest_path, options)
            
        return str(dest_path)

    def export_batch_results(
        self,
        batch_results: List[Dict[str, Any]],
        format_type: str,
        destination: str,
        options: Optional[Dict[str, Any]] = None,
        session_id: Optional[str] = None,
    ) -> str:
        """Export batch results."""
        dest_path = self.output_path / destination

        # If destination is a directory, use session_id or timestamp for filename
        if dest_path.is_dir() or not dest_path.suffix:
            if session_id:
                filename = f"results_{session_id}.{format_type.lower()}"
            else:
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
            self._write_json(batch_data, dest_path, options)
        elif format_type.lower() == "csv":
            self._write_csv(batch_results, dest_path, options)  # For CSV, only results
        elif format_type.lower() == "parquet":
            self._write_parquet(batch_data, dest_path)
        elif format_type.lower() == "pickle":
            self._write_pickle(batch_data, dest_path)
        elif format_type.lower() == "txt":
            self._write_txt(batch_data, dest_path)
        else:
            self._write_json(batch_data, dest_path, options)

        return str(dest_path)

    def get_supported_formats(self) -> List[str]:
        """List supported formats."""
        return ["json", "csv", "parquet", "pickle", "txt"]

    def export_optimization_results(
        self, optimization_data: Dict[str, Any], destination: str, session_id: Optional[str] = None
    ) -> str:
        """Export optimization results."""
        return self.export_results(optimization_data, "json", destination, session_id=session_id)

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

    def _write_json(
        self, data: Any, dest_path: Path, options: Optional[Dict[str, Any]] = None
    ) -> None:
        """Write data in JSON format."""
        # Get JSON format options
        json_options = options.get("json", {}) if options else {}
        indent = json_options.get("indent", 2)
        ensure_ascii = json_options.get("ensure_ascii", False)

        with open(dest_path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=indent, ensure_ascii=ensure_ascii, default=str)

    def _write_csv(
        self, data: Any, dest_path: Path, options: Optional[Dict[str, Any]] = None
    ) -> None:
        """Write data in CSV format."""
        # Get CSV format options
        csv_options = options.get("csv", {}) if options else {}
        separator = csv_options.get("separator", ",")
        encoding = csv_options.get("encoding", "utf-8")
        decimal = csv_options.get("decimal", ".")

        if isinstance(data, list) and data:
            # Assume list of dictionaries
            fieldnames = set()
            for item in data:
                if isinstance(item, dict):
                    fieldnames.update(item.keys())

            fieldnames = sorted(fieldnames)

            with open(dest_path, "w", newline="", encoding=encoding) as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter=separator)
                writer.writeheader()
                for item in data:
                    if isinstance(item, dict):
                        # Convert decimal separator if needed
                        if decimal != ".":
                            formatted_item = {}
                            for k, v in item.items():
                                if isinstance(v, float):
                                    formatted_item[k] = str(v).replace(".", decimal)
                                else:
                                    formatted_item[k] = v
                            writer.writerow(formatted_item)
                        else:
                            writer.writerow(item)
        elif isinstance(data, dict):
            # Try flatten for CSV
            with open(dest_path, "w", newline="", encoding=encoding) as f:
                writer = csv.writer(f, delimiter=separator)
                for key, value in data.items():
                    if isinstance(value, float) and decimal != ".":
                        value = str(value).replace(".", decimal)
                    writer.writerow([key, value])
        else:
            # Fallback: convert to string
            with open(dest_path, "w", encoding=encoding) as f:
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

    def _write_parquet(self, data: Any, dest_path: Path) -> None:
        """Write data in Parquet format."""
        try:
            import pandas as pd

            # Convert data to DataFrame
            if isinstance(data, list) and data:
                df = pd.DataFrame(data)
            elif isinstance(data, dict):
                # Flatten dict or create single-row DataFrame
                if all(
                    isinstance(v, (list, tuple))
                    and len(
                        set(
                            len(v) if hasattr(v, "__len__") else 1
                            for v in data.values()
                        )
                    )
                    == 1
                    for v in data.values()
                ):
                    # All values are lists of same length
                    df = pd.DataFrame(data)
                else:
                    # Mixed data - create single row
                    df = pd.DataFrame([data])
            else:
                # Convert to single-column DataFrame
                df = pd.DataFrame({"data": [data]})

            df.to_parquet(dest_path, index=False)

        except ImportError:
            # Fallback to JSON if pandas/pyarrow not available
            self._write_json(data, dest_path.with_suffix(".json"))
        except Exception:
            # Fallback to JSON on any error
            self._write_json(data, dest_path.with_suffix(".json"))

    def _write_pickle(self, data: Any, dest_path: Path) -> None:
        """Write data in Pickle format."""
        try:
            # First attempt direct serialization
            with open(dest_path, "wb") as f:
                pickle.dump(data, f)
        except (TypeError, AttributeError) as e:
            if "cannot pickle" in str(e):
                # Sanitize data by removing non-serializable objects
                sanitized_data = self._sanitize_for_pickle(data)
                with open(dest_path, "wb") as f:
                    pickle.dump(sanitized_data, f)
            else:
                raise

    def _sanitize_for_pickle(self, data: Any) -> Any:
        """
        Recursively sanitize data for pickle serialization.

        Removes or converts non-serializable objects like module references,
        function objects, and other problematic types.
        """
        import types

        if isinstance(data, dict):
            sanitized = {}
            for key, value in data.items():
                try:
                    # Test if the value can be pickled
                    pickle.dumps(value)
                    sanitized[key] = self._sanitize_for_pickle(value)
                except (TypeError, AttributeError):
                    # Convert non-serializable objects to string representation
                    if isinstance(value, types.ModuleType):
                        sanitized[key] = f"<module '{value.__name__}'>"
                    elif callable(value):
                        sanitized[key] = (
                            f"<function '{getattr(value, '__name__', str(value))}'>"
                        )
                    elif hasattr(value, "__dict__") and hasattr(value, "__class__"):
                        # For complex objects, try to extract serializable attributes
                        sanitized[key] = {
                            "__class__": value.__class__.__name__,
                            "__repr__": str(value),
                        }
                    else:
                        sanitized[key] = str(value)
            return sanitized
        elif isinstance(data, (list, tuple)):
            sanitized_items = []
            for item in data:
                try:
                    pickle.dumps(item)
                    sanitized_items.append(self._sanitize_for_pickle(item))
                except (TypeError, AttributeError):
                    sanitized_items.append(str(item))
            return type(data)(sanitized_items)
        else:
            # For primitive types, return as-is
            return data
