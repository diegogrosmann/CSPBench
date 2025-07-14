"""
Exportador base para arquivos.

Implementa a interface ExportPort e fornece funcionalidade comum para todos os exportadores.
"""

import csv
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from src.application.ports import ExportPort


class FileExporter(ExportPort):
    """Exportador base para arquivos."""

    def __init__(self, output_path: str):
        """Inicializa exportador com caminho de saída."""
        self.output_path = Path(output_path)
        self.output_path.mkdir(parents=True, exist_ok=True)

    def export_results(
        self, results: Dict[str, Any], format_type: str, destination: str
    ) -> str:
        """Exporta resultados em formato específico."""
        dest_path = self.output_path / destination

        # Se destination é um diretório, gerar nome de arquivo
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
            # Default para JSON
            self._write_json(results, dest_path)

        return str(dest_path)

    def export_batch_results(
        self, batch_results: List[Dict[str, Any]], format_type: str, destination: str
    ) -> str:
        """Exporta resultados de batch."""
        dest_path = self.output_path / destination

        # Se destination é um diretório, usar nome de arquivo fixo
        if dest_path.is_dir() or not dest_path.suffix:
            # Usar nome fixo para resultado principal
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
            self._write_csv(batch_results, dest_path)  # Para CSV, só os resultados
        elif format_type.lower() == "txt":
            self._write_txt(batch_data, dest_path)
        else:
            self._write_json(batch_data, dest_path)

        return str(dest_path)

    def get_supported_formats(self) -> List[str]:
        """Lista formatos suportados."""
        return ["json", "csv", "txt"]

    def export_optimization_results(
        self, optimization_data: Dict[str, Any], destination: str
    ) -> str:
        """Exporta resultados de otimização."""
        return self.export_results(optimization_data, "json", destination)

    def export(self, data: Any, filename: Optional[str] = None) -> None:
        """Exporta dados para arquivo (método legacy)."""
        if filename:
            file_path = self.output_path / filename
        else:
            from datetime import datetime

            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            file_path = self.output_path / f"export_{timestamp}.json"

        file_path.parent.mkdir(parents=True, exist_ok=True)
        self._write_json(data, file_path)

    def _write_json(self, data: Any, dest_path: Path) -> None:
        """Escreve dados em formato JSON."""
        with open(dest_path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, ensure_ascii=False, default=str)

    def _write_csv(self, data: Any, dest_path: Path) -> None:
        """Escreve dados em formato CSV."""
        if isinstance(data, list) and data:
            # Assumir lista de dicionários
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
            # Tentar flatten para CSV
            with open(dest_path, "w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                for key, value in data.items():
                    writer.writerow([key, value])
        else:
            # Fallback: converter para string
            with open(dest_path, "w", encoding="utf-8") as f:
                f.write(str(data))

    def _write_txt(self, data: Any, dest_path: Path) -> None:
        """Escreve dados em formato texto."""
        with open(dest_path, "w", encoding="utf-8") as f:
            if isinstance(data, (list, tuple)):
                for item in data:
                    f.write(f"{item}\n")
            elif isinstance(data, dict):
                for key, value in data.items():
                    f.write(f"{key}: {value}\n")
            else:
                f.write(str(data))
