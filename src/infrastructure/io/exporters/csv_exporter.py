"""
Exportador para formato CSV.

Especialização do FileExporter para formato CSV simples.
"""

from .file_exporter import FileExporter


class CsvExporter(FileExporter):
    """Exportador para formato CSV."""

    def get_supported_formats(self) -> list[str]:
        """Lista formatos suportados."""
        return ["csv"]
