"""
Exportador para formato TXT.

Especialização do FileExporter para formato texto simples.
"""

from .file_exporter import FileExporter


class TxtExporter(FileExporter):
    """Exportador para formato TXT."""

    def get_supported_formats(self) -> list[str]:
        """Lista formatos suportados."""
        return ["txt"]
