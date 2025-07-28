"""
CSV Format Exporter.

FileExporter specialization for simple CSV format.
"""

from .file_exporter import FileExporter


class CsvExporter(FileExporter):
    """CSV format exporter."""

    def get_supported_formats(self) -> list[str]:
        """List supported formats."""
        return ["csv"]
