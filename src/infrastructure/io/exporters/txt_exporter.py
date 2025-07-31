"""
TXT Format Exporter.

FileExporter specialization for simple text format.
"""

from .file_exporter import FileExporter


class TxtExporter(FileExporter):
    """TXT format exporter."""

    def get_supported_formats(self) -> list[str]:
        """List supported formats."""
        return ["txt"]
