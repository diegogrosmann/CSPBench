"""
Exportadores de dados para diferentes formatos.

Implementa interfaces para exportar resultados seguindo o padrão Port/Adapter.
"""

from .file_exporter import FileExporter
from .json_exporter import JsonExporter
from .csv_exporter import CsvExporter
from .txt_exporter import TxtExporter

__all__ = [
    "FileExporter",
    "JsonExporter", 
    "CsvExporter",
    "TxtExporter",
]
