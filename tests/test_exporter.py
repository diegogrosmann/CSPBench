"""
Testes unitários para CSPExporter.

Este módulo contém testes para validar o sistema de exportação
de resultados CSP para diferentes formatos.
"""

import json
import tempfile
from pathlib import Path
from unittest import mock

import pytest

from src.core.io.exporter import CSPExporter


class TestCSPExporter:
    """Testes para CSPExporter."""

    def setup_method(self):
        """Setup para cada teste."""
        self.exporter = CSPExporter()
        self.sample_results = {
            "baseline": [
                {"center": "ATCG", "distance": 2.5, "time": 1.2},
                {"center": "TACG", "distance": 2.8, "time": 1.1},
            ],
            "blf_ga": [
                {"center": "ATCG", "distance": 2.0, "time": 0.8},
                {"center": "TACG", "distance": 2.1, "time": 0.9},
            ],
        }

    def test_export_to_csv_basic(self):
        """Test basic CSV export functionality."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            temp_path = f.name

        try:
            self.exporter.export_to_csv(self.sample_results, temp_path)

            # Check file exists
            assert Path(temp_path).exists()

            # Check content
            with open(temp_path) as f:
                content = f.read()
                assert "algoritmo,execucao,melhor_string,distancia,tempo" in content
                assert "baseline" in content
                assert "blf_ga" in content

        finally:
            Path(temp_path).unlink(missing_ok=True)

    def test_export_to_csv_with_extra_info(self):
        """Test CSV export with extra information."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            temp_path = f.name

        try:
            extra_info = {
                "params": {
                    "n": 10,
                    "L": 100,
                    "alphabet": "ATCG",
                }
            }
            self.exporter.export_to_csv(self.sample_results, temp_path, extra_info=extra_info)

            # Check file exists
            assert Path(temp_path).exists()

            # Check content includes extra info
            with open(temp_path) as f:
                content = f.read()
                assert "10" in content  # n parameter
                assert "100" in content  # L parameter
                assert "ATCG" in content  # alphabet parameter

        finally:
            Path(temp_path).unlink(missing_ok=True)

    def test_export_to_csv_empty_results(self):
        """Test CSV export with empty results."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            temp_path = f.name

        try:
            self.exporter.export_to_csv({}, temp_path)

            # Check file exists
            assert Path(temp_path).exists()

            # Check content has at least headers
            with open(temp_path) as f:
                content = f.read()
                assert "algoritmo,execucao,melhor_string,distancia,tempo" in content

        finally:
            Path(temp_path).unlink(missing_ok=True)

    def test_export_to_csv_create_directories(self):
        """Test CSV export creates directories if they don't exist."""
        with tempfile.TemporaryDirectory() as tmpdir:
            temp_path = Path(tmpdir) / "subdir" / "test.csv"

            self.exporter.export_to_csv(self.sample_results, str(temp_path))

            # Check file exists
            assert temp_path.exists()

    def test_export_to_json_basic(self):
        """Test basic JSON export functionality."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            temp_path = f.name

        try:
            self.exporter.export_to_json(self.sample_results, temp_path)

            # Check file exists
            assert Path(temp_path).exists()

            # Check content
            with open(temp_path) as f:
                loaded_data = json.load(f)
                assert "baseline" in loaded_data
                assert "blf_ga" in loaded_data

        finally:
            Path(temp_path).unlink(missing_ok=True)

    def test_export_to_json_without_key_conversion(self):
        """Test JSON export without key conversion."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            temp_path = f.name

        try:
            data = {"1": "test", "2": "test2"}
            self.exporter.export_to_json(data, temp_path, convert_keys=False)

            # Check file exists
            assert Path(temp_path).exists()

            # Check content
            with open(temp_path) as f:
                loaded_data = json.load(f)
                assert loaded_data == {"1": "test", "2": "test2"}

        finally:
            Path(temp_path).unlink(missing_ok=True)

    def test_export_batch_json_to_csv(self):
        """Test batch JSON to CSV export."""
        with tempfile.TemporaryDirectory() as tmpdir:
            json_path = Path(tmpdir) / "batch.json"
            csv_path = Path(tmpdir) / "batch.csv"

            # Create sample batch JSON with expected structure
            batch_data = {
                "execucoes": [
                    {
                        "config_nome": "config1",
                        "bases_info": [
                            {
                                "base_idx": 0,
                                "params": {"n": 10, "L": 100, "alphabet": "ATCG"},
                            }
                        ],
                        "algoritmos_executados": {
                            "baseline": {
                                "base_0": {
                                    "execucoes_detalhadas": [
                                        {"distancia": 2.5, "tempo": 1.2},
                                        {"distancia": 2.8, "tempo": 1.1},
                                    ],
                                    "dist_base": 3.0,
                                }
                            }
                        },
                    }
                ]
            }

            with open(json_path, "w") as f:
                json.dump(batch_data, f)

            self.exporter.export_batch_json_to_csv(str(json_path), str(csv_path))

            # Check file exists
            assert csv_path.exists()

            # Check content
            with open(csv_path) as f:
                content = f.read()
                assert "config_nome" in content
                assert "algoritmo" in content
                assert "config1" in content
                assert "baseline" in content

    def test_export_batch_json_to_csv_missing_file(self):
        """Test batch JSON to CSV export with missing file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            json_path = Path(tmpdir) / "missing.json"
            csv_path = Path(tmpdir) / "batch.csv"

            with pytest.raises(FileNotFoundError):
                self.exporter.export_batch_json_to_csv(str(json_path), str(csv_path))

    def test_create_batch_row(self):
        """Test creating batch row from execution data."""
        exec_data = {"distancia": 2.5, "tempo": 1.2}
        result = {"dist": 2.5, "time": 1.2}
        base_params = {"dist_base": 3.0, "n": 10, "L": 100}

        row = self.exporter._create_batch_row("config1", "baseline", "0", 1, exec_data, result, base_params)

        assert row["config_nome"] == "config1"
        assert row["algoritmo"] == "baseline"
        assert row["base_idx"] == "0"
        assert row["execucao_idx"] == 1
        assert row["dist"] == 2.5
        assert row["tempo"] == 1.2
        assert row["n"] == 10
        assert row["L"] == 100

    def test_create_batch_row_missing_fields(self):
        """Test creating batch row with missing fields."""
        exec_data = {"distancia": 2.5}  # Missing tempo
        result = {"dist": 2.5}
        base_params = {}  # Missing other fields

        row = self.exporter._create_batch_row("config1", "baseline", "0", 1, exec_data, result, base_params)

        assert row["config_nome"] == "config1"
        assert row["algoritmo"] == "baseline"
        assert row["base_idx"] == "0"
        assert row["execucao_idx"] == 1
        assert row["dist"] == 2.5
        assert row["tempo"] == ""  # Default
        assert row["n"] == ""  # Default
        assert row["L"] == ""  # Default

    def test_write_batch_csv(self):
        """Test writing batch data to CSV."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            temp_path = f.name

        try:
            rows = [
                {
                    "config_nome": "config1",
                    "algoritmo": "baseline",
                    "base_idx": "0",
                    "execucao_idx": 1,
                    "dist": 2.5,
                    "tempo": 1.2,
                    "dist_base": 3.0,
                }
            ]

            self.exporter._write_batch_csv(temp_path, rows)

            # Check file exists
            assert Path(temp_path).exists()

            # Check content
            with open(temp_path) as f:
                content = f.read()
                assert "config_nome" in content
                assert "algoritmo" in content
                assert "config1" in content
                assert "baseline" in content

        finally:
            Path(temp_path).unlink(missing_ok=True)

    def test_write_batch_csv_empty_rows(self):
        """Test writing empty batch data to CSV."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            temp_path = f.name

        try:
            self.exporter._write_batch_csv(temp_path, [])

            # Check file exists
            assert Path(temp_path).exists()

            # Check content has at least headers
            with open(temp_path) as f:
                content = f.read()
                assert "config_nome" in content
                assert "algoritmo" in content

        finally:
            Path(temp_path).unlink(missing_ok=True)

    def test_convert_dict_keys_to_str(self):
        """Test converting dictionary keys to strings."""
        data = {
            1: "test",
            2: {"nested": {3: "value"}},
            "string_key": "value",
        }

        result = self.exporter._convert_dict_keys_to_str(data)

        assert "1" in result
        assert "2" in result
        assert "string_key" in result
        assert "3" in result["2"]["nested"]

    def test_convert_dict_keys_to_str_non_dict(self):
        """Test converting non-dict values."""
        data = "not_a_dict"
        result = self.exporter._convert_dict_keys_to_str(data)
        assert result == "not_a_dict"

    def test_convert_dict_keys_to_str_list(self):
        """Test converting list with dict elements."""
        data = [{"1": "value1"}, {"2": "value2"}]
        result = self.exporter._convert_dict_keys_to_str(data)
        assert result == [{"1": "value1"}, {"2": "value2"}]

    @mock.patch("src.core.io.exporter.Path.mkdir")
    def test_export_to_csv_mkdir_error(self, mock_mkdir):
        """Test CSV export with directory creation error."""
        mock_mkdir.side_effect = OSError("Permission denied")

        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            temp_path = f.name

        try:
            with pytest.raises(OSError):
                self.exporter.export_to_csv(self.sample_results, temp_path)
        finally:
            Path(temp_path).unlink(missing_ok=True)

    @mock.patch("builtins.open", side_effect=OSError("Write error"))
    def test_export_to_csv_write_error(self, mock_open):
        """Test CSV export with write error."""
        with pytest.raises(OSError):
            self.exporter.export_to_csv(self.sample_results, "test.csv")

    @mock.patch("builtins.open", side_effect=OSError("Write error"))
    def test_export_to_json_write_error(self, mock_open):
        """Test JSON export with write error."""
        with pytest.raises(OSError):
            self.exporter.export_to_json(self.sample_results, "test.json")

    def test_export_complex_batch_structure(self):
        """Test complex batch structure export."""
        with tempfile.TemporaryDirectory() as tmpdir:
            json_path = Path(tmpdir) / "complex_batch.json"
            csv_path = Path(tmpdir) / "complex_batch.csv"

            # Create complex batch data with correct structure
            complex_batch = {
                "execucoes": [
                    {
                        "config_nome": "config1",
                        "bases_info": [
                            {
                                "base_idx": 0,
                                "params": {"n": 10, "L": 100, "alphabet": "ATCG"},
                            },
                            {
                                "base_idx": 1,
                                "params": {"n": 10, "L": 100, "alphabet": "ATCG"},
                            },
                        ],
                        "algoritmos_executados": {
                            "baseline": {
                                "base_0": {
                                    "execucoes_detalhadas": [
                                        {"distancia": 2.5, "tempo": 1.2},
                                        {"distancia": 2.8, "tempo": 1.1},
                                    ],
                                    "dist_base": 3.0,
                                },
                                "base_1": {
                                    "execucoes_detalhadas": [
                                        {"distancia": 2.3, "tempo": 1.0},
                                    ],
                                    "dist_base": 2.8,
                                },
                            },
                            "blf_ga": {
                                "base_0": {
                                    "execucoes_detalhadas": [
                                        {"distancia": 2.0, "tempo": 0.8},
                                        {"distancia": 2.1, "tempo": 0.9},
                                    ],
                                    "dist_base": 3.0,
                                },
                            },
                        },
                    },
                    {
                        "config_nome": "config2",
                        "bases_info": [
                            {
                                "base_idx": 0,
                                "params": {"n": 10, "L": 100, "alphabet": "ATCG"},
                            },
                        ],
                        "algoritmos_executados": {
                            "baseline": {
                                "base_0": {
                                    "execucoes_detalhadas": [
                                        {"distancia": 2.7, "tempo": 1.3},
                                    ],
                                    "dist_base": 3.2,
                                },
                            },
                        },
                    },
                ]
            }

            with open(json_path, "w") as f:
                json.dump(complex_batch, f)

            self.exporter.export_batch_json_to_csv(str(json_path), str(csv_path))

            # Check file exists
            assert csv_path.exists()

            # Check content for all configs and algorithms
            with open(csv_path) as f:
                content = f.read()
                assert "config1" in content
                assert "config2" in content
                assert "baseline" in content
                assert "blf_ga" in content
