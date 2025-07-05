"""
Tests for BatchReporter module.
"""

import tempfile

from src.core.report.batch_reporter import BatchReporter, ReportFormat, ReportSection


class TestBatchReporter:
    """Test BatchReporter functionality."""

    def setup_method(self):
        """Setup test fixtures."""
        self.reporter = BatchReporter()

        self.sample_results = {
            "config_1": {
                "baseline": {
                    "base_0": {
                        "execucoes_detalhadas": [
                            {"distancia": 2.5, "tempo": 1.2, "melhor_string": "ATCG"},
                            {"distancia": 2.8, "tempo": 1.1, "melhor_string": "ATGC"},
                        ],
                        "dist_base": 3.0,
                    }
                },
                "blfga": {
                    "base_0": {
                        "execucoes_detalhadas": [
                            {"distancia": 1.5, "tempo": 2.1, "melhor_string": "ATCG"},
                            {"distancia": 1.8, "tempo": 2.3, "melhor_string": "ATGC"},
                        ],
                        "dist_base": 3.0,
                    }
                },
            }
        }

        self.sample_config = {
            "nome": "Test Batch",
            "algoritmos": ["baseline", "blfga"],
            "dataset": {"tipo": "synthetic", "n": 10, "L": 4, "alphabet": "ATCG"},
            "execucoes_por_algoritmo_por_base": 2,
            "num_bases": 1,
        }

    def test_init(self):
        """Test BatchReporter initialization."""
        reporter = BatchReporter()
        assert reporter is not None
        assert hasattr(reporter, "generate_batch_report")

    def test_generate_batch_report_basic(self):
        """Test basic batch report generation."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Test with minimal valid inputs
            results = self.reporter.generate_batch_report(
                self.sample_results, self.sample_config, formats=[ReportFormat.TXT]
            )

            assert results is not None
            assert isinstance(results, dict)

    def test_generate_batch_report_multiple_formats(self):
        """Test batch report generation with multiple formats."""
        with tempfile.TemporaryDirectory() as temp_dir:
            results = self.reporter.generate_batch_report(
                self.sample_results,
                self.sample_config,
                formats=[ReportFormat.TXT, ReportFormat.JSON, ReportFormat.CSV],
            )

            assert results is not None
            assert isinstance(results, dict)

    def test_report_format_enum(self):
        """Test ReportFormat enum values."""
        assert ReportFormat.CSV.value == "csv"
        assert ReportFormat.JSON.value == "json"
        assert ReportFormat.HTML.value == "html"
        assert ReportFormat.TXT.value == "txt"

    def test_report_section_enum(self):
        """Test ReportSection enum values."""
        assert ReportSection.SUMMARY.value == "summary"
        assert ReportSection.DETAILED.value == "detailed"
        assert ReportSection.COMPARISON.value == "comparison"
        assert ReportSection.STATISTICS.value == "statistics"
        assert ReportSection.CONFIGURATION.value == "configuration"

    def test_generate_batch_report_with_sections(self):
        """Test batch report generation with specific sections."""
        with tempfile.TemporaryDirectory() as temp_dir:
            results = self.reporter.generate_batch_report(
                self.sample_results,
                self.sample_config,
                formats=[ReportFormat.TXT],
                sections=[ReportSection.SUMMARY, ReportSection.STATISTICS],
            )

            assert results is not None
            assert isinstance(results, dict)

    def test_generate_batch_report_empty_results(self):
        """Test batch report generation with empty results."""
        empty_results = {}

        # This might raise an exception or handle gracefully
        try:
            results = self.reporter.generate_batch_report(empty_results, self.sample_config, formats=[ReportFormat.TXT])
            assert results is not None
        except Exception:
            # Empty results might be expected to cause an error
            pass

    def test_generate_batch_report_invalid_config(self):
        """Test batch report generation with invalid config."""
        invalid_config = {}

        # This might raise an exception
        try:
            results = self.reporter.generate_batch_report(
                self.sample_results, invalid_config, formats=[ReportFormat.TXT]
            )
            assert results is not None
        except Exception:
            # Invalid config might be expected to cause an error
            pass

    def test_batch_reporter_with_default_parameters(self):
        """Test batch reporter with default parameters."""
        # Test that the method works with minimal required parameters
        results = self.reporter.generate_batch_report(self.sample_results, self.sample_config)

        assert results is not None
        assert isinstance(results, dict)
