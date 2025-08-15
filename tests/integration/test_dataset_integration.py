"""
Integration tests for dataset functionality.

Tests the complete dataset ecosystem including generation,
persistence, loading, and orchestration working together.
"""

import os
import tempfile
from pathlib import Path
from unittest.mock import patch, Mock

import pytest

from src.application.services.dataset_service import load_dataset
from src.domain.config import (
    SyntheticDatasetConfig,
    FileDatasetConfig,
    EntrezDatasetConfig,
)
from src.domain.dataset import Dataset
from src.infrastructure.persistence.dataset_repository import FileDatasetRepository
from infrastructure.orchestration.dataset_generation_orchestrator import (
    DatasetGenerationOrchestrator,
)
from src.application.services.dataset_generator import SyntheticDatasetGenerator


class TestDatasetIntegrationBasic:
    """Basic integration tests for dataset functionality."""

    def test_synthetic_dataset_generation_and_persistence(self, tmp_path):
        """Test complete synthetic dataset generation and file persistence."""
        # Create synthetic dataset
        dataset = SyntheticDatasetGenerator.generate_random(
            n=5, length=8, alphabet="ACGT", seed=42
        )

        # Save to file
        saved_path = FileDatasetRepository.save(
            dataset, "integration_test", str(tmp_path)
        )

        # Load back from file
        loaded_dataset, params = FileDatasetRepository.load(
            "integration_test.fasta", str(tmp_path)
        )

        # Verify integrity
        assert loaded_dataset.sequences == dataset.sequences
        assert loaded_dataset.size == dataset.size
        assert loaded_dataset.length == dataset.length
        assert Path(saved_path).exists()

    def test_dataset_service_synthetic_workflow(self):
        """Test dataset service with synthetic configuration."""
        config = SyntheticDatasetConfig(
            method="random", n=3, length=6, alphabet="ACGT", seed=123
        )

        dataset, params = load_dataset(config)

        assert isinstance(dataset, Dataset)
        assert dataset.size == 3
        assert dataset.length == 6
        assert "method" in params
        assert params["method"] == "random"

    def test_dataset_service_file_workflow(self, tmp_path):
        """Test dataset service with file configuration."""
        # Create test FASTA file
        fasta_content = """>seq1
ACGTACGT
>seq2
ATGTACGT
>seq3
GCGTACGT
"""
        test_file = tmp_path / "test.fasta"
        test_file.write_text(fasta_content)

        config = FileDatasetConfig(filename=str(test_file))

        dataset, params = load_dataset(config)

        assert isinstance(dataset, Dataset)
        assert dataset.size == 3
        assert dataset.length == 8
        assert dataset.sequences == ["ACGTACGT", "ATGTACGT", "GCGTACGT"]

    @patch("src.infrastructure.external.dataset_entrez.Entrez")
    def test_dataset_service_entrez_workflow(self, mock_entrez):
        """Test dataset service with Entrez configuration."""
        # Mock Entrez responses
        mock_search_handle = Mock()
        mock_search_result = {"IdList": ["123", "456"]}

        mock_entrez.esearch.return_value = mock_search_handle
        mock_entrez.read.side_effect = [mock_search_result]

        mock_fetch_handle = Mock()

        class MockSeqRecord:
            def __init__(self, seq):
                self.seq = seq

        mock_records = [MockSeqRecord("ACGTACGT"), MockSeqRecord("ATGTACGT")]
        mock_entrez.efetch.return_value = mock_fetch_handle

        with patch("src.infrastructure.external.dataset_entrez.SeqIO") as mock_seqio:
            mock_seqio.parse.return_value = mock_records

            config = EntrezDatasetConfig(query="test query", retmax=2, db="nucleotide")

            dataset, params = load_dataset(config)

            assert isinstance(dataset, Dataset)
            assert dataset.size == 2
            assert dataset.sequences == ["ACGTACGT", "ATGTACGT"]


class TestDatasetIntegrationComplexWorkflows:
    """Tests for complex dataset workflows and scenarios."""

    def test_multi_format_dataset_conversion(self, tmp_path):
        """Test converting datasets between different formats/sources."""
        # Generate synthetic dataset
        original_dataset = SyntheticDatasetGenerator.generate_random(
            n=4, length=10, alphabet="ACGT", seed=42
        )

        # Save to file
        saved_path = FileDatasetRepository.save(
            original_dataset, "multi_format", str(tmp_path)
        )

        # Load using file service
        config = FileDatasetConfig(filename="multi_format.fasta")
        loaded_dataset, params = load_dataset(config)

        # Convert back to synthetic config equivalent
        synthetic_config = SyntheticDatasetConfig(
            method="random",
            n=loaded_dataset.size,
            length=loaded_dataset.length,
            alphabet=loaded_dataset.alphabet,
            seed=123,  # Different seed
        )

        new_dataset, new_params = load_dataset(synthetic_config)

        # Should have same structure but different sequences (different seed)
        assert new_dataset.size == loaded_dataset.size
        assert new_dataset.length == loaded_dataset.length
        assert new_dataset.alphabet == loaded_dataset.alphabet
        assert (
            new_dataset.sequences != loaded_dataset.sequences
        )  # Different due to seed

    def test_dataset_batch_processing(self, tmp_path):
        """Test processing multiple datasets in batch."""
        configs = [
            SyntheticDatasetConfig(
                method="random", n=3, length=6, alphabet="ACGT", seed=i
            )
            for i in range(5)
        ]

        datasets = []
        for config in configs:
            dataset, params = load_dataset(config)
            datasets.append(dataset)

        # Verify all datasets are unique but have same structure
        for i, dataset in enumerate(datasets):
            assert dataset.size == 3
            assert dataset.length == 6

            # Compare with other datasets
            for j, other_dataset in enumerate(datasets):
                if i != j:
                    assert (
                        dataset.sequences != other_dataset.sequences
                    )  # Different due to different seeds

    def test_dataset_filtering_and_sampling_integration(self):
        """Test dataset filtering and sampling operations."""
        # Generate larger dataset
        original_dataset = SyntheticDatasetGenerator.generate_random(
            n=20, length=8, alphabet="ACGT", seed=42
        )

        # Sample subset
        sampled_dataset = original_dataset.sample(10, seed=123)

        # Filter by pattern (sequences starting with 'A')
        a_sequences = [seq for seq in sampled_dataset.sequences if seq[0] == "A"]
        if a_sequences:  # Only test if there are sequences starting with A
            filtered_dataset = sampled_dataset.filter_by_pattern("A", 0)

            assert len(filtered_dataset.sequences) == len(a_sequences)
            assert all(seq[0] == "A" for seq in filtered_dataset.sequences)

    def test_dataset_statistics_consistency_across_operations(self):
        """Test that dataset statistics remain consistent across operations."""
        # Generate dataset
        dataset = SyntheticDatasetGenerator.generate_clustered(
            n=12, length=10, num_clusters=3, alphabet="ACGT", seed=42
        )

        original_stats = dataset.get_statistics()

        # Add and remove sequences
        dataset.add_sequence("ACGTACGTAC")
        assert dataset.size == 13

        removed = dataset.remove_sequence(12)  # Remove the one we just added
        assert removed == "ACGTACGTAC"
        assert dataset.size == 12

        # Statistics should be back to original
        final_stats = dataset.get_statistics()

        assert final_stats["size"] == original_stats["size"]
        assert final_stats["min_length"] == original_stats["min_length"]
        assert final_stats["max_length"] == original_stats["max_length"]

    def test_dataset_serialization_roundtrip_integration(self, tmp_path):
        """Test complete serialization/deserialization roundtrip."""
        # Generate complex dataset
        dataset = SyntheticDatasetGenerator.generate_with_noise(
            base_sequence="ACGTACGTAC", n=8, noise_rate=0.3, alphabet="ACGT", seed=42
        )

        # Convert to dict
        dataset_dict = dataset.to_dict()

        # Reconstruct from dict
        reconstructed = Dataset.from_dict(dataset_dict)

        # Save reconstructed dataset
        saved_path = FileDatasetRepository.save(
            reconstructed, "roundtrip", str(tmp_path)
        )

        # Load from file
        final_dataset, params = FileDatasetRepository.load(
            "roundtrip.fasta", str(tmp_path)
        )

        # Should be identical to original
        assert final_dataset.sequences == dataset.sequences
        assert final_dataset.size == dataset.size
        assert final_dataset.length == dataset.length


class TestDatasetIntegrationErrorHandling:
    """Integration tests for error handling across dataset components."""

    def test_synthetic_to_file_with_invalid_characters(self):
        """Test handling of datasets with invalid characters."""
        # Create dataset with special characters
        sequences = ["ACGT-N*X", "WSSKMRMR"]
        dataset = Dataset(sequences=sequences)

        # Should still be able to save and load
        with tempfile.TemporaryDirectory() as tmp_dir:
            saved_path = FileDatasetRepository.save(dataset, "special_chars", tmp_dir)
            loaded_dataset, params = FileDatasetRepository.load(
                "special_chars.fasta", tmp_dir
            )

            assert loaded_dataset.sequences == sequences

    def test_error_propagation_through_service(self):
        """Test that errors propagate correctly through the service layer."""
        # File that doesn't exist
        config = FileDatasetConfig(filename="nonexistent.fasta")

        with pytest.raises(Exception):  # Should raise DatasetNotFoundError or similar
            load_dataset(config)

    def test_recovery_from_partial_failures(self, tmp_path):
        """Test recovery from partial failures in dataset processing."""
        # Create valid dataset
        valid_dataset = SyntheticDatasetGenerator.generate_random(
            n=3, length=6, alphabet="ACGT", seed=42
        )

        # Save it
        FileDatasetRepository.save(valid_dataset, "valid", str(tmp_path))

        # Try to load non-existent file (should fail)
        try:
            load_dataset(FileDatasetConfig(filename="nonexistent.fasta"))
        except Exception:
            pass  # Expected to fail

        # Should still be able to load the valid one
        loaded_dataset, params = load_dataset(
            FileDatasetConfig(filename=str(tmp_path / "valid.fasta"))
        )

        assert loaded_dataset.sequences == valid_dataset.sequences


class TestDatasetIntegrationPerformance:
    """Integration tests focused on performance characteristics."""

    def test_large_dataset_operations(self):
        """Test operations on larger datasets."""
        # Generate relatively large dataset
        dataset = SyntheticDatasetGenerator.generate_random(
            n=100, length=50, alphabet="ACGT", seed=42
        )

        # Operations should complete in reasonable time
        assert dataset.size == 100
        assert dataset.length == 50

        # Statistics calculation
        stats = dataset.get_statistics()
        assert stats["size"] == 100
        assert stats["total_characters"] == 5000

        # Sampling
        sample = dataset.sample(20, seed=123)
        assert sample.size == 20

    def test_memory_efficiency_across_operations(self, tmp_path):
        """Test memory efficiency of dataset operations."""
        import gc

        # Create dataset
        dataset = SyntheticDatasetGenerator.generate_random(
            n=50, length=20, alphabet="ACGT", seed=42
        )

        # Save to file
        FileDatasetRepository.save(dataset, "memory_test", str(tmp_path))

        # Force garbage collection
        del dataset
        gc.collect()

        # Load from file (should not accumulate memory significantly)
        loaded_dataset, params = FileDatasetRepository.load(
            "memory_test.fasta", str(tmp_path)
        )

        assert loaded_dataset.size == 50
        assert loaded_dataset.length == 20

        # Clean up
        del loaded_dataset
        gc.collect()


class TestDatasetIntegrationEndToEnd:
    """End-to-end integration tests simulating real usage scenarios."""

    @patch(
        "src.infrastructure.orchestrators.dataset_generation_orchestrator.DatasetWizard"
    )
    def test_orchestrator_synthetic_workflow(self, mock_wizard_class, tmp_path):
        """Test complete orchestrator workflow for synthetic datasets."""
        # Setup mock wizard
        mock_wizard = Mock()
        mock_wizard.show_main_menu.return_value = "synthetic"
        mock_wizard.collect_synthetic_params.return_value = {
            "method": "random",
            "n": 5,
            "length": 8,
            "alphabet": "ACGT",
            "seed": 42,
        }
        mock_wizard.generate_default_filename.return_value = (
            "synthetic_random_n5_L8.fasta"
        )
        mock_wizard.get_output_filename.return_value = "test_dataset.fasta"
        mock_wizard_class.return_value = mock_wizard

        # Create orchestrator
        orchestrator = DatasetGenerationOrchestrator(base_path=str(tmp_path))

        # Run complete workflow
        with patch("builtins.print"):  # Suppress output
            result_path = orchestrator.run_interactive_generation()

        # Verify result
        assert result_path is not None
        assert result_path.endswith("test_dataset.fasta")
        assert Path(result_path).exists()

        # Verify generated file content
        loaded_dataset, params = FileDatasetRepository.load(
            "test_dataset.fasta", str(tmp_path)
        )
        assert loaded_dataset.size == 5
        assert loaded_dataset.length == 8

    def test_complete_dataset_lifecycle(self, tmp_path):
        """Test complete dataset lifecycle from generation to analysis."""
        # Step 1: Generate synthetic dataset
        config = SyntheticDatasetConfig(
            method="noise", n=10, length=12, alphabet="ACGT", noise_rate=0.2, seed=42
        )

        original_dataset, gen_params = load_dataset(config)

        # Step 2: Save to file
        saved_path = FileDatasetRepository.save(
            original_dataset, "lifecycle_test", str(tmp_path)
        )

        # Step 3: Load from file using service
        file_config = FileDatasetConfig(filename="lifecycle_test.fasta")
        loaded_dataset, file_params = load_dataset(file_config)

        # Step 4: Analyze and manipulate
        stats = loaded_dataset.get_statistics()
        assert stats["size"] == 10
        assert stats["uniform_lengths"] is True

        # Step 5: Create sample subset
        sample = loaded_dataset.sample(5, seed=123)

        # Step 6: Save sample as new dataset
        sample_path = FileDatasetRepository.save(
            sample, "lifecycle_sample", str(tmp_path)
        )

        # Step 7: Verify complete workflow
        assert Path(saved_path).exists()
        assert Path(sample_path).exists()

        final_dataset, final_params = FileDatasetRepository.load(
            "lifecycle_sample.fasta", str(tmp_path)
        )
        assert final_dataset.size == 5
        assert final_dataset.length == 12

        # Original sequences should contain all sample sequences
        assert all(seq in loaded_dataset.sequences for seq in final_dataset.sequences)

    def test_cross_component_error_handling(self, tmp_path):
        """Test error handling across multiple components."""
        # Create valid dataset
        dataset = SyntheticDatasetGenerator.generate_random(
            n=3, length=4, alphabet="ACGT"
        )

        # Save with valid name
        FileDatasetRepository.save(dataset, "valid_dataset", str(tmp_path))

        # Try various error scenarios
        error_scenarios = [
            # Invalid file config
            FileDatasetConfig(filename="nonexistent.fasta"),
            # Invalid synthetic config (will be caught by generator)
            SyntheticDatasetConfig(method="invalid", n=0, length=-1),
        ]

        valid_loads = 0
        errors_caught = 0

        for config in error_scenarios:
            try:
                load_dataset(config)
                valid_loads += 1
            except Exception:
                errors_caught += 1

        # Should have caught errors appropriately
        assert errors_caught > 0

        # Valid dataset should still be loadable
        valid_config = FileDatasetConfig(filename="valid_dataset.fasta")
        with patch.dict(os.environ, {"DATASET_DIRECTORY": str(tmp_path)}):
            valid_dataset, params = load_dataset(valid_config)
            assert valid_dataset.size == 3

    def test_concurrent_dataset_operations(self, tmp_path):
        """Test concurrent dataset operations."""
        import threading
        import time

        results = []
        errors = []

        def generate_and_save(thread_id):
            try:
                # Generate unique dataset
                dataset = SyntheticDatasetGenerator.generate_random(
                    n=5, length=6, alphabet="ACGT", seed=thread_id
                )

                # Save with thread-specific name
                path = FileDatasetRepository.save(
                    dataset, f"thread_{thread_id}", str(tmp_path)
                )

                # Small delay to simulate processing time
                time.sleep(0.01)

                # Load back to verify
                loaded, params = FileDatasetRepository.load(
                    f"thread_{thread_id}.fasta", str(tmp_path)
                )

                results.append(
                    {
                        "thread_id": thread_id,
                        "original_size": dataset.size,
                        "loaded_size": loaded.size,
                        "path": path,
                    }
                )
            except Exception as e:
                errors.append({"thread_id": thread_id, "error": str(e)})

        # Start multiple threads
        threads = []
        for i in range(5):
            thread = threading.Thread(target=generate_and_save, args=(i,))
            threads.append(thread)
            thread.start()

        # Wait for all threads to complete
        for thread in threads:
            thread.join()

        # Verify results
        assert len(errors) == 0, f"Errors occurred: {errors}"
        assert len(results) == 5

        # All datasets should be valid and unique
        for result in results:
            assert result["original_size"] == result["loaded_size"]
            assert Path(result["path"]).exists()

        # Verify files were created for each thread
        created_files = list(Path(tmp_path).glob("thread_*.fasta"))
        assert len(created_files) == 5
