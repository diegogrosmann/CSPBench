from unittest.mock import patch
from pathlib import Path

from src.infrastructure.orchestration.dataset_generation_orchestrator import (
    DatasetGenerationOrchestrator,
)
from src.domain.dataset import Dataset


def test_orchestrator_synthetic_flow(tmp_path):
    orch = DatasetGenerationOrchestrator(base_path=str(tmp_path))
    with patch(
        "src.presentation.cli.dataset_wizard.DatasetWizard.show_main_menu",
        return_value="synthetic",
    ):
        with patch(
            "src.presentation.cli.dataset_wizard.DatasetWizard.collect_synthetic_params",
            return_value={"method": "random", "n": 3, "length": 5, "alphabet": "ACGT"},
        ):
            with patch(
                "src.presentation.cli.dataset_wizard.DatasetWizard.get_output_filename",
                return_value="out.fasta",
            ):
                with patch(
                    "src.application.services.dataset_generator.SyntheticDatasetGenerator.generate_random",
                    return_value=(
                        Dataset(id="test_dataset", name="test", sequences=["ACGTA", "ACGTT", "ACGTC"], alphabet="ACGT"),
                        {},
                    ),
                ):
                    path = orch.run_interactive_generation()
                    assert path is not None
                    assert Path(path).exists()
