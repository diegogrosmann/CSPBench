"""
Dataset Generation Orchestrator - High-Level Orchestration

Refactored to use new structures:
 - Domain: Dataset
 - Services: SyntheticDatasetGenerator (synthetic generation)
 - External: EntrezDatasetDownloader (NCBI download)
 - Persistence: FileDatasetRepository (save/load .fasta)

Orchestrates interactive flow via DatasetWizard.
"""

from pathlib import Path
from typing import Any, Dict, Optional, Tuple

from src.application.services.dataset_generator import SyntheticDatasetGenerator
from src.domain.dataset import Dataset
from src.infrastructure.external.dataset_entrez import EntrezDatasetDownloader
from src.infrastructure.persistence.dataset_repository import FileDatasetRepository

try:
    # Nova implementação CLI substitui antigo módulo TUI removido
    from src.presentation.cli.dataset_wizard import DatasetWizard  # type: ignore
except ModuleNotFoundError as _e:  # pragma: no cover
    raise RuntimeError(
        "DatasetWizard implementation not found. Ensure src/presentation/cli/dataset_wizard.py exists."
    ) from _e


class DatasetGenerationOrchestrator:
    """Orchestrator for interactive dataset generation."""

    def __init__(self, base_path: str = "datasets"):
        """Initialize orchestrator.

        Args:
            base_path: Base directory to save datasets
        """
        self.base_path = base_path
        self.wizard = DatasetWizard()

        # Repository (file-based)
        self.repository = FileDatasetRepository

    def run_interactive_generation(self) -> Optional[str]:
        """Execute complete interactive dataset generation process.

        Returns:
            Generated file path or None if cancelled
        """
        try:
            # Show main menu
            choice = self.wizard.show_main_menu()

            if choice == "exit":
                return None

            # Collect parameters based on choice
            if choice == "synthetic":
                return self._handle_synthetic_generation()
            elif choice == "real":
                return self._handle_real_generation()

        except (KeyboardInterrupt, EOFError):
            print("\n👋 Operation cancelled by user!")
            return None
        except (FileNotFoundError, PermissionError, OSError) as e:
            print(f"\n❌ File error: {e}")
            return None

    def _handle_synthetic_generation(self) -> Optional[str]:
        """Handle synthetic dataset generation.
        
        Returns:
            Generated file path or None if cancelled
        """
        # Collect parameters
        params = self.wizard.collect_synthetic_params()

        # Generate default name
        default_filename = self.wizard.generate_default_filename("synthetic", params)

        # Request final name
        filename = self.wizard.get_output_filename(default_filename)

        # Check if file already exists
        output_path = Path(self.base_path) / filename
        if output_path.exists():
            overwrite = (
                input(f"\n⚠️  File '{filename}' already exists. Overwrite? (y/N): ")
                .strip()
                .lower()
            )
            if overwrite not in ["s", "sim", "y", "yes"]:
                print("📄 Operation cancelled.")
                return None

        # Generate dataset (dispatch by method)
        print("\n🧪 Generating synthetic dataset...")
        method = params.get("method", "random")
        ds, used = self._generate_synthetic_dispatch(params)

        # Save dataset
        saved_path = self._save_dataset(ds, filename)

        # Show summary
        summary = self._build_summary(ds, used)
        self._show_generation_summary("Synthetic", saved_path, summary)

        return saved_path

    def _handle_real_generation(self) -> Optional[str]:
        """Handle real dataset generation.
        
        Returns:
            Generated file path or None if cancelled
        """
        # Collect parameters
        params = self.wizard.collect_real_params()

        # Generate default name
        default_filename = self.wizard.generate_default_filename("real", params)

        # Request final name
        filename = self.wizard.get_output_filename(default_filename)

        # Check if file already exists
        output_path = Path(self.base_path) / filename
        if output_path.exists():
            overwrite = (
                input(f"\n⚠️  File '{filename}' already exists. Overwrite? (y/N): ")
                .strip()
                .lower()
            )
            if overwrite not in ["s", "sim", "y", "yes"]:
                print("📄 Operation cancelled.")
                return None

        # Download/import dataset
        print("\n🌐 Processing real dataset...")
        # Only NCBI supported
        entrez_params: Dict[str, Any] = {
            "query": params.get("query", "*"),
            "db": params.get("db"),
            "max_sequences": params.get("max_sequences"),
            "min_length": params.get("min_length"),
            "max_length": params.get("max_length"),
        }
        ds, used = EntrezDatasetDownloader.download(entrez_params)

        # Save dataset
        saved_path = self._save_dataset(ds, filename)

        # Show summary
        summary = self._build_summary(ds, used)
        self._show_generation_summary("Real", saved_path, summary)

        return saved_path

    def _show_generation_summary(
        self, dataset_type: str, saved_path: str, summary: Dict[str, Any]
    ) -> None:
        """Show generation summary.
        
        Args:
            dataset_type: Type of dataset generated
            saved_path: Path where dataset was saved
            summary: Summary statistics dictionary
        """
        print(f"\n✅ {dataset_type} dataset generated successfully!")
        print("=" * 50)
        print(f"📁 File saved: {saved_path}")
        print(f"📊 Sequences: {summary['total_sequences']}")
        print(f"📏 Length: {summary['sequence_length']}")
        print(f"🔤 Alphabet: {summary['alphabet']}")
        print(f"💾 Estimated size: {summary['estimated_size_kb']:.1f} KB")
        print("=" * 50)

    # -------- Helpers --------
    def _save_dataset(self, dataset: Dataset, filename: str) -> str:
        """Persist dataset using FileDatasetRepository and return saved path.
        
        Args:
            dataset: Dataset object to save
            filename: Desired filename
            
        Returns:
            Actual saved file path
        """
        name_no_ext = Path(filename).stem
        return self.repository.save(dataset, name_no_ext, base_path=self.base_path)

    def _build_summary(
        self, ds: Dataset, used_params: Dict[str, Any] | None = None
    ) -> Dict[str, Any]:
        """Build summary statistics for dataset.
        
        Args:
            ds: Dataset object
            used_params: Parameters used for generation
            
        Returns:
            Summary statistics dictionary
        """
        stats = ds.get_statistics()
        est_kb = stats["total_characters"] / 1024.0
        return {
            "total_sequences": stats["size"],
            "sequence_length": stats["L"],
            "alphabet": stats["alphabet"],
            "estimated_size_kb": est_kb,
            "params": used_params or {},
        }

    def _generate_synthetic_dispatch(
        self, params: Dict[str, Any]
    ) -> Tuple[Dataset, Dict[str, Any]]:
        """Dispatch synthetic generation to appropriate method using new defaults.
        
        Args:
            params: Generation parameters dictionary
            
        Returns:
            Tuple of (generated_dataset, used_parameters)
        """
        method = params.get("method", "random")

        # Extract common parameters - let generator methods handle defaults
        common_params = {
            "n": params.get("n"),
            "alphabet": params.get("alphabet"),
            "seed": params.get("seed"),
        }

        if method == "random":
            return SyntheticDatasetGenerator.generate_random(
                length=params.get("length"), **common_params
            )

        elif method == "noise":
            return SyntheticDatasetGenerator.generate_with_noise(
                base_sequence=params.get("center_sequence")
                or params.get("base_sequence"),
                noise_rate=params.get("noise_rate"),
                **common_params,
            )

        elif method == "clustered":
            return SyntheticDatasetGenerator.generate_clustered(
                length=params.get("length"),
                num_clusters=params.get("num_clusters"),
                cluster_distance=params.get("cluster_noise")
                or params.get("cluster_distance"),
                **common_params,
            )

        elif method == "mutations":
            return SyntheticDatasetGenerator.generate_with_mutations(
                base_sequence=params.get("base_sequence"),
                mutation_rate=params.get("mutation_rate"),
                **common_params,
            )

        # fallback
        return SyntheticDatasetGenerator.generate_random(**common_params)
