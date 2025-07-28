"""
Dataset Generation Orchestrator - High-Level Orchestration

Coordinates the interactive wizard with generation services,
maintaining the responsibility separation of hexagonal architecture.
"""

from pathlib import Path
from typing import Any, Dict, Optional

from src.application.services.dataset_generation_service import DatasetGenerationService
from src.infrastructure.persistence.dataset_repository import FileDatasetRepository
from src.presentation.tui.dataset_wizard import DatasetWizard


class DatasetGenerationOrchestrator:
    """Orchestrator for interactive dataset generation."""

    def __init__(self, base_path: str = "datasets"):
        """
        Initialize orchestrator.

        Args:
            base_path: Base directory to save datasets
        """
        self.base_path = base_path
        self.wizard = DatasetWizard()

        # Create repository and service
        self.repository = FileDatasetRepository(base_path)
        self.service = DatasetGenerationService(self.repository)

    def run_interactive_generation(self) -> Optional[str]:
        """
        Execute complete interactive dataset generation process.

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
        """Handle synthetic dataset generation."""
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

        # Generate dataset
        print("\n🧪 Generating synthetic dataset...")
        print(f"   Parameters: {params}")

        dataset = self.service.generate_synthetic_dataset(params)

        # Save dataset
        saved_path = self.service.save_dataset(dataset, filename, self.base_path)

        # Show summary
        summary = self.service.get_generation_summary(dataset, params)
        self._show_generation_summary("Synthetic", saved_path, summary)

        return saved_path

    def _handle_real_generation(self) -> Optional[str]:
        """Handle real dataset generation."""
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
        print(f"   Parameters: {params}")

        dataset = self.service.download_real_dataset(params)

        # Save dataset
        saved_path = self.service.save_dataset(dataset, filename, self.base_path)

        # Show summary
        summary = self.service.get_generation_summary(dataset, params)
        self._show_generation_summary("Real", saved_path, summary)

        return saved_path

    def _show_generation_summary(
        self, dataset_type: str, saved_path: str, summary: Dict[str, Any]
    ) -> None:
        """Show generation summary."""
        print(f"\n✅ {dataset_type} dataset generated successfully!")
        print("=" * 50)
        print(f"📁 File saved: {saved_path}")
        print(f"📊 Sequences: {summary['total_sequences']}")
        print(f"📏 Length: {summary['sequence_length']}")
        print(f"🔤 Alphabet: {summary['alphabet']}")
        print(f"💾 Estimated size: {summary['estimated_size_kb']:.1f} KB")
        print("=" * 50)
