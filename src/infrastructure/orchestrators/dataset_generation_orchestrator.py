"""
Dataset Generation Orchestrator - Orquestração de Alto Nível

Coordena o wizard interativo com os serviços de geração,
mantendo a separação de responsabilidades da arquitetura hexagonal.
"""

from pathlib import Path
from typing import Any, Dict, Optional

from src.application.services.dataset_generation_service import DatasetGenerationService
from src.infrastructure.persistence.dataset_repository import FileDatasetRepository
from src.presentation.tui.dataset_wizard import DatasetWizard


class DatasetGenerationOrchestrator:
    """Orquestrador para geração interativa de datasets."""

    def __init__(self, base_path: str = "datasets"):
        """
        Inicializa o orquestrador.

        Args:
            base_path: Diretório base para salvar datasets
        """
        self.base_path = base_path
        self.wizard = DatasetWizard()

        # Criar repositório e serviço
        self.repository = FileDatasetRepository(base_path)
        self.service = DatasetGenerationService(self.repository)

    def run_interactive_generation(self) -> Optional[str]:
        """
        Executa o processo interativo completo de geração de dataset.

        Returns:
            Caminho do arquivo gerado ou None se cancelado
        """
        try:
            # Mostrar menu principal
            choice = self.wizard.show_main_menu()

            if choice == "exit":
                return None

            # Coletar parâmetros baseado na escolha
            if choice == "synthetic":
                return self._handle_synthetic_generation()
            elif choice == "real":
                return self._handle_real_generation()

        except (KeyboardInterrupt, EOFError):
            print("\n👋 Operação cancelada pelo usuário!")
            return None
        except (FileNotFoundError, PermissionError, OSError) as e:
            print(f"\n❌ Erro de arquivo: {e}")
            return None

    def _handle_synthetic_generation(self) -> Optional[str]:
        """Manipula geração de dataset sintético."""
        # Coletar parâmetros
        params = self.wizard.collect_synthetic_params()

        # Gerar nome padrão
        default_filename = self.wizard.generate_default_filename("synthetic", params)

        # Solicitar nome final
        filename = self.wizard.get_output_filename(default_filename)

        # Verificar se arquivo já existe
        output_path = Path(self.base_path) / filename
        if output_path.exists():
            overwrite = (
                input(f"\n⚠️  Arquivo '{filename}' já existe. Sobrescrever? (s/N): ")
                .strip()
                .lower()
            )
            if overwrite not in ["s", "sim", "y", "yes"]:
                print("📄 Operação cancelada.")
                return None

        # Gerar dataset
        print("\n🧪 Gerando dataset sintético...")
        print(f"   Parâmetros: {params}")

        dataset = self.service.generate_synthetic_dataset(params)

        # Salvar dataset
        saved_path = self.service.save_dataset(dataset, filename, self.base_path)

        # Mostrar resumo
        summary = self.service.get_generation_summary(dataset, params)
        self._show_generation_summary("Sintético", saved_path, summary)

        return saved_path

    def _handle_real_generation(self) -> Optional[str]:
        """Manipula geração de dataset real."""
        # Coletar parâmetros
        params = self.wizard.collect_real_params()

        # Gerar nome padrão
        default_filename = self.wizard.generate_default_filename("real", params)

        # Solicitar nome final
        filename = self.wizard.get_output_filename(default_filename)

        # Verificar se arquivo já existe
        output_path = Path(self.base_path) / filename
        if output_path.exists():
            overwrite = (
                input(f"\n⚠️  Arquivo '{filename}' já existe. Sobrescrever? (s/N): ")
                .strip()
                .lower()
            )
            if overwrite not in ["s", "sim", "y", "yes"]:
                print("📄 Operação cancelada.")
                return None

        # Baixar/importar dataset
        print("\n🌐 Processando dataset real...")
        print(f"   Parâmetros: {params}")

        dataset = self.service.download_real_dataset(params)

        # Salvar dataset
        saved_path = self.service.save_dataset(dataset, filename, self.base_path)

        # Mostrar resumo
        summary = self.service.get_generation_summary(dataset, params)
        self._show_generation_summary("Real", saved_path, summary)

        return saved_path

    def _show_generation_summary(
        self, dataset_type: str, saved_path: str, summary: Dict[str, Any]
    ) -> None:
        """Mostra resumo da geração."""
        print(f"\n✅ Dataset {dataset_type} gerado com sucesso!")
        print("=" * 50)
        print(f"📁 Arquivo salvo: {saved_path}")
        print(f"📊 Sequências: {summary['total_sequences']}")
        print(f"📏 Comprimento: {summary['sequence_length']}")
        print(f"🔤 Alfabeto: {summary['alphabet']}")
        print(f"💾 Tamanho estimado: {summary['estimated_size_kb']:.1f} KB")
        print("=" * 50)
