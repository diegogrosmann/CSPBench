"""
Dataset Generation Wizard - Interactive Interface for Dataset Generation

Implements interactive wizards for parameter collection and generation of
synthetic and real datasets, following hexagonal architecture.
"""

import re
from pathlib import Path
from typing import Any, Dict, Optional


class DatasetWizard:
    """Interactive wizard for dataset generation."""

    def __init__(self):
        """Initialize the dataset wizard."""
        pass

    def show_main_menu(self) -> str:
        """
        Show main menu and return selected option.

        Returns:
            str: 'synthetic', 'real' or 'exit'
        """
        print("\n" + "=" * 60)
        print("🧬 Dataset Generation Wizard")
        print("=" * 60)
        print("\nChoose dataset type:")
        print("  1. Synthetic Dataset  - Generate random sequences")
        print("  2. Real Dataset       - Download from databases")
        print("  0. Exit")
        print()

        while True:
            choice = input("Select option (1/2/0): ").strip()
            if choice == "1":
                return "synthetic"
            elif choice == "2":
                return "real"
            elif choice == "0":
                return "exit"
            else:
                print("❌ Invalid option. Please choose 1, 2 or 0.")

    def collect_synthetic_params(self) -> Dict[str, Any]:
        """
        Collect parameters for synthetic dataset generation.

        Returns:
            Dict with parameters: n, length, alphabet, noise, seed
        """
        print("\n" + "-" * 40)
        print("🧪 Synthetic Dataset Configuration")
        print("-" * 40)

        params = {}

        # Number of sequences
        params["n"] = self._get_int_input(
            "Number of sequences", default=20, min_val=2, max_val=1000
        )

        # Sequence length
        params["length"] = self._get_int_input(
            "Sequence length", default=50, min_val=4, max_val=10000
        )

        # Alphabet
        params["alphabet"] = self._get_alphabet_input()

        # Noise level
        params["noise"] = self._get_float_input(
            "Noise level (0.0-1.0)", default=0.1, min_val=0.0, max_val=1.0
        )

        # Seed for reproducibility
        params["seed"] = self._get_optional_int_input(
            "Seed for reproducibility (Enter for random)", default=None
        )

        return params

    def collect_real_params(self) -> Dict[str, Any]:
        """
        Collect parameters for real dataset download.

        Returns:
            Dict with parameters: source, query, max_sequences, min_length, max_length
        """
        print("\n" + "-" * 40)
        print("🌐 Real Dataset Configuration")
        print("-" * 40)

        params = {}

        # Fonte dos dados
        params["source"] = self._get_source_input()

        if params["source"] == "ncbi":
            params.update(self._collect_ncbi_params())
        elif params["source"] == "file":
            params.update(self._collect_file_params())

        return params

    def _collect_ncbi_params(self) -> Dict[str, Any]:
        """Coleta parâmetros específicos do NCBI."""
        params = {}

        print("\n📚 Parâmetros do NCBI:")

        # Query de busca
        params["query"] = self._get_text_input(
            "Termo de busca", default="COI[Gene]", required=True
        )

        # Número máximo de sequências
        params["max_sequences"] = self._get_int_input(
            "Número máximo de sequências", default=50, min_val=1, max_val=10000
        )

        # Filtros de comprimento
        params["min_length"] = self._get_int_input(
            "Comprimento mínimo das sequências", default=100, min_val=1, max_val=50000
        )

        params["max_length"] = self._get_int_input(
            "Comprimento máximo das sequências",
            default=2000,
            min_val=params["min_length"],
            max_val=50000,
        )

        return params

    def _collect_file_params(self) -> Dict[str, Any]:
        """Coleta parâmetros para upload de arquivo."""
        params = {}

        print("\n📁 Parâmetros do Arquivo:")

        # Caminho do arquivo
        while True:
            file_path = self._get_text_input("Caminho do arquivo FASTA", required=True)

            path = Path(file_path)
            if path.exists() and path.is_file():
                params["file_path"] = str(path.absolute())
                break
            else:
                print("❌ Arquivo não encontrado! Tente novamente.")

        return params

    def generate_default_filename(
        self, dataset_type: str, params: Dict[str, Any]
    ) -> str:
        """
        Gera nome padrão para o arquivo baseado nos parâmetros.

        Args:
            dataset_type: 'synthetic' ou 'real'
            params: Parâmetros do dataset

        Returns:
            Nome do arquivo padrão
        """
        if dataset_type == "synthetic":
            n = params.get("n", "n")
            length = params.get("length", "L")
            alphabet = params.get("alphabet", "ACTG")
            noise = params.get("noise", 0.1)

            # Formato: synthetic_n20_L50_noise0.1_ACTG.fasta
            return f"synthetic_n{n}_L{length}_noise{noise}_{alphabet}.fasta"

        elif dataset_type == "real":
            source = params.get("source", "unknown")

            if source == "ncbi":
                query = params.get("query", "query")
                max_seq = params.get("max_sequences", "n")
                min_len = params.get("min_length", "min")
                max_len = params.get("max_length", "max")

                # Sanitizar query para nome de arquivo
                clean_query = re.sub(r"[^\w\s-]", "", query).strip()
                clean_query = re.sub(r"[-\s]+", "_", clean_query)

                # Formato: ncbi_COI_Gene_n50_100-2000bp.fasta
                return f"ncbi_{clean_query}_n{max_seq}_{min_len}-{max_len}bp.fasta"

            elif source == "file":
                original_path = Path(params.get("file_path", "unknown"))
                return f"imported_{original_path.stem}.fasta"

        return "dataset.fasta"

    def get_output_filename(self, default_name: str) -> str:
        """
        Solicita nome do arquivo de saída.

        Args:
            default_name: Nome padrão gerado

        Returns:
            Nome final do arquivo
        """
        print("\n📁 Nome do arquivo de saída:")
        print(f"   Padrão sugerido: {default_name}")

        filename = input("💾 Nome do arquivo (Enter para usar padrão): ").strip()

        if not filename:
            filename = default_name

        # Garantir extensão .fasta
        if not filename.lower().endswith(".fasta"):
            filename += ".fasta"

        return filename

    def _get_int_input(
        self, prompt: str, default: int, min_val: int, max_val: int
    ) -> int:
        """Obtém entrada inteira com validação."""
        while True:
            try:
                user_input = input(f"🔢 {prompt} (padrão: {default}): ").strip()

                if not user_input:
                    return default

                value = int(user_input)

                if min_val <= value <= max_val:
                    return value
                else:
                    print(f"❌ Valor deve estar entre {min_val} e {max_val}")

            except ValueError:
                print("❌ Por favor, insira um número válido")
            except (KeyboardInterrupt, EOFError):
                print("\n👋 Operação cancelada!")
                raise

    def _get_float_input(
        self, prompt: str, default: float, min_val: float, max_val: float
    ) -> float:
        """Obtém entrada float com validação."""
        while True:
            try:
                user_input = input(f"🔢 {prompt} (padrão: {default}): ").strip()

                if not user_input:
                    return default

                value = float(user_input)

                if min_val <= value <= max_val:
                    return value
                else:
                    print(f"❌ Valor deve estar entre {min_val} e {max_val}")

            except ValueError:
                print("❌ Por favor, insira um número válido")
            except (KeyboardInterrupt, EOFError):
                print("\n👋 Operação cancelada!")
                raise

    def _get_optional_int_input(
        self, prompt: str, default: Optional[int]
    ) -> Optional[int]:
        """Obtém entrada inteira opcional."""
        try:
            user_input = input(f"🔢 {prompt}: ").strip()

            if not user_input:
                return default

            return int(user_input)

        except ValueError:
            print("❌ Valor inválido, usando padrão")
            return default
        except (KeyboardInterrupt, EOFError):
            print("\n👋 Operação cancelada!")
            raise

    def _get_text_input(
        self, prompt: str, default: str = "", required: bool = False
    ) -> str:
        """Obtém entrada de texto."""
        while True:
            try:
                if default:
                    user_input = input(f"📝 {prompt} (padrão: {default}): ").strip()
                else:
                    user_input = input(f"📝 {prompt}: ").strip()

                if not user_input and default:
                    return default
                elif not user_input and required:
                    print("❌ Este campo é obrigatório")
                    continue
                else:
                    return user_input

            except (KeyboardInterrupt, EOFError):
                print("\n👋 Operação cancelada!")
                raise

    def _get_alphabet_input(self) -> str:
        """Obtém entrada para alfabeto."""
        print("\n🔤 Alfabetos disponíveis:")
        print("  1. DNA (ACTG)")
        print("  2. RNA (ACUG)")
        print("  3. Proteína (20 aminoácidos)")
        print("  4. Binário (01)")
        print("  5. Personalizado")

        alphabets = {"1": "ACTG", "2": "ACUG", "3": "ACDEFGHIKLMNPQRSTVWY", "4": "01"}

        while True:
            try:
                choice = input("🔤 Selecione o alfabeto (1-5): ").strip()

                if choice in alphabets:
                    return alphabets[choice]
                elif choice == "5":
                    custom = (
                        input("🔤 Digite o alfabeto personalizado: ").strip().upper()
                    )
                    if len(custom) >= 2:
                        return custom
                    else:
                        print("❌ Alfabeto deve ter pelo menos 2 caracteres")
                else:
                    print("❌ Opção inválida")

            except (KeyboardInterrupt, EOFError):
                print("\n👋 Operação cancelada!")
                raise

    def _get_source_input(self) -> str:
        """Obtém entrada para fonte dos dados."""
        print("\n🌐 Fontes disponíveis:")
        print("  1. NCBI - Download do banco de dados")
        print("  2. Arquivo - Importar arquivo local")

        sources = {"1": "ncbi", "2": "file"}

        while True:
            try:
                choice = input("🌐 Selecione a fonte (1-2): ").strip()

                if choice in sources:
                    return sources[choice]
                else:
                    print("❌ Opção inválida")

            except (KeyboardInterrupt, EOFError):
                print("\n👋 Operação cancelada!")
                raise
