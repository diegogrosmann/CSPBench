"""
Dataset Generation Service - Camada de Aplicação

Orquestra a geração de datasets sintéticos e reais, coordenando
entre diferentes geradores e repositórios.
"""

from pathlib import Path
from typing import Any, Dict, Optional

from src.domain import Dataset, SyntheticDatasetGenerator
from src.infrastructure.persistence.dataset_repository import FileDatasetRepository


class DatasetGenerationService:
    """Serviço de geração de datasets."""

    def __init__(self, dataset_repository: FileDatasetRepository):
        """
        Inicializa o serviço.

        Args:
            dataset_repository: Repositório para persistência de datasets
        """
        self.dataset_repository = dataset_repository
        self.synthetic_generator = SyntheticDatasetGenerator()

    def generate_synthetic_dataset(self, params: Dict[str, Any]) -> Dataset:
        """
        Gera dataset sintético baseado nos parâmetros.

        Args:
            params: Parâmetros de geração (n, length, alphabet, noise, seed)

        Returns:
            Dataset gerado
        """
        return self.synthetic_generator.generate_random(
            n=params["n"],
            length=params["length"],
            alphabet=params["alphabet"],
            seed=params.get("seed"),
        )

    def download_real_dataset(self, params: Dict[str, Any]) -> Dataset:
        """
        Baixa dataset real baseado nos parâmetros.

        Args:
            params: Parâmetros de download (source, query, etc.)

        Returns:
            Dataset baixado

        Raises:
            NotImplementedError: Para fontes não implementadas
        """
        source = params.get("source")

        if source == "ncbi":
            return self._download_from_ncbi(params)
        elif source == "file":
            return self._import_from_file(params)
        else:
            raise NotImplementedError(f"Fonte '{source}' não implementada")

    def save_dataset(
        self, dataset: Dataset, filename: str, base_path: str = "datasets"
    ) -> str:
        """
        Salva dataset em arquivo.

        Args:
            dataset: Dataset para salvar
            filename: Nome do arquivo
            base_path: Diretório base

        Returns:
            Caminho completo do arquivo salvo
        """
        output_path = Path(base_path) / filename

        # Criar diretório se não existir
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Salvar em formato FASTA
        with open(output_path, "w", encoding="utf-8") as f:
            for i, seq in enumerate(dataset.sequences):
                f.write(f">seq_{i+1}\n{seq}\n")

        return str(output_path.absolute())

    def _download_from_ncbi(self, params: Dict[str, Any]) -> Dataset:
        """
        Baixa dataset do NCBI.

        Args:
            params: Parâmetros NCBI (query, max_sequences, min_length, max_length)

        Returns:
            Dataset baixado
        """
        # Por enquanto, simulamos o download com dados sintéticos
        # Na implementação real, usaria Bio.Entrez ou similar

        print(f"🌐 Simulando download do NCBI...")
        print(f"   Query: {params['query']}")
        print(f"   Max sequences: {params['max_sequences']}")
        print(f"   Length range: {params['min_length']}-{params['max_length']}")

        # Simular com dataset sintético baseado nos parâmetros
        synthetic_params = {
            "n": params["max_sequences"],
            "length": (params["min_length"] + params["max_length"]) // 2,
            "alphabet": "ACTG",  # DNA padrão
            "noise": 0.05,  # Baixo ruído para dados "reais"
            "seed": 42,  # Determinístico para simulação
        }

        return self.generate_synthetic_dataset(synthetic_params)

    def _import_from_file(self, params: Dict[str, Any]) -> Dataset:
        """
        Importa dataset de arquivo.

        Args:
            params: Parâmetros do arquivo (file_path)

        Returns:
            Dataset importado
        """
        file_path = params["file_path"]

        print(f"📁 Importando arquivo: {file_path}")

        # Usar o repositório para carregar o arquivo
        dataset_name = Path(file_path).stem
        return self.dataset_repository.load(dataset_name)

    def get_generation_summary(
        self, dataset: Dataset, params: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Gera resumo da geração do dataset.

        Args:
            dataset: Dataset gerado
            params: Parâmetros usados

        Returns:
            Resumo da geração
        """
        return {
            "total_sequences": len(dataset.sequences),
            "sequence_length": len(dataset.sequences[0]) if dataset.sequences else 0,
            "alphabet": dataset.alphabet,
            "parameters_used": params,
            "estimated_size_kb": sum(len(seq) for seq in dataset.sequences) / 1024,
        }
