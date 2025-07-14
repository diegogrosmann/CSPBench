"""
Domínio: Dataset

Este módulo contém a entidade principal Dataset e suas operações,
representando sequências de strings e metadados associados.
Implementação pura sem dependências externas.
"""

import random
from typing import Any, Dict, List, Optional


class Dataset:
    """
    Entidade Dataset representando um conjunto de strings para CSP.

    Attributes:
        sequences: Lista de strings do dataset
        metadata: Metadados do dataset (tamanho, origem, etc.)
    """

    def __init__(self, sequences: List[str], metadata: Optional[Dict[str, Any]] = None):
        """
        Inicializa um dataset com sequências e metadados.

        Args:
            sequences: Lista de strings
            metadata: Dicionário com metadados opcionais
        """
        if not sequences:
            raise ValueError("Dataset não pode estar vazio")

        # Validar que todas as strings têm mesmo comprimento
        length = len(sequences[0])
        if not all(len(seq) == length for seq in sequences):
            raise ValueError("Todas as strings devem ter mesmo comprimento")

        self.sequences = sequences
        self.metadata = metadata or {}

        # Inferir metadados básicos
        self.metadata.update(
            {
                "n": len(sequences),
                "L": length,
                "alphabet": self._infer_alphabet(),
                "diversity": self._calculate_diversity(),
            }
        )

    def _infer_alphabet(self) -> str:
        """Infere alfabeto a partir das sequências."""
        alphabet_set = set()
        for seq in self.sequences:
            alphabet_set.update(seq)
        return "".join(sorted(alphabet_set))

    def _calculate_diversity(self) -> float:
        """Calcula diversidade média do dataset."""
        if len(self.sequences) < 2:
            return 0.0

        total_distance = 0
        total_pairs = 0

        for i in range(len(self.sequences)):
            for j in range(i + 1, len(self.sequences)):
                distance = sum(
                    c1 != c2 for c1, c2 in zip(self.sequences[i], self.sequences[j])
                )
                total_distance += distance
                total_pairs += 1

        avg_distance = total_distance / total_pairs if total_pairs > 0 else 0
        max_possible = len(self.sequences[0])
        return avg_distance / max_possible if max_possible > 0 else 0

    @property
    def size(self) -> int:
        """Retorna número de sequências."""
        return len(self.sequences)

    @property
    def length(self) -> int:
        """Retorna comprimento das sequências."""
        return len(self.sequences[0]) if self.sequences else 0

    @property
    def alphabet(self) -> str:
        """Retorna alfabeto do dataset."""
        return self.metadata.get("alphabet", "")

    def validate(self) -> bool:
        """
        Valida consistência do dataset.

        Returns:
            bool: True se dataset está válido
        """
        if not self.sequences:
            return False

        # Verificar comprimentos uniformes
        length = len(self.sequences[0])
        if not all(len(seq) == length for seq in self.sequences):
            return False

        # Verificar caracteres válidos
        alphabet_set = set(self.alphabet)
        for seq in self.sequences:
            if not all(c in alphabet_set for c in seq):
                return False

        return True

    def get_statistics(self) -> Dict[str, Any]:
        """
        Retorna estatísticas detalhadas do dataset.

        Returns:
            dict: Estatísticas do dataset
        """
        return {
            "size": self.size,
            "length": self.length,
            "alphabet": self.alphabet,
            "alphabet_size": len(self.alphabet),
            "diversity": self.metadata.get("diversity", 0),
            "total_characters": self.size * self.length,
            "metadata": self.metadata.copy(),
        }

    def to_dict(self) -> Dict[str, Any]:
        """
        Converte dataset para dicionário.

        Returns:
            dict: Representação em dicionário
        """
        return {"sequences": self.sequences.copy(), "metadata": self.metadata.copy()}

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Dataset":
        """
        Cria dataset a partir de dicionário.

        Args:
            data: Dicionário com sequences e metadata

        Returns:
            Dataset: Instância criada
        """
        return cls(data["sequences"], data.get("metadata"))

    def add_sequence(self, sequence: str) -> None:
        """
        Adiciona nova sequência ao dataset.

        Args:
            sequence: String a ser adicionada

        Raises:
            ValueError: Se sequência tem comprimento diferente
        """
        if len(sequence) != self.length:
            raise ValueError(f"Sequência deve ter comprimento {self.length}")

        self.sequences.append(sequence)

        # Atualizar metadados
        self.metadata["n"] = len(self.sequences)
        self.metadata["alphabet"] = self._infer_alphabet()
        self.metadata["diversity"] = self._calculate_diversity()

    def remove_sequence(self, index: int) -> str:
        """
        Remove sequência por índice.

        Args:
            index: Índice da sequência a ser removida

        Returns:
            str: Sequência removida

        Raises:
            IndexError: Se índice inválido
        """
        if not 0 <= index < len(self.sequences):
            raise IndexError("Índice fora do intervalo")

        removed = self.sequences.pop(index)

        # Atualizar metadados
        self.metadata["n"] = len(self.sequences)
        if self.sequences:
            self.metadata["diversity"] = self._calculate_diversity()

        return removed

    def filter_by_pattern(self, pattern: str, position: int) -> "Dataset":
        """
        Filtra sequências que têm padrão específico em posição.

        Args:
            pattern: Caractere ou padrão a buscar
            position: Posição para verificar

        Returns:
            Dataset: Novo dataset com sequências filtradas
        """
        filtered_sequences = [seq for seq in self.sequences if seq[position] == pattern]

        new_metadata = self.metadata.copy()
        new_metadata["filter_applied"] = f"position_{position}={pattern}"

        return Dataset(filtered_sequences, new_metadata)

    def sample(self, n: int, seed: Optional[int] = None) -> "Dataset":
        """
        Retorna amostra aleatória do dataset.

        Args:
            n: Número de sequências na amostra
            seed: Semente para reprodutibilidade

        Returns:
            Dataset: Novo dataset com amostra
        """
        if n > len(self.sequences):
            raise ValueError("Tamanho da amostra maior que dataset")

        rng = random.Random(seed)
        sampled_sequences = rng.sample(self.sequences, n)

        new_metadata = self.metadata.copy()
        new_metadata["sampled_from"] = len(self.sequences)
        new_metadata["sample_seed"] = seed

        return Dataset(sampled_sequences, new_metadata)


class SyntheticDatasetGenerator:
    """Gerador de datasets sintéticos para teste de algoritmos CSP."""

    @staticmethod
    def generate_from_center(
        center: str,
        n: int,
        noise_rate: float,
        alphabet: str,
        seed: Optional[int] = None,
    ) -> Dataset:
        """
        Gera dataset baseado em string central com ruído.

        Args:
            center: String central
            n: Número de sequências a gerar
            noise_rate: Taxa de ruído (0-1)
            alphabet: Alfabeto válido
            seed: Semente para reprodutibilidade

        Returns:
            Dataset: Dataset sintético gerado
        """
        rng = random.Random(seed)
        sequences = []

        for _ in range(n):
            sequence = list(center)

            # Aplicar ruído
            for i in range(len(sequence)):
                if rng.random() < noise_rate:
                    # Trocar por caractere diferente do alfabeto
                    available = [c for c in alphabet if c != sequence[i]]
                    if available:
                        sequence[i] = rng.choice(available)

            sequences.append("".join(sequence))

        metadata = {
            "type": "synthetic",
            "center_string": center,
            "noise_rate": noise_rate,
            "generation_seed": seed,
            "alphabet_used": alphabet,
        }

        return Dataset(sequences, metadata)

    @staticmethod
    def generate_random(
        n: int, length: int, alphabet: str, seed: Optional[int] = None
    ) -> Dataset:
        """
        Gera dataset completamente aleatório.

        Args:
            n: Número de sequências
            length: Comprimento das sequências
            alphabet: Alfabeto válido
            seed: Semente para reprodutibilidade

        Returns:
            Dataset: Dataset aleatório gerado
        """
        rng = random.Random(seed)
        sequences = []

        for _ in range(n):
            sequence = "".join(rng.choice(alphabet) for _ in range(length))
            sequences.append(sequence)

        metadata = {
            "type": "random",
            "generation_seed": seed,
            "alphabet_used": alphabet,
        }

        return Dataset(sequences, metadata)

    @staticmethod
    def generate_clustered(
        n_clusters: int,
        sequences_per_cluster: int,
        length: int,
        alphabet: str,
        noise_rate: float = 0.1,
        seed: Optional[int] = None,
    ) -> Dataset:
        """
        Gera dataset com clusters de sequências similares.

        Args:
            n_clusters: Número de clusters
            sequences_per_cluster: Sequências por cluster
            length: Comprimento das sequências
            alphabet: Alfabeto válido
            noise_rate: Taxa de ruído dentro dos clusters
            seed: Semente para reprodutibilidade

        Returns:
            Dataset: Dataset com clusters
        """
        rng = random.Random(seed)
        all_sequences = []
        cluster_centers = []

        # Gerar centros dos clusters
        for _ in range(n_clusters):
            center = "".join(rng.choice(alphabet) for _ in range(length))
            cluster_centers.append(center)

            # Gerar sequências para este cluster
            cluster_dataset = SyntheticDatasetGenerator.generate_from_center(
                center,
                sequences_per_cluster,
                noise_rate,
                alphabet,
                rng.randint(0, 1000000),
            )
            all_sequences.extend(cluster_dataset.sequences)

        metadata = {
            "type": "clustered",
            "n_clusters": n_clusters,
            "sequences_per_cluster": sequences_per_cluster,
            "noise_rate": noise_rate,
            "cluster_centers": cluster_centers,
            "generation_seed": seed,
            "alphabet_used": alphabet,
        }

        return Dataset(all_sequences, metadata)
