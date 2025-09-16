"""
Synthetic Dataset Generator.

Provides synthetic dataset generation capabilities for testing and
benchmarking purposes. Supports multiple generation modes including
random sequences, noise-based variations, mutation-based sequences,
and clustered data.

Features:
    - Random sequence generation
    - Noise-based sequence variation
    - Mutation-based sequence evolution (substitution, insertion, deletion)
    - Clustered sequence generation
    - Configurable parameters for all generation modes
    - Deterministic generation with seed support
    - Compatibility with existing test infrastructure

Default Values:
    - n: 20 sequences
    - length: 50 characters
    - alphabet: "ACGT" (DNA)
    - noise_rate: 0.1 (10%)
    - mutation_rate: 0.1 (10%)
    - num_clusters: 2
    - cluster_distance: 0.2 (20%)
    - pad_char: "N"
"""

import random
from typing import List, Optional

from src.domain.config import SyntheticDatasetConfig
from src.domain.dataset import Dataset

# Default values for synthetic generation
SYNTHETIC_DEFAULTS = {
    "n": 20,
    "length": 50,
    "alphabet": "ACGT",
    "seed": None,
    "noise_rate": 0.1,
    "mutation_rate": 0.1,
    "mutation_types": None,  # will be ["substitution"] if None
    "num_clusters": 2,
    "cluster_distance": 0.2,
    "base_sequence": None,  # will be generated automatically if None
    "pad_char": "N",
}


class SyntheticDatasetGenerator:
    """
    Synthetic dataset generator with compatibility for existing tests.

    Provides static methods for generating various types of synthetic
    datasets with configurable parameters. All methods return a tuple
    of (Dataset, parameters) for consistency.

    Methods:
        generate_random: Random sequence generation
        generate_with_noise: Noise-based variation
        generate_with_mutations: Mutation-based evolution
        generate_clustered: Clustered sequence generation
        generate_from_config: Configuration-based generation
    """

    @staticmethod
    def _ensure_rng(
        seed: Optional[int] = None, rng: Optional[random.Random] = None
    ) -> random.Random:
        """
        Return an RNG from an existing one or create a new one based on seed.

        Prioritizes the injected rng. If none provided, creates a new one using the seed.

        Args:
            seed (Optional[int]): Random seed for reproducibility.
            rng (Optional[random.Random]): Existing random number generator.

        Returns:
            random.Random: Random number generator instance.
        """
        return rng if rng is not None else random.Random(seed)

    @staticmethod
    def generate_random(
        n: Optional[int] = None,
        length: Optional[int] = None,
        alphabet: Optional[str] = None,
        seed: Optional[int] = None,
        rng: Optional[random.Random] = None,
    ) -> tuple[Dataset, dict]:
        """
        Generate completely random sequences.

        Creates n sequences of specified length using characters from
        the given alphabet. Each position is independently random.

        Args:
            n (Optional[int]): Number of sequences to generate. Defaults to 20.
            length (Optional[int]): Length of each sequence. Defaults to 50.
            alphabet (Optional[str]): Character alphabet to use. Defaults to "ACGT".
            seed (Optional[int]): Random seed for reproducibility.
            rng (Optional[random.Random]): Existing random number generator.

        Returns:
            Tuple[Dataset, dict]: A tuple containing:
                - Dataset: Dataset with random sequences
                - dict: Generation parameters

        Raises:
            ValueError: If n or length <= 0, or alphabet is empty.
        """
        # Apply default values
        n = n if n is not None else SYNTHETIC_DEFAULTS["n"]
        length = length if length is not None else SYNTHETIC_DEFAULTS["length"]
        alphabet = alphabet if alphabet is not None else SYNTHETIC_DEFAULTS["alphabet"]

        if n <= 0 or length <= 0:
            raise ValueError("n and length must be > 0")
        if not alphabet:
            raise ValueError("alphabet cannot be empty")
        _rng = SyntheticDatasetGenerator._ensure_rng(seed, rng)
        seqs = ["".join(_rng.choice(alphabet) for _ in range(length)) for _ in range(n)]
        ds = Dataset(
            id="synthetic_random",
            name="synthetic_random",
            sequences=seqs,
            alphabet=alphabet,
        )
        params = {
            "generator": "generate_random",
            "n": n,
            "length": length,
            "alphabet": alphabet,
            "seed_arg": seed,
            "rng_provided": rng is not None,
        }
        return ds, params

    @staticmethod
    def generate_with_noise(
        base_sequence: Optional[str] = None,
        n: Optional[int] = None,
        noise_rate: Optional[float] = None,
        alphabet: Optional[str] = None,
        seed: Optional[int] = None,
        rng: Optional[random.Random] = None,
    ) -> tuple[Dataset, dict]:
        """
        Generate sequences with noise variation from a base sequence.

        Creates n sequences by applying random noise to a base sequence.
        Each position has a probability of noise_rate to be replaced
        with a random character from the alphabet.

        Args:
            base_sequence (Optional[str]): Template sequence for variation.
                Auto-generated if None.
            n (Optional[int]): Number of sequences to generate. Defaults to 20.
            noise_rate (Optional[float]): Probability of noise per position (0.0-1.0).
                Defaults to 0.1.
            alphabet (Optional[str]): Character alphabet to use. Defaults to "ACGT".
            seed (Optional[int]): Random seed for reproducibility.
            rng (Optional[random.Random]): Existing random number generator.

        Returns:
            Tuple[Dataset, dict]: A tuple containing:
                - Dataset: Dataset with noisy sequences
                - dict: Generation parameters

        Raises:
            ValueError: If n <= 0, noise_rate not in [0,1], or base_sequence empty.
        """
        # Apply default values
        n = n if n is not None else SYNTHETIC_DEFAULTS["n"]
        noise_rate = (
            noise_rate if noise_rate is not None else SYNTHETIC_DEFAULTS["noise_rate"]
        )
        alphabet = alphabet if alphabet is not None else SYNTHETIC_DEFAULTS["alphabet"]

        # If base_sequence not provided, generate one automatically
        if base_sequence is None:
            length = SYNTHETIC_DEFAULTS["length"]
            _temp_rng = SyntheticDatasetGenerator._ensure_rng(seed, rng)
            base_sequence = "".join(_temp_rng.choice(alphabet) for _ in range(length))

        if n <= 0:
            raise ValueError("n must be > 0")
        if not 0.0 <= noise_rate <= 1.0:
            raise ValueError("noise_rate must be in [0,1]")
        if not base_sequence:
            raise ValueError("base_sequence empty")
        _rng = SyntheticDatasetGenerator._ensure_rng(seed, rng)

        def mutate_char(ch: str) -> str:
            if _rng.random() < noise_rate:
                choices = [c for c in alphabet if c != ch] or [ch]
                return _rng.choice(choices)
            return ch

        seqs = ["".join(mutate_char(ch) for ch in base_sequence) for _ in range(n)]
        ds = Dataset(
            id="synthetic_noise",
            name="synthetic_noise",
            sequences=seqs,
            alphabet=alphabet,
        )
        params = {
            "generator": "generate_with_noise",
            "n": n,
            "length": len(base_sequence),
            "alphabet": alphabet,
            "base_sequence": base_sequence,
            "noise_rate": noise_rate,
            "seed_arg": seed,
            "rng_provided": rng is not None,
        }
        return ds, params

    @staticmethod
    def generate_with_mutations(
        base_sequence: Optional[str] = None,
        n: Optional[int] = None,
        mutation_types: Optional[List[str]] = None,
        mutation_rate: Optional[float] = None,
        alphabet: Optional[str] = None,
        seed: Optional[int] = None,
        rng: Optional[random.Random] = None,
        pad_char: Optional[str] = None,
    ) -> tuple[Dataset, dict]:
        """
        Generate sequences with mutation operations from a base sequence.

        Creates n sequences by applying mutation operations (substitution,
        insertion, deletion) to a base sequence. Supports different
        mutation types and rates.

        Args:
            base_sequence (Optional[str]): Template sequence for mutation.
                Auto-generated if None.
            n (Optional[int]): Number of sequences to generate. Defaults to 20.
            mutation_types (Optional[List[str]]): Types of mutations to apply.
                Defaults to ["substitution"]. Options: ["substitution", "insertion", "deletion"].
            mutation_rate (Optional[float]): Probability of mutation per position (0.0-1.0).
                Defaults to 0.1.
            alphabet (Optional[str]): Character alphabet to use. Defaults to "ACGT".
            seed (Optional[int]): Random seed for reproducibility.
            rng (Optional[random.Random]): Existing random number generator.
            pad_char (Optional[str]): Character for padding sequences to uniform length.
                Defaults to "N".

        Returns:
            Tuple[Dataset, dict]: A tuple containing:
                - Dataset: Dataset with mutated sequences
                - dict: Generation parameters

        Raises:
            ValueError: If n <= 0, mutation_rate not in [0,1], or base_sequence empty.
        """
        # Apply default values
        n = n if n is not None else SYNTHETIC_DEFAULTS["n"]
        mutation_rate = (
            mutation_rate
            if mutation_rate is not None
            else SYNTHETIC_DEFAULTS["mutation_rate"]
        )
        alphabet = alphabet if alphabet is not None else SYNTHETIC_DEFAULTS["alphabet"]
        mutation_types = (
            mutation_types if mutation_types is not None else ["substitution"]
        )
        pad_char = pad_char if pad_char is not None else SYNTHETIC_DEFAULTS["pad_char"]

        # If base_sequence not provided, generate one automatically
        if base_sequence is None:
            length = SYNTHETIC_DEFAULTS["length"]
            _temp_rng = SyntheticDatasetGenerator._ensure_rng(seed, rng)
            base_sequence = "".join(_temp_rng.choice(alphabet) for _ in range(length))

        if n <= 0:
            raise ValueError("n must be > 0")
        if not base_sequence:
            raise ValueError("base_sequence empty")
        if mutation_rate < 0 or mutation_rate > 1:
            raise ValueError("mutation_rate must be in [0,1]")
        _rng = SyntheticDatasetGenerator._ensure_rng(seed, rng)

        def mutate(seq: str) -> str:
            s = list(seq)
            # Iterate with controlled index to support insertions/deletions safely
            i = 0
            ops = 0
            max_ops = max(10 * len(s) + 100, 200)  # safety limit for termination
            while i < len(s) and ops < max_ops:
                ops += 1
                if _rng.random() < mutation_rate:
                    m = _rng.choice(mutation_types)
                    if m == "substitution":
                        choices = [c for c in alphabet if c != s[i]] or [s[i]]
                        s[i] = _rng.choice(choices)
                        i += 1
                    elif m == "insertion":
                        s.insert(i, _rng.choice(alphabet))
                        # Advance two positions to ensure progress even with successive insertions
                        i += 2
                    elif m == "deletion":
                        if len(s) > 1:
                            s.pop(i)
                            # don't increment i to evaluate the new item that took the position
                        else:
                            i += 1
                    else:
                        # Unknown type: just advance
                        i += 1
                else:
                    i += 1
            return "".join(s)

        seqs = [mutate(base_sequence) for _ in range(n)]
        # Optionally normalize to largest length via pad to maintain valid Dataset
        L = max(len(s) for s in seqs)
        # Configurable pad_char; if used and not in alphabet, include it
        effective_pad = pad_char if pad_char else "N"
        seqs = [s + (effective_pad * (L - len(s))) for s in seqs]

        out_alphabet = alphabet
        if effective_pad and effective_pad not in out_alphabet:
            out_alphabet = out_alphabet + effective_pad

        ds = Dataset(
            id="synthetic_mutations",
            name="synthetic_mutations",
            sequences=seqs,
            alphabet=out_alphabet,
        )
        params = {
            "generator": "generate_with_mutations",
            "n": n,
            "length_base": len(base_sequence),
            "length_max": L,
            "alphabet": out_alphabet,
            "base_sequence": base_sequence,
            "mutation_types": mutation_types,
            "mutation_rate": mutation_rate,
            "pad_char": effective_pad,
            "seed_arg": seed,
            "rng_provided": rng is not None,
        }
        return ds, params

    @staticmethod
    def generate_clustered(
        n: Optional[int] = None,
        length: Optional[int] = None,
        num_clusters: Optional[int] = None,
        alphabet: Optional[str] = None,
        cluster_distance: Optional[float] = None,
        seed: Optional[int] = None,
        rng: Optional[random.Random] = None,
    ) -> tuple[Dataset, dict]:
        """
        Generate clustered sequences around multiple centers.

        Creates n sequences organized around num_clusters center sequences.
        Each sequence is generated by varying the assigned center sequence
        with a specified cluster_distance probability.

        Args:
            n (Optional[int]): Number of sequences to generate. Defaults to 20.
            length (Optional[int]): Length of each sequence. Defaults to 50.
            num_clusters (Optional[int]): Number of cluster centers. Defaults to 2.
            alphabet (Optional[str]): Character alphabet to use. Defaults to "ACGT".
            cluster_distance (Optional[float]): Probability of variation from cluster
                center (0.0-1.0). Defaults to 0.2.
            seed (Optional[int]): Random seed for reproducibility.
            rng (Optional[random.Random]): Existing random number generator.

        Returns:
            Tuple[Dataset, dict]: A tuple containing:
                - Dataset: Dataset with clustered sequences
                - dict: Generation parameters

        Raises:
            ValueError: If parameters are invalid or num_clusters > n.
        """
        # Apply default values
        n = n if n is not None else SYNTHETIC_DEFAULTS["n"]
        length = length if length is not None else SYNTHETIC_DEFAULTS["length"]
        alphabet = alphabet if alphabet is not None else SYNTHETIC_DEFAULTS["alphabet"]
        num_clusters = (
            num_clusters
            if num_clusters is not None
            else SYNTHETIC_DEFAULTS["num_clusters"]
        )
        cluster_distance = (
            cluster_distance
            if cluster_distance is not None
            else SYNTHETIC_DEFAULTS["cluster_distance"]
        )

        if n <= 0 or length <= 0:
            raise ValueError("n and length must be > 0")
        if num_clusters <= 0 or num_clusters > n:
            raise ValueError("num_clusters invalid")
        if not 0.0 <= cluster_distance <= 1.0:
            raise ValueError("cluster_distance must be in [0,1]")
        _rng = SyntheticDatasetGenerator._ensure_rng(seed, rng)
        # Create centers
        centers = [
            "".join(_rng.choice(alphabet) for _ in range(length))
            for _ in range(num_clusters)
        ]
        # Allocate by round-robin generating variations around centers
        seqs: List[str] = []
        for i in range(n):
            c = centers[i % num_clusters]
            p = cluster_distance
            s = "".join(_rng.choice(alphabet) if _rng.random() < p else ch for ch in c)
            seqs.append(s)

        ds = Dataset(
            id="synthetic_clustered",
            name="synthetic_clustered",
            sequences=seqs,
            alphabet=alphabet,
        )
        params = {
            "generator": "generate_clustered",
            "n": n,
            "length": length,
            "alphabet": alphabet,
            "num_clusters": num_clusters,
            "cluster_distance": cluster_distance,
            "seed_arg": seed,
            "rng_provided": rng is not None,
        }
        return ds, params

    @staticmethod
    def generate_from_config(
        cfg: SyntheticDatasetConfig, rng: Optional[random.Random] = None
    ) -> tuple[Dataset, dict]:
        """
        Generate Dataset from SyntheticDatasetConfig with support for multiple modes.

        Main entry point for configuration-based dataset generation.
        Supports multiple generation modes as specified in the configuration.

        Args:
            cfg (SyntheticDatasetConfig): Synthetic dataset configuration.
            rng (Optional[random.Random]): Optional random number generator
                for deterministic generation.

        Returns:
            Tuple[Dataset, dict]: A tuple containing:
                - Dataset: Generated dataset
                - dict: Generation parameters

        Raises:
            ValueError: If generation mode is invalid.

        Supported Modes:
            - random: Completely random sequences
            - noise: Noise-based variation from base sequence
            - mutations: Mutation-based evolution
            - clustered: Clustered sequence generation
        """
        mode = getattr(cfg, "mode", "random")
        alphabet = cfg.alphabet or "ACGT"
        seed = cfg.seed
        params_mode = getattr(cfg, "parameters_mode", {})
        base_rng = SyntheticDatasetGenerator._ensure_rng(seed, rng)
        dataset_name = getattr(cfg, "name", None)
        dataset_id = getattr(cfg, "id", f"synthetic_{mode}")

        if mode == "random":
            ds, p = SyntheticDatasetGenerator.generate_random(
                cfg.n, cfg.L, alphabet, seed=None, rng=base_rng
            )
            if dataset_name:
                ds.name = dataset_name
            ds.id = dataset_id
            p.update({"mode": mode, "config_seed": seed, "from_config": True})
            return ds, p

        if mode == "noise":
            # Extract parameters from parameters_mode
            base_sequence = params_mode.get("base_sequence")
            noise = params_mode.get("noise", 0.1)

            # No noise => behaves like random
            if not noise or noise <= 0:
                ds, p = SyntheticDatasetGenerator.generate_random(
                    cfg.n, cfg.L, alphabet, seed=None, rng=base_rng
                )
                if dataset_name:
                    ds.name = dataset_name
                ds.id = dataset_id
                p.update({"mode": "random", "config_seed": seed, "from_config": True})
                return ds, p

            # If base_sequence is None, generate random base sequence
            if base_sequence is None:
                ds_tmp, _ = SyntheticDatasetGenerator.generate_random(
                    1, cfg.L, alphabet, seed=None, rng=base_rng
                )
                base_sequence = ds_tmp.sequences[0]

            # Generate n sequences with noise using the same base RNG (unified)
            ds, p = SyntheticDatasetGenerator.generate_with_noise(
                base_sequence=base_sequence,
                n=int(cfg.n),
                noise_rate=float(noise),
                alphabet=alphabet,
                seed=None,
                rng=base_rng,
            )
            if dataset_name:
                ds.name = dataset_name
            ds.id = dataset_id
            p.update({"mode": mode, "config_seed": seed, "from_config": True})
            return ds, p

        if mode == "mutations":
            # Extract parameters from parameters_mode
            base_sequence = params_mode.get("base_sequence")
            mutation_types = params_mode.get("mutation_types", ["substitution"])
            mutation_rate = params_mode.get("mutation_rate", 0.1)
            pad_char = params_mode.get("pad_char", "N")

            if not base_sequence:
                # Generate random base sequence of size L
                ds_tmp, _ = SyntheticDatasetGenerator.generate_random(
                    1, cfg.L, alphabet, seed=None, rng=base_rng
                )
                base_sequence = ds_tmp.sequences[0]

            ds, p = SyntheticDatasetGenerator.generate_with_mutations(
                base_sequence=base_sequence,
                n=int(cfg.n),
                mutation_types=mutation_types,
                mutation_rate=float(mutation_rate),
                alphabet=alphabet,
                seed=None,
                rng=base_rng,
                pad_char=pad_char,
            )
            if dataset_name:
                ds.name = dataset_name
            ds.id = dataset_id
            p.update({"mode": mode, "config_seed": seed, "from_config": True})
            return ds, p

        if mode == "clustered":
            # Extract parameters from parameters_mode
            num_clusters = params_mode.get("num_clusters", 2)
            cluster_distance = params_mode.get("cluster_distance", 0.2)
            ds, p = SyntheticDatasetGenerator.generate_clustered(
                n=int(cfg.n),
                length=int(cfg.L),
                num_clusters=int(num_clusters),
                alphabet=alphabet,
                cluster_distance=float(cluster_distance),
                seed=None,
                rng=base_rng,
            )
            if dataset_name:
                ds.name = dataset_name
            ds.id = dataset_id
            p.update({"mode": mode, "config_seed": seed, "from_config": True})
            return ds, p

        raise ValueError(f"Invalid generation mode: {mode!r}")
