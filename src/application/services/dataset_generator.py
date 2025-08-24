import random
from typing import List, Optional
from src.domain.config import SyntheticDatasetConfig
from src.domain.dataset import Dataset

# Valores padrão para geração sintética
SYNTHETIC_DEFAULTS = {
    "n": 20,
    "length": 50,
    "alphabet": "ACGT",
    "seed": None,
    "noise_rate": 0.1,
    "mutation_rate": 0.1,
    "mutation_types": None,  # será ["substitution"] se None
    "num_clusters": 2,
    "cluster_distance": 0.2,
    "base_sequence": None,  # será gerada automaticamente se None
    "pad_char": "N",
}


class SyntheticDatasetGenerator:
    """Gerador de datasets sintéticos (compatibilidade com testes existentes)."""

    @staticmethod
    def _ensure_rng(
        seed: Optional[int] = None, rng: Optional[random.Random] = None
    ) -> random.Random:
        """Retorna um RNG a partir de um existente ou criando um novo com base na seed.

        Prioriza o rng injetado. Se não houver, cria um novo usando a seed.
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
        # Aplicar valores padrão
        n = n if n is not None else SYNTHETIC_DEFAULTS["n"]
        length = length if length is not None else SYNTHETIC_DEFAULTS["length"]
        alphabet = alphabet if alphabet is not None else SYNTHETIC_DEFAULTS["alphabet"]

        if n <= 0 or length <= 0:
            raise ValueError("n e length devem ser > 0")
        if not alphabet:
            raise ValueError("alphabet não pode ser vazio")
        _rng = SyntheticDatasetGenerator._ensure_rng(seed, rng)
        seqs = ["".join(_rng.choice(alphabet) for _ in range(length)) for _ in range(n)]
        ds = Dataset(id="synthetic_random", name="synthetic_random", sequences=seqs, alphabet=alphabet)
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
        # Aplicar valores padrão
        n = n if n is not None else SYNTHETIC_DEFAULTS["n"]
        noise_rate = (
            noise_rate if noise_rate is not None else SYNTHETIC_DEFAULTS["noise_rate"]
        )
        alphabet = alphabet if alphabet is not None else SYNTHETIC_DEFAULTS["alphabet"]

        # Se base_sequence não foi fornecida, gerar uma automaticamente
        if base_sequence is None:
            length = SYNTHETIC_DEFAULTS["length"]
            _temp_rng = SyntheticDatasetGenerator._ensure_rng(seed, rng)
            base_sequence = "".join(_temp_rng.choice(alphabet) for _ in range(length))

        if n <= 0:
            raise ValueError("n deve ser > 0")
        if not 0.0 <= noise_rate <= 1.0:
            raise ValueError("noise_rate deve estar em [0,1]")
        if not base_sequence:
            raise ValueError("base_sequence vazia")
        _rng = SyntheticDatasetGenerator._ensure_rng(seed, rng)

        def mutate_char(ch: str) -> str:
            if _rng.random() < noise_rate:
                choices = [c for c in alphabet if c != ch] or [ch]
                return _rng.choice(choices)
            return ch

        seqs = ["".join(mutate_char(ch) for ch in base_sequence) for _ in range(n)]
        ds = Dataset(id="synthetic_noise", name="synthetic_noise", sequences=seqs, alphabet=alphabet)
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
        # Aplicar valores padrão
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

        # Se base_sequence não foi fornecida, gerar uma automaticamente
        if base_sequence is None:
            length = SYNTHETIC_DEFAULTS["length"]
            _temp_rng = SyntheticDatasetGenerator._ensure_rng(seed, rng)
            base_sequence = "".join(_temp_rng.choice(alphabet) for _ in range(length))

        if n <= 0:
            raise ValueError("n deve ser > 0")
        if not base_sequence:
            raise ValueError("base_sequence vazia")
        if mutation_rate < 0 or mutation_rate > 1:
            raise ValueError("mutation_rate deve estar em [0,1]")
        _rng = SyntheticDatasetGenerator._ensure_rng(seed, rng)

        def mutate(seq: str) -> str:
            s = list(seq)
            # Itera com índice controlado para suportar inserções/remoções com segurança
            i = 0
            ops = 0
            max_ops = max(10 * len(s) + 100, 200)  # limite de segurança para término
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
                        # Avança duas posições para garantir progresso mesmo com inserções sucessivas
                        i += 2
                    elif m == "deletion":
                        if len(s) > 1:
                            s.pop(i)
                            # não incrementa i para avaliar o novo item que ocupou a posição
                        else:
                            i += 1
                    else:
                        # Tipo desconhecido: apenas avança
                        i += 1
                else:
                    i += 1
            return "".join(s)

        seqs = [mutate(base_sequence) for _ in range(n)]
        # Uniformiza opcionalmente ao maior comprimento via pad para manter Dataset válido
        L = max(len(s) for s in seqs)
        # pad_char configurável; se usado e não estiver no alphabet, incluí-lo
        effective_pad = pad_char if pad_char else "N"
        seqs = [s + (effective_pad * (L - len(s))) for s in seqs]

        out_alphabet = alphabet
        if effective_pad and effective_pad not in out_alphabet:
            out_alphabet = out_alphabet + effective_pad

        ds = Dataset(id="synthetic_mutations", name="synthetic_mutations", sequences=seqs, alphabet=out_alphabet)
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
        # Aplicar valores padrão
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
            raise ValueError("n e length devem ser > 0")
        if num_clusters <= 0 or num_clusters > n:
            raise ValueError("num_clusters inválido")
        if not 0.0 <= cluster_distance <= 1.0:
            raise ValueError("cluster_distance deve estar em [0,1]")
        _rng = SyntheticDatasetGenerator._ensure_rng(seed, rng)
        # cria centros
        centers = [
            "".join(_rng.choice(alphabet) for _ in range(length))
            for _ in range(num_clusters)
        ]
        # aloca por round-robin gerando variações ao redor dos centros
        seqs: List[str] = []
        for i in range(n):
            c = centers[i % num_clusters]
            p = cluster_distance
            s = "".join(_rng.choice(alphabet) if _rng.random() < p else ch for ch in c)
            seqs.append(s)

        ds = Dataset(id="synthetic_clustered", name="synthetic_clustered", sequences=seqs, alphabet=alphabet)
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
        """Gera Dataset a partir de SyntheticDatasetConfig (suporta modos)."""
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
            # Extrai parâmetros do parameters_mode
            base_sequence = params_mode.get("base_sequence")
            noise = params_mode.get("noise", 0.1)

            # Sem ruído => comporta-se como random
            if not noise or noise <= 0:
                ds, p = SyntheticDatasetGenerator.generate_random(
                    cfg.n, cfg.L, alphabet, seed=None, rng=base_rng
                )
                if dataset_name:
                    ds.name = dataset_name
                ds.id = dataset_id
                p.update({"mode": "random", "config_seed": seed, "from_config": True})
                return ds, p

            # Se base_sequence for None, gera sequência aleatória base
            if base_sequence is None:
                ds_tmp, _ = SyntheticDatasetGenerator.generate_random(
                    1, cfg.L, alphabet, seed=None, rng=base_rng
                )
                base_sequence = ds_tmp.sequences[0]

            # Gera n sequências com ruído usando o mesmo RNG base (uniformizado)
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
            # Extrai parâmetros do parameters_mode
            base_sequence = params_mode.get("base_sequence")
            mutation_types = params_mode.get("mutation_types", ["substitution"])
            mutation_rate = params_mode.get("mutation_rate", 0.1)
            pad_char = params_mode.get("pad_char", "N")

            if not base_sequence:
                # Gera sequência base aleatória de tamanho L
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
            # Extrai parâmetros do parameters_mode
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

        raise ValueError(f"Modo de geração inválido: {mode!r}")
