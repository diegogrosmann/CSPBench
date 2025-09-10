import random

import pytest

from src.application.services.dataset_generator import SyntheticDatasetGenerator
from src.domain.config import SyntheticDatasetConfig


def _same_dataset(d1, d2):
    # d1 e d2 são tuplas (Dataset, dict) - extrair apenas o Dataset
    dataset1 = d1[0] if isinstance(d1, tuple) else d1
    dataset2 = d2[0] if isinstance(d2, tuple) else d2
    return (
        dataset1.sequences == dataset2.sequences
        and dataset1.alphabet == dataset2.alphabet
    )


def test_generate_random_reproducible_by_seed():
    ds1 = SyntheticDatasetGenerator.generate_random(
        n=5, length=8, alphabet="ACGT", seed=42
    )
    ds2 = SyntheticDatasetGenerator.generate_random(
        n=5, length=8, alphabet="ACGT", seed=42
    )
    assert _same_dataset(ds1, ds2)


def test_generate_random_diff_with_diff_seed():
    ds1 = SyntheticDatasetGenerator.generate_random(
        n=5, length=8, alphabet="ACGT", seed=1
    )
    ds2 = SyntheticDatasetGenerator.generate_random(
        n=5, length=8, alphabet="ACGT", seed=2
    )
    assert ds1[0].sequences != ds2[0].sequences


@pytest.mark.parametrize("noise_rate", [0.0, 0.3, 1.0])
def test_generate_with_noise_reproducibility_varied_noise(noise_rate):
    base = "ACGTACGT"
    ds1 = SyntheticDatasetGenerator.generate_with_noise(
        base_sequence=base, n=6, noise_rate=noise_rate, alphabet="ACGT", seed=123
    )
    ds2 = SyntheticDatasetGenerator.generate_with_noise(
        base_sequence=base, n=6, noise_rate=noise_rate, alphabet="ACGT", seed=123
    )
    assert _same_dataset(ds1, ds2)
    if noise_rate == 0.0:
        assert all(seq == base for seq in ds1[0].sequences)


@pytest.mark.parametrize("seed", [1, 42])
@pytest.mark.parametrize("n,length", [(1, 4), (5, 8)])
@pytest.mark.parametrize("alphabet", ["ACGT", "ACGTN"])
def test_generate_random_parametrized(seed, n, length, alphabet):
    # seed -> determinístico
    a = SyntheticDatasetGenerator.generate_random(
        n=n, length=length, alphabet=alphabet, seed=seed
    )
    b = SyntheticDatasetGenerator.generate_random(
        n=n, length=length, alphabet=alphabet, seed=seed
    )
    assert _same_dataset(a, b)
    # rng injetado -> determinístico
    rng1 = random.Random(seed)
    rng2 = random.Random(seed)
    c = SyntheticDatasetGenerator.generate_random(
        n=n, length=length, alphabet=alphabet, seed=None, rng=rng1
    )
    d = SyntheticDatasetGenerator.generate_random(
        n=n, length=length, alphabet=alphabet, seed=None, rng=rng2
    )
    assert _same_dataset(c, d)


@pytest.mark.parametrize(
    "noise_rate", [0.0, 0.1, 0.5, 1.0], ids=["nr0", "nr01", "nr05", "nr1"]
)
@pytest.mark.parametrize("base_len", [8, 12])
@pytest.mark.parametrize("n", [1, 6])
def test_generate_with_noise_parametrized(noise_rate, base_len, n):
    base = SyntheticDatasetGenerator.generate_random(
        1, base_len, alphabet="ACGT", seed=77
    )[0].sequences[0]
    a = SyntheticDatasetGenerator.generate_with_noise(
        base_sequence=base, n=n, noise_rate=noise_rate, alphabet="ACGT", seed=123
    )
    b = SyntheticDatasetGenerator.generate_with_noise(
        base_sequence=base, n=n, noise_rate=noise_rate, alphabet="ACGT", seed=123
    )
    assert _same_dataset(a, b)
    if noise_rate == 0.0:
        assert all(seq == base for seq in a[0].sequences)


def test_generate_with_mutations_reproducibility_and_pad_alphabet():
    base = "AAAAA"
    ds1 = SyntheticDatasetGenerator.generate_with_mutations(
        base_sequence=base,
        n=8,
        mutation_types=["substitution", "insertion", "deletion"],
        mutation_rate=0.5,
        alphabet="ACGT",
        seed=999,
        pad_char="Z",
    )
    ds2 = SyntheticDatasetGenerator.generate_with_mutations(
        base_sequence=base,
        n=8,
        mutation_types=["substitution", "insertion", "deletion"],
        mutation_rate=0.5,
        alphabet="ACGT",
        seed=999,
        pad_char="Z",
    )
    assert _same_dataset(ds1, ds2)
    assert "Z" in ds1[0].alphabet

    # Reprodutibilidade com RNG injetado (dois RNGs com a mesma seed inicial)
    rng1 = random.Random(321)
    rng2 = random.Random(321)
    ds3 = SyntheticDatasetGenerator.generate_with_mutations(
        base_sequence=base,
        n=8,
        mutation_types=["substitution", "insertion", "deletion"],
        mutation_rate=0.5,
        alphabet="ACGT",
        seed=None,
        rng=rng1,
        pad_char="Z",
    )
    ds4 = SyntheticDatasetGenerator.generate_with_mutations(
        base_sequence=base,
        n=8,
        mutation_types=["substitution", "insertion", "deletion"],
        mutation_rate=0.5,
        alphabet="ACGT",
        seed=None,
        rng=rng2,
        pad_char="Z",
    )
    assert _same_dataset(ds3, ds4)


@pytest.mark.parametrize(
    "mutation_types",
    [
        ["substitution"],
        ["insertion"],
        ["deletion"],
        ["substitution", "insertion", "deletion"],
    ],
    ids=["subs", "ins", "del", "all"],
)
@pytest.mark.parametrize("mutation_rate", [0.0, 0.2, 0.5, 1.0])
@pytest.mark.parametrize("pad_char", [None, "N", "Z"], ids=["padNone", "padN", "padZ"])
def test_generate_with_mutations_parametrized(mutation_types, mutation_rate, pad_char):
    base = "A" * 7
    ds1 = SyntheticDatasetGenerator.generate_with_mutations(
        base_sequence=base,
        n=6,
        mutation_types=mutation_types,
        mutation_rate=mutation_rate,
        alphabet="ACGT",
        seed=555,
        pad_char=pad_char,
    )
    ds2 = SyntheticDatasetGenerator.generate_with_mutations(
        base_sequence=base,
        n=6,
        mutation_types=mutation_types,
        mutation_rate=mutation_rate,
        alphabet="ACGT",
        seed=555,
        pad_char=pad_char,
    )
    assert _same_dataset(ds1, ds2)
    effective_pad = pad_char if pad_char else "N"
    assert effective_pad in ds1[0].alphabet


def test_generate_clustered_reproducibility_and_bounds():
    for p in [0.0, 0.4, 1.0]:
        ds1 = SyntheticDatasetGenerator.generate_clustered(
            n=12, length=10, num_clusters=3, alphabet="ACGT", cluster_distance=p, seed=7
        )
        ds2 = SyntheticDatasetGenerator.generate_clustered(
            n=12, length=10, num_clusters=3, alphabet="ACGT", cluster_distance=p, seed=7
        )
        assert _same_dataset(ds1, ds2)

    with pytest.raises(ValueError):
        SyntheticDatasetGenerator.generate_clustered(
            n=10,
            length=5,
            num_clusters=2,
            alphabet="ACGT",
            cluster_distance=-0.1,
            seed=1,
        )
    with pytest.raises(ValueError):
        SyntheticDatasetGenerator.generate_clustered(
            n=10,
            length=5,
            num_clusters=2,
            alphabet="ACGT",
            cluster_distance=1.1,
            seed=1,
        )


@pytest.mark.parametrize("cluster_distance", [0.0, 0.2, 0.5, 1.0])
@pytest.mark.parametrize("num_clusters", [1, 2, 4])
@pytest.mark.parametrize("length", [5, 10])
def test_generate_clustered_parametrized(cluster_distance, num_clusters, length):
    n = 12
    ds1 = SyntheticDatasetGenerator.generate_clustered(
        n=n,
        length=length,
        num_clusters=num_clusters,
        alphabet="ACGT",
        cluster_distance=cluster_distance,
        seed=77,
    )
    ds2 = SyntheticDatasetGenerator.generate_clustered(
        n=n,
        length=length,
        num_clusters=num_clusters,
        alphabet="ACGT",
        cluster_distance=cluster_distance,
        seed=77,
    )
    assert _same_dataset(ds1, ds2)


def test_generate_from_config_random_reproducibility():
    cfg = SyntheticDatasetConfig(
        id="ds1", name="RandomDS", mode="random", n=5, L=8, alphabet="ACGT", seed=77
    )
    ds1 = SyntheticDatasetGenerator.generate_from_config(cfg)
    ds2 = SyntheticDatasetGenerator.generate_from_config(cfg)
    assert _same_dataset(ds1, ds2)


def test_generate_from_config_noise_reproducibility_base_provided_and_auto():
    # base_sequence fornecida
    cfg1 = SyntheticDatasetConfig(
        id="ds2",
        name="NoiseProvided",
        mode="noise",
        n=6,
        L=8,
        alphabet="ACGT",
        seed=101,
        parameters_mode={"base_sequence": "ACGTACGT", "noise": 0.2},
    )
    a1 = SyntheticDatasetGenerator.generate_from_config(cfg1)
    a2 = SyntheticDatasetGenerator.generate_from_config(cfg1)
    assert _same_dataset(a1, a2)

    # base_sequence gerada a partir do RNG do config
    cfg2 = SyntheticDatasetConfig(
        id="ds3",
        name="NoiseAutoBase",
        mode="noise",
        n=6,
        L=8,
        alphabet="ACGT",
        seed=202,
        parameters_mode={"noise": 0.2},
    )
    b1 = SyntheticDatasetGenerator.generate_from_config(cfg2)
    b2 = SyntheticDatasetGenerator.generate_from_config(cfg2)
    assert _same_dataset(b1, b2)


def test_generate_from_config_mutations_reproducibilidade_and_pad():
    cfg = SyntheticDatasetConfig(
        id="ds4",
        name="Mutations",
        mode="mutations",
        n=7,
        L=10,
        alphabet="ACGT",
        seed=303,
        parameters_mode={
            "base_sequence": "AAAAAAA",
            "mutation_types": ["substitution", "insertion", "deletion"],
            "mutation_rate": 0.3,
            "pad_char": "Z",
        },
    )
    d1 = SyntheticDatasetGenerator.generate_from_config(cfg)
    d2 = SyntheticDatasetGenerator.generate_from_config(cfg)
    assert _same_dataset(d1, d2)
    assert "Z" in d1[0].alphabet


def test_generate_from_config_clustered_reproducibility():
    cfg = SyntheticDatasetConfig(
        id="ds5",
        name="Clustered",
        mode="clustered",
        n=15,
        L=12,
        alphabet="ACGT",
        seed=404,
        parameters_mode={"num_clusters": 3, "cluster_distance": 0.35},
    )
    e1 = SyntheticDatasetGenerator.generate_from_config(cfg)
    e2 = SyntheticDatasetGenerator.generate_from_config(cfg)
    assert _same_dataset(e1, e2)


@pytest.mark.parametrize(
    "mode,params_mode,n,L,alphabet,seed",
    [
        ("random", {}, 5, 8, "ACGT", 11),
        ("noise", {"base_sequence": "ACGTACGT", "noise": 0.2}, 6, 8, "ACGT", 22),
        ("noise", {"noise": 0.3}, 6, 10, "ACGT", 33),  # base auto
        (
            "mutations",
            {
                "base_sequence": "AAAAAA",
                "mutation_types": ["substitution"],
                "mutation_rate": 0.4,
                "pad_char": "N",
            },
            7,
            10,
            "ACGT",
            44,
        ),
        (
            "mutations",
            {
                "base_sequence": "AAAAAA",
                "mutation_types": ["substitution", "insertion", "deletion"],
                "mutation_rate": 0.5,
                "pad_char": "Z",
            },
            7,
            10,
            "ACGT",
            55,
        ),
        (
            "clustered",
            {"num_clusters": 3, "cluster_distance": 0.35},
            12,
            12,
            "ACGT",
            66,
        ),
    ],
)
def test_generate_from_config_parametrized(mode, params_mode, n, L, alphabet, seed):
    cfg = SyntheticDatasetConfig(
        id=f"cfg_{mode}",
        name="Cfg",
        mode=mode,
        n=n,
        L=L,
        alphabet=alphabet,
        seed=seed,
        parameters_mode=params_mode,
    )
    r1 = SyntheticDatasetGenerator.generate_from_config(cfg)
    r2 = SyntheticDatasetGenerator.generate_from_config(cfg)
    assert _same_dataset(r1, r2)


def test_rng_injection_uniform_across_modes():
    # Usar rng injetado e seed=None no config para todos os modos
    for mode, pm in [
        ("random", {}),
        ("noise", {"base_sequence": "ACGTACGT", "noise": 0.25}),
        (
            "mutations",
            {
                "base_sequence": "AAAAAA",
                "mutation_types": ["substitution"],
                "mutation_rate": 0.4,
                "pad_char": "N",
            },
        ),
        ("clustered", {"num_clusters": 2, "cluster_distance": 0.2}),
    ]:
        cfg = SyntheticDatasetConfig(
            id=f"rng_{mode}",
            name=f"RNG {mode}",
            mode=mode,
            n=8,
            L=8,
            alphabet="ACGT",
            seed=None,
            parameters_mode=pm,
        )
        rng1 = random.Random(555)
        rng2 = random.Random(555)
        d1 = SyntheticDatasetGenerator.generate_from_config(cfg, rng=rng1)
        d2 = SyntheticDatasetGenerator.generate_from_config(cfg, rng=rng2)
        assert _same_dataset(d1, d2)


@pytest.mark.parametrize(
    "mode,params_mode",
    [
        ("random", {}),
        ("noise", {"base_sequence": "ACGTACGT", "noise": 0.25}),
        (
            "mutations",
            {
                "base_sequence": "AAAAAA",
                "mutation_types": ["substitution"],
                "mutation_rate": 0.4,
                "pad_char": "N",
            },
        ),
        ("clustered", {"num_clusters": 2, "cluster_distance": 0.2}),
    ],
)
@pytest.mark.parametrize("seed_rng", [111, 555])
def test_rng_injection_parametrized_across_modes(mode, params_mode, seed_rng):
    cfg = SyntheticDatasetConfig(
        id=f"rng_{mode}",
        name=f"RNG {mode}",
        mode=mode,
        n=8,
        L=8,
        alphabet="ACGT",
        seed=None,
        parameters_mode=params_mode,
    )
    rng1 = random.Random(seed_rng)
    rng2 = random.Random(seed_rng)
    d1 = SyntheticDatasetGenerator.generate_from_config(cfg, rng=rng1)
    d2 = SyntheticDatasetGenerator.generate_from_config(cfg, rng=rng2)
    assert _same_dataset(d1, d2)
