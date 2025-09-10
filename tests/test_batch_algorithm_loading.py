"""
Testes para carregamento de batch com algoritmos comentados.

Este módulo testa se o sistema corretamente carrega apenas algoritmos
não comentados nos arquivos YAML de configuração de batch.
"""

from pathlib import Path
from typing import Any, Dict, List

import pytest

from src.domain.config import CSPBenchConfig
from src.domain.errors import BatchConfigurationError


class TestBatchAlgorithmLoading:
    """Testes para carregamento de algoritmos em configurações de batch."""

    def get_test_batch_path(self, filename: str) -> Path:
        """Retorna o caminho para um arquivo de teste de batch."""
        return Path(__file__).parent / "batches" / filename

    def load_batch_config(self, filename: str) -> CSPBenchConfig:
        """Carrega uma configuração de batch de teste."""
        import tempfile

        from src.domain.config import load_cspbench_config

        # Criar um settings.yaml básico temporário
        settings_content = """
system:
  global_seed: 123
  max_workers: 1
  memory_limit: 1024

output:
  base_path: './'
  format: 'json'

resources:
  cpu_limit: 1
  memory_limit: 1024
"""

        batch_path = self.get_test_batch_path(filename)
        assert batch_path.exists(), f"Arquivo de teste não encontrado: {batch_path}"

        # Criar arquivo temporário para settings
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            f.write(settings_content)
            settings_path = f.name

        try:
            config = load_cspbench_config(
                batch_path=str(batch_path), settings_path=settings_path
            )
            return config
        finally:
            import os

            os.unlink(settings_path)

    def get_algorithms_from_preset(
        self, config: CSPBenchConfig, preset_id: str
    ) -> List[str]:
        """Extrai a lista de algoritmos de um preset de configuração."""
        preset = config.algorithms.get(preset_id)
        assert preset is not None, f"Preset '{preset_id}' não encontrado"
        return [item.name for item in preset.items]

    def test_commented_algorithms_are_ignored(self):
        """Testa se algoritmos comentados são ignorados no carregamento."""
        config = self.load_batch_config("batch_test_commented_algorithms.yaml")

        # Obter a lista de algoritmos do preset
        algorithms = self.get_algorithms_from_preset(config, "test_commented")

        # Verificar que apenas algoritmos não comentados foram carregados
        # Baseado no YAML: apenas "BLF-GA" e "H2-CSP" devem estar presentes
        expected_algorithms = ["BLF-GA", "H2-CSP"]

        assert set(algorithms) == set(expected_algorithms), (
            f"Algoritmos carregados: {algorithms}, " f"Esperados: {expected_algorithms}"
        )  # Verificar que algoritmos comentados foram ignorados
        commented_algorithms = ["Baseline", "CSC", "DP-CSP"]
        for alg in commented_algorithms:
            assert (
                alg not in algorithms
            ), f"Algoritmo comentado '{alg}' foi carregado incorretamente"

    def test_all_enabled_algorithms_are_loaded(self):
        """Testa se todos os algoritmos são carregados quando não comentados."""
        config = self.load_batch_config("batch_test_all_enabled.yaml")

        algorithms = self.get_algorithms_from_preset(config, "test_all_enabled")

        # Todos os algoritmos devem estar presentes
        expected_algorithms = ["Baseline", "BLF-GA", "CSC", "H2-CSP", "DP-CSP"]

        assert set(algorithms) == set(expected_algorithms), (
            f"Algoritmos carregados: {algorithms}, " f"Esperados: {expected_algorithms}"
        )

    def test_all_commented_algorithms_result_in_empty_list(self):
        """Testa se todos os algoritmos comentados resultam em lista vazia."""
        config = self.load_batch_config("batch_test_all_commented.yaml")

        algorithms = self.get_algorithms_from_preset(config, "test_all_commented")

        # A lista deve estar vazia
        assert len(algorithms) == 0, f"Esperava lista vazia, mas obteve: {algorithms}"

    def test_mixed_algorithms_configuration(self):
        """Testa configurações mistas com múltiplos presets."""
        config = self.load_batch_config("batch_test_mixed_algorithms.yaml")

        # Testar primeiro preset
        algorithms_1 = self.get_algorithms_from_preset(config, "config_mixed_1")
        expected_1 = ["Baseline", "CSC", "DP-CSP"]  # BLF-GA e H²-CSP comentados

        assert set(algorithms_1) == set(expected_1), (
            f"Preset 1 - Algoritmos carregados: {algorithms_1}, "
            f"Esperados: {expected_1}"
        )

        # Testar segundo preset
        algorithms_2 = self.get_algorithms_from_preset(config, "config_mixed_2")
        expected_2 = ["BLF-GA", "H2-CSP"]  # Baseline, CSC e DP-CSP comentados

        assert set(algorithms_2) == set(expected_2), (
            f"Preset 2 - Algoritmos carregados: {algorithms_2}, "
            f"Esperados: {expected_2}"
        )

    def test_no_algorithm_list_loads_all_params(self):
        """Testa se ausência da lista algorithms carrega todos os algorithm_params."""
        config = self.load_batch_config("batch_test_no_algorithm_list.yaml")

        algorithms = self.get_algorithms_from_preset(config, "test_no_algorithm_list")

        # Deve carregar todos os algoritmos definidos em algorithm_params
        expected_algorithms = ["Baseline", "BLF-GA", "CSC"]

        assert set(algorithms) == set(expected_algorithms), (
            f"Algoritmos carregados: {algorithms}, " f"Esperados: {expected_algorithms}"
        )

    def test_empty_algorithm_list_results_in_empty_preset(self):
        """Testa se lista explicitamente vazia resulta em preset vazio."""
        config = self.load_batch_config("batch_test_empty_algorithm_list.yaml")

        algorithms = self.get_algorithms_from_preset(
            config, "test_empty_algorithm_list"
        )

        # A lista deve estar vazia
        assert len(algorithms) == 0, f"Esperava lista vazia, mas obteve: {algorithms}"

    def test_algorithm_params_are_preserved(self):
        """Testa se os parâmetros dos algoritmos são preservados corretamente."""
        config = self.load_batch_config("batch_test_commented_algorithms.yaml")

        preset = config.algorithms["test_commented"]

        # Verificar que apenas algoritmos não comentados têm parâmetros carregados
        loaded_alg_names = [item.name for item in preset.items]

        for item in preset.items:
            assert item.name in [
                "BLF-GA",
                "H2-CSP",
            ], f"Algoritmo inesperado carregado: {item.name}"
            assert isinstance(
                item.params, dict
            ), f"Parâmetros do algoritmo {item.name} devem ser um dicionário"
            assert (
                len(item.params) > 0
            ), f"Algoritmo {item.name} deve ter parâmetros definidos"

    def test_preset_metadata_is_preserved(self):
        """Testa se metadados do preset são preservados durante o carregamento."""
        config = self.load_batch_config("batch_test_commented_algorithms.yaml")

        preset = config.algorithms["test_commented"]

        assert preset.id == "test_commented"
        assert preset.name == "Configuration with Commented Algorithms"
        assert "Testing algorithms" in preset.description
        assert isinstance(preset.items, list)

    def test_batch_execution_with_commented_algorithms(self):
        """Testa se a execução do batch funciona com algoritmos comentados."""
        config = self.load_batch_config("batch_test_commented_algorithms.yaml")

        # Verificar que a configuração de experimento está correta
        assert config.tasks.type == "experiment"

        # Verificar que apenas algoritmos não comentados serão executados
        experiment_tasks = config.tasks.items
        assert len(experiment_tasks) == 1

        task = experiment_tasks[0]
        assert "test_commented" in task.algorithms

        # Verificar que o preset existe e tem algoritmos válidos
        preset = config.algorithms["test_commented"]
        assert len(preset.items) > 0, "Preset deve ter pelo menos um algoritmo válido"

    def test_mixed_datasets_with_filtered_algorithms(self):
        """Testa execução com múltiplos datasets e algoritmos filtrados."""
        config = self.load_batch_config("batch_test_mixed_algorithms.yaml")

        # Verificar múltiplos datasets
        dataset_ids = list(config.datasets.keys())
        assert "dataset_test_1" in dataset_ids
        assert "dataset_test_2" in dataset_ids

        # Verificar que ambos os presets têm algoritmos válidos
        preset_1 = config.algorithms["config_mixed_1"]
        preset_2 = config.algorithms["config_mixed_2"]

        assert len(preset_1.items) > 0, "Preset 1 deve ter algoritmos válidos"
        assert len(preset_2.items) > 0, "Preset 2 deve ter algoritmos válidos"

        # Verificar que algoritmos são diferentes entre presets
        algs_1 = {item.name for item in preset_1.items}
        algs_2 = {item.name for item in preset_2.items}

        assert algs_1 != algs_2, "Presets devem ter conjuntos diferentes de algoritmos"

    @pytest.mark.parametrize(
        "filename,expected_error",
        [
            # Adicionar casos de erro se necessário
        ],
    )
    def test_batch_loading_error_cases(self, filename: str, expected_error: str):
        """Testa casos de erro no carregamento de batch."""
        # Implementar se houver casos de erro específicos
        pass

    def test_yaml_comment_parsing_integrity(self):
        """Testa se o parser YAML não interpreta comentários como valores."""
        config = self.load_batch_config("batch_test_commented_algorithms.yaml")

        # Carregar o YAML diretamente para verificar parsing
        import yaml

        yaml_path = self.get_test_batch_path("batch_test_commented_algorithms.yaml")

        with open(yaml_path) as f:
            raw_data = yaml.safe_load(f)

        # Verificar que a lista algorithms no YAML não contém elementos comentados
        preset_data = None
        for preset in raw_data.get("algorithms", []):
            if preset.get("id") == "test_commented":
                preset_data = preset
                break

        assert preset_data is not None, "Preset de teste não encontrado no YAML"

        # A lista algorithms deve conter apenas elementos não comentados
        yaml_algorithms = preset_data.get("algorithms", [])

        # Verificar que comentários não aparecem como valores
        for alg in yaml_algorithms:
            assert not alg.startswith("#"), f"Comentário encontrado como valor: {alg}"
            assert alg in ["BLF-GA", "H2-CSP"], f"Algoritmo inesperado no YAML: {alg}"

    def test_performance_with_large_configuration(self):
        """Testa performance do carregamento com configurações maiores."""
        # Usar o arquivo mixed que tem mais configurações
        import time

        start_time = time.time()
        config = self.load_batch_config("batch_test_mixed_algorithms.yaml")
        load_time = time.time() - start_time

        # Carregamento deve ser rápido (menos de 1 segundo)
        assert load_time < 1.0, f"Carregamento demorou muito: {load_time:.2f}s"

        # Verificar que a configuração foi carregada corretamente
        assert len(config.algorithms) == 2, "Deve ter 2 presets carregados"
        assert len(config.datasets) == 2, "Deve ter 2 datasets carregados"

    def test_algorithm_parameter_consistency(self):
        """Testa se parâmetros são consistentes entre diferentes configurações."""
        # Testar múltiplos arquivos para verificar consistência
        configs = [
            self.load_batch_config("batch_test_commented_algorithms.yaml"),
            self.load_batch_config("batch_test_all_enabled.yaml"),
            self.load_batch_config("batch_test_mixed_algorithms.yaml"),
        ]

        # Verificar que parâmetros do BLF-GA são consistentes onde definidos
        blfga_params = []
        for config in configs:
            for preset in config.algorithms.values():
                for item in preset.items:
                    if item.name == "BLF-GA":
                        blfga_params.append(item.params)

        assert (
            len(blfga_params) > 0
        ), "Deve encontrar pelo menos um conjunto de parâmetros BLF-GA"

        # Verificar que todos têm campos obrigatórios
        for params in blfga_params:
            assert "pop_size" in params, "BLF-GA deve ter parâmetro pop_size"
            assert "max_gens" in params, "BLF-GA deve ter parâmetro max_gens"
            assert isinstance(params["pop_size"], int), "pop_size deve ser inteiro"
            assert isinstance(params["max_gens"], int), "max_gens deve ser inteiro"
