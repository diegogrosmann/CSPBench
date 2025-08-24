"""Testes para submit_combinations."""

import pytest
import tempfile
from pathlib import Path
from src.infrastructure.persistence.work_state.core import WorkStatePersistence


class TestSubmitCombinations:
    """Testes para a funcionalidade submit_combinations."""

    @pytest.fixture
    def temp_db(self):
        """Cria um banco temporário para testes."""
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = Path(tmpdir) / "test.db"
            persistence = WorkStatePersistence(db_path)
            yield persistence

    def test_submit_combinations_inserts_all(self, temp_db):
        """Testa que submit_combinations insere todas as combinações fornecidas."""
        work_id = "test_work_123"

        combinations = [
            {
                "task_id": "task1",
                "dataset_id": "dataset1",
                "preset_id": "preset1",
                "algorithm_id": "algo1",
                "mode": "experiment",
                "total_sequences": 5,
            },
            {
                "task_id": "task1",
                "dataset_id": "dataset2",
                "preset_id": "preset1",
                "algorithm_id": "algo1",
                "mode": "experiment",
                "total_sequences": 3,
            },
            {
                "task_id": "task2",
                "dataset_id": "dataset1",
                "preset_id": "preset1",
                "algorithm_id": "algo2",
                "mode": "optimization",
                "total_sequences": 10,
            },
        ]

        # Inserir combinações
        inserted_count = temp_db.submit_combinations(work_id, combinations)

        # Verificar quantas foram inseridas
        assert inserted_count == 3

        # Verificar se todas estão no banco
        all_combinations = temp_db.get_combinations(work_id=work_id)
        assert len(all_combinations) == 3

        # Verificar dados específicos
        task1_combos = [c for c in all_combinations if c["task_id"] == "task1"]
        assert len(task1_combos) == 2

        task2_combos = [c for c in all_combinations if c["task_id"] == "task2"]
        assert len(task2_combos) == 1
        assert task2_combos[0]["mode"] == "optimization"

    def test_submit_combinations_ignores_duplicates(self, temp_db):
        """Testa que submit_combinations ignora duplicatas (INSERT OR IGNORE)."""
        work_id = "test_work_456"

        combination = {
            "task_id": "task1",
            "dataset_id": "dataset1",
            "preset_id": "preset1",
            "algorithm_id": "algo1",
            "mode": "experiment",
            "total_sequences": 5,
        }

        # Inserir primeira vez
        inserted_count1 = temp_db.submit_combinations(work_id, [combination])
        assert inserted_count1 == 1

        # Inserir novamente (duplicata)
        inserted_count2 = temp_db.submit_combinations(work_id, [combination])
        assert inserted_count2 == 0  # Nenhuma nova inserção

        # Verificar que há apenas uma no banco
        all_combinations = temp_db.get_combinations(work_id=work_id)
        assert len(all_combinations) == 1

    def test_submit_combinations_empty_list(self, temp_db):
        """Testa comportamento com lista vazia."""
        work_id = "test_work_empty"

        inserted_count = temp_db.submit_combinations(work_id, [])
        assert inserted_count == 0

        all_combinations = temp_db.get_combinations(work_id=work_id)
        assert len(all_combinations) == 0

    def test_submit_combinations_mixed_new_and_existing(self, temp_db):
        """Testa inserção mista: algumas novas, algumas existentes."""
        work_id = "test_work_mixed"

        # Inserir algumas combinações iniciais
        initial_combinations = [
            {
                "task_id": "task1",
                "dataset_id": "dataset1",
                "preset_id": "preset1",
                "algorithm_id": "algo1",
                "mode": "experiment",
                "total_sequences": 5,
            },
            {
                "task_id": "task1",
                "dataset_id": "dataset2",
                "preset_id": "preset1",
                "algorithm_id": "algo1",
                "mode": "experiment",
                "total_sequences": 3,
            },
        ]

        inserted_count1 = temp_db.submit_combinations(work_id, initial_combinations)
        assert inserted_count1 == 2

        # Inserir misto: uma existente, duas novas
        mixed_combinations = [
            # Esta já existe
            {
                "task_id": "task1",
                "dataset_id": "dataset1",
                "preset_id": "preset1",
                "algorithm_id": "algo1",
                "mode": "experiment",
                "total_sequences": 5,
            },
            # Estas são novas
            {
                "task_id": "task2",
                "dataset_id": "dataset1",
                "preset_id": "preset1",
                "algorithm_id": "algo1",
                "mode": "experiment",
                "total_sequences": 8,
            },
            {
                "task_id": "task1",
                "dataset_id": "dataset3",
                "preset_id": "preset1",
                "algorithm_id": "algo1",
                "mode": "experiment",
                "total_sequences": 2,
            },
        ]

        inserted_count2 = temp_db.submit_combinations(work_id, mixed_combinations)
        assert inserted_count2 == 2  # Apenas as duas novas

        # Verificar total no banco
        all_combinations = temp_db.get_combinations(work_id=work_id)
        assert len(all_combinations) == 4  # 2 iniciais + 2 novas

    def test_submit_combinations_default_values(self, temp_db):
        """Testa valores padrão para campos opcionais."""
        work_id = "test_work_defaults"

        # Combinação sem mode e total_sequences
        combination = {
            "task_id": "task1",
            "dataset_id": "dataset1",
            "preset_id": "preset1",
            "algorithm_id": "algo1",
        }

        inserted_count = temp_db.submit_combinations(work_id, [combination])
        assert inserted_count == 1

        all_combinations = temp_db.get_combinations(work_id=work_id)
        assert len(all_combinations) == 1

        combo = all_combinations[0]
        assert combo["mode"] == "experiment"  # valor padrão
        assert combo["total_sequences"] == 0  # valor padrão
        assert combo["status"] == "queued"  # status inicial
