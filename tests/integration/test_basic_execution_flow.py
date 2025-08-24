"""Teste de fluxo básico de execução."""

import pytest
import tempfile
import time
from pathlib import Path
from unittest.mock import Mock, patch

from src.application.services.execution_manager import ExecutionManager
from src.domain.config import (
    CSPBenchConfig,
    MetadataConfig,
    SyntheticDatasetConfig,
    AlgParamsConfig,
    AlgorithmsPresetConfig,
    ExperimentTaskConfig,
    TasksGroupConfig,
    SystemConfig,
)
from src.domain.status import BaseStatus
from src.infrastructure.persistence.work_state import WorkStatePersistence


class TestBasicExecutionFlow:
    """Testes para fluxo básico de execução do pipeline."""

    @pytest.fixture
    def minimal_config(self):
        """Cria uma configuração mínima para testes."""
        return CSPBenchConfig(
            metadata=MetadataConfig(
                name="Test Config",
                description="Configuração mínima para testes",
                author="Test Suite",
                version="1.0.0",
                creation_date="2025-01-01",
                tags=["test"],
            ),
            datasets={
                "test_dataset": SyntheticDatasetConfig(
                    id="test_dataset",
                    name="Test Synthetic Dataset",
                    type="synthetic",
                    mode="random",
                    n=3,  # 3 sequências
                    L=10,  # comprimento 10
                    alphabet="ACGT",
                    seed=42,
                )
            },
            algorithms={
                "baseline_preset": AlgorithmsPresetConfig(
                    id="baseline_preset",
                    name="Baseline Algorithms",
                    description="Preset com algoritmo baseline",
                    items=[AlgParamsConfig(name="baseline", params={})],
                )
            },
            tasks=TasksGroupConfig(
                type="experiment",
                items=[
                    ExperimentTaskConfig(
                        id="test_task",
                        name="Test Experiment",
                        type="experiment",
                        datasets=["test_dataset"],
                        algorithms=["baseline_preset"],
                        repetitions=2,
                    )
                ],
            ),
            system=SystemConfig(
                global_seed=42, distance_method="hamming", enable_distance_cache=True
            ),
        )

    @pytest.fixture
    def temp_work_dir(self):
        """Cria diretório temporário para work."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_execution_manager_validates_config(self, minimal_config):
        """Testa que ExecutionManager valida configuração antes da execução."""
        execution_manager = ExecutionManager()

        # Configuração válida não deve levantar exceção
        try:
            execution_manager._validate_config(minimal_config)
        except Exception as e:
            pytest.fail(f"Configuração válida foi rejeitada: {e}")

        # Configuração sem datasets deve falhar
        invalid_config = CSPBenchConfig(
            metadata=minimal_config.metadata,
            datasets={},  # vazio
            algorithms=minimal_config.algorithms,
            tasks=minimal_config.tasks,
            system=minimal_config.system,
        )

        with pytest.raises(ValueError, match="At least one dataset"):
            execution_manager._validate_config(invalid_config)

    def test_pipeline_runner_generates_combinations(
        self, minimal_config, temp_work_dir
    ):
        """Testa que PipelineRunner gera combinações corretamente."""
        from src.infrastructure.orchestration.pipeline_runner import PipelineRunner
        from src.infrastructure.persistence.work_state.wrappers.work_scoped import (
            WorkScopedPersistence,
        )
        from src.domain.work import WorkItem

        # Criar work item e persistence
        work_id = "test_work_123"
        work_item = WorkItem(
            id=work_id,
            config=minimal_config,
            status=BaseStatus.QUEUED,
            created_at=time.time(),
            updated_at=time.time(),
            output_path=str(temp_work_dir / "results"),
        )

        # Inicializar persistence
        base_store = WorkStatePersistence(temp_work_dir / "state.db")
        base_store.submit_work(work_item)
        work_store = WorkScopedPersistence(base_store, work_id)

        # Criar runner e gerar combinações
        runner = PipelineRunner(work_store=work_store)
        runner._generate_pipeline_combinations(minimal_config)

        # Verificar combinações geradas
        combinations = work_store.get_combinations()

        # Devemos ter 1 combinação: 1 task × 1 dataset × 1 preset × 1 algorithm
        assert len(combinations) == 1

        combo = combinations[0]
        assert combo["task_id"] == "test_task"
        assert combo["dataset_id"] == "test_dataset"
        assert combo["preset_id"] == "baseline_preset"
        assert combo["algorithm_id"] == "baseline"
        assert combo["mode"] == "experiment"
        assert combo["status"] == "queued"

    def test_pipeline_runner_handles_reset_combinations(
        self, minimal_config, temp_work_dir
    ):
        """Testa que PipelineRunner adiciona novas combinações após reset."""
        from src.infrastructure.orchestration.pipeline_runner import PipelineRunner
        from src.infrastructure.persistence.work_state.wrappers.work_scoped import (
            WorkScopedPersistence,
        )
        from src.domain.work import WorkItem

        work_id = "test_work_reset"
        work_item = WorkItem(
            id=work_id,
            config=minimal_config,
            status=BaseStatus.QUEUED,
            created_at=time.time(),
            updated_at=time.time(),
            output_path=str(temp_work_dir / "results"),
        )

        base_store = WorkStatePersistence(temp_work_dir / "state.db")
        base_store.submit_work(work_item)
        work_store = WorkScopedPersistence(base_store, work_id)
        runner = PipelineRunner(work_store=work_store)

        # Primeira geração
        runner._generate_pipeline_combinations(minimal_config)
        combinations_first = work_store.get_combinations()
        assert len(combinations_first) == 1

        # Simular execução parcial: marcar como running
        combo = combinations_first[0]
        work_store.update_combination_status(
            combo["task_id"],
            combo["dataset_id"],
            combo["preset_id"],
            combo["algorithm_id"],
            "running",
        )

        # Criar configuração estendida com novo algoritmo
        extended_config = CSPBenchConfig(
            metadata=minimal_config.metadata,
            datasets=minimal_config.datasets,
            algorithms={
                "baseline_preset": AlgorithmsPresetConfig(
                    id="baseline_preset",
                    name="Baseline Algorithms",
                    description="Preset estendido",
                    items=[
                        AlgParamsConfig(name="baseline", params={}),
                        AlgParamsConfig(name="csc", params={}),  # novo algoritmo
                    ],
                )
            },
            tasks=minimal_config.tasks,
            system=minimal_config.system,
        )

        # Segunda geração com configuração estendida
        runner._generate_pipeline_combinations(extended_config)
        combinations_after = work_store.get_combinations()

        # Devemos ter 2 combinações agora: original + nova
        assert len(combinations_after) == 2

        # Verificar que a original foi reiniciada para 'queued'
        queued_combos = [c for c in combinations_after if c["status"] == "queued"]
        assert len(queued_combos) == 2  # ambas devem estar em 'queued' após reset

    @patch("src.infrastructure.orchestration.pipeline_runner.PipelineRunner.run")
    def test_execution_manager_workflow(self, mock_runner_run, minimal_config):
        """Testa fluxo completo do ExecutionManager (com mock do pipeline)."""
        from src.application.services.work_service import get_work_service

        # Mock do work service
        mock_work_service = Mock()
        mock_work_service.submit.return_value = "test_work_123"
        mock_work_service.get.return_value = {
            "id": "test_work_123",
            "config_json": minimal_config.to_dict(),
            "status": "queued",
            "created_at": time.time(),
            "updated_at": time.time(),
            "output_path": "/tmp/test_work",
            "error": None,
            "extra_json": "{}",
        }
        mock_work_service.mark_running.return_value = True
        mock_work_service.mark_finished.return_value = True

        # Criar ExecutionManager com mock
        execution_manager = ExecutionManager(work_service=mock_work_service)

        # Executar
        work_id = execution_manager.execute(minimal_config)

        # Verificar chamadas
        assert work_id == "test_work_123"
        mock_work_service.submit.assert_called_once()

        # Aguardar thread (simplificado - na prática seria mais robusto)
        time.sleep(0.1)

        # Verificar que o pipeline seria executado
        # (mock_runner_run será chamado pela thread em background)

    def test_status_normalization_in_pipeline(self, minimal_config, temp_work_dir):
        """Testa que normalização de status funciona no contexto do pipeline."""
        from src.domain.status import normalize_status
        from src.infrastructure.persistence.work_state.utils.validation import (
            validate_status,
        )

        # Testar vários formatos de status
        test_statuses = [
            (BaseStatus.RUNNING, "running"),
            ("COMPLETED", "completed"),
            ("  Failed  ", "failed"),
            ("CanceleD", "canceled"),
        ]

        for input_status, expected in test_statuses:
            # Normalização deve funcionar
            normalized = normalize_status(input_status)
            assert normalized == expected

            # Validação não deve levantar exceção
            validate_status(input_status)

    def test_combinations_submission_integration(self, temp_work_dir, minimal_config):
        """Testa integração da submissão de combinações no contexto real."""
        from src.infrastructure.persistence.work_state.wrappers.work_scoped import (
            WorkScopedPersistence,
        )
        from src.domain.work import WorkItem

        work_id = "integration_test"
        work_item = WorkItem(
            id=work_id,
            config=minimal_config,
            status=BaseStatus.QUEUED,
            created_at=time.time(),
            updated_at=time.time(),
            output_path=str(temp_work_dir),
        )

        base_store = WorkStatePersistence(temp_work_dir / "state.db")
        base_store.submit_work(work_item)
        work_store = WorkScopedPersistence(base_store, work_id)

        # Submeter múltiplas combinações
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
                "dataset_id": "dataset1",
                "preset_id": "preset1",
                "algorithm_id": "algo2",
                "mode": "experiment",
                "total_sequences": 5,
            },
        ]

        # Submeter e verificar
        count = work_store.submit_combinations(combinations)
        assert count == 2

        retrieved_combinations = work_store.get_combinations()
        assert len(retrieved_combinations) == 2

        # Verificar que todas têm status 'queued'
        for combo in retrieved_combinations:
            assert combo["status"] == "queued"
            assert combo["work_id"] == work_id
