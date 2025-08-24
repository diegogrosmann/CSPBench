"""
ExecutionManager - Unified Pipeline Execution Service
Centralizes execution orchestration for both CLI and Web interfaces.
"""

import logging
from pathlib import Path
import threading
from typing import Any, Dict, Optional
import time

from src.infrastructure.orchestration.pipeline_runner import PipelineRunner
from src.infrastructure.persistence.work_state import (
    WorkStatePersistence,
    WorkScopedPersistence,
)
from src.domain.config import CSPBenchConfig
from src.domain.work import WorkItem, WorkStatus
from src.domain.status import BaseStatus

logger = logging.getLogger(__name__)


class ExecutionManager:
    """
    Unified execution manager for both CLI and Web interfaces.

    Handles work submission, monitoring, and status management through WorkService.

    Features:
    - Configuration validation before execution
    - Automatic configuration backup for audit purposes
    - Enhanced logging with configuration details
    - Configuration summary generation for monitoring
    - Standardized error handling and status management
    """

    def __init__(self, work_service=None):
        """
        Initialize ExecutionManager with dependencies.

        Args:
            work_service: Work service instance for dependency injection
        """

        if work_service is None:
            # Lazy load work service se não injetado
            from src.application.services.work_service import get_work_service

            self._work_service = get_work_service()
        else:
            self._work_service = work_service

    def execute(
        self, config: CSPBenchConfig, extra: Optional[Dict[str, Any]] = None
    ) -> str:
        """
        Execute pipeline configuration with unified workflow.

        Args:
            config: Pipeline configuration to execute
            extra: Additional metadata to store with work item

        Returns:
            work_id
        """
        if extra is None:
            extra = {}

        try:
            # Validate configuration before execution
            self._validate_config(config)

            # Submit work to persistent storage
            work_id = self._work_service.submit(config=config, extra=extra)

            logger.info(f"Work submitted: {work_id}")

            thread = threading.Thread(
                target=self._execute_work, args=(work_id, config), daemon=True
            )
            thread.start()

            return work_id

        except Exception as e:
            logger.error(f"Failed to execute work: {e}")
            raise

    def restart(self, work_id: str) -> bool:
        """
        Restart an existing work item.

        Args:
            work_id: ID of the work item to restart

        Returns:
            True if restart was successful, False otherwise
        """
        try:
            # Get work item (may be dict from repository wrapper)
            work_data = self._work_service.get(work_id)
            if not work_data:
                logger.error(f"Work item {work_id} not found")
                return False
            if isinstance(work_data, dict):
                work_item = WorkItem.from_dict(work_data)
            else:
                work_item = work_data
            config = work_item.config
            if not config:
                logger.error(f"No config found for work item {work_id}")
                return False

            # Validate configuration before restart
            try:
                self._validate_config(config)
            except ValueError as e:
                logger.error(f"Configuration validation failed for work {work_id}: {e}")
                return False

            # Restart the work item (resets status to QUEUED)
            if not self._work_service.restart(work_id):
                logger.error(f"Failed to restart work item {work_id}")
                return False

            logger.info(f"Restarting work: {work_id} (config: {config.metadata.name})")

            thread = threading.Thread(
                target=self._execute_work, args=(work_id, config), daemon=True
            )
            thread.start()

            return True

        except Exception as e:
            logger.error(f"Failed to restart work {work_id}: {e}")
            return False

    def _execute_work(self, work_id: str, config: CSPBenchConfig) -> Dict[str, Any]:
        """
        Execute work item with monitoring and status updates using standardized BaseStatus.
        """
        try:

            # Get work item for pipeline execution (convert dict if needed)
            work_data = self._work_service.get(work_id)
            if not work_data:
                raise ValueError(f"Work item {work_id} not found")
            if isinstance(work_data, dict):
                work_item = WorkItem.from_dict(work_data)
            else:
                work_item = work_data

            # Execute pipeline directly
            self._run_pipeline(work_item, config)

            # Mark as completed using standardized status
            self._work_service.mark_finished(work_id)

            # Return standardized status
            return {
                "status": BaseStatus.COMPLETED.value,
                "output_path": str(work_item.output_path),
            }

        except Exception as e:
            # Mark as failed using standardized status
            self._work_service.mark_error(work_id, str(e))
            logger.error(f"Work {work_id} failed: {e}")
            return {"status": BaseStatus.FAILED.value, "error": str(e)}

    def _run_pipeline(self, work: WorkItem, config: CSPBenchConfig) -> None:
        logger.info(f"Iniciando execução do pipeline para work_id: {work.id}")
        logger.debug(f"Configuração: {config.metadata.name} v{config.metadata.version}")

        # Log detalhes da configuração para auditoria
        logger.debug(f"Autor: {config.metadata.author}")
        logger.debug(f"Descrição: {config.metadata.description}")
        logger.debug(f"Data de criação: {config.metadata.creation_date}")
        if config.metadata.tags:
            logger.debug(f"Tags: {', '.join(config.metadata.tags)}")

        try:
            work_dir = work.output_path

            # Inicializa store SQLite por work
            base_store = WorkStatePersistence(Path(work_dir) / "state.db")
            # Inicia o store
            base_store.submit_work(work)
            # Criar store scoped para facilitar uso
            work_store = WorkScopedPersistence(base_store, work.id)

            logger.debug("Item de trabalho salvo no store")

            logger.info("Criando PipelineRunner")
            runner = PipelineRunner(work_store=work_store)
            try:
                logger.info("Iniciando execução do runner")
                runner.run(config)
                logger.info("Pipeline executado com sucesso")
            except Exception as e:  # noqa: BLE001
                logger.error(f"Erro durante execução do pipeline: {e}", exc_info=True)

                # Atualizar status do work com erro usando nova estrutura
                work_store.update_work_status(BaseStatus.FAILED.value, error=str(e))
                logger.debug("Item de erro salvo no store")
                raise
            finally:
                try:
                    # Fechar conexão do store adequadamente
                    if hasattr(base_store, "_conn"):
                        base_store._conn.close()
                    logger.debug("WorkStatePersistence fechada")
                except Exception:  # noqa: BLE001
                    logger.warning("Erro ao fechar WorkStatePersistence")
                    pass
        except Exception as e:  # noqa: BLE001
            logger.error(f"Falha geral na execução do work: {e}")
            self._work_service.mark_error(work.id, str(e))
            raise

    def _validate_config(self, config: CSPBenchConfig) -> None:
        """
        Validate configuration before execution.

        Args:
            config: Configuration to validate

        Raises:
            ValueError: If configuration is invalid
        """
        # Validate metadata
        if not config.metadata:
            raise ValueError("Metadata configuration is required")

        if not config.metadata.name:
            raise ValueError("Configuration name is required")

        # Validate that we have datasets
        if not config.datasets:
            raise ValueError("At least one dataset must be configured")

        # Validate that we have algorithms
        if not config.algorithms:
            raise ValueError("At least one algorithm preset must be configured")

        # Validate that we have tasks
        if not config.tasks or not config.tasks.items:
            raise ValueError("At least one task must be configured")

        # Validate that task references exist
        for task in config.tasks.items:
            # Check dataset references
            for dataset_id in task.datasets:
                if dataset_id not in config.datasets:
                    raise ValueError(
                        f"Task '{task.name}' references unknown dataset '{dataset_id}'"
                    )

            # Check algorithm references
            for algorithm_id in task.algorithms:
                if algorithm_id not in config.algorithms:
                    raise ValueError(
                        f"Task '{task.name}' references unknown algorithm preset '{algorithm_id}'"
                    )

        logger.debug("Configuration validation completed successfully")

    def get_config_summary(self, config: CSPBenchConfig) -> Dict[str, Any]:
        """
        Get a summary of the configuration for logging or display purposes.

        Args:
            config: Configuration to summarize

        Returns:
            Dictionary with configuration summary
        """
        summary = {
            "metadata": {
                "name": config.metadata.name,
                "version": config.metadata.version,
                "author": config.metadata.author,
                "description": config.metadata.description,
                "creation_date": config.metadata.creation_date,
                "tags": config.metadata.tags,
            },
            "datasets": {
                "count": len(config.datasets),
                "types": list(set(ds.type for ds in config.datasets.values())),
                "ids": list(config.datasets.keys()),
            },
            "algorithms": {
                "count": len(config.algorithms),
                "presets": [
                    {
                        "id": preset.id,
                        "name": preset.name,
                        "algorithm_count": len(preset.items),
                    }
                    for preset in config.algorithms.values()
                ],
            },
            "tasks": {
                "type": config.tasks.type,
                "count": len(config.tasks.items),
                "items": [
                    {
                        "id": task.id,
                        "name": task.name,
                        "type": task.type,
                        "datasets_count": len(task.datasets),
                        "algorithms_count": len(task.algorithms),
                    }
                    for task in config.tasks.items
                ],
            },
        }

        # Add optional sections if present
        if config.output:
            summary["output"] = {
                "logging": config.output.logging,
                "formats": {
                    "csv": config.output.results.formats.csv,
                    "json": config.output.results.formats.json,
                    "parquet": config.output.results.formats.parquet,
                    "pickle": config.output.results.formats.pickle,
                },
                "partial_results": config.output.results.partial_results,
            }

        if config.resources:
            summary["resources"] = {
                "cpu": {
                    "exclusive_cores": config.resources.cpu.exclusive_cores,
                    "max_workers": config.resources.cpu.max_workers,
                    "internal_jobs": config.resources.cpu.internal_jobs,
                },
                "memory": {"max_memory_gb": config.resources.memory.max_memory_gb},
                "timeouts": {
                    "timeout_per_item": config.resources.timeouts.timeout_per_item,
                    "timeout_total_batch": config.resources.timeouts.timeout_total_batch,
                },
            }

        if config.system:
            summary["system"] = {"global_seed": config.system.global_seed}

        return summary
