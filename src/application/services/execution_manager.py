"""
ExecutionManager - Unified Pipeline Execution Service
Centralizes execution orchestration for both CLI and Web interfaces.
"""

import logging
import threading
from typing import Any, Dict, Optional
import time

from infrastructure.orchestration.pipeline_runner import PipelineRunner
from infrastructure.persistence.work_state_persistence import WorkStatePersistence
from src.domain.config import CSPBenchConfig
from src.domain.work import WorkItem, WorkStatus

logger = logging.getLogger(__name__)

class ExecutionManager:
    """
    Unified execution manager for both CLI and Web interfaces.
    Handles work submission, monitoring, and status management through WorkService.
    """
    
    def __init__(self, work_service=None):
        """
        Initialize ExecutionManager with dependencies.
        
        Args:
            work_service: Work service instance for dependency injection
        """
        self._work_service = work_service
    
    @property
    def work_service(self):
        """Lazy load work service if not injected."""
        if self._work_service is None:
            from application.services.work_service import get_work_service
            self._work_service = get_work_service()
        return self._work_service
    
    def execute(
        self,
        config: CSPBenchConfig,
        extra: Optional[Dict[str, Any]] = None
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
            # Submit work to persistent storage
            work_id = self.work_service.submit(config=config, extra=extra)
            
            logger.info(f"Work submitted: {work_id}")
            
            thread = threading.Thread(
                target=self._execute_work,
                args=(work_id, config),
                daemon=True
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
            # Get work item to extract config
            work_item = self.work_service.get(work_id)
            if not work_item:
                logger.error(f"Work item {work_id} not found")
                return False
            
            # Restart the work item (resets status to QUEUED)
            if not self.work_service.restart(work_id):
                logger.error(f"Failed to restart work item {work_id}")
                return False
            
            # Extract config from work item
            config = work_item.config
            if not config:
                logger.error(f"No config found for work item {work_id}")
                return False
            
            logger.info(f"Restarting work: {work_id}")
            
            thread = threading.Thread(
                target=self._execute_work,
                args=(work_id, config),
                daemon=True
            )
            thread.start()
            
            return True
            
        except Exception as e:
            logger.error(f"Failed to restart work {work_id}: {e}")
            return False
        
    def _execute_work(self, work_id: str, config: CSPBenchConfig) -> Dict[str, Any]:
        """
        Execute work item with monitoring and status updates using standardized WorkStatus.
        """
        try:
            # Mark as running
            self.work_service.mark_running(work_id)
            
            # Get work item for pipeline execution
            work_item = self.work_service.get(work_id)
            if not work_item:
                raise ValueError(f"Work item {work_id} not found")
            
            # Execute pipeline directly
            self._run_pipeline(work_item, config)
            
            # Mark as completed using standardized status
            self.work_service.mark_finished(work_id)
            
            # Return standardized status
            return {
                "status": WorkStatus.COMPLETED.value,
                "output_path": str(work_item.output_path)
            }
            
        except Exception as e:
            # Mark as failed using standardized status
            self.work_service.mark_error(work_id, str(e))
            logger.error(f"Work {work_id} failed: {e}")
            return {
                "status": WorkStatus.FAILED.value,
                "error": str(e)
            }

    def _run_pipeline(
        self,
        work: WorkItem,
        config: CSPBenchConfig
    ) -> None:
        logger.info(f"Iniciando execução do pipeline para work_id: {work.id}")
        logger.debug(f"Configuração: {config.metadata.name}")

        try:
            work_dir = work.output_path
            
            # Inicializa store SQLite por work
            store = WorkStatePersistence(work_dir / "state.db")
          
            store.upsert_work(work_dict=work.to_dict())
            logger.debug("Item de trabalho salvo no store")

            store.save_config(config)
            logger.debug("Configuração salva no store")

            logger.info("Criando PipelineRunner")
            runner = PipelineRunner(work_id=work.id, store=store)
            try:
                logger.info("Iniciando execução do runner")
                runner.run(config)
                logger.info("Pipeline executado com sucesso")
            except Exception as e:  # noqa: BLE001
                logger.error(f"Erro durante execução do pipeline: {e}", exc_info=True)
                
                # Criar dicionário do work com erro para persistir
                work_with_error = work.to_dict()
                work_with_error["status"] = "FAILED"
                work_with_error["error"] = str(e)
                work_with_error["updated_at"] = time.time()
                
                store.upsert_work(work_with_error)
                logger.debug("Item de erro salvo no store")
                raise
            finally:
                try:
                    store.close()
                    logger.debug("WorkStateStore fechada")
                except Exception:  # noqa: BLE001
                    logger.warning("Erro ao fechar WorkStateStore")
                    pass
        except Exception as e:  # já marcado error acima; fallback silencioso
            logger.error("Falha geral na execução do pipeline")
            pass
