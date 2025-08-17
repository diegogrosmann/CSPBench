"""High-level facade to execute a CSPBenchConfig through the new pipeline."""

from __future__ import annotations

from pathlib import Path
import threading

from src.application.work.manager import get_work_manager
from src.domain.config import CSPBenchConfig
from src.infrastructure.monitoring.monitor_interface import Monitor, NoOpMonitor
from src.infrastructure.orchestration.pipeline_runner import PipelineRunner
from src.infrastructure.persistence.work_state_store import (
    WorkStateStore,
    register_work_store,
)
from src.infrastructure.logging_config import get_logger
import os

# Create module logger
logger = get_logger("CSPBench.PipelineService")


class PipelineService:
    """High-level service for pipeline execution."""
    
    def __init__(self):
        logger.info("PipelineService inicializado")

    @staticmethod
    def _execute(
        work_id: str,
        config: CSPBenchConfig,
        monitor: Monitor,
    ) -> None:
        logger.info(f"Iniciando execução do pipeline para work_id: {work_id}")
        logger.debug(f"Configuração: {config.metadata.name}")
        
        wm = get_work_manager()
        item_dict = wm.get(work_id)
        if not item_dict:
            error_msg = f"Work: {work_id} não está registrado"
            logger.error(error_msg)
            raise Exception(error_msg)
            
        try:
            logger.info(f"Marcando work_id {work_id} como 'running'")
            wm.mark_running(work_id)
            work_dir = Path(item_dict["output_path"])
            logger.debug(f"Diretório de trabalho: {work_dir}")

            # Inicializa store SQLite por work
            store = WorkStateStore(work_dir / "state.db")
            register_work_store(work_id, store)
            logger.info("WorkStateStore inicializada e registrada")

            # Salva config e estado inicial
            if item_dict:
                store.upsert_work(item_dict)
                logger.debug("Item de trabalho salvo no store")

            store.save_config(config)
            logger.debug("Configuração salva no store")

            logger.info("Criando PipelineRunner")
            runner = PipelineRunner(monitor=monitor, work_id=work_id, store=store)
            try:
                logger.info("Iniciando execução do runner")
                runner.run(config)
                logger.info("Pipeline executado com sucesso")
                
                # Considera o caminho de saída como o próprio state.db desta execução
                wm.mark_finished(work_id)
                logger.info(f"Work_id {work_id} marcado como finalizado")
                
                item_final = wm.get(work_id)
                if item_final:
                    store.upsert_work(item_final)
                    logger.debug("Item final atualizado no store")
            except Exception as e:  # noqa: BLE001
                logger.error(f"Erro durante execução do pipeline: {e}", exc_info=True)
                wm.mark_error(work_id, str(e))
                item_err = wm.get(work_id)
                if item_err:
                    store.upsert_work(item_err)
                    logger.debug("Item de erro salvo no store")
                raise
            finally:
                try:
                    store.close()
                    logger.debug("WorkStateStore fechada")
                except Exception:  # noqa: BLE001
                    logger.warning("Erro ao fechar WorkStateStore")
                    pass
        except Exception:  # já marcado error acima; fallback silencioso
            logger.error("Falha geral na execução do pipeline")
            pass

    @classmethod
    def run(
        cls,
        config: CSPBenchConfig,
        monitor: Monitor | None = None,
    ) -> str:
        logger.info("PipelineService.run() iniciado")
        logger.info(f"Configuração: {config.metadata.name}")

        mon = monitor if monitor is not None else NoOpMonitor()
        results_path = Path(os.environ.get("OUTPUT_BASE_DIRECTORY", "output"))
        results_path.mkdir(parents=True, exist_ok=True)
        logger.debug(f"Diretório de resultados: {results_path}")

        wm = get_work_manager()

        logger.info("Submetendo trabalho ao WorkManager")
        work_id, item = wm.submit(config=config, extra=None)
        logger.info(f"Trabalho submetido com work_id: {work_id}")

        # Cria diretório da execução e registra estado 'queued' no state.db da execução
        work_dir = item["output_path"]
        Path(work_dir).mkdir(parents=True, exist_ok=True)
        logger.debug(f"Diretório de trabalho criado: {work_dir}")

        try:
            store = WorkStateStore(Path(work_dir) / "state.db")
            store.init_schema()
            logger.debug("Schema do banco de dados inicializado")

            if item:
                store.upsert_work(item)
                logger.debug("Item de trabalho inicial salvo")
            store.save_config(config)
            logger.debug("Configuração salva no banco")
            # Não registra globalmente aqui para evitar manter conexão aberta antes da thread
        finally:
            try:
                store.close()
            except Exception:  # noqa: BLE001
                pass

        logger.info("Iniciando thread de execução")
        thread = threading.Thread(
            target=cls._execute,
            args=(work_id, config, mon),
            name=f"pipeline-{work_id}",
            daemon=True,
        )
        thread.start()
        logger.info(f"Thread iniciada para work_id: {work_id}")
        return work_id
