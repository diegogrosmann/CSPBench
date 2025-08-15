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
import os


class PipelineService:
    """Facade simplificada para execução do pipeline.

    Agora não bloqueante: run(config) agenda execução em uma thread e retorna imediatamente o work_id.
    Status inicial: queued -> running (quando thread inicia) -> finished|error.
    """

    @staticmethod
    def _execute(
        work_id: str,
        config: CSPBenchConfig,
        monitor: Monitor,
    ) -> None:
        wm = get_work_manager()
        item_dict = wm.get(work_id)
        if not item_dict:
            raise Exception(f"Work: {work_id} não está registrado")
        try:
            wm.mark_running(work_id)
            work_dir = Path(item_dict["output_path"])

            # Inicializa store SQLite por work
            store = WorkStateStore(work_dir / "state.db")
            register_work_store(work_id, store)

            # Salva config e estado inicial
            if item_dict:
                store.upsert_work(item_dict)

            store.save_config(config)

            runner = PipelineRunner(monitor=monitor, work_id=work_id, store=store)
            try:
                runner.run(config)
                # Considera o caminho de saída como o próprio state.db desta execução
                wm.mark_finished(work_id)
                item_final = wm.get(work_id)
                if item_final:
                    store.upsert_work(item_final)
            except Exception as e:  # noqa: BLE001
                wm.mark_error(work_id, str(e))
                item_err = wm.get(work_id)
                if item_err:
                    store.upsert_work(item_err)
                raise
            finally:
                try:
                    store.close()
                except Exception:  # noqa: BLE001
                    pass
        except Exception:  # já marcado error acima; fallback silencioso
            pass

    @classmethod
    def run(
        cls,
        config: CSPBenchConfig,
        monitor: Monitor = NoOpMonitor(),
    ) -> str:

        mon = monitor
        results_path = Path(os.environ.get("OUTPUT_BASE_DIRECTORY", "output"))
        results_path.mkdir(parents=True, exist_ok=True)

        wm = get_work_manager()

        work_id, item = wm.submit(config=config, extra=None)

        # Cria diretório da execução e registra estado 'queued' no state.db da execução
        work_dir = item["output_path"]
        Path(work_dir).mkdir(parents=True, exist_ok=True)

        wm.attach_directory(work_id, work_dir)
        try:
            store = WorkStateStore(Path(work_dir) / "state.db")
            store.init_schema()

            if item:
                store.upsert_work(item)
            store.save_config(config)
            # Não registra globalmente aqui para evitar manter conexão aberta antes da thread
        finally:
            try:
                store.close()
            except Exception:  # noqa: BLE001
                pass

        thread = threading.Thread(
            target=cls._execute,
            args=(work_id, config, mon),
            name=f"pipeline-{work_id}",
            daemon=True,
        )
        thread.start()
        return work_id
