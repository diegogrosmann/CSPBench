"""
Work Monitor Manager for WebSocket real-time monitoring.
"""

import asyncio
import hashlib
import json
import logging
import time
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

from src.infrastructure.persistence.work_state.core import WorkPersistence

from .schemas import (
    EventType,
    ExecutionChanges,
    MessageType,
    ProgressMessage,
    ProgressSnapshot,
    ProgressUpdate,
    serialize_error_summary,
    serialize_execution_detail,
    serialize_progress_summary,
)

logger = logging.getLogger(__name__)


class WorkMonitorSession:
    """Individual monitoring session for a work."""

    def __init__(self, work_id: str, update_interval: float = 0.5):
        self.work_id = work_id
        self.update_interval = update_interval
        self.running = False

        # State tracking
        self.last_progress_hash: Optional[str] = None
        self.last_executions_state: Dict[str, Dict] = {}
        # Cache de execuções por combination_id (para servir full list sob demanda)
        self._executions_cache: Dict[int, List] = {}
        self._executions_state_cache: Dict[int, Dict[str, Dict]] = {}
        # Combinações já "flushadas" (estado final enviado) para evitar envios duplicados
        self._flushed_combinations: Set[int] = set()
        self.last_combination_id: Optional[int] = None
        self.last_work_status: Optional[str] = None
        self.last_error_timestamp: float = 0
        self.last_warning_timestamp: float = 0
        self.last_event_timestamp: float = 0

        # Subscribers
        self.subscribers: Set[asyncio.Queue] = set()

        # Internal task
        self._monitor_task: Optional[asyncio.Task] = None
        self._lock = asyncio.Lock()
        # Flag para indicar se snapshot inicial já foi enviado
        self._initial_snapshot_sent: bool = False

        # Create persistence instance following console pattern
        self._persistence = WorkPersistence()

    async def add_subscriber(self, queue: asyncio.Queue) -> None:
        """Add a subscriber queue."""
        async with self._lock:
            self.subscribers.add(queue)
            logger.debug(
                f"Added subscriber to work {self.work_id}. Total: {len(self.subscribers)}"
            )

            # Se o monitor já está rodando e snapshot inicial já foi enviado, enviamos
            # um snapshot atualizado somente para o novo subscriber (late joiner).
            if self.running and self._initial_snapshot_sent:
                try:
                    snapshot_message = await self._build_snapshot_message()
                    if snapshot_message:
                        ws_message = snapshot_message.to_websocket_message().to_json()
                        await queue.put(ws_message)
                        logger.debug(
                            f"Sent on-demand snapshot to late subscriber of {self.work_id}"
                        )
                except Exception as e:
                    logger.warning(
                        f"Failed to send on-demand snapshot to new subscriber for {self.work_id}: {e}"
                    )

    async def remove_subscriber(self, queue: asyncio.Queue) -> None:
        """Remove a subscriber queue."""
        async with self._lock:
            self.subscribers.discard(queue)
            logger.debug(
                f"Removed subscriber from work {self.work_id}. Total: {len(self.subscribers)}"
            )

            # Stop monitoring if no subscribers
            if not self.subscribers and self.running:
                logger.info(
                    f"No more subscribers for work {self.work_id}, stopping monitor"
                )
                await self.stop()

    async def start(self) -> None:
        """Start monitoring loop."""
        if self.running:
            return

        self.running = True
        self._monitor_task = asyncio.create_task(self._monitor_loop())
        logger.info(f"[WEBSOCKET] Iniciado monitoramento para work {self.work_id}")

    async def stop(self) -> None:
        """Stop monitoring loop."""
        if not self.running:
            return

        logger.info(f"[WEBSOCKET] Parando monitoramento para work {self.work_id}")
        self.running = False
        if self._monitor_task:
            self._monitor_task.cancel()
            try:
                await self._monitor_task
            except asyncio.CancelledError:
                pass

    async def _monitor_loop(self) -> None:
        """Main monitoring loop."""
        loop_count = 0
        try:
            # Send initial snapshot
            await self._send_initial_snapshot()

            while self.running:
                try:
                    await self._check_for_updates()
                    await asyncio.sleep(self.update_interval)
                except Exception as e:
                    logger.error(f"[WEBSOCKET] Erro no loop de monitoramento para work {self.work_id}: {e}")
                    await asyncio.sleep(self.update_interval)

        except asyncio.CancelledError:
            pass
        except Exception as e:
            logger.error(f"[WEBSOCKET] Loop de monitoramento falhou para work {self.work_id}: {e}")

    async def _send_initial_snapshot(self) -> None:
        """Send initial complete snapshot."""
        try:
            message = await self._build_snapshot_message(update_tracking=True)
            if message:
                await self._broadcast_message(message)
                self._initial_snapshot_sent = True

        except Exception as e:
            logger.error(
                f"Failed to send initial snapshot for work {self.work_id}: {e}"
            )
            error_msg = ProgressMessage(
                type=MessageType.ERROR,
                work_id=self.work_id,
                timestamp=time.time(),
                error={"code": "SNAPSHOT_FAILED", "message": str(e)},
            )
            await self._broadcast_message(error_msg)

    async def _build_snapshot_message(
        self, update_tracking: bool = False
    ) -> Optional[ProgressMessage]:
        """Constroi mensagem de snapshot atual. Se update_tracking, atualiza estado interno."""
        try:
            progress_summary = self._persistence.get_work_progress_summary(self.work_id)
            if not progress_summary:
                logger.debug(
                    f"Snapshot build: work {self.work_id} progress_summary not found"
                )
                return ProgressMessage(
                    type=MessageType.ERROR,
                    work_id=self.work_id,
                    timestamp=time.time(),
                    error={
                        "code": "WORK_NOT_FOUND",
                        "message": f"Work {self.work_id} not found",
                    },
                )

            # Buscar TODAS as execuções de TODAS as combinações
            all_executions = []
            try:
                combinations_list = self._persistence.list_combinations(self.work_id)
                for combo in combinations_list:
                    combo_id = combo.get("combination_id")
                    if combo_id:
                        combo_executions = self._persistence.get_combination_executions_detail(
                            combo_id
                        )
                        all_executions.extend(combo_executions)
                        # Atualiza cache para cada combinação
                        self._executions_cache[combo_id] = combo_executions
                        self._executions_state_cache[combo_id] = {
                            e.unit_id: {"status": e.status, "progress": e.progress}
                            for e in combo_executions
                        }
            except Exception as e:
                logger.warning(f"Failed to load all executions: {e}")
                # Fallback: apenas combinação atual
                if progress_summary.current_combination_details:
                    combination_id = progress_summary.current_combination_details.get(
                        "combination_id"
                    )
                    if combination_id:
                        all_executions = self._persistence.get_combination_executions_detail(
                            combination_id
                        )
                        self._executions_cache[combination_id] = all_executions
                        self._executions_state_cache[combination_id] = {
                            e.unit_id: {"status": e.status, "progress": e.progress}
                            for e in all_executions
                        }

            # Lista de combinações para permitir seleção arbitrária
            try:
                combinations_list = self._persistence.list_combinations(self.work_id)
            except Exception:
                combinations_list = []

            errors = self._persistence.get_error_summary(self.work_id, limit=10)
            warnings = self._persistence.get_execution_warnings(self.work_id, limit=10)

            # Buscar todos os eventos (se método existir)
            try:
                events = self._persistence.event_list_by_work(self.work_id)[-50:]  # últimos 50
            except AttributeError:
                events = []

            snapshot = ProgressSnapshot(
                progress=serialize_progress_summary(progress_summary),
                executions=[serialize_execution_detail(e) for e in all_executions],
                logs={
                    "errors": [serialize_error_summary(e) for e in errors],
                    "warnings": warnings,
                },
                events=events,
                executions_full=[serialize_execution_detail(e) for e in all_executions],
                combinations=combinations_list,
            )

            if update_tracking:
                self._update_tracking_state(progress_summary, all_executions)
                logger.debug(
                    f"Snapshot build: tracking updated for work {self.work_id}"
                )

            return ProgressMessage(
                type=MessageType.SNAPSHOT,
                work_id=self.work_id,
                timestamp=time.time(),
                snapshot=snapshot,
            )
        except Exception as e:
            logger.error(f"Error building snapshot for work {self.work_id}: {e}")
            return ProgressMessage(
                type=MessageType.ERROR,
                work_id=self.work_id,
                timestamp=time.time(),
                error={
                    "code": "SNAPSHOT_ERROR",
                    "message": str(e),
                },
            )

    async def _check_for_updates(self) -> None:
        """Check for updates and send incremental changes."""
        try:
            # Verificar status do work INDEPENDENTE do progress_summary
            work_data = self._persistence.work_get(self.work_id)
            current_work_status = work_data.get("status") if work_data else None
            
            final_statuses = ["completed", "failed", "cancelled", "error"]
            is_final_status = current_work_status in final_statuses
            
            if is_final_status:
                logger.info(f"[WEBSOCKET] Work {self.work_id} finalizou com status {current_work_status}")
                asyncio.create_task(self._handle_work_completion())

            # Obter progress_summary para atualizações de progresso
            progress_summary = self._persistence.get_work_progress_summary(self.work_id)

            # Check for work status change
            if self.last_work_status != current_work_status:
                await self._send_work_status_event(
                    self.last_work_status, current_work_status
                )
                self.last_work_status = current_work_status

            # Se não há progress_summary, interromper processamento de atualizações
            if not progress_summary:
                return

            # Check for progress changes
            progress_hash = self._hash_progress(progress_summary)
            progress_changed = progress_hash != self.last_progress_hash

            # Execuções atuais
            executions = []
            combination_changed = False
            current_combination_id = None

            previous_combination_id = self.last_combination_id
            if progress_summary.current_combination_details:
                current_combination_id = (
                    progress_summary.current_combination_details.get(
                        "combination_id"
                    )
                )

                # Caso de mudança de combinação (existe nova e difere da anterior)
                if (
                    current_combination_id
                    and previous_combination_id
                    and current_combination_id != previous_combination_id
                ):
                    try:
                        # Flush final da combinação anterior antes de trocar
                        prev_execs = self._persistence.get_combination_executions_detail(
                            previous_combination_id
                        )
                        self._executions_cache[previous_combination_id] = prev_execs
                        self._executions_state_cache[previous_combination_id] = {
                            e.unit_id: {"status": e.status, "progress": e.progress}
                            for e in prev_execs
                        }
                        # Envia update final dessa combinação (força cliente a ver completions finais)
                        await self._send_update(
                            progress_summary=None,
                            execution_changes=None,
                            new_logs=None,
                            executions=[
                                serialize_execution_detail(e) for e in prev_execs
                            ],
                            executions_combination_id=previous_combination_id,
                        )
                        self._flushed_combinations.add(previous_combination_id)
                    except Exception as e:
                        logger.debug(
                            f"Flush previous combination {previous_combination_id} failed: {e}"
                        )

                # Caso nenhuma nova combinação mas havia anterior e ainda não flushado (fim do trabalho ou pausa)
                if (
                    not current_combination_id
                    and previous_combination_id is not None
                    and previous_combination_id not in self._flushed_combinations
                ):
                    try:
                        prev_execs = self._persistence.get_combination_executions_detail(
                            previous_combination_id
                        )
                        await self._send_update(
                            progress_summary=None,
                            execution_changes=None,
                            new_logs=None,
                            executions=[
                                serialize_execution_detail(e) for e in prev_execs
                            ],
                            executions_combination_id=previous_combination_id,
                        )
                        self._flushed_combinations.add(previous_combination_id)
                    except Exception as e:
                        logger.debug(
                            f"Final flush (no new combination) {previous_combination_id} failed: {e}"
                        )

                # Agora processa execuções da combinação corrente (se houver)
                executions = []
                if current_combination_id:
                    executions = self._persistence.get_combination_executions_detail(
                        current_combination_id
                    )
                    if previous_combination_id != current_combination_id:
                        combination_changed = True
                        await self._send_combination_changed_event(
                            previous_combination_id, current_combination_id
                        )
                        # Reset estado de execuções para nova combinação para detectar 'new'
                        self.last_executions_state = {}
                        self.last_combination_id = current_combination_id

                # Check for execution changes
                execution_changes = self._detect_execution_changes(executions)
                executions_changed = execution_changes is not None

                # Check for new logs
                new_logs = await self._check_new_logs(self._persistence)
                logs_changed = bool(new_logs["errors"] or new_logs["warnings"])

                # Check for new events
                new_events = await self._check_new_events(self._persistence)
                events_changed = bool(new_events)

                # Send updates if anything changed
                if (
                    progress_changed
                    or executions_changed
                    or logs_changed
                    or events_changed
                    or combination_changed
                ):
                    full_execs_serialized = None
                    if combination_changed or executions_changed:
                        # Buscar TODAS as execuções para atualização incremental também
                        try:
                            all_executions_for_update = []
                            combinations_list = self._persistence.list_combinations(self.work_id)
                            for combo in combinations_list:
                                combo_id = combo.get("combination_id")
                                if combo_id:
                                    combo_executions = (
                                        self._persistence.get_combination_executions_detail(
                                            combo_id
                                        )
                                    )
                                    all_executions_for_update.extend(combo_executions)
                            full_execs_serialized = [
                                serialize_execution_detail(e)
                                for e in all_executions_for_update
                            ]
                        except Exception as e:
                            logger.warning(
                                f"Failed to load all executions for update: {e}"
                            )
                            # Fallback: apenas combinação atual
                            full_execs_serialized = [
                                serialize_execution_detail(e) for e in executions
                            ]

                    await self._send_update(
                        progress_summary if progress_changed else None,
                        execution_changes if executions_changed else None,
                        new_logs if logs_changed else None,
                        new_events if events_changed else None,
                        executions=full_execs_serialized,
                        executions_combination_id=(
                            current_combination_id
                            if full_execs_serialized is not None
                            else None
                        ),
                    )

                    if progress_changed:
                        self.last_progress_hash = progress_hash
                    if executions_changed:
                        self._update_executions_state(executions)

        except Exception as e:
            logger.error(f"[WEBSOCKET] Exception em _check_for_updates para work {self.work_id}: {e}")

    async def _send_work_status_event(
        self, old_status: Optional[str], new_status: str
    ) -> None:
        """Send work status change event."""
        logger.info(f"[WEBSOCKET] 🔄 Enviando evento de mudança de status para work {self.work_id}: {old_status} -> {new_status}")
        
        event = ProgressMessage(
            type=MessageType.EVENT,
            work_id=self.work_id,
            timestamp=time.time(),
            event={
                "event_type": EventType.WORK_STATUS_CHANGED,
                "old_status": old_status,
                "new_status": new_status,
            },
        )
        await self._broadcast_message(event)

    async def _send_combination_changed_event(
        self, old_id: Optional[int], new_id: int
    ) -> None:
        """Send combination change event."""
        event = ProgressMessage(
            type=MessageType.EVENT,
            work_id=self.work_id,
            timestamp=time.time(),
            event={
                "event_type": EventType.COMBINATION_CHANGED,
                "old_combination_id": old_id,
                "new_combination_id": new_id,
            },
        )
        await self._broadcast_message(event)

    async def _send_update(
        self,
        progress_summary=None,
        execution_changes=None,
        new_logs=None,
        new_events=None,
        executions=None,
        executions_combination_id=None,
    ) -> None:
        """Send incremental update."""
        update_data = {}

        if progress_summary:
            update_data["progress"] = serialize_progress_summary(progress_summary)

        if execution_changes:
            update_data["executions_changed"] = asdict(execution_changes)

        if new_logs:
            update_data["logs_appended"] = new_logs

        if new_events:
            update_data["events_appended"] = new_events

        if executions is not None:
            update_data["executions"] = executions
            if executions_combination_id is not None:
                update_data["executions_combination_id"] = executions_combination_id
        if update_data:
            message = ProgressMessage(
                type=MessageType.UPDATE,
                work_id=self.work_id,
                timestamp=time.time(),
                update=ProgressUpdate(**update_data),
            )
            await self._broadcast_message(message)

    async def _check_new_logs(self, persistence: WorkPersistence) -> Dict[str, List]:
        """Check for new error and warning logs."""
        new_logs = {"errors": [], "warnings": []}

        try:
            # Get recent errors
            errors = persistence.get_error_summary(self.work_id, limit=5)
            for error in errors:
                if error.timestamp > self.last_error_timestamp:
                    new_logs["errors"].append(serialize_error_summary(error))
                    self.last_error_timestamp = max(
                        self.last_error_timestamp, error.timestamp
                    )

            # Get recent warnings
            warnings = persistence.get_execution_warnings(self.work_id, limit=5)
            for warning in warnings:
                warning_timestamp = warning.get("timestamp", 0)
                if warning_timestamp > self.last_warning_timestamp:
                    new_logs["warnings"].append(warning)
                    self.last_warning_timestamp = max(
                        self.last_warning_timestamp, warning_timestamp
                    )

        except Exception as e:
            logger.error(f"Error checking new logs for work {self.work_id}: {e}")

        return new_logs

    async def _check_new_events(
        self, persistence: WorkPersistence
    ) -> List[Dict[str, Any]]:
        """Check for new events from events table."""
        new_events = []

        try:
            events = persistence.get_events(
                self.work_id, limit=10
            )
            for event in events:
                event_timestamp = event.get("timestamp", 0)
                if event_timestamp > self.last_event_timestamp:
                    new_events.append(event)
                    self.last_event_timestamp = max(
                        self.last_event_timestamp, event_timestamp
                    )

        except Exception as e:
            logger.error(f"Error checking new events for work {self.work_id}: {e}")

        return new_events

    def _hash_progress(self, progress_summary) -> str:
        """Generate hash for progress summary."""
        # Create simplified dict for hashing (exclude timestamps)
        hash_data = {
            "global_progress": round(progress_summary.global_progress, 3),
            "global_execution": progress_summary.global_execution,
            "tasks": progress_summary.tasks,
            "datasets": progress_summary.datasets,
            "configs": progress_summary.configs,
            "algorithms": progress_summary.algorithms,
            "execution": progress_summary.execution,
        }
        return hashlib.md5(json.dumps(hash_data, sort_keys=True).encode()).hexdigest()

    def _detect_execution_changes(self, executions: List) -> Optional[ExecutionChanges]:
        """Detect changes in executions list."""
        current_state = {}
        for exec_detail in executions:
            current_state[exec_detail.unit_id] = {
                "status": exec_detail.status,
                # Mantém precisão integral para refletir pequenos avanços
                "progress": exec_detail.progress,
            }

        if not self.last_executions_state:
            # First time, consider all as new
            self.last_executions_state = current_state
            return (
                ExecutionChanges(new=[uid for uid in current_state.keys()])
                if current_state
                else None
            )

        changes = ExecutionChanges()

        # Check for updates and completions
        for unit_id, current in current_state.items():
            if unit_id in self.last_executions_state:
                last = self.last_executions_state[unit_id]
                if current != last:
                    if not changes.updated:
                        changes.updated = []
                    changes.updated.append(unit_id)

                    # Check if completed
                    if current["status"] in ["completed", "failed", "error"] and last[
                        "status"
                    ] not in ["completed", "failed", "error"]:
                        if not changes.completed:
                            changes.completed = []
                        changes.completed.append(unit_id)
            else:
                # New execution
                if not changes.new:
                    changes.new = []
                changes.new.append(unit_id)

        # Check for removed executions
        for unit_id in self.last_executions_state:
            if unit_id not in current_state:
                if not changes.removed:
                    changes.removed = []
                changes.removed.append(unit_id)

        # Return changes only if there are actual changes
        return (
            changes
            if any([changes.updated, changes.completed, changes.new, changes.removed])
            else None
        )

    def _update_tracking_state(self, progress_summary, executions):
        """Update internal tracking state."""
        self.last_progress_hash = self._hash_progress(progress_summary)
        self._update_executions_state(executions)

        # Inicializar last_work_status na primeira vez
        work_data = self._persistence.work_get(self.work_id)
        self.last_work_status = work_data.get("status") if work_data else None

        # Atualizar timestamp de eventos para o mais recente
        try:
            events = self._persistence.get_events(self.work_id, limit=1)
            if events:
                self.last_event_timestamp = max(
                    self.last_event_timestamp, events[0]["timestamp"]
                )
        except AttributeError:
            # Method may not exist in all versions
            pass

        if progress_summary.current_combination_details:
            self.last_combination_id = progress_summary.current_combination_details.get(
                "combination_id"
            )

    def _update_executions_state(self, executions):
        """Update executions tracking state."""
        self.last_executions_state = {}
        for exec_detail in executions:
            self.last_executions_state[exec_detail.unit_id] = {
                "status": exec_detail.status,
                "progress": exec_detail.progress,
            }
        # Mantém cache coerente para combination corrente
        if self.last_combination_id is not None:
            self._executions_cache[self.last_combination_id] = executions
            self._executions_state_cache[self.last_combination_id] = (
                self.last_executions_state.copy()
            )

    async def get_full_executions_for_combination(
        self, combination_id: int
    ) -> Optional[List[Dict]]:
        """Retorna lista completa serializada para uma combination específica."""
        from .schemas import serialize_execution_detail

        try:
            if combination_id not in self._executions_cache:
                # Carregar do banco
                try:
                    executions = self._persistence.get_combination_executions_detail(
                        combination_id
                    )
                    self._executions_cache[combination_id] = executions
                    self._executions_state_cache[combination_id] = {
                        e.unit_id: {"status": e.status, "progress": e.progress}
                        for e in executions
                    }
                except Exception as e:
                    logger.error(f"Error loading executions for combination {combination_id}: {e}")
                    return None
            return [
                serialize_execution_detail(e)
                for e in self._executions_cache.get(combination_id, [])
            ]
        except Exception as e:
            logger.warning(
                f"Failed to load executions for combination {combination_id}: {e}"
            )
            return None

    async def _handle_work_completion(self) -> None:
        """Handle work completion - send final snapshot and cleanup."""
        try:
            logger.info(f"[WEBSOCKET] Finalizando work {self.work_id} - enviando snapshot final")
            
            # Enviar snapshot final
            try:
                final_message = await self._build_snapshot_message(update_tracking=False)
                if final_message:
                    await self._broadcast_message(final_message)
            except Exception as snapshot_error:
                logger.error(f"[WEBSOCKET] Erro ao enviar snapshot final para work {self.work_id}: {snapshot_error}")
            
            # Período de graça antes de parar o monitoramento
            await asyncio.sleep(30)
            
            if self.running:
                logger.info(f"[WEBSOCKET] Parando monitoramento para work {self.work_id}")
                await self.stop()
            
        except asyncio.CancelledError:
            logger.debug(f"[WEBSOCKET] _handle_work_completion cancelado para work {self.work_id}")
            raise
        except Exception as e:
            logger.error(f"[WEBSOCKET] Erro em _handle_work_completion para work {self.work_id}: {e}")

    async def _broadcast_message(self, message: ProgressMessage) -> None:
        """Broadcast message to all subscribers."""
        if not self.subscribers:
            return

        # Log important events only
        if message.type in [MessageType.ERROR]:
            logger.info(f"[WEBSOCKET] Broadcasting {message.type.value} para {len(self.subscribers)} subscribers do work {self.work_id}")

        # Convert to WebSocket message
        ws_message = message.to_websocket_message()
        message_json = ws_message.to_json()

        # Send to all subscribers
        failed_queues = set()
        for queue in self.subscribers.copy():
            try:
                await queue.put(message_json)
            except Exception as e:
                logger.warning(f"Failed to send message to subscriber: {e}")
                failed_queues.add(queue)

        # Remove failed queues
        if failed_queues:
            async with self._lock:
                self.subscribers -= failed_queues


class WorkMonitorManager:
    """Manager for all work monitoring sessions."""

    def __init__(self):
        self.sessions: Dict[str, WorkMonitorSession] = {}
        self._lock = asyncio.Lock()

    async def get_or_create_session(
        self, work_id: str
    ) -> WorkMonitorSession:
        """Get existing session or create new one."""
        async with self._lock:
            if work_id not in self.sessions:
                session = WorkMonitorSession(work_id)
                self.sessions[work_id] = session
                logger.info(f"Created new monitoring session for work {work_id}")
            return self.sessions[work_id]

    async def remove_session(self, work_id: str) -> None:
        """Remove session when no longer needed."""
        async with self._lock:
            if work_id in self.sessions:
                session = self.sessions[work_id]
                await session.stop()
                del self.sessions[work_id]
                logger.info(f"Removed monitoring session for work {work_id}")

    async def cleanup(self) -> None:
        """Cleanup all sessions."""
        async with self._lock:
            for session in self.sessions.values():
                await session.stop()
            self.sessions.clear()
            logger.info("Cleaned up all monitoring sessions")


# Global manager instance
work_monitor_manager = WorkMonitorManager()
