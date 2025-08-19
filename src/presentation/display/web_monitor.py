"""Monitor para uso com interface web.

Armazena de forma organizada as informações de progresso recebidas pelo sistema
e fornece métodos para obter um snapshot serializável (pronto para JSON) que
será usado pela API do frontend.

Esta implementação foca em manter um estado coerente e thread-safe. Não
possui dependências de saída (não escreve no terminal) — apenas armazena o
estado em memória.
"""
from __future__ import annotations

import copy
import logging
import threading
from collections import defaultdict
from datetime import datetime
from typing import Any, Dict, Optional


_LOGGER = logging.getLogger(__name__)


def _now_iso() -> str:
    return datetime.utcnow().isoformat() + "Z"


def _serialize_val(v: Any) -> Any:
    """Serializa valores simples para JSON (converte datetimes)."""
    if isinstance(v, datetime):
        return v.isoformat() + "Z"
    return v


class WebMonitor:
    """Monitor que guarda estado para API web.

    Contrato mínimo:
    - handle_event(event): aceita objetos de evento com atributos acessíveis via
      getattr(). A implementação tenta extrair campos conhecidos.
    - get_state(): retorna um dicionário pronto para JSON com snapshot do estado.
    """

    def __init__(self, max_events: int = 500) -> None:
        self._lock = threading.RLock()
        self._start_time: Optional[datetime] = None

        # Informações gerais da tarefa/pipeline
        self.task: Dict[str, Any] = {}

        # Execuções em progresso ou finalizadas
        # chave -> { name, status, start_time, end_time, total_items, completed_items, failed_items, current_item, progress, metadata, datasets, algorithms }
        self.executions: Dict[str, Dict[str, Any]] = {}
        self.current_execution: Optional[str] = None

        # Lista de eventos recentes (raw-ish) para inspeção/debug — limitado
        self.recent_events: list[Dict[str, Any]] = []
        self.max_events = max_events

    # ------------------------- Interface principal -------------------------
    def handle_event(self, event: Any) -> None:
        """Recebe um evento e atualiza o estado interno.

        O evento pode ser qualquer objeto; atributos são obtidos via getattr.
        """
        with self._lock:
            try:
                etype = type(event).__name__
                handler_map = {
                    "TaskStartedEvent": self._handle_task_started,
                    "TaskFinishedEvent": self._handle_task_finished,
                    "ExecutionStartedEvent": self._handle_execution_started,
                    "ExecutionProgressEvent": self._handle_execution_progress,
                    "ExecutionFinishedEvent": self._handle_execution_finished,
                    "ExperimentStartedEvent": self._handle_execution_started,
                    "ExperimentProgressEvent": self._handle_execution_progress,
                    "ExperimentFinishedEvent": self._handle_execution_finished,
                    "AlgorithmProgressEvent": self._handle_algorithm_progress,
                    "ErrorEvent": self._handle_error,
                    "WarningEvent": self._handle_warning,
                    "DisplayEvent": self._handle_display_event,
                }

                # Append a small raw event summary for debugging
                self._append_recent_event(event)

                handler = handler_map.get(etype, self._handle_generic)
                handler(event)
            except Exception as exc:  # pragma: no cover - defensive
                _LOGGER.exception("Erro ao processar evento: %s", exc)

    # ------------------------- Handlers de eventos -------------------------
    def _handle_task_started(self, event: Any) -> None:
        self._start_time = datetime.utcnow()
        name = getattr(event, "task_name", getattr(event, "name", "task"))
        ttype = getattr(event, "task_type", None)
        metadata = getattr(event, "metadata", {}) or {}

        self.task = {
            "name": name,
            "type": getattr(ttype, "value", str(ttype)) if ttype is not None else None,
            "metadata": dict(metadata),
            "status": "running",
            "start_time": self._start_time,
        }

    def _handle_task_finished(self, event: Any) -> None:
        success = getattr(event, "success", True)
        error_message = getattr(event, "error_message", None)
        results = getattr(event, "results", None) or {}

        if self.task:
            self.task["status"] = "completed" if success else "failed"
            self.task["end_time"] = datetime.utcnow()
            if error_message:
                self.task.setdefault("errors", []).append(str(error_message))
            if results:
                self.task.setdefault("results", {}).update(results)

    def _handle_execution_started(self, event: Any) -> None:
        name = getattr(event, "execution_name", getattr(event, "experiment_name", "execution"))
        metadata = getattr(event, "metadata", {}) or {}
        total_items = int(getattr(event, "total_items", 0) or 0)

        exec_state = {
            "name": name,
            "status": "running",
            "start_time": datetime.utcnow(),
            "end_time": None,
            "total_items": total_items,
            "completed_items": 0,
            "failed_items": 0,
            "current_item": None,
            "progress": 0.0,
            "metadata": dict(metadata),
            "datasets": {},
            "algorithms": {},
        }

        self.executions[name] = exec_state
        self.current_execution = name

    def _handle_execution_progress(self, event: Any) -> None:
        if not self.current_execution:
            return
        execution = self.executions.get(self.current_execution)
        if not execution:
            return

        progress = float(getattr(event, "progress_percent", getattr(event, "progress", 0.0) or 0.0))
        item_name = getattr(event, "item_name", None) or getattr(event, "item", None)
        context = getattr(event, "context", {}) or {}

        execution["progress"] = progress
        execution["current_item"] = item_name

        dataset_name = context.get("dataset_name") or item_name
        algorithm_name = context.get("algorithm_name")
        config_name = context.get("algorithm_config_name")

        if dataset_name:
            ds = execution["datasets"].setdefault(dataset_name, {"status": "running", "progress": 0.0})
            ds["status"] = "running"
            ds["progress"] = progress

        if algorithm_name:
            alg = execution["algorithms"].setdefault(algorithm_name, {})
            alg.update({"status": "running", "progress": progress, "dataset": dataset_name, "config": config_name})

    def _handle_execution_finished(self, event: Any) -> None:
        name = getattr(event, "execution_name", getattr(event, "experiment_name", self.current_execution))
        if not name or name not in self.executions:
            return

        execution = self.executions[name]
        success = getattr(event, "success", True)
        execution["status"] = "completed" if success else "failed"
        execution["end_time"] = datetime.utcnow()

    def _handle_algorithm_progress(self, event: Any) -> None:
        if not self.current_execution:
            return
        execution = self.executions.get(self.current_execution)
        if not execution:
            return

        name = getattr(event, "algorithm_name", "algorithm")
        progress = float(getattr(event, "progress_percent", getattr(event, "progress", 0.0) or 0.0))
        message = getattr(event, "message", None)

        alg = execution["algorithms"].setdefault(name, {})
        alg.update({"progress": progress, "message": message, "status": "completed" if progress >= 100 else "running"})

    def _handle_error(self, event: Any) -> None:
        # Guarda erro associado à execução atual ou globalmente
        error_type = getattr(event, "error_type", None) or type(event).__name__
        error_message = getattr(event, "error_message", getattr(event, "message", None))

        target = self.current_execution or "__global__"
        bucket = self.executions.setdefault(target, {"name": target, "errors": []})
        bucket.setdefault("errors", []).append({"type": str(error_type), "message": str(error_message), "time": _now_iso()})

    def _handle_warning(self, event: Any) -> None:
        msg = getattr(event, "message", None)
        target = self.current_execution or "__global__"
        bucket = self.executions.setdefault(target, {"name": target, "warnings": []})
        bucket.setdefault("warnings", []).append({"message": str(msg), "time": _now_iso()})

    def _handle_display_event(self, event: Any) -> None:
        # Eventos Display podem trazer informações ricas (phase, payload, trial,etc.)
        phase = getattr(event, "phase", None)
        dataset = getattr(event, "dataset_id", None) or getattr(event, "dataset_name", None)
        algorithm = getattr(event, "algorithm_name", None)
        message = getattr(event, "message", None)
        payload = getattr(event, "payload", {}) or {}

        target = self.current_execution or "__global__"
        bucket = self.executions.setdefault(target, {"name": target})
        bucket.setdefault("display_events", []).append({
            "phase": getattr(phase, "value", str(phase)) if phase is not None else None,
            "dataset": dataset,
            "algorithm": algorithm,
            "message": message,
            "payload": payload,
            "time": _now_iso(),
        })

    def _handle_generic(self, event: Any) -> None:
        # Guarda uma representação mínima
        target = self.current_execution or "__global__"
        bucket = self.executions.setdefault(target, {"name": target})
        bucket.setdefault("events", []).append({"type": type(event).__name__, "time": _now_iso()})

    # ------------------------- Utilitários -------------------------
    def _append_recent_event(self, event: Any) -> None:
        try:
            summary = {"type": type(event).__name__}
            # tenta extrair campos comuns sem falhar
            for key in ("item_name", "execution_name", "algorithm_name", "message"):
                val = getattr(event, key, None)
                if val is not None:
                    summary[key] = val
            summary["time"] = _now_iso()
            self.recent_events.append(summary)
            if len(self.recent_events) > self.max_events:
                self.recent_events.pop(0)
        except Exception:
            pass

    def get_state(self) -> Dict[str, Any]:
        """Retorna um snapshot pronto para serialização JSON.

        Converte datetimes para strings ISO e faz deep copy para evitar races.
        """
        with self._lock:
            state = {
                "task": copy.deepcopy(self.task),
                "executions": copy.deepcopy(self.executions),
                "current_execution": self.current_execution,
                "recent_events": list(self.recent_events),
                "start_time": _serialize_val(self._start_time) if self._start_time else None,
            }

        # Serialize datetimes inside nested structures
        def _walk(obj: Any):
            if isinstance(obj, dict):
                return {k: _walk(v) for k, v in obj.items()}
            if isinstance(obj, list):
                return [_walk(v) for v in obj]
            return _serialize_val(obj)

        return _walk(state)

    # Pequenos helpers para APIs específicas
    def get_execution(self, name: str) -> Optional[Dict[str, Any]]:
        with self._lock:
            e = self.executions.get(name)
            return copy.deepcopy(e) if e is not None else None

    def clear(self) -> None:
        with self._lock:
            self.task = {}
            self.executions = {}
            self.current_execution = None
            self.recent_events = []


__all__ = ["WebMonitor"]
import threading
from collections import defaultdict, deque
from datetime import datetime
from typing import Any, Dict, List, Optional


def _now_iso() -> str:
    return datetime.utcnow().isoformat() + "Z"


class WebMonitor:
    """
    Monitor para serviço web.

    - Armazena estado de tarefas, execuções, datasets e algoritmos em estruturas
      serializáveis (dicts/lists).
    - Seguro para acesso concorrente via RLock.
    - Fornece snapshots prontos para serem retornados por uma API.
    """

    def __init__(self, history_size: int = 200):
        self._lock = threading.RLock()
        self._start_time: Optional[str] = None

        # Estruturas principais
        self.tasks: Dict[str, Dict[str, Any]] = {}
        self.executions: Dict[str, Dict[str, Any]] = {}
        self.events: deque = deque(maxlen=history_size)

        # Estado auxiliar
        self.current_execution: Optional[str] = None
        self._last_update: Optional[str] = None

    # ---------- Helpers ----------
    def _push_event(self, evt: Dict[str, Any]) -> None:
        with self._lock:
            evt["_ts"] = _now_iso()
            self.events.append(evt)
            self._last_update = evt["_ts"]

    def _ensure_execution(self, name: str) -> Dict[str, Any]:
        if name not in self.executions:
            self.executions[name] = {
                "name": name,
                "status": "unknown",
                "start_time": None,
                "end_time": None,
                "progress": 0.0,
                "current_item": None,
                "metadata": {},
                "datasets": {},
                "algorithms": {},
                "completed_items": 0,
                "failed_items": 0,
            }
        return self.executions[name]

    # ---------- Public API methods ----------
    def get_state(self) -> Dict[str, Any]:
        """Retorna snapshot serializável do estado do monitor."""
        with self._lock:
            return {
                "start_time": self._start_time,
                "last_update": self._last_update,
                "current_execution": self.current_execution,
                "tasks": dict(self.tasks),
                "executions": dict(self.executions),
                "events": list(self.events),
            }

    def get_execution(self, name: str) -> Optional[Dict[str, Any]]:
        with self._lock:
            return dict(self.executions.get(name)) if name in self.executions else None

    def list_active_executions(self) -> List[str]:
        with self._lock:
            return [n for n, e in self.executions.items() if e.get("status") == "running"]

    def reset(self) -> None:
        with self._lock:
            self.__init__(history_size=self.events.maxlen)

    # ---------- Event handling entrypoint ----------
    def handle_event(self, event: Any) -> None:
        """
        Recebe um evento (objeto qualquer) e atualiza o estado interno.
        O monitor tenta extrair atributos comuns via getattr.
        """
        with self._lock:
            try:
                event_type = type(event).__name__
                handler_map = {
                    "TaskStartedEvent": self._handle_task_started,
                    "TaskFinishedEvent": self._handle_task_finished,
                    "ExecutionStartedEvent": self._handle_execution_started,
                    "ExecutionProgressEvent": self._handle_execution_progress,
                    "ExecutionFinishedEvent": self._handle_execution_finished,
                    "ExperimentStartedEvent": self._handle_execution_started,
                    "ExperimentProgressEvent": self._handle_execution_progress,
                    "ExperimentFinishedEvent": self._handle_execution_finished,
                    "AlgorithmProgressEvent": self._handle_algorithm_progress,
                    "ErrorEvent": self._handle_error,
                    "WarningEvent": self._handle_warning,
                    "DisplayEvent": self._handle_display_event,
                }

                handler = handler_map.get(event_type, self._handle_generic)
                handler(event)
            except Exception:
                # Não levantar exceções para não quebrar o serviço que publica eventos
                self._push_event({"type": "MonitorError", "message": "failed to process event"})

    # ---------- Handlers ----------
    def _handle_task_started(self, event: Any) -> None:
        task_id = getattr(event, "task_id", getattr(event, "task_name", None)) or "task"
        name = getattr(event, "task_name", task_id)
        task_type = getattr(event, "task_type", None)
        meta = getattr(event, "metadata", {}) or {}

        self._start_time = self._start_time or _now_iso()
        self.tasks[task_id] = {
            "id": task_id,
            "name": name,
            "type": getattr(task_type, "value", str(task_type)) if task_type else "task",
            "metadata": dict(meta),
            "status": "running",
            "start_time": _now_iso(),
            "end_time": None,
        }

        self._push_event({"type": "TaskStarted", "task_id": task_id, "name": name})

    def _handle_task_finished(self, event: Any) -> None:
        task_id = getattr(event, "task_id", getattr(event, "task_name", None)) or "task"
        success = getattr(event, "success", True)
        results = getattr(event, "results", {}) or {}
        error_message = getattr(event, "error_message", None)

        task = self.tasks.get(task_id)
        if task:
            task["status"] = "completed" if success else "failed"
            task["end_time"] = _now_iso()
            task["results"] = dict(results)
            if error_message:
                task["error"] = str(error_message)

        self._push_event({"type": "TaskFinished", "task_id": task_id, "success": bool(success)})

    def _handle_execution_started(self, event: Any) -> None:
        name = getattr(event, "execution_name", getattr(event, "experiment_name", "execution"))
        metadata = getattr(event, "metadata", {}) or {}
        total_items = int(getattr(event, "total_items", 0) or 0)

        exec_obj = self._ensure_execution(name)
        exec_obj.update(
            {
                "status": "running",
                "start_time": _now_iso(),
                "end_time": None,
                "progress": 0.0,
                "current_item": None,
                "metadata": dict(metadata),
                "total_items": total_items,
            }
        )

        self.current_execution = name
        self._push_event({"type": "ExecutionStarted", "execution": name})

    def _handle_execution_progress(self, event: Any) -> None:
        exec_name = getattr(event, "execution_name", getattr(event, "experiment_name", self.current_execution))
        if not exec_name:
            return

        exec_obj = self._ensure_execution(exec_name)
        progress = float(getattr(event, "progress_percent", getattr(event, "progress", exec_obj.get("progress", 0.0))) or 0.0)
        item_name = getattr(event, "item_name", None)
        context = getattr(event, "context", {}) or {}

        exec_obj["progress"] = progress
        exec_obj["current_item"] = item_name or exec_obj.get("current_item")
        exec_obj["last_update"] = _now_iso()

        # Atualiza dataset e algorithm se presentes no contexto
        dataset_name = context.get("dataset_name")
        algorithm_name = context.get("algorithm_name")
        alg_config = context.get("algorithm_config_name")

        if dataset_name:
            ds = exec_obj["datasets"].setdefault(dataset_name, {"name": dataset_name, "progress": 0.0, "status": "running"})
            ds["progress"] = progress

        if algorithm_name:
            ag = exec_obj["algorithms"].setdefault(algorithm_name, {"name": algorithm_name, "progress": 0.0, "status": "running", "config": alg_config})
            ag["progress"] = progress
            if alg_config:
                ag["config"] = alg_config

        self._push_event({"type": "ExecutionProgress", "execution": exec_name, "progress": progress})

    def _handle_execution_finished(self, event: Any) -> None:
        exec_name = getattr(event, "execution_name", getattr(event, "experiment_name", self.current_execution))
        success = getattr(event, "success", True)

        if not exec_name or exec_name not in self.executions:
            return

        exec_obj = self.executions[exec_name]
        exec_obj["status"] = "completed" if success else "failed"
        exec_obj["end_time"] = _now_iso()
        exec_obj["progress"] = 100.0 if success else exec_obj.get("progress", 0.0)

        # Estatísticas opcionais
        exec_obj["summary"] = getattr(event, "summary", {}) or {}

        self._push_event({"type": "ExecutionFinished", "execution": exec_name, "success": bool(success)})

    def _handle_algorithm_progress(self, event: Any) -> None:
        algorithm_name = getattr(event, "algorithm_name", "unknown")
        progress = float(getattr(event, "progress_percent", getattr(event, "progress", 0.0)) or 0.0)
        message = getattr(event, "message", "")

        if not self.current_execution:
            # tenta localizar execução no evento
            exec_name = getattr(event, "execution_name", None)
        else:
            exec_name = getattr(event, "execution_name", self.current_execution)

        if not exec_name:
            return

        exec_obj = self._ensure_execution(exec_name)
        alg = exec_obj["algorithms"].setdefault(algorithm_name, {"name": algorithm_name, "progress": 0.0, "status": "running", "message": ""})
        alg["progress"] = progress
        alg["message"] = str(message)
        alg["status"] = "completed" if progress >= 100.0 else "running"

        self._push_event({"type": "AlgorithmProgress", "execution": exec_name, "algorithm": algorithm_name, "progress": progress})

    def _handle_error(self, event: Any) -> None:
        error_type = getattr(event, "error_type", None) or type(event).__name__
        error_message = getattr(event, "error_message", getattr(event, "message", ""))
        ctx = getattr(event, "context", {}) or {}

        self._push_event({"type": "Error", "error_type": str(error_type), "message": str(error_message), "context": dict(ctx)})

    def _handle_warning(self, event: Any) -> None:
        message = getattr(event, "message", "warning")
        self._push_event({"type": "Warning", "message": str(message)})

    def _handle_display_event(self, event: Any) -> None:
        # Eventos de Display podem carregar fase/contexto; armazenar o payload
        phase = getattr(event, "phase", None)
        dataset = getattr(event, "dataset_id", None)
        algorithm = getattr(event, "algorithm_name", None)
        message = getattr(event, "message", None)
        payload = getattr(event, "payload", {}) or {}

        record = {
            "type": "Display",
            "phase": getattr(phase, "value", str(phase)) if phase is not None else None,
            "dataset": dataset,
            "algorithm": algorithm,
            "message": message,
            "payload": dict(payload),
        }

        # Para eventos de progresso embutido, atualizar algoritmo se aplicável
        try:
            self._push_event(record)
        except Exception:
            self._push_event({"type": "Display", "error": "failed to store display event"})

    def _handle_generic(self, event: Any) -> None:
        # Guarda representação simples do evento para inspeção posterior
        data = {"type": "Generic", "repr": repr(event)}
        self._push_event(data)
