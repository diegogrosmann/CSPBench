"""PersistenceMonitor

Monitor responsável por persistir eventos de execução (progresso, avisos, erros)
utilizando `ExecutionScopedPersistence`.

Camada: infraestrutura.
Depende apenas do protocolo de domínio `AlgorithmMonitor` e do wrapper de
persistência com escopo de execução.

Funcionalidades:
- on_progress: registra progresso (com throttling opcional)
- on_warning: registra warning associado à unidade
- on_error: registra erro (ou warning se só mensagem)
- is_cancelled: permite cancelamento externo via ExecutionController ou flag interna

Throttling: evita gravação excessiva de progresso com base em delta mínimo de
percentual e/ou intervalo mínimo de tempo. Sempre registra 0% inicial (se
invocado) e 100% final.
"""

from __future__ import annotations

from typing import Any, Optional, TYPE_CHECKING
import time
from src.domain.monitoring import AlgorithmMonitor
from src.infrastructure.persistence.work_state.wrappers.execution_scoped import (
    ExecutionScopedPersistence,
)

if TYPE_CHECKING:
    from src.infrastructure.execution_control import ExecutionController


class PersistenceMonitor(AlgorithmMonitor):
    """Monitor que persiste eventos no store via `ExecutionScopedPersistence`.

    Args:
        exec_store: Wrapper de persistência com escopo de execução.
        min_delta_pct: Delta mínimo de progresso para persistir (default 1.0).
        min_interval_s: Intervalo mínimo entre persistências (default 0.5s).
    """

    __slots__ = (
        "_exec_store",
        "_min_delta_pct",
        "_min_interval_s",
        "_last_progress",
        "_last_time",
        "_cancelled",
        "_execution_controller",
    )

    def __init__(
        self,
        exec_store: ExecutionScopedPersistence,
        *,
        min_delta_pct: float = 0.01,
        min_interval_s: float = 0.00,
        execution_controller: Optional["ExecutionController"] = None,
    ) -> None:
        self._exec_store = exec_store
        self._min_delta_pct = min_delta_pct
        self._min_interval_s = min_interval_s
        self._last_progress = -float("inf")
        self._last_time = 0.0
        self._cancelled = False
        self._execution_controller = execution_controller

    @classmethod
    def with_execution_controller(
        cls,
        exec_store: ExecutionScopedPersistence,
        execution_controller: "ExecutionController",
        *,
        min_delta_pct: float = 0.01,
        min_interval_s: float = 0.00,
    ) -> "PersistenceMonitor":
        """
        Convenience factory method to create PersistenceMonitor with ExecutionController.
        
        Args:
            exec_store: Execution-scoped persistence store
            execution_controller: Controller for status checks and resource management
            min_delta_pct: Minimum progress delta to persist
            min_interval_s: Minimum time interval between persistencies
            
        Returns:
            Configured PersistenceMonitor instance
        """
        return cls(
            exec_store=exec_store,
            min_delta_pct=min_delta_pct,
            min_interval_s=min_interval_s,
            execution_controller=execution_controller,
        )

    # ---------------------------------------------------------------------
    # API AlgorithmMonitor
    # ---------------------------------------------------------------------
    def on_progress(
        self, progress: float, message: str, /, **data: Any
    ) -> None:  # noqa: D401
        try:
            now = time.time()
            pct_ok = (progress - self._last_progress) >= self._min_delta_pct
            time_ok = (now - self._last_time) >= self._min_interval_s

            if (pct_ok and time_ok):
                self._exec_store.add_progress(progress, message)
                self._last_progress = progress
                self._last_time = now
        except Exception:  # pragma: no cover - não deve interromper o algoritmo
            pass

    def on_warning(self, message: str, /, **data: Any) -> None:  # noqa: D401
        try:
            context = data if data else None
            self._exec_store.unit_warning(message, context)
        except Exception:  # pragma: no cover
            pass

    def on_error(
        self, message: str, exc: Exception | None = None, /, **data: Any
    ) -> None:  # noqa: D401,E501
        try:
            if exc is not None:
                # Registra erro completo
                self._exec_store.unit_error(exc)
            else:
                # Sem exceção concreta -> tratar como warning categorizado
                context = (
                    {"as_error_message": message, **data}
                    if data
                    else {"as_error_message": message}
                )
                self._exec_store.unit_warning(message, context)
        except Exception:  # pragma: no cover
            pass

    def is_cancelled(self) -> bool:  # noqa: D401
        """Check if execution is cancelled via ExecutionController or internal flag."""
        # Check ExecutionController first if available
        if self._execution_controller:
            try:
                from src.domain.status import BaseStatus
                status = self._execution_controller.check_status()
                # Only CANCELED status means cancelled, not PAUSED
                if status == BaseStatus.CANCELED:
                    return True
            except Exception:  # pragma: no cover
                pass
        
        # Fallback to internal flag
        return self._cancelled

    # ------------------------------------------------------------------
    # Extensão
    # ------------------------------------------------------------------
    def cancel(self) -> None:
        """Sinaliza cancelamento externo."""
        self._cancelled = True

    # Conveniência para algoritmos que desejem acesso direto (evitar import extra)
    @property
    def execution(self) -> ExecutionScopedPersistence:  # noqa: D401
        return self._exec_store
    def cancel(self) -> None:
        """Sinaliza cancelamento externo."""
        self._cancelled = True