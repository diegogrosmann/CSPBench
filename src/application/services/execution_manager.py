"""
DEPRECATED: ExecutionManager

Esta classe foi movida para WorkManager como parte da refatoração de simplificação
da arquitetura. O ExecutionManager foi incorporado ao WorkManager para unificar
a gestão de trabalhos e execução de pipelines.

Migração:
- Substitua ExecutionManager() por get_work_service()
- O método execute() está disponível no WorkManager
- O método restart() foi renomeado para restart_execution()

Exemplo de migração:
    # Antes:
    from src.application.services.execution_manager import ExecutionManager
    manager = ExecutionManager()
    work_id = manager.execute(config)
    
    # Depois:
    from src.application.services.work_service import get_work_service
    manager = get_work_service()
    work_id = manager.execute(config)

Este arquivo será removido em versões futuras.
"""

import warnings
from typing import Any, Dict, Optional

from src.application.services.work_service import get_work_service
from src.domain.config import CSPBenchConfig


class ExecutionManager:
    """
    DEPRECATED: Use WorkManager via get_work_service() instead.
    
    This class is a compatibility wrapper that delegates to WorkManager.
    It will be removed in future versions.
    """
    
    def __init__(self, work_service=None):
        warnings.warn(
            "ExecutionManager is deprecated. Use get_work_service() instead.",
            DeprecationWarning,
            stacklevel=2
        )
        if work_service is None:
            self._work_service = get_work_service()
        else:
            self._work_service = work_service
    
    def execute(self, config: CSPBenchConfig, extra: Optional[Dict[str, Any]] = None) -> str:
        """Delegate to WorkManager.execute()"""
        return self._work_service.execute(config, extra)
    
    def restart(self, work_id: str) -> bool:
        """Delegate to WorkManager.restart_execution()"""
        return self._work_service.restart_execution(work_id)