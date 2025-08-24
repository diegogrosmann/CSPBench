from enum import Enum


class BaseStatus(str, Enum):
    """Standardized work execution statuses across the entire system."""

    QUEUED = "queued"  # Trabalho na fila, aguardando para começar
    RUNNING = "running"  # Trabalho em execução atualmente
    PAUSED = "paused"  # Trabalho pausado, pode ser retomado
    CANCELED = "canceled"  # Trabalho cancelado antes de terminar
    COMPLETED = "completed"  # Trabalho finalizado com sucesso
    FAILED = "failed"  # Falha impediu a execução de prosseguir
    ERROR = "error"  # Erro durante execução, mas prosseguiu até o fim

    @property
    def is_final(self) -> bool:
        """Indica se o status representa um estado finalizado (terminal)."""
        return self in FINAL_STATUSES

    @property
    def is_incomplete(self) -> bool:
        """Indica se o status representa um estado incompleto (não finalizado)."""
        return self in INCOMPLETE_STATUSES


# Status finalizados (terminais) - trabalho chegou ao fim
# Estes status indicam que o trabalho foi concluído de alguma forma
FINAL_STATUSES: set[BaseStatus] = {
    BaseStatus.COMPLETED,  # Sucesso - trabalho executado completamente
    BaseStatus.FAILED,  # Falha - execução impedida de prosseguir
    BaseStatus.ERROR,  # Erro - execução prosseguiu mas com problemas
}

# Status incompletos (não finalizados) - trabalho não chegou ao fim
# Estes status indicam que o trabalho ainda não foi concluído
INCOMPLETE_STATUSES: set[BaseStatus] = {
    BaseStatus.QUEUED,  # Aguardando - trabalho na fila para execução
    BaseStatus.RUNNING,  # Ativo - trabalho sendo executado no momento
    BaseStatus.PAUSED,  # Suspenso - trabalho pausado, pode ser retomado
    BaseStatus.CANCELED,  # Cancelado - trabalho interrompido antes do fim
}

ALLOWEDSTATUS: dict[BaseStatus, set[BaseStatus]] = {
    BaseStatus.QUEUED: {BaseStatus.RUNNING, BaseStatus.CANCELED},
    BaseStatus.RUNNING: {
        BaseStatus.PAUSED,
        BaseStatus.COMPLETED,
        BaseStatus.FAILED,
        BaseStatus.CANCELED,
        BaseStatus.ERROR,
    },
    BaseStatus.PAUSED: {
        BaseStatus.RUNNING,
        BaseStatus.CANCELED,
    },
}


def normalize_status(value) -> str:
    """
    Normaliza um valor de status para string lowercase.

    Aceita:
    - Instâncias de BaseStatus (retorna .value)
    - Strings (converte para lowercase e strip)
    - Outros tipos (converte para string e aplica lowercase)

    Args:
        value: Valor a ser normalizado (BaseStatus, str, ou outro)

    Returns:
        String normalizada em lowercase

    Raises:
        ValueError: Se o valor normalizado não corresponder a um status válido
    """
    if isinstance(value, BaseStatus):
        return value.value

    if isinstance(value, str):
        normalized = value.strip().lower()
    else:
        normalized = str(value).strip().lower()

    # Validar se é um status conhecido
    valid_values = {s.value for s in BaseStatus}
    if normalized not in valid_values:
        raise ValueError(f"Status inválido: {value}. Válidos: {sorted(valid_values)}")

    return normalized


def normalize_status(value: str | BaseStatus) -> str:
    """Normaliza um status para string lowercase validada.

    Args:
        value: Status como BaseStatus ou string (case-insensitive)
    Returns:
        str normalizada (um dos valores de BaseStatus)
    Raises:
        ValueError se o valor não representar um status válido
    """
    if isinstance(value, BaseStatus):
        return value.value
    if not isinstance(value, str):  # fallback
        value = str(value)
    norm = value.strip().lower()
    if norm not in {s.value for s in BaseStatus}:
        raise ValueError(f"Status inválido: {value}")
    return norm
