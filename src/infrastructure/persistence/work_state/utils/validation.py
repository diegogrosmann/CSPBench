"""Validation utilities for work state persistence."""

from typing import Any

from src.domain.status import BaseStatus, normalize_status

# Conjunto de valores de status válidos (strings) para validação rápida
VALID_WORK_STATUSES = {s.value for s in BaseStatus}
VALID_EVENT_TYPES = {"error", "warning"}
VALID_EVENT_CATEGORIES = {
    "work",
    "task",
    "dataset",
    "preset",
    "combination",
    "unit",
    "other",
}


def validate_required_field(value: Any, field_name: str) -> None:
    """Validate that a required field is present and not empty.

    Args:
        value: Field value to validate
        field_name: Name of the field (for error messages)
    Raises:
        ValueError: If value is None or empty string
    """
    if value is None or (isinstance(value, str) and value.strip() == ""):
        raise ValueError(f"Campo obrigatório '{field_name}' ausente ou vazio")


def validate_status(status: str | BaseStatus) -> None:
    """Valida status aceitando enum ou string.

    Utiliza o helper normalize_status para padronização e validação.
    """
    try:
        normalize_status(status)  # Vai levantar ValueError se inválido
    except ValueError:
        # Re-raise com mensagem mais específica para contexto de validação
        raise


def validate_event_type(event_type: str) -> None:
    """Validate event type."""
    if event_type not in VALID_EVENT_TYPES:
        raise ValueError(
            f"Tipo de evento inválido: {event_type}. Válidos: {VALID_EVENT_TYPES}"
        )


def validate_event_category(event_category: str) -> None:
    """Validate event category."""
    if event_category not in VALID_EVENT_CATEGORIES:
        raise ValueError(
            f"Categoria de evento inválida: {event_category}. Válidas: {VALID_EVENT_CATEGORIES}"
        )
