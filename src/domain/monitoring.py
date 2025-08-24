"""Monitoring interfaces for CSP algorithms.

Módulo de domínio puro contendo APENAS o protocolo de monitoramento
(`AlgorithmMonitor`). Implementações concretas (Null, Composite, Throttled,
adapters para persistência, etc.) devem viver fora desta unidade para manter
o domínio mínimo e sem dependências desnecessárias.

Extensões futuras podem adicionar novos métodos ao protocolo; para manter
retrocompatibilidade, forneça defaults aqui ou use verificações com `hasattr`.
"""

from __future__ import annotations

from typing import Protocol, runtime_checkable, Any


@runtime_checkable
class AlgorithmMonitor(Protocol):
    """Contrato mínimo para monitoramento de execução de algoritmos.

    Todos os métodos devem ser best-effort: implementações **não** devem
    levantar exceções para não interromper o fluxo do algoritmo.

    Parâmetros comuns:
        progress (float): 0.0 a 100.0 representando percentual (ou pode extrapolar >100 se fizer sentido em casos especiais).
        message (str): descrição curta do evento.
        **data: payload adicional serializável / simples.
    """

    # Métodos principais --------------------------------------------------
    def on_progress(self, progress: float, message: str, /, **data: Any) -> None:
        """Evento de progresso.

        Recomendação: Implementações podem aplicar throttling (tempo ou delta
        mínimo de progresso) para evitar flooding.
        """
        ...  # pragma: no cover

    def on_warning(self, message: str, /, **data: Any) -> None:
        """Evento de aviso não fatal."""
        ...  # pragma: no cover

    # Métodos opcionais ---------------------------------------------------
    def on_error(
        self, message: str, exc: Exception | None = None, /, **data: Any
    ) -> None:  # noqa: D401,E501
        """Evento de erro (opcional)."""
        ...  # pragma: no cover

    def is_cancelled(self) -> bool:
        """Indica se a execução foi cancelada externamente.

        Algoritmos podem consultar em loops longos e abortar graciosamente.
        Default: False.
        """
        return False  # pragma: no cover
