"""Deprecated exporters package.

All previous exporter implementations were removed in favor of the unified
`FinalizationService` located at `src.infrastructure.export.finalization_service`.

This module is intentionally left minimal to avoid import errors on stale
references while signaling deprecation.
"""

__all__: list[str] = []  # no public symbols
