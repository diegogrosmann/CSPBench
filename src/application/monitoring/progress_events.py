"""Stub progress events (legacy compatibility).
Provide TaskType enum for existing ExperimentService imports.
"""

from __future__ import annotations
from enum import Enum


class TaskType(str, Enum):
    EXPERIMENT = "experiment"
    OPTIMIZATION = "optimization"
    SENSITIVITY = "sensitivity"
