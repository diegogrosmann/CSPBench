"""
Web application configuration and initialization.
"""

import logging
from pathlib import Path
from typing import Any, Dict, Optional

import yaml

from src.application.services.experiment_service import ExperimentService
from src.domain.algorithms import global_registry
from src.infrastructure.orchestrators.session_manager import SessionManager

logger = logging.getLogger(__name__)


class WebConfig:
    """Web application configuration manager."""

    def __init__(self):
        self.config: Optional[Dict[str, Any]] = None
        self.experiment_service: Optional[ExperimentService] = None
        self.session_manager: Optional[SessionManager] = None
        self._initialized = False

    def load_config(self) -> Dict[str, Any]:
        """Load configuration from settings.yaml file."""
        if self.config is not None:
            return self.config

        try:
            config_path = Path("config/settings.yaml")

            if not config_path.exists():
                logger.warning(
                    f"Configuration file not found: {config_path}, using defaults"
                )
                self.config = self._get_default_config()
            else:
                with open(config_path, encoding="utf-8") as f:
                    self.config = yaml.safe_load(f)

            return self.config

        except Exception as e:
            logger.error(f"Error loading configuration: {e}")
            self.config = self._get_default_config()
            return self.config

    def _get_default_config(self) -> Dict[str, Any]:
        """Get default configuration."""
        return {
            "application": {"name": "CSPBench", "version": "0.1.0"},
            "infrastructure": {"executor": {"config": {"timeout_default": 300}}},
            "web": {
                "security": {
                    "max_dataset_size": 50 * 1024 * 1024,  # 50MB
                    "allowed_extensions": [".fasta", ".fa", ".txt", ".seq"],
                    "max_sessions_per_ip": 10,
                    "session_timeout": 7200,  # 2 hours
                    "max_filename_length": 255,
                },
                "rate_limiting": {"enabled": True, "requests_per_minute": 60},
            },
        }

    def initialize_services(self) -> bool:
        """Initialize infrastructure services."""
        if self._initialized:
            return True

        try:
            # Import algorithms to populate global_registry
            import algorithms

            logger.info(f"Loaded {len(global_registry)} algorithms")

            # Load configuration
            config = self.load_config()

            # Initialize infrastructure components
            from src.infrastructure.io.exporters.json_exporter import JsonExporter
            from src.infrastructure.orchestrators.executors import Executor
            from src.infrastructure.persistence.algorithm_registry import (
                DomainAlgorithmRegistry,
            )
            from src.infrastructure.persistence.dataset_repository import (
                FileDatasetRepository,
            )

            dataset_repository = FileDatasetRepository("./datasets")
            algorithm_registry = DomainAlgorithmRegistry()
            executor = Executor()
            exporter = JsonExporter("./outputs")

            self.experiment_service = ExperimentService(
                dataset_repo=dataset_repository,
                exporter=exporter,
                executor=executor,
                algo_registry=algorithm_registry,
            )

            self.session_manager = SessionManager(config)

            self._initialized = True
            logger.info("All services initialized successfully")
            return True

        except Exception as e:
            logger.error(f"Failed to initialize services: {e}")
            return False

    def get_experiment_service(self) -> Optional[ExperimentService]:
        """Get experiment service instance."""
        if not self._initialized:
            self.initialize_services()
        return self.experiment_service

    def get_session_manager(self) -> Optional[SessionManager]:
        """Get session manager instance."""
        if not self._initialized:
            self.initialize_services()
        return self.session_manager


# Global configuration instance
web_config = WebConfig()
