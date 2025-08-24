"""Web application configuration module."""

import os
import threading
from pathlib import Path
from typing import Optional

# from src.application.services.pipeline_service import PipelineService  # Module doesn't exist
from src.domain.config import CSPBenchConfig
from src.infrastructure.logging_config import get_logger

logger = get_logger("CSPBench.WebConfig")


class WebConfig:
    """Configuration management for CSPBench web interface."""

    _instance: Optional["WebConfig"] = None
    _lock = threading.Lock()

    def __init__(self):
        """Initialize the configuration."""
        self.config: Optional[CSPBenchConfig] = None
        # self.pipeline_service: Optional[PipelineService] = None  # Service doesn't exist

        # Web-specific configuration
        self.web_host = os.getenv("WEB_HOST", "0.0.0.0")
        self.web_port = int(os.getenv("WEB_PORT", "8000"))
        self.debug = bool(os.getenv("DEBUG", "false").lower() in ["true", "1", "yes"])
        self.log_level = os.getenv("LOG_LEVEL", "info").upper()
        self.access_log = bool(
            os.getenv("ACCESS_LOG", "true").lower() in ["true", "1", "yes"]
        )

        # Paths
        self.datasets_path = Path(os.getenv("DATASET_PATH", "datasets"))
        self.batches_path = Path(os.getenv("BATCH_PATH", "batches"))
        self.outputs_path = Path(os.getenv("OUTPUT_PATH", "outputs"))

        # Ensure directories exist
        self.datasets_path.mkdir(exist_ok=True)
        self.batches_path.mkdir(exist_ok=True)
        self.outputs_path.mkdir(exist_ok=True)

        logger.info("Web configuration initialized")

    @classmethod
    def get_instance(cls) -> "WebConfig":
        """Get singleton instance."""
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = cls()
        return cls._instance

    def initialize_services(self):
        """Initialize CSPBench configuration and services."""
        try:
            logger.info("Loading CSPBench configuration...")
            from src.domain.config import load_cspbench_config

            self.config = load_cspbench_config()

            logger.info("Initializing pipeline service...")
            # self.pipeline_service = PipelineService(self.config)  # Service doesn't exist

            logger.info("Web services initialized successfully")
            return True

        except Exception as e:
            logger.error(f"Failed to initialize web services: {e}")
            return False

    def get_pipeline_service(self):  # -> PipelineService: Service doesn't exist
        """Get pipeline service instance."""
        # if self.pipeline_service is None:
        #     self.initialize_services()
        # return self.pipeline_service
        return None  # Service doesn't exist

    def get_config(self) -> CSPBenchConfig:
        """Get CSPBench configuration."""
        if self.config is None:
            self.initialize_services()
        return self.config


# Global instance
web_config = WebConfig.get_instance()
