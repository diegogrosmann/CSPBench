"""Web application configuration module.

This module provides configuration management for the CSPBench web interface,
including environment variables, path management, and service initialization.
"""

import os
import threading
from typing import Optional

# from src.application.services.pipeline_service import PipelineService  # Module doesn't exist
from src.domain.config import CSPBenchConfig
from src.infrastructure.logging_config import get_logger
from src.infrastructure.utils.path_utils import (
    get_batch_directory,
    get_dataset_directory,
    get_output_base_directory,
)

logger = get_logger("CSPBench.WebConfig")


class WebConfig:
    """Configuration management for CSPBench web interface.
    
    This class manages web application configuration including environment
    variables, paths, and service initialization. It implements the
    singleton pattern to ensure consistent configuration across the application.
    
    Attributes:
        config (Optional[CSPBenchConfig]): Core CSPBench configuration.
        web_host (str): Web server bind address.
        web_port (int): Web server port number.
        debug (bool): Debug mode flag.
        log_level (str): Logging level.
        access_log (bool): Enable access logging.
        datasets_path (Path): Directory for dataset storage.
        batches_path (Path): Directory for batch configurations.
        outputs_path (Path): Directory for output files.
    """

    _instance: Optional["WebConfig"] = None
    _lock = threading.Lock()

    def __init__(self):
        """Initialize the configuration.
        
        Loads configuration from environment variables and sets up
        default paths using the standardized path utilities.
        """
        self.config: Optional[CSPBenchConfig] = None
        # self.pipeline_service: Optional[PipelineService] = None  # Service doesn't exist

        # Web-specific configuration
        self.web_host = os.getenv("WEB_HOST", "0.0.0.0")
        self.web_port = int(os.getenv("PORT", "8080"))
        self.debug = bool(os.getenv("DEBUG", "false").lower() in ["true", "1", "yes"])
        self.log_level = os.getenv("LOG_LEVEL", "info").upper()
        self.access_log = bool(
            os.getenv("ACCESS_LOG", "true").lower() in ["true", "1", "yes"]
        )

        # Paths - use standardized utility functions
        self.datasets_path = get_dataset_directory()
        self.batches_path = get_batch_directory()
        self.outputs_path = get_output_base_directory()

        # Note: path_utils functions already handle directory creation
        # self.datasets_path.mkdir(exist_ok=True)
        # self.batches_path.mkdir(exist_ok=True)
        # self.outputs_path.mkdir(exist_ok=True)

        logger.info("Web configuration initialized")

    @classmethod
    def get_instance(cls) -> "WebConfig":
        """Get singleton instance.
        
        Returns the singleton WebConfig instance, creating it if necessary.
        This method is thread-safe using double-checked locking pattern.
        
        Returns:
            WebConfig: The singleton configuration instance.
        """
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = cls()
        return cls._instance

    def initialize_services(self):
        """Initialize CSPBench configuration and services.
        
        Loads the core CSPBench configuration and initializes required
        services for the web application.
        
        Returns:
            bool: True if initialization succeeded, False otherwise.
            
        Note:
            This method handles exceptions gracefully and logs errors
            to help with debugging configuration issues.
        """
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
        """Get pipeline service instance.
        
        Returns the pipeline service instance, initializing services
        if they haven't been initialized yet.
        
        Returns:
            PipelineService or None: Pipeline service instance or None if not available.
            
        Note:
            Currently returns None as the PipelineService doesn't exist yet.
        """
        # if self.pipeline_service is None:
        #     self.initialize_services()
        # return self.pipeline_service
        return None  # Service doesn't exist

    def get_config(self) -> CSPBenchConfig:
        """Get CSPBench configuration.
        
        Returns the core CSPBench configuration, initializing it
        if it hasn't been loaded yet.
        
        Returns:
            CSPBenchConfig: The core configuration object.
        """
        if self.config is None:
            self.initialize_services()
        return self.config


# Global instance - lazy initialization
# web_config = WebConfig.get_instance()


def get_web_config() -> WebConfig:
    """Get the global web config instance.
    
    Convenience function to access the singleton WebConfig instance.
    
    Returns:
        WebConfig: The global configuration instance.
    """
    return WebConfig.get_instance()
