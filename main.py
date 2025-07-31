#!/usr/bin/env python3
"""
CSPBench Main Entry Point - Hexagonal Architecture

This is the main entry point for the refactored CSPBench with hexagonal architecture.
Implements dependency injection and centralized configuration.

Functionality:
    - Loading configurations from config/settings.yaml
    - Configuration-based dependency injection
    - Support for override via environment variables
    - Modern CLI using Typer

Main Usage (according to guidelines):
    ```bash
    # Interactive menu (no arguments)
    python main.py

    # General help
    python main.py --help

    # List available algorithms
    python main.py --algorithms

    # Generate/save synthetic datasets
    python main.py --datasetsave

    # Execute batch directly (.yaml/.yml file)
    python main.p    view-report <n>               Open report in browser
    web                              Start web interface (with options)

Examples:atches/example.yaml
    python main.py configuration.yml
    ```

Specific Commands:
    ```bash
    # Basic system test
    python main.py test

    # Execute algorithm on dataset
    python main.py run <algorithm> <dataset>

    # Execute batch via command
    python main.py batch <file.yaml>

    # Session management
    python main.py sessions
    python main.py show-session <name>
    python main.py view-report <name>
    python main.py cleanup

    # System information
    python main.py config-info
    python main.py algorithms
    ```

Environment Variables:
    - EXECUTOR_IMPL: Override executor (Executor)
    - EXPORT_FORMAT: Override export format (json, csv, txt)
    - DATASET_PATH: Override dataset base path
    - NCBI_EMAIL: Email for NCBI API
    - NCBI_API_KEY: NCBI API key
"""

import os
import sys
from pathlib import Path
from typing import Any, Dict, Optional

import typer
import yaml

# Add the root directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

# IMPORTANT: Import algorithms first to load the global_registry
import algorithms
from src.application.services.experiment_service import ExperimentService
from src.domain import SyntheticDatasetGenerator
from src.infrastructure import (
    CsvExporter,
    DomainAlgorithmRegistry,
    Executor,
    FileDatasetRepository,
    JsonExporter,
    SessionManager,
    TxtExporter,
)
from src.infrastructure.logging_config import LoggerConfig, get_logger
from src.presentation.cli.commands import register_commands

app = typer.Typer(
    name="cspbench",
    help="CSPBench - Framework for Closest String Problem",
    add_completion=False,
    no_args_is_help=False,  # Don't show help when no arguments
    rich_markup_mode=None,  # Disable rich to avoid errors
)

# Global variables for DI
experiment_service: Optional[ExperimentService] = None
config: Optional[Dict[str, Any]] = None
session_manager: Optional[SessionManager] = None


def load_config() -> Dict[str, Any]:
    """Load configuration from settings.yaml file."""
    config_path = Path("config/settings.yaml")

    if not config_path.exists():
        typer.echo(f"❌ Configuration file not found: {config_path}")
        raise typer.Exit(1)

    try:
        with open(config_path, "r", encoding="utf-8") as f:
            return yaml.safe_load(f)
    except Exception as e:
        typer.echo(f"❌ Error loading configuration: {e}")
        raise typer.Exit(1)


def create_dataset_repository(config: Dict[str, Any]) -> FileDatasetRepository:
    """Create dataset repository based on configuration."""
    repo_config = config["infrastructure"]["dataset_repository"]["config"]
    base_path = os.getenv("DATASET_PATH", repo_config["base_path"])

    return FileDatasetRepository(base_path)


def create_algorithm_registry(config: Dict[str, Any]) -> DomainAlgorithmRegistry:
    """Create algorithm registry based on configuration."""
    return DomainAlgorithmRegistry()


def create_entrez_repository(config: Dict[str, Any]) -> Optional[Any]:
    """Create Entrez dataset repository based on configuration."""
    try:
        from src.infrastructure.persistence.entrez_dataset_repository import NCBIEntrezDatasetRepository
        
        # Check if required environment variables are set
        email = os.getenv("NCBI_EMAIL")
        if not email:
            typer.echo("⚠️  NCBI_EMAIL not set. Entrez datasets will not be available.", err=True)
            return None
        
        # Create repository
        entrez_repo = NCBIEntrezDatasetRepository()
        
        # Test availability
        if entrez_repo.is_available():
            typer.echo("✅ Entrez datasets enabled")
            return entrez_repo
        else:
            typer.echo("⚠️  Entrez service not available. Check network connection.", err=True)
            return None
            
    except Exception as e:
        typer.echo(f"⚠️  Failed to initialize Entrez repository: {str(e)}", err=True)
        return None


def create_executor(config: Dict[str, Any]) -> Executor:
    """Create executor based on configuration."""
    executor_impl = os.getenv("EXECUTOR_IMPL", "Executor")

    if executor_impl == "Executor":
        return Executor()
    else:
        typer.echo(f"⚠️  Executor '{executor_impl}' not implemented, using Executor")
        return Executor()


def create_exporter(config: Dict[str, Any]):
    """Create exporter based on configuration."""
    global session_manager

    # Configuração do exportador (se existir na configuração legada)
    exporter_config = (
        config.get("infrastructure", {}).get("exporter", {}).get("config", {})
    )
    export_fmt = os.getenv("EXPORT_FORMAT", exporter_config.get("format", "json"))

    # Use session directory if SessionManager exists
    if session_manager:
        output_path = session_manager.get_result_dir()
    else:
        # Use unified output configuration
        output_config = config.get("output", {})
        if output_config:
            base_directory = output_config.get("base_directory", "outputs/{session}")
            if "{session}" in base_directory:
                output_path = base_directory.replace("{session}", "")
            else:
                output_path = base_directory
        else:
            # Fallback to old configuration
            result_config = config.get("infrastructure", {}).get("result", {})
            output_path = os.getenv(
                "OUTPUT_PATH", result_config.get("base_result_dir", "./outputs/results")
            )

    if export_fmt.lower() == "json":
        return JsonExporter(output_path, config)  # Pass complete configuration
    elif export_fmt.lower() == "csv":
        return CsvExporter(output_path)
    elif export_fmt.lower() == "txt":
        return TxtExporter(output_path)
    else:
        typer.echo(f"⚠️  Format '{export_fmt}' not supported, using JSON")
        return JsonExporter(output_path, config)  # Pass complete configuration


def create_monitoring_service(config: Dict[str, Any]):
    """Create monitoring service based on configuration."""
    from src.application.services.basic_monitoring_service import BasicMonitoringService

    return BasicMonitoringService(config)


def _initialize_logging(config: Dict[str, Any]) -> None:
    """Initialize logging system based on configuration."""
    global session_manager

    # Default configurations
    default_level = "INFO"
    default_log_dir = "outputs/logs"

    try:
        # Create SessionManager if it doesn't exist
        if session_manager is None:
            session_manager = SessionManager(config)

        # Create new session
        session_folder = session_manager.create_session()

        # Try to extract new unified output configuration
        output_config = config.get("output", {})
        log_config = output_config.get("logging", {})

        # Check if logs are enabled
        if not log_config.get("enabled", True):
            return

        # Extract configurations
        log_level = log_config.get("level", default_level)

        # Use session path for logs
        log_file_path = session_manager.get_log_path()

        # Initialize LoggerConfig with specific session path
        LoggerConfig.initialize(
            level=log_level,
            log_file_path=log_file_path,
            base_name="cspbench",
        )

        # Create logger and register initialization
        logger = get_logger(__name__)
        logger.info("=" * 60)
        logger.info("CSPBench - Logging system initialized")
        logger.info(f"Log level: {log_level}")
        logger.info(f"Session: {session_folder}")
        logger.info(f"Log file: {log_file_path}")
        logger.info(f"Timestamp: {__import__('datetime').datetime.now()}")
        logger.info("=" * 60)

    except Exception as e:
        # Fallback to basic configuration in case of error
        LoggerConfig.initialize(level=default_level, log_dir=default_log_dir)
        logger = get_logger(__name__)
        logger.warning(f"Error configuring logging, using defaults: {e}")


def _update_logging_from_batch_config(batch_config: Dict[str, Any]) -> None:
    """Update logging configuration based on specific batch file."""
    try:
        # Try new unified output configuration
        output_config = batch_config.get("output", {})
        log_config = output_config.get("logging", {})

        # Check if logs are enabled
        if not log_config.get("enabled", True):
            return

        # Update level if specified
        if "level" in log_config:
            new_level = log_config["level"]
            LoggerConfig.set_level(new_level)
            logger = get_logger(__name__)
            logger.info(f"Log level updated to: {new_level}")
    except Exception as e:
        logger = get_logger(__name__)
        logger.warning(f"Error updating logging configuration from batch: {e}")


def initialize_service() -> ExperimentService:
    """Initialize experiment service with dependency injection."""
    global config, experiment_service

    if experiment_service is None:
        if config is None:
            config = load_config()

        # Initialize logging system
        _initialize_logging(config)

        # Create infrastructure components
        dataset_repo = create_dataset_repository(config)
        algo_registry = create_algorithm_registry(config)
        executor = create_executor(config)
        exporter = create_exporter(config)
        entrez_repo = create_entrez_repository(config)
        monitoring_service = create_monitoring_service(config)

        # Create service with DI
        experiment_service = ExperimentService(
            dataset_repo=dataset_repo,
            exporter=exporter,
            executor=executor,
            algo_registry=algo_registry,
            entrez_repo=entrez_repo,
            monitoring_service=monitoring_service,
            session_manager=session_manager,
        )

        typer.echo("✅ CSPBench initialized successfully!")

    return experiment_service


def show_interactive_menu() -> None:
    """Show interactive interface for batch file selection."""
    print("CSPBench v0.1.0 - Framework for Closest String Problem")
    print("=" * 60)
    print()

    # List available batch files
    batches_dir = Path("batches")
    if not batches_dir.exists():
        _show_commands_help("📁 Directory 'batches/' not found")
        return

    batch_files = list(batches_dir.glob("*.yaml")) + list(batches_dir.glob("*.yml"))

    if not batch_files:
        _show_commands_help("📋 No batch files found in 'batches/'")
        return

    _display_batch_files(batch_files)
    _display_manual_commands()

    selected_file = _get_user_selection(batch_files)
    if selected_file:
        _execute_selected_batch(selected_file)


def _show_commands_help(message: str) -> None:
    """Show help with available commands."""
    print(message)
    print()
    print("Available commands:")
    print("  test         - Basic system test")
    print("  algorithms   - List available algorithms")
    print("  config-info  - Show configuration")
    print("  run          - Execute single algorithm")
    print()
    print("Usage: python main.py <command> [arguments]")
    print("Example: python main.py test")


def _display_batch_files(batch_files: list) -> None:
    """Display list of batch files with descriptions."""
    print("📋 Available batch files:")
    print()

    for i, batch_file in enumerate(batch_files, 1):
        description = _extract_file_description(batch_file)
        print(f"  {i}. {batch_file.name}")
        print(f"     {description}")
        print()


def _extract_file_description(batch_file: Path) -> str:
    """Extract description from batch file from comments."""
    try:
        with open(batch_file, "r", encoding="utf-8") as f:
            first_lines = f.read(300)  # First 300 chars

        for line in first_lines.split("\n"):
            line = line.strip()
            if line.startswith("#") and len(line) > 1:
                desc_text = line[1:].strip()
                if desc_text and not desc_text.startswith("="):
                    return desc_text

        return "Batch configuration file"
    except:
        return "Batch configuration file"


def _display_manual_commands() -> None:
    """Display list of available manual commands."""
    print("📋 Available manual commands:")
    print("  test         - Basic system test")
    print("  algorithms   - List available algorithms")
    print("  config-info  - Show configuration")
    print("  run          - Execute single algorithm")
    print()


def _get_user_selection(batch_files: list) -> Optional[Path]:
    """Get user selection and validate."""
    try:
        choice = input(
            "💡 Select a file (number) or press Enter to exit: "
        ).strip()

        if choice == "":
            print("👋 Goodbye!")
            return None

        if choice.isdigit():
            choice_num = int(choice)
            if 1 <= choice_num <= len(batch_files):
                return batch_files[choice_num - 1]
            else:
                print("❌ Invalid number!")
                return None
        else:
            print("❌ Invalid input!")
            return None

    except KeyboardInterrupt:
        print("\n👋 Goodbye!")
        return None
    except EOFError:
        print("\n👋 Goodbye!")
        return None


def _execute_selected_batch(selected_file: Path) -> None:
    """Execute the selected batch file."""
    print(f"\n🚀 Executing: {selected_file.name}")
    print("-" * 40)

    # Execute selected batch
    import sys

    sys.argv = ["main.py", "batch", str(selected_file)]
    app()


def show_algorithms_and_exit():
    """Show available algorithms and exit."""
    try:
        # Import from algorithms module to activate auto-discovery
        import algorithms
        from src.domain.algorithms import global_registry

        print("🧠 Available algorithms:")
        for name, cls in global_registry.items():
            print(f"  • {name}: {cls.__doc__ or 'No description'}")

        if not global_registry:
            print("  (No algorithms registered)")

    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1)

    sys.exit(0)


def show_datasetsave_and_exit():
    """Execute interactive dataset generation wizard and exit."""
    try:
        from src.infrastructure.orchestrators.dataset_generation_orchestrator import (
            DatasetGenerationOrchestrator,
        )

        # Create and execute orchestrator
        orchestrator = DatasetGenerationOrchestrator()
        result_path = orchestrator.run_interactive_generation()

        if result_path:
            print(f"\n🎉 Dataset saved successfully!")
        else:
            print("\n❌ Operation cancelled.")

    except Exception as e:
        print(f"❌ Error: {e}")


def start_web_interface():
    """Start the web interface with default settings and exit."""
    try:
        print("🌐 Starting CSPBench Web Interface...")
        print("🖥️  Host: 0.0.0.0")
        print("🔌 Port: 8000")
        print("🛠️  Mode: Production")
        
        # Check if web dependencies are available
        try:
            import uvicorn
            from src.presentation.web.app import app as web_app
        except ImportError as e:
            print(f"❌ Web dependencies not installed: {e}")
            print("💡 Install with: pip install -r requirements.web.txt")
            sys.exit(1)
        
        print(f"\n🚀 Web interface starting at http://localhost:8000")
        print("🔗 Click the link above or paste it into your browser")
        print("⏹️  Press Ctrl+C to stop the server")
        
        # Start uvicorn server with default settings
        uvicorn.run(
            "src.presentation.web.app:app",
            host="0.0.0.0",
            port=8000,
            reload=False,
            log_level="info",
            access_log=False
        )
        
    except KeyboardInterrupt:
        print("\n🛑 Web server stopped")
        sys.exit(0)
    except Exception as e:
        print(f"❌ Error starting web server: {e}")
        sys.exit(1)


def _display_manual_commands() -> None:
    """Display list of available manual commands."""
    print("📋 Available manual commands:")
    print("  test         - Basic system test")
    print("  algorithms   - List available algorithms")
    print("  config-info  - Show configuration")
    print("  run          - Execute single algorithm")
    print()


def _get_user_selection(batch_files: list) -> Optional[Path]:
    """Get user selection and validate."""
    try:
        choice = input(
            "💡 Select a file (number) or press Enter to exit: "
        ).strip()

        if choice == "":
            print("👋 Goodbye!")
            return None

        if choice.isdigit():
            choice_num = int(choice)
            if 1 <= choice_num <= len(batch_files):
                return batch_files[choice_num - 1]
            else:
                print("❌ Invalid number!")
                return None
        else:
            print("❌ Invalid input!")
            return None

    except KeyboardInterrupt:
        print("\n👋 Goodbye!")
        return None
    except EOFError:
        print("\n👋 Goodbye!")
        return None


def _execute_selected_batch(selected_file: Path) -> None:
    """Execute the selected batch file."""
    print(f"\n🚀 Executing: {selected_file.name}")
    print("-" * 40)

    # Execute the selected batch
    import sys

    sys.argv = ["main.py", "batch", str(selected_file)]
    app()


def show_algorithms_and_exit():
    """Show available algorithms and exit."""
    try:
        # Import from algorithms module to activate auto-discovery
        import algorithms
        from src.domain.algorithms import global_registry

        print("🧠 Available algorithms:")
        for name, cls in global_registry.items():
            print(f"  • {name}: {cls.__doc__ or 'No description'}")

        if not global_registry:
            print("  (No algorithms registered)")

    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1)

    sys.exit(0)


def show_datasetsave_and_exit():
    """Execute interactive dataset generation wizard and exit."""
    try:
        from src.infrastructure.orchestrators.dataset_generation_orchestrator import (
            DatasetGenerationOrchestrator,
        )

        # Create and execute orchestrator
        orchestrator = DatasetGenerationOrchestrator()
        result_path = orchestrator.run_interactive_generation()

        if result_path:
            print(f"\n🎉 Dataset saved successfully!")
        else:
            print("\n❌ Operation cancelled.")

    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1)

    sys.exit(0)


def initialize_service_with_batch_config(batch_config: Dict[str, Any]) -> ExperimentService:
    """Initialize experiment service with batch-specific configuration."""
    global config, experiment_service, session_manager

    if config is None:
        config = load_config()

    # Merge batch output configuration with system config for session management
    merged_config = config.copy()
    if 'output' in batch_config:
        merged_config['output'] = batch_config['output']

    # Initialize logging system with batch configuration
    _initialize_logging(merged_config)

    # Create infrastructure components
    dataset_repo = create_dataset_repository(config)
    algo_registry = create_algorithm_registry(config)
    executor = create_executor(config)
    exporter = create_exporter(config)
    entrez_repo = create_entrez_repository(config)
    monitoring_service = create_monitoring_service(config)

    # Create service with DI
    experiment_service = ExperimentService(
        dataset_repo=dataset_repo,
        exporter=exporter,
        executor=executor,
        algo_registry=algo_registry,
        entrez_repo=entrez_repo,
        monitoring_service=monitoring_service,
        session_manager=session_manager,
    )

    typer.echo("✅ CSPBench initialized successfully!")
    return experiment_service


def execute_batch_file(batch_file: str):
    """Execute batch file directly."""
    try:
        batch_path = Path(batch_file)

        if not batch_path.exists():
            print(f"❌ File not found: {batch_file}")
            sys.exit(1)

        if not batch_path.suffix.lower() in [".yaml", ".yml"]:
            print(f"❌ File must be .yaml or .yml: {batch_file}")
            sys.exit(1)

        # Load batch configuration FIRST to use its output configuration for logging
        print(f"📋 Executing batch: {batch_file}...")
        
        # Load batch configuration to extract output config for proper session management
        from src.application.services.config_parser import ConfigParser
        batch_config = ConfigParser.parse_config(str(batch_path))
        batch_dict = {
            'output': {
                'base_directory': batch_config.export.destination if batch_config.export else "outputs/{session}",
                'logging': {
                    'enabled': True,  # Default to enabled
                    'level': batch_config.logging.level if batch_config.logging else "INFO"
                }
            }
        }
        
        # Initialize service with batch-specific configuration for session management
        service = initialize_service_with_batch_config(batch_dict)

        result = service.run_batch(str(batch_path))
        print(f"✅ Batch completed: {result.get('summary', 'Completed')}")

    except KeyboardInterrupt:
        print("\n🚫 Operation cancelled by user (Ctrl+C)")
        print("📋 Batch interrupted during execution")
        sys.exit(0)  # Exit with code 0 to indicate intentional cancellation
    except Exception as e:
        print(f"❌ Batch error: {e}")
        sys.exit(1)


# Register all CLI commands
register_commands(app, initialize_service)


def main(args: Optional[list] = None):
    """Main function for programmatic execution."""
    import sys

    if args is None:
        args = sys.argv[1:]

    try:
        # If it's a batch file, execute directly
        if len(args) == 1 and args[0].endswith((".yaml", ".yml")):
            execute_batch_file(args[0])
            return

        # Otherwise, use Typer command system
        original_argv = sys.argv[:]
        try:
            sys.argv = ["main.py"] + args
            app()
        finally:
            sys.argv = original_argv

    except KeyboardInterrupt:
        print("\n🚫 Operation cancelled by user (Ctrl+C)")
        sys.exit(0)
    except SystemExit:
        # Allow SystemExit (including sys.exit) to pass without handling
        raise
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    import sys

    try:
        # If executed without arguments, show interactive interface
        if len(sys.argv) == 1:
            show_interactive_menu()
        elif len(sys.argv) == 2:
            arg = sys.argv[1]

            # Implement special flags
            if arg == "--algorithms":
                show_algorithms_and_exit()
            elif arg == "--datasetsave":
                show_datasetsave_and_exit()
            elif arg == "--web":
                # Start web interface directly
                start_web_interface()
            elif arg == "--help":
                print(
                    """
CSPBench v0.1.0 - Framework for Closest String Problem

Usage:
    python main.py                    Interactive menu
    python main.py --help            This help
    python main.py --algorithms      List available algorithms
    python main.py --datasetsave     Generate/save datasets
    python main.py --web             Start web interface
    python main.py <file.yaml>       Execute batch
    python main.py <command>         Execute specific command

Available commands:
    test                             Basic system test
    run <algorithm> <dataset>        Execute algorithm on dataset
    batch <file.yaml>                Execute batch
    algorithms                       List algorithms
    config-info                      Show configuration
    sessions                         List sessions
    cleanup                          Remove old sessions
    show-session <name>              Show session details
    view-report <name>               Open report in browser

Examples:
    python main.py test
    python main.py run Baseline test.fasta
    python main.py batch batches/example.yaml
    python main.py batches/example.yaml
    python main.py --web
    python main.py web --dev --port 8080
"""
                )
                sys.exit(0)
            elif not arg.startswith("-") and (
                arg.endswith(".yaml") or arg.endswith(".yml")
            ):
                # It's a batch file - execute directly
                execute_batch_file(arg)
                sys.exit(0)
            else:
                # Use normal Typer CLI
                app()
        else:
            # Use normal Typer CLI
            app()
    except KeyboardInterrupt:
        print("\n🚫 Operation cancelled by user (Ctrl+C)")
        sys.exit(0)
    except SystemExit:
        # Allow SystemExit (including sys.exit) to pass without handling
        raise
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        sys.exit(1)
