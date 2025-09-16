"""CLI Commands Registration Module

This module centralizes the registration of CLI commands using the unified WorkManager.
It provides command implementations for batch execution, work management, algorithm
listing, web interface startup, and dataset generation.

The module handles:
- Batch file execution with optional progress monitoring
- Work item management (restart, list, status)
- Algorithm registry display
- Web interface startup with configurable parameters
- Interactive dataset generation wizard
"""

import os
from pathlib import Path

import typer

from src.domain.config import load_cspbench_config
from src.infrastructure.logging_config import get_logger

# Create module logger
command_logger = get_logger("CSPBench.CLI.Commands")


def _run_work_monitor(
    work_id: str, show_final_status: bool = True, fallback_message: bool = True
) -> bool:
    """Run progress monitor for a work item.

    This function starts the curses-based progress monitor for the specified work.
    It handles initialization errors gracefully and provides fallback messaging
    when the monitor cannot be started.

    Args:
        work_id: Work identifier to monitor
        show_final_status: Whether to show final status message after monitoring
        fallback_message: Whether to show fallback message if monitor fails

    Returns:
        True if monitoring completed successfully, False if failed
    """
    try:
        from src.presentation.cli.progress_monitor import ProgressMonitor

        monitor = ProgressMonitor(work_id)
        monitor.start()

        if show_final_status:
            _show_final_status_message(work_id)

        return True

    except Exception as e:
        if fallback_message:
            typer.echo(f"⚠️  Could not start monitor: {e}")
            typer.echo("🔄 Work continues running in background...")
            typer.echo(
                f"💡 Use 'cspbench monitor {work_id}' to try monitoring again later"
            )

        return False


def _show_final_status_message(work_id: str) -> None:
    """Show appropriate exit message based on final work status.

    Retrieves the final status of the work and displays an appropriate message
    with relevant information such as output paths or error details.

    Args:
        work_id: Work identifier to check status for
    """
    try:
        from src.application.services.work_service import get_work_service
        from src.domain.status import BaseStatus

        work_service = get_work_service()
        details = work_service.get(work_id)

        if not details:
            typer.echo(f"⚠️  Could not get final status for work {work_id}")
            return

        status = details.status.value

        # Show message based on final status
        if status == BaseStatus.COMPLETED.value:
            typer.echo("\n🎉 Work completed successfully!")
            typer.echo(f"✅ Final status: {status}")
            if details.output_path:
                typer.echo(f"📁 Results saved at: {details.output_path}")

        elif status == BaseStatus.FAILED.value:
            typer.echo("\n❌ Work failed during execution")
            typer.echo(f"💥 Final status: {status}")
            if details.error:
                typer.echo(f"🔍 Error: {details.error}")

        elif status == BaseStatus.ERROR.value:
            typer.echo("\n⚠️  Work finished with errors")
            typer.echo(f"🟡 Final status: {status}")
            typer.echo("🔍 Check logs for more details")
            if details.error:
                typer.echo(f"💬 Last error: {details.error}")

        elif status == BaseStatus.CANCELED.value:
            typer.echo("\n🛑 Work was canceled")
            typer.echo(f"⏹️  Final status: {status}")
            typer.echo("📝 Execution interrupted by user")

        else:
            # Handle other statuses (paused, running, queued, etc.)
            typer.echo(f"\n🔄 Work finished with status: {status}")
            if status in [
                BaseStatus.PAUSED.value,
                BaseStatus.RUNNING.value,
                BaseStatus.QUEUED.value,
            ]:
                typer.echo("💡 Work can be resumed at any time.")

    except Exception as e:
        typer.echo(f"⚠️  Error checking final status: {e}")


def register_commands(app: typer.Typer) -> None:
    """Register all CLI commands in the Typer application.

    This function sets up all available CLI commands including batch execution,
    work management, algorithm listing, web interface, and dataset generation.

    Args:
        app: Typer application instance to register commands with
    """
    command_logger.info("Registering CLI commands in Typer")

    @app.command()
    def batch(
        batch: Path = typer.Argument(
            ..., exists=True, readable=True, help="Batch YAML file"
        ),
        no_monitor: bool = typer.Option(
            False, "--no-monitor", help="Disable progress monitoring interface"
        ),
    ):
        """Execute pipeline using unified WorkManager.

        Uses WorkManager for consistent execution workflow between CLI and Web.
        Automatically shows progress monitoring interface unless --no-monitor is used.

        Args:
            batch: Path to batch configuration YAML file
            no_monitor: Flag to disable real-time progress monitoring
        """
        try:
            typer.echo("")
            typer.echo(f"🚀 Executing: batch={batch}")
            config = load_cspbench_config(batch)

            extra = {"origin": "cli", "batch_file": str(batch)}

            from src.application.services.work_service import get_work_service

            work_manager = get_work_service()

            work_id = work_manager.execute(config=config, extra=extra)

            typer.echo(f"✅ Work submitted successfully. work_id={work_id}")

            if not no_monitor:
                # Use the reusable monitor function
                _run_work_monitor(
                    work_id, show_final_status=True, fallback_message=True
                )
            else:
                # Wait for completion without monitoring
                from src.application.services.work_service import get_work_service

                work_service = get_work_service()
                work_service.wait_until_terminal(work_id)
                _show_final_status_message(work_id)

        except Exception as e:  # noqa: BLE001
            typer.echo(f"❌ Error: {e}")
            raise typer.Exit(1)

    # --- Work management commands ---
    work_app = typer.Typer(help="Manage WorkItems using WorkService")

    @work_app.command("restart")
    def work_restart(
        work_id: str = typer.Argument(..., help="Work ID to be restarted"),
        no_monitor: bool = typer.Option(
            False, "--no-monitor", help="Disable progress monitoring interface"
        ),
    ):
        """Restart a work item using WorkService.

        Args:
            work_id: Identifier of the work item to restart
            no_monitor: Flag to disable real-time progress monitoring
        """
        try:
            from src.application.services.work_service import get_work_service

            work_service = get_work_service()
            if not work_service.get(work_id):
                typer.echo("❌ Work not found")
                return

            work_manager = get_work_service()

            work_manager.restart(work_id)

            if not no_monitor:
                # Use the reusable monitor function
                _run_work_monitor(
                    work_id, show_final_status=True, fallback_message=True
                )
            else:
                # Wait for completion without monitoring
                from src.application.services.work_service import get_work_service

                work_service = get_work_service()
                work_service.wait_until_terminal(work_id)
                _show_final_status_message(work_id)

        except Exception as e:  # noqa: BLE001
            typer.echo(f"❌ Error: {e}")
            raise typer.Exit(1)

    @work_app.command("list")
    def work_list():
        """List all work items using WorkService.

        Displays a summary of all work items including their IDs and current status.
        """
        from src.application.services.work_service import get_work_service

        work_service = get_work_service()
        items = work_service.list()

        if not items:
            typer.echo("📭 No work items found")
            return

        typer.echo("📋 Work items list:")
        for item in items:
            status = item.status.value
            wid = item.id
            typer.echo(f"  • {wid}: {status}")

    @work_app.command("status")
    def work_status(work_id: str):
        """Get work item status using WorkService.

        Displays detailed status information for a specific work item including
        creation time, last update, and any error messages.

        Args:
            work_id: Identifier of the work item to check
        """
        from src.application.services.work_service import get_work_service

        work_service = get_work_service()
        work_item = work_service.get(work_id)

        if not work_item:
            typer.echo(f"❌ Work {work_id} not found")
            raise typer.Exit(1)

        typer.echo(f"📊 Status for work {work_id}:")
        typer.echo(f"  Status: {work_item.status.value}")
        typer.echo(f"  Created: {work_item.created_at}")
        typer.echo(f"  Updated: {work_item.updated_at}")
        if work_item.error:
            typer.echo(f"  Error: {work_item.error}")

    app.add_typer(work_app, name="work")

    @app.command()
    def algorithms():
        """List available algorithms.

        Displays all registered algorithms in the global algorithm registry
        along with their descriptions.
        """
        try:
            from src.domain.algorithms import global_registry  # lazy import

            typer.echo("🧠 Available algorithms:")
            for name, cls in global_registry.items():
                typer.echo("  • %s: %s" % (name, cls.__doc__ or "No description"))
            if not global_registry:
                typer.echo("  (No algorithms registered)")
            return
        except Exception as e:  # noqa: BLE001
            typer.echo("❌ Error listing algorithms: %s" % e)
            raise typer.Exit(1)

    @app.command()
    def web(
        host: str = typer.Option(None, "--host", "-h", help="Host to bind to"),
        port: int = typer.Option(None, "--port", "-p", help="Port to bind to"),
        dev: bool = typer.Option(None, "--dev", help="Run in development mode"),
    ) -> None:
        """Start the web interface.

        Launches the FastAPI web server with configurable host, port, and development
        mode settings. Uses environment variables as defaults when options are not provided.

        Args:
            host: Network interface to bind to (default: 0.0.0.0)
            port: Port number to listen on (default: 8080)
            dev: Enable development mode with auto-reload and debug logging
        """
        try:
            host = host if host is not None else os.getenv("WEB_HOST", "0.0.0.0")
            # Use PORT environment variable as primary, fallback to WEB_PORT for compatibility
            default_port = int(os.getenv("PORT", os.getenv("WEB_PORT", "8080")))
            port = int(port) if port is not None else default_port

            # Use dev flag or WEB_DEBUG env var
            debug_env = os.getenv("WEB_DEBUG", "false").lower() == "true"
            debug = dev if dev is not None else debug_env

            log_level = os.getenv("WEB_LOG_LEVEL", "info" if debug else "warning")
            access_log = os.getenv("WEB_ACCESS_LOG", str(debug)).lower() == "true"

            typer.echo("🌐 Starting CSPBench Web Interface...")
            typer.echo("🖥️  Host: %s" % host)
            typer.echo("🔌 Port: %s" % port)
            typer.echo("🛠️  Mode: %s" % ("Development" if debug else "Production"))
            typer.echo("📜 Log Level: %s" % log_level)
            typer.echo("📈 Access Log: %s" % ("Enabled" if access_log else "Disabled"))

            try:
                import uvicorn  # type: ignore

                # Import application to register routes
                from src.presentation.web.app import app as web_app  # noqa: F401
            except ImportError as e:  # noqa: BLE001
                typer.echo("❌ Web dependencies not installed: %s" % e)
                typer.echo("💡 Install with: pip install -r requirements.web.txt")
                raise typer.Exit(1)

            typer.echo("\n🚀 Web interface starting at http://%s:%s" % (host, port))
            typer.echo("🔗 Open the link in your browser")
            typer.echo("⏹️  Press Ctrl+C to stop the server")

            uvicorn.run(
                "src.presentation.web.app:app",
                host=host,
                port=port,
                reload=debug,
                log_level=log_level,
                access_log=access_log,
            )
            return
        except KeyboardInterrupt:
            typer.echo("\n🛑 Web server stopped")
            return
        except Exception as e:  # noqa: BLE001
            typer.echo("❌ Error starting web server: %s" % e)
            raise typer.Exit(1)

    @app.command(name="datasetsave")
    def dataset_save() -> None:
        """Interactive synthetic dataset generation wizard.

        Launches an interactive command-line wizard that guides users through
        the process of generating synthetic datasets or downloading real datasets
        from NCBI databases.
        """
        try:
            from src.infrastructure.orchestration.dataset_generation_orchestrator import (
                DatasetGenerationOrchestrator,
            )
            from src.infrastructure.utils.path_utils import get_dataset_directory

            dataset_path = get_dataset_directory()
            typer.echo("📁 Using dataset path: %s" % dataset_path)
            orchestrator = DatasetGenerationOrchestrator(base_path=str(dataset_path))
            result_path = orchestrator.run_interactive_generation()
            if result_path:
                typer.echo("\n🎉 Dataset saved successfully!")
            else:
                typer.echo("\n❌ Operation cancelled.")
            return
        except KeyboardInterrupt:
            typer.echo("\n🚫 Operation cancelled by user")
            raise typer.Exit(1)
        except Exception as e:  # noqa: BLE001
            typer.echo("❌ Error in dataset wizard: %s" % e)
            raise typer.Exit(1)
