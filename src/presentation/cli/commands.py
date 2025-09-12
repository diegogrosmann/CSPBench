"""CLI Commands Registration Module

Centraliza registro de comandos usando WorkManager unificado.
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
    """
    Run progress monitor for a work item.

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
            typer.echo(f"⚠️  Não foi possível iniciar o monitor: {e}")
            typer.echo("🔄 O trabalho continua executando em segundo plano...")
            typer.echo(
                f"💡 Use 'cspbench monitor {work_id}' para tentar monitorar novamente mais tarde"
            )

        return False


def _show_final_status_message(work_id: str) -> None:
    """
    Show appropriate exit message based on final work status.

    Args:
        work_id: Work identifier to check status for
    """
    try:
        from src.application.services.work_service import get_work_service
        from src.domain.status import BaseStatus

        work_service = get_work_service()
        details = work_service.get(work_id)

        if not details:
            typer.echo(f"⚠️  Não foi possível obter status final do trabalho {work_id}")
            return

        status = details.status.value

        # Show message based on final status
        if status == BaseStatus.COMPLETED.value:
            typer.echo("\n🎉 Trabalho concluído com sucesso!")
            typer.echo(f"✅ Status final: {status}")
            if details.output_path:
                typer.echo(f"📁 Resultados salvos em: {details.output_path}")

        elif status == BaseStatus.FAILED.value:
            typer.echo("\n❌ Trabalho falhou durante a execução")
            typer.echo(f"💥 Status final: {status}")
            if details.error:
                typer.echo(f"🔍 Erro: {details.error}")

        elif status == BaseStatus.ERROR.value:
            typer.echo("\n⚠️  Trabalho finalizado com erros")
            typer.echo(f"🟡 Status final: {status}")
            typer.echo("🔍 Verifique os logs para mais detalhes")
            if details.error:
                typer.echo(f"💬 Último erro: {details.error}")

        elif status == BaseStatus.CANCELED.value:
            typer.echo("\n🛑 Trabalho foi cancelado")
            typer.echo(f"⏹️  Status final: {status}")
            typer.echo("📝 Execução interrompida pelo usuário")

        else:
            # Handle other statuses (paused, running, queued, etc.)
            typer.echo(f"\n🔄 Trabalho terminou com status: {status}")
            if status in [
                BaseStatus.PAUSED.value,
                BaseStatus.RUNNING.value,
                BaseStatus.QUEUED.value,
            ]:
                typer.echo("💡 O trabalho pode ser retomado a qualquer momento.")

    except Exception as e:
        typer.echo(f"⚠️  Erro ao verificar status final: {e}")


def register_commands(app: typer.Typer) -> None:
    """
    Register all CLI commands in the Typer application.
    """
    command_logger.info("Registrando comandos CLI no Typer")

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
        """
        try:
            typer.echo("")
            typer.echo(f"🚀 Executando: batch={batch}")
            config = load_cspbench_config(batch)

            extra = {"origin": "cli", "batch_file": str(batch)}
            
            from src.application.services.work_service import get_work_service
            
            work_manager = get_work_service()

            work_id = work_manager.execute(config=config, extra=extra)

            typer.echo("✅ Trabalho submetido com sucesso. work_id={work_id}")

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
            typer.echo(f"❌ Erro: {e}")
            raise typer.Exit(1)

    # --- Work management commands ---
    work_app = typer.Typer(help="Gerencia WorkItems usando WorkService")

    @work_app.command("restart")
    def work_restart(
        work_id: str = typer.Argument(
            ..., exists=True, readable=True, help="Work ID a ser Reiniciado"
        ),
        no_monitor: bool = typer.Option(
            False, "--no-monitor", help="Disable progress monitoring interface"
        ),
    ):
        """Restart a work item using WorkService."""
        try:
            from src.application.services.work_service import get_work_service

            work_service = get_work_service()
            if not work_service.get(work_id):
                typer.echo("❌ Work não encontrado")
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
            typer.echo(f"❌ Erro: {e}")
            raise typer.Exit(1)

    @work_app.command("list")
    def work_list():
        """List all work items using WorkService."""
        from src.application.services.work_service import get_work_service

        work_service = get_work_service()
        items = work_service.list()

        if not items:
            typer.echo("📭 Nenhum trabalho encontrado")
            return

        typer.echo("📋 Lista de trabalhos:")
        for item in items:
            status = item.status.value
            wid = item.id
            typer.echo(f"  • {wid}: {status}")

    @work_app.command("status")
    def work_status(work_id: str):
        """Get work item status using WorkService."""
        from src.application.services.work_service import get_work_service

        work_service = get_work_service()
        work_item = work_service.get(work_id)

        if not work_item:
            typer.echo(f"❌ Trabalho {work_id} não encontrado")
            raise typer.Exit(1)

        typer.echo(f"📊 Status do trabalho {work_id}:")
        typer.echo(f"  Status: {work_item.status.value}")
        typer.echo(f"  Criado: {work_item.created_at}")
        typer.echo(f"  Atualizado: {work_item.updated_at}")
        if work_item.error:
            typer.echo(f"  Erro: {work_item.error}")

    app.add_typer(work_app, name="work")

    @app.command()
    def algorithms():
        """List available algorithms."""
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
        """Start the web interface."""
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
        """Interactive synthetic dataset generation wizard."""
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
