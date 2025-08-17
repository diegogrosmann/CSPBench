"""CLI Commands Registration Module

Centraliza registro de comandos. Adicionado comando `pipeline run` para novo orquestrador.
"""

import os
import json
from pathlib import Path

import typer

from src.domain.config import load_cspbench_config
from src.application.services.pipeline_service import PipelineService
from src.application.work.manager import get_work_manager
from src.presentation.display.terminal_monitor import TerminalMonitor
from src.infrastructure.monitoring.monitor_interface import LoggingMonitor, NoOpMonitor
from src.infrastructure.logging_config import setup_logging_from_env, get_logger
import logging

# Create module logger
command_logger = get_logger("CSPBench.CLI.Commands")


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
        monitor: str = typer.Option(
            "terminal", 
            help="Monitor type: 'none' (no monitoring), 'terminal' (visual terminal monitor), 'log' (logging monitor)"
        ),
    ):
        """Submete pipeline (experiment|optimization|sensitivity) ao WorkManager.

        Usa PipelineService para orquestrar execução assíncrona e registra WorkItem via WorkManager.
        """
        try:
            typer.echo(f"🚀 Submetendo pipeline: batch={batch}")
            config = load_cspbench_config(batch)
            wm = get_work_manager()

            # Configurar monitor baseado na opção escolhida
            display_monitor = None
            if monitor == "terminal":
                display_monitor = TerminalMonitor()
                typer.echo("📊 Monitor de terminal ativado")
            elif monitor == "log":
                # Configurar logging para arquivo usando variáveis de ambiente
                setup_logging_from_env()
                logger = get_logger("CSPBench.Pipeline")
                display_monitor = LoggingMonitor(logger)
                typer.echo("📝 Monitor de log ativado")
            elif monitor == "none":
                display_monitor = NoOpMonitor()
                typer.echo("🔇 Monitor desabilitado")
                command_logger.info("Monitor desabilitado")
            else:
                error_msg = f"❌ Tipo de monitor inválido: {monitor}"
                typer.echo(error_msg)
                typer.echo("💡 Opções válidas: 'none', 'terminal', 'log'")
                command_logger.error(f"Tipo de monitor inválido: {monitor}")
                raise typer.Exit(1)

            wid = PipelineService.run(config, monitor=display_monitor)

            item = wm.get(wid)
            if not item:
                typer.echo("❌ Falha ao registrar WorkItem")
                raise typer.Exit(1)
            typer.echo(f"🆔 work_id={wid} status={item['status']}")

            typer.echo(f"⏳ Aguardando o término...")
            final_status = wm.wait_until_terminal(wid)
            typer.echo(f"✅ Status final: {final_status}")

        except Exception as e:  # noqa: BLE001
            typer.echo(f"❌ Erro: {e}")
            raise typer.Exit(1)

    # --- Subcomandos para WorkManager ---
    work_app = typer.Typer(help="Gerencia WorkItems em execução")

    @work_app.command("restart")
    def work_restart(work_id: str):
        wm = get_work_manager()
        if wm.restart(work_id):
            typer.echo("🔁 Restarted (queued)")
        else:
            typer.echo("❌ Não foi possível reiniciar (estado inválido?)")
            raise typer.Exit(1)

    app.add_typer(work_app, name="work")

    @app.command()
    def algorithms():  # thin wrapper
        """List available algorithms (delegated)."""
        """Lista algoritmos registrados retornando exit code lógico."""
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
        """Start the web interface (delegated)."""
        try:
            host = host if host is not None else os.getenv("WEB_HOST", "0.0.0.0")
            port = int(port) if port is not None else int(os.getenv("WEB_PORT", "8000"))

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

                # Importa aplicação para registrar rotas
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

            dataset_path = os.getenv("DATASET_PATH", "datasets")
            typer.echo("📁 Using dataset path: %s" % dataset_path)
            orchestrator = DatasetGenerationOrchestrator(base_path=dataset_path)
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
