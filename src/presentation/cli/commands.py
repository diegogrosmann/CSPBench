"""
CLI Commands Registration Module

Centralizes the registration of all CLI commands for modularity.
"""

import json
import os
from pathlib import Path
from typing import Optional

import typer

from src.application.services.experiment_service import ExperimentService
from src.domain import SyntheticDatasetGenerator
from src.infrastructure.orchestrators.session_manager import SessionManager


def load_config():
    """Load configuration from settings.yaml file."""
    import yaml

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


def list_sessions() -> None:
    """List all available sessions."""
    try:
        config = load_config()
        session_mgr = SessionManager(config)
        sessions = session_mgr.list_sessions()

        if not sessions:
            typer.echo("📂 No sessions found.")
            return

        typer.echo("📂 Available sessions:")
        typer.echo("-" * 60)

        # Sort by creation date (most recent first)
        sorted_sessions = sorted(
            sessions.items(),
            key=lambda x: x[1]["created"],
            reverse=True,
        )

        for session_name, info in sorted_sessions:
            created_str = info["created"].strftime("%Y-%m-%d %H:%M:%S")
            logs_status = "✅" if info["logs"] else "❌"
            results_status = "✅" if info["results"] else "❌"

            typer.echo(f"  {session_name}")
            typer.echo(f"    Created: {created_str}")
            typer.echo(f"    Logs: {logs_status}  Results: {results_status}")
            typer.echo()

    except Exception as e:
        typer.echo(f"❌ Error listing sessions: {e}")


def cleanup_old_sessions(keep_last: int = 10) -> None:
    """Remove old sessions, keeping only the most recent ones."""
    try:
        config = load_config()
        session_mgr = SessionManager(config)
        session_mgr.cleanup_old_sessions(keep_last)
        typer.echo(
            f"🧹 Cleanup completed. Kept the {keep_last} most recent sessions."
        )
    except Exception as e:
        typer.echo(f"❌ Cleanup error: {e}")


def register_commands(app: typer.Typer, experiment_service_getter) -> None:
    """
    Register all CLI commands in the Typer application.

    Args:
        app: Typer instance where to register commands
        experiment_service_getter: Function that returns initialized ExperimentService
    """

    @app.command()
    def test():
        """Basic system test."""
        try:
            service = experiment_service_getter()
            assert service is not None

            # Create synthetic dataset for testing
            generator = SyntheticDatasetGenerator()
            dataset = generator.generate_random(n=10, length=20, alphabet="ACTG")

            typer.echo(
                f"📊 Dataset generated: {len(dataset.sequences)} strings of size {len(dataset.sequences[0])}"
            )

            # Test available algorithm
            import algorithms
            from algorithms import global_registry

            available_algorithms = list(global_registry.keys())
            if not available_algorithms:
                typer.echo("⚠️ No algorithms available", color=True)
                return

            # Use the first available algorithm
            algorithm_name = available_algorithms[0]
            algorithm_class = global_registry[algorithm_name]

            algorithm = algorithm_class(
                strings=dataset.sequences, alphabet=dataset.alphabet
            )

            result_string, max_distance, metadata = algorithm.run()

            typer.echo(f"🎯 Result: {result_string}")
            typer.echo(f"📏 Maximum distance: {max_distance}")
            typer.echo(f"📋 Metadata: {metadata}")
            typer.echo("✅ Test completed successfully!")

        except Exception as e:
            typer.echo(f"❌ Test error: {e}")
            raise typer.Exit(1)

    @app.command()
    def run(
        algorithm: str = typer.Argument(..., help="Algorithm name"),
        dataset: str = typer.Argument(..., help="Dataset path"),
        params: Optional[str] = typer.Option(
            None, "--params", "-p", help="JSON with parameters"
        ),
        timeout: Optional[int] = typer.Option(
            None, "--timeout", "-t", help="Timeout in seconds"
        ),
        output: Optional[str] = typer.Option(
            None, "--output", "-o", help="Output file"
        ),
    ):
        """Execute an algorithm on a dataset."""
        try:
            service = experiment_service_getter()
            assert service is not None

            # Parse JSON parameters if provided
            params_dict = json.loads(params) if params else {}

            typer.echo(f"🚀 Executing {algorithm} on {dataset}...")

            result = service.run_single_experiment(
                algorithm, dataset, params=params_dict, timeout=timeout
            )

            typer.echo(f"🎯 Result: {result}")

        except Exception as e:
            typer.echo(f"❌ Error: {e}")
            raise typer.Exit(1)

    @app.command()
    def batch(
        cfg: Path = typer.Argument(
            ..., exists=True, readable=True, help="Batch YAML file"
        ),
        verbose: bool = typer.Option(
            False, "--verbose", "-v", help="Show result details"
        ),
    ):
        """Execute a batch file (runs, optimizations or sensitivity analysis)."""
        try:
            service = experiment_service_getter()
            assert service is not None

            typer.echo(f"📋 Executing batch: {cfg}...")

            result = service.run_batch(str(cfg))

            if verbose:
                typer.echo(f"📊 Detailed results:")
                for i, res in enumerate(
                    result["results"][:5]
                ):  # First 5 results
                    typer.echo(f"  Result {i+1}: {res}")
                if len(result["results"]) > 5:
                    typer.echo(f"  ... and {len(result['results']) - 5} more results")

            typer.echo(f"✅ Batch completed: {result['summary']}")

        except Exception as e:
            typer.echo(f"❌ Batch error: {e}")
            raise typer.Exit(1)

    @app.command()
    def algorithms():
        """List available algorithms."""
        try:
            service = experiment_service_getter()  # To ensure initialization

            # Import from algorithms module to activate auto-discovery
            import algorithms
            from algorithms import global_registry

            typer.echo("🧠 Available algorithms:")
            for name, cls in global_registry.items():
                typer.echo(f"  • {name}: {cls.__doc__ or 'No description'}")

            if not global_registry:
                typer.echo("  (No algorithms registered)")

        except Exception as e:
            typer.echo(f"❌ Error: {e}")
            raise typer.Exit(1)

    @app.command()
    def config_info():
        """Show configuration information."""
        try:
            import yaml

            config_path = Path("config/settings.yaml")
            if not config_path.exists():
                typer.echo(f"❌ Configuration file not found: {config_path}")
                raise typer.Exit(1)

            with open(config_path, "r", encoding="utf-8") as f:
                config = yaml.safe_load(f)

            app_info = config["application"]
            typer.echo(f"📋 {app_info['name']} v{app_info['version']}")
            typer.echo(f"📝 {app_info['description']}")

            typer.echo("\\n🔧 Infrastructure configuration:")
            for component, info in config["infrastructure"].items():
                typer.echo(f"  • {component}: {info['type']}")

            typer.echo("\\n🌍 Environment variables:")
            typer.echo(
                f"  • NCBI_EMAIL: {'defined' if os.getenv('NCBI_EMAIL') else 'not defined'}"
            )
            typer.echo(
                f"  • NCBI_API_KEY: {'defined' if os.getenv('NCBI_API_KEY') else 'not defined'}"
            )
            typer.echo(
                f"  • EXECUTOR_IMPL: {os.getenv('EXECUTOR_IMPL', 'not defined')}"
            )
            typer.echo(f"  • EXPORT_FMT: {os.getenv('EXPORT_FMT', 'not defined')}")
            typer.echo(f"  • DATASET_PATH: {os.getenv('DATASET_PATH', 'not defined')}")

        except Exception as e:
            typer.echo(f"❌ Error: {e}")
            raise typer.Exit(1)

    @app.command()
    def sessions() -> None:
        """
        List all available sessions with their information.
        """
        list_sessions()

    @app.command()
    def cleanup(
        keep: int = typer.Option(
            10, "--keep", "-k", help="Number of most recent sessions to keep"
        )
    ) -> None:
        """
        Remove old sessions, keeping only the most recent ones.
        """
        cleanup_old_sessions(keep)

    @app.command()
    def show_session(
        session_name: str = typer.Argument(
            ..., help="Session name (format: YYYYMMDD_HHMMSS)"
        )
    ) -> None:
        """
        Show details of a specific session.
        """
        try:
            config = load_config()
            session_mgr = SessionManager(config)
            sessions = session_mgr.list_sessions()

            if session_name not in sessions:
                typer.echo(f"❌ Session '{session_name}' not found.")
                typer.echo("\n📂 Available sessions:")
                for name in sorted(sessions.keys(), reverse=True):
                    typer.echo(f"  • {name}")
                return

            info = sessions[session_name]
            created = info["created"].strftime("%Y-%m-%d %H:%M:%S")

            typer.echo(f"🗂️  Session: {session_name}")
            typer.echo(f"📅 Created: {created}")
            typer.echo()

            # Show logs if they exist
            if info["logs"]:
                log_path = session_mgr.get_log_path(session_name)
                typer.echo(f"📄 Log: {log_path}")
                if log_path.exists():
                    stat = log_path.stat()
                    size_kb = stat.st_size / 1024
                    typer.echo(f"   📊 Size: {size_kb:.1f} KB")

            # Show results if they exist
            if info["results"]:
                result_path = session_mgr.get_result_path(session_name)
                typer.echo(f"🗃️  Result: {result_path}")
                if result_path.exists():
                    stat = result_path.stat()
                    size_kb = stat.st_size / 1024
                    typer.echo(f"   📊 Size: {size_kb:.1f} KB")

                    # Try to show result summary
                    try:
                        import json

                        with open(result_path, "r") as f:
                            result_data = json.load(f)

                        if "summary" in result_data:
                            summary = result_data["summary"]
                            typer.echo(f"   📈 Summary: {summary}")

                    except Exception:
                        pass  # Ignore errors when reading result

        except Exception as e:
            typer.echo(f"❌ Error showing session: {e}")

    @app.command()
    def view_report(
        session_name: str = typer.Argument(
            ..., help="Session name (format: YYYYMMDD_HHMMSS)"
        )
    ) -> None:
        """
        Open the HTML report of a session in the browser.
        """
        try:
            config = load_config()
            session_mgr = SessionManager(config)
            sessions = session_mgr.list_sessions()

            if session_name not in sessions:
                typer.echo(f"❌ Session '{session_name}' not found.")
                return

            # Build report path
            result_base_dir = Path(
                config["infrastructure"]["result"]["base_result_dir"]
            )
            report_path = result_base_dir / session_name / "report" / "report.html"

            if not report_path.exists():
                typer.echo(f"❌ Report not found for session '{session_name}'.")
                typer.echo(f"   Expected at: {report_path}")
                return

            # Open in browser
            import webbrowser

            webbrowser.open(f"file://{report_path.absolute()}")
            typer.echo(f"🌐 Opening report in browser: {report_path}")

        except Exception as e:
            typer.echo(f"❌ Error opening report: {e}")
