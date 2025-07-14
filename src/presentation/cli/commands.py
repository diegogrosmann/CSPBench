"""
Módulo de Registro de Comandos CLI

Centraliza o registro de todos os comandos da CLI para modularidade.
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
    """Carrega configuração do arquivo settings.yaml."""
    import yaml

    config_path = Path("config/settings.yaml")

    if not config_path.exists():
        typer.echo(f"❌ Arquivo de configuração não encontrado: {config_path}")
        raise typer.Exit(1)

    try:
        with open(config_path, "r", encoding="utf-8") as f:
            return yaml.safe_load(f)
    except Exception as e:
        typer.echo(f"❌ Erro ao carregar configuração: {e}")
        raise typer.Exit(1)


def list_sessions() -> None:
    """Lista todas as sessões disponíveis."""
    try:
        config = load_config()
        session_mgr = SessionManager(config)
        sessions = session_mgr.list_sessions()

        if not sessions:
            typer.echo("📂 Nenhuma sessão encontrada.")
            return

        typer.echo("📂 Sessões disponíveis:")
        typer.echo("-" * 60)

        # Ordenar por data de criação (mais recentes primeiro)
        sorted_sessions = sorted(
            sessions.items(),
            key=lambda x: x[1]["created"],
            reverse=True,
        )

        for session_name, info in sorted_sessions:
            created = info["created"].strftime("%Y-%m-%d %H:%M:%S")
            logs_status = "✅" if info["logs"] else "❌"
            results_status = "✅" if info["results"] else "❌"

            typer.echo(f"🗂️  {session_name}")
            typer.echo(f"   📅 Criado: {created}")
            typer.echo(f"   📄 Logs: {logs_status}  🗃️  Resultados: {results_status}")
            typer.echo()

    except Exception as e:
        typer.echo(f"❌ Erro ao listar sessões: {e}")


def cleanup_old_sessions(keep_last: int = 10) -> None:
    """Remove sessões antigas, mantendo apenas as mais recentes."""
    try:
        config = load_config()
        session_mgr = SessionManager(config)
        session_mgr.cleanup_old_sessions(keep_last)
        typer.echo(
            f"🧹 Limpeza concluída. Mantidas as {keep_last} sessões mais recentes."
        )
    except Exception as e:
        typer.echo(f"❌ Erro na limpeza: {e}")


def register_commands(app: typer.Typer, experiment_service_getter) -> None:
    """
    Registra todos os comandos da CLI na aplicação Typer.

    Args:
        app: Instância do Typer onde registrar comandos
        experiment_service_getter: Função que retorna o ExperimentService inicializado
    """

    @app.command()
    def test():
        """Teste básico do sistema."""
        try:
            service = experiment_service_getter()
            assert service is not None

            # Cria dataset sintético para teste
            generator = SyntheticDatasetGenerator()
            dataset = generator.generate_random(n=10, length=20, alphabet="ACTG")

            typer.echo(
                f"📊 Dataset gerado: {len(dataset.sequences)} strings de tamanho {len(dataset.sequences[0])}"
            )

            # Testa algoritmo disponível
            import algorithms
            from algorithms import global_registry

            available_algorithms = list(global_registry.keys())
            if not available_algorithms:
                typer.echo("⚠️ Nenhum algoritmo disponível", color=True)
                return

            # Usa o primeiro algoritmo disponível
            algorithm_name = available_algorithms[0]
            algorithm_class = global_registry[algorithm_name]

            algorithm = algorithm_class(
                strings=dataset.sequences, alphabet=dataset.alphabet
            )

            result_string, max_distance, metadata = algorithm.run()

            typer.echo(f"🎯 Resultado: {result_string}")
            typer.echo(f"📏 Distância máxima: {max_distance}")
            typer.echo(f"📋 Metadados: {metadata}")
            typer.echo("✅ Teste concluído com sucesso!")

        except Exception as e:
            typer.echo(f"❌ Erro no teste: {e}")
            raise typer.Exit(1)

    @app.command()
    def run(
        algorithm: str = typer.Argument(..., help="Nome do algoritmo"),
        dataset: str = typer.Argument(..., help="Caminho do dataset"),
        params: Optional[str] = typer.Option(
            None, "--params", "-p", help="JSON com parâmetros"
        ),
        timeout: Optional[int] = typer.Option(
            None, "--timeout", "-t", help="Timeout em segundos"
        ),
        output: Optional[str] = typer.Option(
            None, "--output", "-o", help="Arquivo de saída"
        ),
    ):
        """Executa um algoritmo em um dataset."""
        try:
            service = experiment_service_getter()
            assert service is not None

            # Parse dos parâmetros JSON se fornecidos
            params_dict = json.loads(params) if params else {}

            typer.echo(f"🚀 Executando {algorithm} em {dataset}...")

            result = service.run_single_experiment(
                algorithm, dataset, params=params_dict, timeout=timeout
            )

            typer.echo(f"🎯 Resultado: {result}")

        except Exception as e:
            typer.echo(f"❌ Erro: {e}")
            raise typer.Exit(1)

    @app.command()
    def batch(
        cfg: Path = typer.Argument(
            ..., exists=True, readable=True, help="YAML do batch"
        ),
        verbose: bool = typer.Option(
            False, "--verbose", "-v", help="Mostrar detalhes dos resultados"
        ),
    ):
        """Executa um arquivo de batch (runs, otimizações ou sensibilidade)."""
        try:
            service = experiment_service_getter()
            assert service is not None

            typer.echo(f"📋 Executando batch: {cfg}...")

            result = service.run_batch(str(cfg))

            if verbose:
                typer.echo(f"📊 Resultados detalhados:")
                for i, res in enumerate(
                    result["results"][:5]
                ):  # Primeiros 5 resultados
                    typer.echo(f"  Resultado {i+1}: {res}")
                if len(result["results"]) > 5:
                    typer.echo(f"  ... e mais {len(result['results']) - 5} resultados")

            typer.echo(f"✅ Batch concluído: {result['summary']}")

        except Exception as e:
            typer.echo(f"❌ Erro no batch: {e}")
            raise typer.Exit(1)

    @app.command()
    def algorithms():
        """Lista algoritmos disponíveis."""
        try:
            service = experiment_service_getter()  # Para garantir inicialização

            # Importa do módulo algorithms para ativar auto-descoberta
            import algorithms
            from algorithms import global_registry

            typer.echo("🧠 Algoritmos disponíveis:")
            for name, cls in global_registry.items():
                typer.echo(f"  • {name}: {cls.__doc__ or 'Sem descrição'}")

            if not global_registry:
                typer.echo("  (Nenhum algoritmo registrado)")

        except Exception as e:
            typer.echo(f"❌ Erro: {e}")
            raise typer.Exit(1)

    @app.command()
    def config_info():
        """Mostra informações de configuração."""
        try:
            import yaml

            config_path = Path("config/settings.yaml")
            if not config_path.exists():
                typer.echo(f"❌ Arquivo de configuração não encontrado: {config_path}")
                raise typer.Exit(1)

            with open(config_path, "r", encoding="utf-8") as f:
                config = yaml.safe_load(f)

            app_info = config["application"]
            typer.echo(f"📋 {app_info['name']} v{app_info['version']}")
            typer.echo(f"📝 {app_info['description']}")

            typer.echo("\\n🔧 Configuração de infraestrutura:")
            for component, info in config["infrastructure"].items():
                typer.echo(f"  • {component}: {info['type']}")

            typer.echo("\\n🌍 Variáveis de ambiente:")
            typer.echo(
                f"  • NCBI_EMAIL: {'definido' if os.getenv('NCBI_EMAIL') else 'não definido'}"
            )
            typer.echo(
                f"  • NCBI_API_KEY: {'definido' if os.getenv('NCBI_API_KEY') else 'não definido'}"
            )
            typer.echo(
                f"  • EXECUTOR_IMPL: {os.getenv('EXECUTOR_IMPL', 'não definida')}"
            )
            typer.echo(f"  • EXPORT_FMT: {os.getenv('EXPORT_FMT', 'não definida')}")
            typer.echo(f"  • DATASET_PATH: {os.getenv('DATASET_PATH', 'não definida')}")

        except Exception as e:
            typer.echo(f"❌ Erro: {e}")
            raise typer.Exit(1)

    @app.command()
    def sessions() -> None:
        """
        Lista todas as sessões disponíveis com suas informações.
        """
        list_sessions()

    @app.command()
    def cleanup(
        keep: int = typer.Option(
            10, "--keep", "-k", help="Número de sessões mais recentes para manter"
        )
    ) -> None:
        """
        Remove sessões antigas, mantendo apenas as mais recentes.
        """
        cleanup_old_sessions(keep)

    @app.command()
    def show_session(
        session_name: str = typer.Argument(
            ..., help="Nome da sessão (formato: YYYYMMDD_HHMMSS)"
        )
    ) -> None:
        """
        Mostra detalhes de uma sessão específica.
        """
        try:
            config = load_config()
            session_mgr = SessionManager(config)
            sessions = session_mgr.list_sessions()

            if session_name not in sessions:
                typer.echo(f"❌ Sessão '{session_name}' não encontrada.")
                typer.echo("\n📂 Sessões disponíveis:")
                for name in sorted(sessions.keys(), reverse=True):
                    typer.echo(f"  • {name}")
                return

            info = sessions[session_name]
            created = info["created"].strftime("%Y-%m-%d %H:%M:%S")

            typer.echo(f"🗂️  Sessão: {session_name}")
            typer.echo(f"📅 Criado: {created}")
            typer.echo()

            # Mostrar logs se existirem
            if info["logs"]:
                log_path = session_mgr.get_log_path(session_name)
                typer.echo(f"📄 Log: {log_path}")
                if log_path.exists():
                    stat = log_path.stat()
                    size_kb = stat.st_size / 1024
                    typer.echo(f"   📊 Tamanho: {size_kb:.1f} KB")

            # Mostrar resultados se existirem
            if info["results"]:
                result_path = session_mgr.get_result_path(session_name)
                typer.echo(f"🗃️  Resultado: {result_path}")
                if result_path.exists():
                    stat = result_path.stat()
                    size_kb = stat.st_size / 1024
                    typer.echo(f"   📊 Tamanho: {size_kb:.1f} KB")

                    # Tentar mostrar resumo do resultado
                    try:
                        import json

                        with open(result_path, "r") as f:
                            result_data = json.load(f)

                        if "summary" in result_data:
                            summary = result_data["summary"]
                            typer.echo(f"   📈 Resumo: {summary}")

                    except Exception:
                        pass  # Ignora erros ao ler resultado

        except Exception as e:
            typer.echo(f"❌ Erro ao mostrar sessão: {e}")

    @app.command()
    def view_report(
        session_name: str = typer.Argument(
            ..., help="Nome da sessão (formato: YYYYMMDD_HHMMSS)"
        )
    ) -> None:
        """
        Abre o relatório HTML de uma sessão no navegador.
        """
        try:
            config = load_config()
            session_mgr = SessionManager(config)
            sessions = session_mgr.list_sessions()

            if session_name not in sessions:
                typer.echo(f"❌ Sessão '{session_name}' não encontrada.")
                return

            # Construir caminho do relatório
            result_base_dir = Path(
                config["infrastructure"]["result"]["base_result_dir"]
            )
            report_path = result_base_dir / session_name / "report" / "report.html"

            if not report_path.exists():
                typer.echo(f"❌ Relatório não encontrado para sessão '{session_name}'.")
                typer.echo(f"   Esperado em: {report_path}")
                return

            # Abrir no navegador
            import webbrowser

            webbrowser.open(f"file://{report_path.absolute()}")
            typer.echo(f"🌐 Abrindo relatório no navegador: {report_path}")

        except Exception as e:
            typer.echo(f"❌ Erro ao abrir relatório: {e}")
