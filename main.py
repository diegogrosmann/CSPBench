#!/usr/bin/env python3
"""
Ponto de Entrada Principal do CSPBench - Arquitetura Hexagonal

Este é o ponto de entrada principal do CSPBench refatorado com arquitetura hexagonal.
Implementa injeção de dependências e configuração centralizada.

Funcionalidade:
    - Carregamento de configurações do config/settings.yaml
    - Injeção de dependências baseada em configuração
    - Suporte a override via variáveis de ambiente
    - CLI moderna usando Typer

Uso Principal (conforme diretrizes):
    ```bash
    # Menu interativo (sem argumentos)
    python main.py

    # Ajuda geral
    python main.py --help

    # Lista algoritmos disponíveis
    python main.py --algorithms

    # Gera/salva datasets sintéticos
    python main.py --datasetsave

    # Executa batch diretamente (arquivo .yaml/.yml)
    python main.py batches/exemplo.yaml
    python main.py configuracao.yml
    ```

Comandos Específicos:
    ```bash
    # Teste básico do sistema
    python main.py test

    # Executa algoritmo em dataset
    python main.py run <algoritmo> <dataset>

    # Executa batch via comando
    python main.py batch <arquivo.yaml>

    # Gerenciamento de sessões
    python main.py sessions
    python main.py show-session <nome>
    python main.py view-report <nome>
    python main.py cleanup

    # Informações do sistema
    python main.py config-info
    python main.py algorithms
    ```

Variáveis de Ambiente:
    - EXECUTOR_IMPL: Override do executor (Executor)
    - EXPORT_FORMAT: Override do formato de exportação (json, csv, txt)
    - DATASET_PATH: Override do caminho base de datasets
    - NCBI_EMAIL: Email para API do NCBI
    - NCBI_API_KEY: Chave da API do NCBI
"""

import os
import sys
from pathlib import Path
from typing import Any, Dict, Optional

import typer
import yaml

# Adiciona o diretório raiz ao path para imports
sys.path.insert(0, str(Path(__file__).parent))

# IMPORTANTE: Importar algorithms primeiro para carregar o global_registry
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
    help="CSPBench - Framework para Closest String Problem",
    add_completion=False,
    no_args_is_help=False,  # Não mostrar help quando sem argumentos
    rich_markup_mode=None,  # Desabilita rich para evitar erros
)

# Variáveis globais para DI
experiment_service: Optional[ExperimentService] = None
config: Optional[Dict[str, Any]] = None
session_manager: Optional[SessionManager] = None


def load_config() -> Dict[str, Any]:
    """Carrega configuração do arquivo settings.yaml."""
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


def create_dataset_repository(config: Dict[str, Any]) -> FileDatasetRepository:
    """Cria repositório de datasets baseado na configuração."""
    repo_config = config["infrastructure"]["dataset_repository"]["config"]
    base_path = os.getenv("DATASET_PATH", repo_config["base_path"])

    return FileDatasetRepository(base_path)


def create_algorithm_registry(config: Dict[str, Any]) -> DomainAlgorithmRegistry:
    """Cria registry de algoritmos baseado na configuração."""
    return DomainAlgorithmRegistry()


def create_executor(config: Dict[str, Any]) -> Executor:
    """Cria executor baseado na configuração."""
    executor_impl = os.getenv("EXECUTOR_IMPL", "Executor")

    if executor_impl == "Executor":
        return Executor()
    else:
        typer.echo(f"⚠️  Executor '{executor_impl}' não implementado, usando Executor")
        return Executor()


def create_exporter(config: Dict[str, Any]):
    """Cria exportador baseado na configuração."""
    global session_manager

    # Configuração do exportador (se existir na configuração legada)
    exporter_config = (
        config.get("infrastructure", {}).get("exporter", {}).get("config", {})
    )
    export_fmt = os.getenv("EXPORT_FORMAT", exporter_config.get("format", "json"))

    # Usar diretório da sessão se SessionManager existe
    if session_manager:
        output_path = session_manager.get_result_dir()
    else:
        # Usar configuração da nova estrutura
        result_config = config.get("infrastructure", {}).get("result", {})
        output_path = os.getenv(
            "OUTPUT_PATH", result_config.get("base_result_dir", "./outputs/results")
        )

    if export_fmt.lower() == "json":
        return JsonExporter(output_path, config)  # Passar configuração completa
    elif export_fmt.lower() == "csv":
        return CsvExporter(output_path)
    elif export_fmt.lower() == "txt":
        return TxtExporter(output_path)
    else:
        typer.echo(f"⚠️  Formato '{export_fmt}' não suportado, usando JSON")
        return JsonExporter(output_path, config)  # Passar configuração completa


def create_monitoring_service(config: Dict[str, Any]):
    """Cria serviço de monitoramento baseado na configuração."""
    from src.application.services.basic_monitoring_service import BasicMonitoringService

    return BasicMonitoringService(config)


def _initialize_logging(config: Dict[str, Any]) -> None:
    """Inicializa sistema de logging baseado na configuração."""
    global session_manager

    # Configurações padrão
    default_level = "INFO"
    default_log_dir = "outputs/logs"

    try:
        # Criar SessionManager se não existir
        if session_manager is None:
            session_manager = SessionManager(config)

        # Criar nova sessão
        session_folder = session_manager.create_session()

        # Tentar extrair configurações do batch (se existir)
        log_config = {}
        if "advanced" in config and "logs" in config["advanced"]:
            log_config = config["advanced"]["logs"]

        # Verificar se logs estão habilitados
        if not log_config.get("enable", True):
            return

        # Extrair configurações
        log_level = log_config.get("log_level", default_level)

        # Usar caminho da sessão para logs
        log_file_path = session_manager.get_log_path()

        # Inicializar LoggerConfig com caminho específico da sessão
        LoggerConfig.initialize(
            level=log_level,
            log_file_path=log_file_path,
            base_name="cspbench",
        )

        # Criar logger e registrar inicialização
        logger = get_logger(__name__)
        logger.info("=" * 60)
        logger.info("CSPBench - Sistema de logging inicializado")
        logger.info(f"Nível de log: {log_level}")
        logger.info(f"Sessão: {session_folder}")
        logger.info(f"Arquivo de log: {log_file_path}")
        logger.info(f"Timestamp: {__import__('datetime').datetime.now()}")
        logger.info("=" * 60)

    except Exception as e:
        # Fallback para configuração básica em caso de erro
        LoggerConfig.initialize(level=default_level, log_dir=default_log_dir)
        logger = get_logger(__name__)
        logger.warning(f"Erro ao configurar logging, usando padrões: {e}")


def _update_logging_from_batch_config(batch_config: Dict[str, Any]) -> None:
    """Atualiza configuração de logging baseado no arquivo de batch específico."""
    try:
        if "advanced" in batch_config and "logs" in batch_config["advanced"]:
            log_config = batch_config["advanced"]["logs"]

            # Verificar se logs estão habilitados
            if not log_config.get("enable", True):
                return

            # Atualizar nível se especificado
            if "log_level" in log_config:
                new_level = log_config["log_level"]
                LoggerConfig.set_level(new_level)
                logger = get_logger(__name__)
                logger.info(f"Nível de log atualizado para: {new_level}")
    except Exception as e:
        logger = get_logger(__name__)
        logger.warning(f"Erro ao atualizar configuração de logging do batch: {e}")


def initialize_service() -> ExperimentService:
    """Inicializa o serviço de experimentos com injeção de dependências."""
    global config, experiment_service

    if experiment_service is None:
        if config is None:
            config = load_config()

        # Inicializar sistema de logging
        _initialize_logging(config)

        # Cria componentes da infraestrutura
        dataset_repo = create_dataset_repository(config)
        algo_registry = create_algorithm_registry(config)
        executor = create_executor(config)
        exporter = create_exporter(config)
        monitoring_service = create_monitoring_service(config)

        # Cria serviço com DI
        experiment_service = ExperimentService(
            dataset_repo=dataset_repo,
            exporter=exporter,
            executor=executor,
            algo_registry=algo_registry,
            monitoring_service=monitoring_service,
        )

        typer.echo("✅ CSPBench inicializado com sucesso!")

    return experiment_service


def show_interactive_menu() -> None:
    """Mostra interface interativa para seleção de arquivos de batch."""
    print("CSPBench v0.1.0 - Framework para Closest String Problem")
    print("=" * 60)
    print()

    # Listar arquivos de batch disponíveis
    batches_dir = Path("batches")
    if not batches_dir.exists():
        _show_commands_help("📁 Diretório 'batches/' não encontrado")
        return

    batch_files = list(batches_dir.glob("*.yaml")) + list(batches_dir.glob("*.yml"))

    if not batch_files:
        _show_commands_help("📋 Nenhum arquivo de batch encontrado em 'batches/'")
        return

    _display_batch_files(batch_files)
    _display_manual_commands()

    selected_file = _get_user_selection(batch_files)
    if selected_file:
        _execute_selected_batch(selected_file)


def _show_commands_help(message: str) -> None:
    """Mostra ajuda com comandos disponíveis."""
    print(message)
    print()
    print("Comandos disponíveis:")
    print("  test         - Teste básico do sistema")
    print("  algorithms   - Lista algoritmos disponíveis")
    print("  config-info  - Mostra configuração")
    print("  run          - Executa algoritmo único")
    print()
    print("Uso: python main.py <comando> [argumentos]")
    print("Exemplo: python main.py test")


def _display_batch_files(batch_files: list) -> None:
    """Exibe lista de arquivos de batch com descrições."""
    print("📋 Arquivos de batch disponíveis:")
    print()

    for i, batch_file in enumerate(batch_files, 1):
        description = _extract_file_description(batch_file)
        print(f"  {i}. {batch_file.name}")
        print(f"     {description}")
        print()


def _extract_file_description(batch_file: Path) -> str:
    """Extrai descrição do arquivo de batch a partir de comentários."""
    try:
        with open(batch_file, "r", encoding="utf-8") as f:
            first_lines = f.read(300)  # Primeiros 300 chars

        for line in first_lines.split("\n"):
            line = line.strip()
            if line.startswith("#") and len(line) > 1:
                desc_text = line[1:].strip()
                if desc_text and not desc_text.startswith("="):
                    return desc_text

        return "Arquivo de configuração de batch"
    except:
        return "Arquivo de configuração de batch"


def _display_manual_commands() -> None:
    """Exibe lista de comandos manuais disponíveis."""
    print("📋 Comandos manuais disponíveis:")
    print("  test         - Teste básico do sistema")
    print("  algorithms   - Lista algoritmos disponíveis")
    print("  config-info  - Mostra configuração")
    print("  run          - Executa algoritmo único")
    print()


def _get_user_selection(batch_files: list) -> Optional[Path]:
    """Obtém seleção do usuário e valida."""
    try:
        choice = input(
            "💡 Selecione um arquivo (número) ou pressione Enter para sair: "
        ).strip()

        if choice == "":
            print("👋 Até logo!")
            return None

        if choice.isdigit():
            choice_num = int(choice)
            if 1 <= choice_num <= len(batch_files):
                return batch_files[choice_num - 1]
            else:
                print("❌ Número inválido!")
                return None
        else:
            print("❌ Entrada inválida!")
            return None

    except KeyboardInterrupt:
        print("\n👋 Até logo!")
        return None
    except EOFError:
        print("\n👋 Até logo!")
        return None


def _execute_selected_batch(selected_file: Path) -> None:
    """Executa o arquivo de batch selecionado."""
    print(f"\n🚀 Executando: {selected_file.name}")
    print("-" * 40)

    # Executa o batch selecionado
    import sys

    sys.argv = ["main.py", "batch", str(selected_file)]
    app()


def show_algorithms_and_exit():
    """Mostra algoritmos disponíveis e sai."""
    try:
        # Importa do módulo algorithms para ativar auto-descoberta
        import algorithms
        from algorithms import global_registry

        print("🧠 Algoritmos disponíveis:")
        for name, cls in global_registry.items():
            print(f"  • {name}: {cls.__doc__ or 'Sem descrição'}")

        if not global_registry:
            print("  (Nenhum algoritmo registrado)")

    except Exception as e:
        print(f"❌ Erro: {e}")
        sys.exit(1)

    sys.exit(0)


def show_datasetsave_and_exit():
    """Executa wizard interativo de geração de datasets e sai."""
    try:
        from src.infrastructure.orchestrators.dataset_generation_orchestrator import (
            DatasetGenerationOrchestrator,
        )

        # Criar e executar orquestrador
        orchestrator = DatasetGenerationOrchestrator()
        result_path = orchestrator.run_interactive_generation()

        if result_path:
            print(f"\n🎉 Dataset salvo com sucesso!")
        else:
            print("\n� Operação cancelada.")

    except Exception as e:
        print(f"❌ Erro: {e}")
        sys.exit(1)

    sys.exit(0)


def execute_batch_file(batch_file: str):
    """Executa arquivo de batch diretamente."""
    try:
        batch_path = Path(batch_file)

        if not batch_path.exists():
            print(f"❌ Arquivo não encontrado: {batch_file}")
            raise typer.Exit(1)

        if not batch_path.suffix.lower() in [".yaml", ".yml"]:
            print(f"❌ Arquivo deve ser .yaml ou .yml: {batch_file}")
            raise typer.Exit(1)

        # Executar batch
        service = initialize_service()
        print(f"📋 Executando batch: {batch_file}...")

        result = service.run_batch(str(batch_path))
        print(f"✅ Batch concluído: {result.get('summary', 'Concluído')}")

    except Exception as e:
        print(f"❌ Erro no batch: {e}")
        raise typer.Exit(1)


# Registrar todos os comandos da CLI
register_commands(app, initialize_service)


def main(args: Optional[list] = None):
    """Função main para execução programática."""
    import sys

    if args is None:
        args = sys.argv[1:]

    # Se for um arquivo de batch, executar diretamente
    if len(args) == 1 and args[0].endswith((".yaml", ".yml")):
        execute_batch_file(args[0])
        return

    # Caso contrário, usar o sistema de comandos do Typer
    original_argv = sys.argv[:]
    try:
        sys.argv = ["main.py"] + args
        app()
    finally:
        sys.argv = original_argv


if __name__ == "__main__":
    import sys

    # Se executado sem argumentos, mostra interface interativa
    if len(sys.argv) == 1:
        show_interactive_menu()
    elif len(sys.argv) == 2:
        arg = sys.argv[1]

        # Implementar flags especiais
        if arg == "--algorithms":
            show_algorithms_and_exit()
        elif arg == "--datasetsave":
            show_datasetsave_and_exit()
        elif arg == "--help":
            print(
                """
CSPBench v0.1.0 - Framework para Closest String Problem

Uso:
    python main.py                    Menu interativo
    python main.py --help            Esta ajuda
    python main.py --algorithms      Lista algoritmos disponíveis
    python main.py --datasetsave     Gera/salva datasets
    python main.py <arquivo.yaml>    Executa batch
    python main.py <comando>         Executa comando específico

Comandos disponíveis:
    test                             Teste básico do sistema
    run <algoritmo> <dataset>        Executa algoritmo em dataset
    batch <arquivo.yaml>             Executa batch
    algorithms                       Lista algoritmos
    config-info                      Mostra configuração
    sessions                         Lista sessões
    cleanup                          Remove sessões antigas
    show-session <nome>              Mostra detalhes da sessão
    view-report <nome>               Abre relatório no navegador

Exemplos:
    python main.py test
    python main.py run Baseline teste.fasta
    python main.py batch batches/exemplo.yaml
    python main.py batches/exemplo.yaml
"""
            )
            sys.exit(0)
        elif not arg.startswith("-") and (
            arg.endswith(".yaml") or arg.endswith(".yml")
        ):
            # É um arquivo de batch - executar diretamente
            execute_batch_file(arg)
            sys.exit(0)
        else:
            # Usar CLI normal do Typer
            app()
    else:
        # Usar CLI normal do Typer
        app()
