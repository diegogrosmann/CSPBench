"""
Extensão da interface CLI para incluir modo curses.

Funções:
    add_curses_batch_option: Adiciona opção de execução curses ao menu
    run_batch_with_curses: Executa batch com interface curses
"""

import logging
from pathlib import Path

from src.core.exec.batch_executor import BatchExecutor
from src.ui.cli.console_manager import console
from src.ui.curses_interface import run_batch_with_curses
from src.utils.config import safe_input

logger = logging.getLogger(__name__)


def add_curses_batch_option_to_menu():
    """Adiciona opção de execução curses ao menu principal."""

    def execute_batch_curses():
        """Executa batch com interface curses."""
        console.print("\n🖥️  Execução em Lote com Interface Curses")
        console.print("=" * 50)

        # Verificar se há arquivos de configuração
        batch_dir = Path("batch_configs")
        if not batch_dir.exists():
            console.print("❌ Diretório batch_configs não encontrado!")
            console.print("💡 Crie arquivos de configuração em batch_configs/")
            return

        # Listar arquivos de configuração
        config_files = []
        for pattern in ["*.yaml", "*.yml", "*.json", "*.xml"]:
            config_files.extend(batch_dir.glob(pattern))

        if not config_files:
            console.print("❌ Nenhum arquivo de configuração encontrado!")
            console.print("💡 Crie arquivos .yaml, .json ou .xml em batch_configs/")
            return

        # Mostrar opções
        console.print("\n📁 Arquivos de configuração disponíveis:")
        for i, config_file in enumerate(config_files, 1):
            console.print(f"  {i}. {config_file.name}")

        # Solicitar seleção
        try:
            choice = safe_input("\n🔢 Selecione o arquivo (número): ")
            choice_idx = int(choice) - 1

            if 0 <= choice_idx < len(config_files):
                config_path = config_files[choice_idx]
                console.print(f"\n🗂️  Configuração selecionada: {config_path.name}")

                # Executar com interface curses
                run_curses_batch_execution(config_path)

            else:
                console.print("❌ Seleção inválida!")

        except (ValueError, KeyboardInterrupt):
            console.print("❌ Operação cancelada!")

    return execute_batch_curses


def run_curses_batch_execution(config_path: Path):
    """Executa batch com interface curses."""

    try:
        # Criar executor com arquivo de configuração
        executor = BatchExecutor(str(config_path))

        # As configurações já foram carregadas no construtor
        execucoes = executor.execucoes

        if not execucoes:
            console.print("❌ Nenhuma configuração válida encontrada!")
            return

        console.print(f"📊 Carregadas {len(execucoes)} configurações")

        # Perguntar se deseja usar interface curses
        use_curses = safe_input("\n🖥️  Usar interface curses? (s/N): ").lower().strip()

        if use_curses in ["s", "sim", "y", "yes"]:
            console.print("\n🚀 Iniciando execução com interface curses...")
            console.print("💡 Pressione 'q' para sair da interface")

            # Pequeno delay para mostrar mensagem
            import time

            time.sleep(1)

            # Executar com interface curses
            run_batch_with_curses(executor, execucoes)

            console.print("\n✅ Execução concluída!")

        else:
            console.print("\n🔄 Executando em modo tradicional...")

            # Executar em modo tradicional
            results = executor.execute_batch()

            console.print("\n✅ Execução concluída!")

    except Exception as e:
        console.print(f"❌ Erro na execução: {e}")
        logger.error(f"Erro na execução curses: {e}")


def show_curses_info():
    """Mostra informações sobre a interface curses."""

    console.print("\n🖥️  Interface Curses - Informações")
    console.print("=" * 40)
    console.print("📌 A interface curses permite acompanhar:")
    console.print("  • Progresso em tempo real dos algoritmos")
    console.print("  • Informações do dataset atual")
    console.print("  • Status de execução paralela")
    console.print("  • Melhor distância encontrada")
    console.print("  • Tempo decorrido por algoritmo")
    console.print("")
    console.print("⌨️  Controles:")
    console.print("  • q/Q: Sair da interface")
    console.print("  • r/R: Forçar atualização")
    console.print("")
    console.print("⚠️  Requisitos:")
    console.print("  • Terminal com suporte a cores")
    console.print("  • Tamanho mínimo: 80x24")
    console.print("  • Biblioteca curses disponível")


def check_curses_support() -> bool:
    """Verifica se o sistema suporta curses."""

    try:
        import curses

        # Testar se curses funciona
        curses.setupterm()
        return True

    except ImportError:
        console.print("❌ Biblioteca curses não disponível!")
        return False

    except Exception as e:
        console.print(f"❌ Erro ao verificar suporte curses: {e}")
        return False


def create_curses_menu_options():
    """Cria opções do menu para interface curses."""

    options = {}

    if check_curses_support():
        options["curses_batch"] = {
            "label": "🖥️  Executar Batch com Interface Curses",
            "function": add_curses_batch_option_to_menu(),
            "description": "Executa configurações em lote com interface visual",
        }

        options["curses_info"] = {
            "label": "ℹ️  Informações sobre Interface Curses",
            "function": show_curses_info,
            "description": "Mostra informações sobre a interface curses",
        }

    else:
        options["curses_unavailable"] = {
            "label": "❌ Interface Curses Indisponível",
            "function": lambda: console.print("❌ Curses não suportado neste sistema"),
            "description": "Interface curses não disponível",
        }

    return options


def integrate_curses_with_app(app_instance):
    """Integra opções curses com a aplicação principal."""

    # Adicionar opções curses ao menu
    curses_options = create_curses_menu_options()

    # Se o app tem método para adicionar opções, usar
    if hasattr(app_instance, "add_menu_options"):
        app_instance.add_menu_options(curses_options)

    # Caso contrário, modificar diretamente o menu
    elif hasattr(app_instance, "menu_options"):
        app_instance.menu_options.update(curses_options)

    return curses_options
