"""
Módulo de Interface de Usuário (UI) - CSPBench

Este módulo fornece um sistema completo de interfaces de usuário para o CSPBench,
oferecendo múltiplas formas de interação com o framework através de CLI interativa,
interfaces visuais e sistemas de automação.

Arquitetura:
    O módulo implementa uma arquitetura modular com:
    - CLI interativa com menus guiados
    - Interface curses para monitoramento visual
    - Sistema de automação para execução silenciosa
    - Wizards especializados para configuração
    - Gerenciamento de console unificado

Interfaces Disponíveis:
    - **CLI Interativa**: Menus amigáveis para usuários
    - **CLI Silenciosa**: Automação e scripts
    - **Interface Curses**: Monitoramento visual em tempo real
    - **Execução em Lote**: Configurações YAML
    - **Wizards**: Assistentes para configuração complexa

Submódulos:
    - cli: Interface de linha de comando principal
    - curses_interface: Interface visual para monitoramento
    - curses_integration: Integração com sistema de execução

Funcionalidades:
    - Navegação intuitiva por menus
    - Validação automática de entrada
    - Feedback visual de progresso
    - Tratamento robusto de erros
    - Suporte a múltiplos modos de execução
    - Configuração flexível via parâmetros

Exemplo de Uso:
    ```python
    from src.ui.cli.app import main
    from src.ui.curses_interface import CursesInterface

    # Execução da aplicação principal
    main()

    # Usar interface curses diretamente
    interface = CursesInterface()
    interface.run_monitoring_session()
    ```

Modos de Operação:
    - **Interativo**: Usuário navega por menus
    - **Silencioso**: Execução automatizada
    - **Visual**: Monitoramento em tempo real
    - **Batch**: Configuração via arquivos YAML
    - **Híbrido**: Combinação de modos

Design de UX:
    - Interface intuitiva e amigável
    - Feedback claro e informativo
    - Tratamento gracioso de erros
    - Progressão lógica de tarefas
    - Documentação contextual

Integração:
    - Sistema de execução do core
    - Algoritmos via registry
    - Datasets através de módulos especializados
    - Configuração centralizada
    - Logging estruturado

Autor: CSPBench Development Team
Data: 2024
"""
