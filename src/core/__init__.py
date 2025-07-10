"""
Módulo Core do CSPBench

Este módulo constitui o núcleo central do CSPBench, fornecendo a infraestrutura
fundamental para execução de algoritmos, gerenciamento de recursos, processamento
de dados e geração de relatórios.

Submódulos:
    - interfaces: Contratos e abstrações do sistema
    - scheduler: Sistema de agendamento e execução paralela
    - io: Entrada/saída e processamento de dados
    - report: Geração e formatação de relatórios
    - config: Gerenciamento de configurações

Arquitetura:
    O módulo implementa uma arquitetura em camadas com:
    - Camada de Interface: Contratos padronizados
    - Camada de Execução: Schedulers e executores
    - Camada de Dados: I/O e persistência
    - Camada de Relatórios: Análise e visualização

Funcionalidades Principais:
    - Execução paralela de algoritmos CSP
    - Gerenciamento inteligente de recursos
    - Sistema de filas e agendamento
    - Processamento robusto de dados
    - Geração automática de relatórios
    - Monitoramento em tempo real

Padrões de Design:
    - Factory Pattern: Criação de objetos especializados
    - Strategy Pattern: Algoritmos intercambiáveis
    - Observer Pattern: Monitoramento de progresso
    - Command Pattern: Encapsulamento de operações

Exemplo de Uso:
    ```python
    from src.core.scheduler.scheduler import ExecutionScheduler
    from src.core.interfaces.factory import create_executor

    # Criar sistema de execução
    scheduler = ExecutionScheduler()
    executor = create_executor('parallel', num_workers=4)

    # Configurar e executar experimentos
    results = scheduler.execute_batch(tasks, executor)
    ```

Integração:
    - Algoritmos: Via interfaces padronizadas
    - Datasets: Sistema unificado de carregamento
    - Configuração: Parâmetros centralizados
    - Logging: Sistema estruturado de logs
    - UI: Interfaces múltiplas (CLI, curses)

Autor: CSPBench Development Team
Data: 2024
"""
