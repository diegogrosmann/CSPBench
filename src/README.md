# Arquitetura do Código-Fonte (`src`)

Este diretório contém o núcleo da lógica de aplicação do projeto CSP-BLFGA. A arquitetura foi projetada para ser modular, extensível e robusta, separando claramente as responsabilidades entre os diferentes componentes.

## Estrutura de Diretórios

-   `core/`: Contém a lógica central de execução, agendamento e gerenciamento de tarefas.
    -   `interfaces/`: Define os contratos (protocolos e classes base) que governam a interação entre os componentes do sistema, como `IExecutor` para executores de algoritmos e `IAlgorithm` para os próprios algoritmos.
    -   `scheduler/`: Implementa o `ExecutionScheduler`, um agendador avançado que gerencia uma fila de tarefas, controla o uso de recursos (CPU, memória) e executa os algoritmos de forma estável.
    -   `data/`: Estruturas de dados centrais, como `TaskResult`.
    -   `report/`: Lógica para geração de relatórios.
-   `datasets/`: Módulos para carregar e gerenciar os conjuntos de dados (datasets) usados nas execuções.
-   `optimization/`: Componentes relacionados a processos de otimização de hiperparâmetros.
-   `ui/`: Contém a lógica da interface do usuário.
    -   `cli/`: Implementação da interface de linha de comando (CLI), usando `argparse` para argumentos e `PyInquirer` para menus interativos.
    -   `curses_integration.py`: Módulo que fornece a visualização em tempo real do progresso das execuções no terminal, usando a biblioteca `curses`.
-   `utils/`: Funções utilitárias e auxiliares usadas em todo o projeto.
-   `main.py`: O ponto de entrada principal da aplicação, responsável por orquestrar a inicialização da CLI e o fluxo geral.

## Fluxo de Execução Principal

1.  **Inicialização**: O `main.py` é executado, que por sua vez invoca a aplicação CLI em `src/ui/cli/app.py`.
2.  **Parsing de Argumentos**: A CLI processa os argumentos da linha de comando para determinar o modo de operação (execução, otimização, etc.), o algoritmo a ser usado, o dataset e outros parâmetros.
3.  **Criação do Executor**: Um `SchedulerExecutor` é instanciado. Este executor é uma implementação da interface `IExecutor` que utiliza o `ExecutionScheduler` para gerenciar a execução das tarefas.
4.  **Submissão de Tarefas**: Para cada execução de algoritmo solicitada, uma instância do algoritmo correspondente (ex: `BLFGAAlgorithm`) é criada e submetida ao `SchedulerExecutor`.
5.  **Agendamento e Execução**:
    -   O `ExecutionScheduler` adiciona a tarefa a uma fila FIFO.
    -   Um loop de agendamento monitora continuamente os recursos do sistema (CPU/memória) e o número de tarefas ativas.
    -   Quando as condições são favoráveis (recursos disponíveis, delay entre tarefas respeitado), o agendador retira uma tarefa da fila e a executa em um `ThreadPoolExecutor`.
6.  **Monitoramento (UI)**: Se o modo visual estiver ativo, a `CursesApp` (`curses_integration.py`) é iniciada. Ela consulta periodicamente o `SchedulerExecutor` para obter o status de todas as tarefas (em fila, em execução, concluídas) e atualiza a tela do terminal em tempo real.
7.  **Coleta de Resultados**: Após a conclusão de todas as tarefas, a aplicação principal coleta os `TaskResult` de cada execução.
8.  **Geração de Relatórios**: Os resultados são processados e salvos em arquivos de relatório (JSON, CSV, etc.) no diretório `outputs/`.

## Design e Extensibilidade

-   **Inversão de Dependência**: O uso de interfaces como `IExecutor` e `IAlgorithm` (definida em `algorithms/base.py`) desacopla a lógica principal da implementação concreta dos algoritmos e dos executores.
-   **Extensibilidade de Algoritmos**: Para adicionar um novo algoritmo, basta criar uma nova classe que herde de `CSPAlgorithm` e usar o decorador `@register_algorithm`. O novo algoritmo será automaticamente descoberto e disponibilizado na CLI.
-   **Robustez**: O `ExecutionScheduler` foi projetado para ser robusto, evitando sobrecarregar o sistema ao executar múltiplos experimentos em paralelo.
