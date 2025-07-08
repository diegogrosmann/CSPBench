# Arquitetura do Código-Fonte (`src/`)

Este diretório contém o **núcleo da lógica de aplicação** do projeto CSP-BLFGA. A arquitetura foi projetada para ser **modular**, **extensível** e **robusta**, separando claramente as responsabilidades entre os diferentes componentes seguindo princípios de **Clean Architecture** e **SOLID**.

## 🏗️ Princípios Arquiteturais

### **Separação de Responsabilidades**
- **UI**: Interface de usuário isolada da lógica de negócio
- **Core**: Lógica central independente de frameworks
- **Utils**: Utilitários reutilizáveis e independentes
- **Interfaces**: Contratos bem definidos entre componentes

### **Inversão de Dependências**
- Uso de **protocolos** e **interfaces abstratas**
- Implementações concretas dependem de abstrações
- Facilita testes unitários e mocking

### **Extensibilidade**
- Sistema de **registro automático** de algoritmos
- **Factory patterns** para criação de objetos
- **Plugin architecture** para novos componentes

## 📂 Estrutura de Diretórios Detalhada

### **`core/` - Lógica Central do Sistema**

#### **`core/interfaces/` - Contratos e Protocolos**
```python
# Interfaces principais que definem os contratos do sistema
├── algorithm.py           # IAlgorithm - Interface para algoritmos CSP
├── executor.py           # IExecutor - Interface para executores
├── console.py            # IConsole - Interface para console/UI
├── task_result.py        # TaskResult - Estrutura padronizada de resultados
└── factory.py            # Factories para criação de objetos
```

**Características**:
- **Protocolos Python**: Uso de `typing.Protocol` para duck typing
- **Contratos Rígidos**: Métodos obrigatórios bem definidos
- **Backwards Compatibility**: Suporte a algoritmos legados

#### **`core/scheduler/` - Agendamento e Execução**
```python
# Sistema avançado de agendamento com controle de recursos
├── scheduler.py          # ExecutionScheduler - Fila FIFO + controle recursos
├── executor.py           # SchedulerExecutor - Wrapper para interface IExecutor
├── resource_monitor.py   # Monitoramento de CPU/Memória em tempo real
└── __init__.py          # Exports padronizados
```

**Características**:
- **Fila FIFO Absoluta**: Ordem rigorosa de execução
- **Controle de Recursos**: Monitoramento automático de CPU/memória
- **Timeout Configurável**: Prevenção de execuções infinitas
- **Process Watching**: Monitoramento de processos filhos
- **Thread Safety**: Operações thread-safe com locks

#### **`core/data/` - Estruturas de Dados**
```python
# Estruturas de dados centrais do sistema
├── task_result.py        # TaskResult - Resultado padronizado de execução
├── task_handle.py        # TaskHandle - Controle de tarefas em execução
└── execution_stats.py    # Estatísticas de execução
```

#### **`core/config/` - Configuração do Sistema**
```python
# Configurações centralizadas e validação
├── scheduler_config.py   # Configurações do scheduler
├── resource_config.py    # Limites e configurações de recursos
└── validation.py         # Validação de parâmetros
```

#### **`core/io/` - Entrada/Saída e Relatórios**
```python
# Sistema de I/O e geração de relatórios
├── results_formatter.py  # Formatação de resultados para diferentes formatos
├── exporter.py          # Exportação para CSV/JSON/outros formatos
└── report_generator.py   # Geração de relatórios detalhados
```

#### **`core/report/` - Geração de Relatórios**
```python
# Utilitários especializados para relatórios
├── report_utils.py       # Funções utilitárias para relatórios
├── statistics.py         # Cálculos estatísticos
└── visualization.py      # Geração de gráficos (quando disponível)
```

### **`ui/` - Interface de Usuário**

#### **`ui/cli/` - Interface de Linha de Comando**
```python
# Interface CLI completa e interativa
├── app.py               # Orquestrador principal da aplicação
├── menu.py              # Sistema de menus interativos
├── console_manager.py   # Gerenciamento thread-safe do console
└── save_wizard.py       # Wizard para salvamento de datasets
```

**Características**:
- **Modo Interativo**: Menus guiados com PyInquirer
- **Modo Silencioso**: Execução automatizada para scripts
- **Validação de Entrada**: Verificação de parâmetros
- **Progress Feedback**: Feedback visual de progresso

#### **`ui/curses_integration.py` - Interface Visual**
```python
# Sistema de monitoramento visual em tempo real
- Monitoramento de múltiplas tarefas simultâneas
- Atualização em tempo real de progresso
- Informações de recursos (CPU/memória)
- Controle de keyboard para interação
```

### **`datasets/` - Gestão de Datasets**
```python
# Módulos especializados para diferentes tipos de datasets
├── dataset_synthetic.py  # Geração de datasets sintéticos parametrizáveis
├── dataset_file.py      # Carregamento de arquivos FASTA/texto
├── dataset_entrez.py    # Download automático do NCBI/Entrez
└── dataset_validation.py # Validação e normalização de dados
```

**Características**:
- **Geração Sintética**: Datasets parametrizáveis com controle de ruído
- **Formatos Múltiplos**: Suporte a FASTA, texto, CSV
- **Validação Automática**: Verificação de integridade dos dados
- **Cache Inteligente**: Reutilização de datasets já processados

### **`optimization/` - Otimização de Hiperparâmetros**
```python
# Sistema de otimização baseado em Optuna
├── optimizer.py         # Otimizador principal com Optuna
├── objective_functions.py # Funções objetivo para otimização
├── parameter_spaces.py   # Definição de espaços de parâmetros
└── sensitivity_analysis.py # Análise de sensibilidade
```

### **`utils/` - Utilitários e Helpers**
```python
# Utilitários reutilizáveis em todo o projeto
├── config.py            # Configurações globais e constantes
├── logging.py           # Sistema de logging padronizado
├── distance.py          # Funções de distância (Hamming, etc.)
└── resource_monitor.py  # Monitoramento de recursos do sistema
```

## 🔄 Fluxo de Execução Detalhado

### **1. Inicialização (main.py → app.py)**
```python
# Sequência de inicialização
main.py → src.ui.cli.app.main() → ArgumentParser → Setup Logging
```

### **2. Processamento de Argumentos**
```python
# Análise de argumentos da CLI
ArgumentParser → Validation → Configuration → Mode Selection
```

### **3. Preparação do Dataset**
```python
# Preparação e validação de dados
Dataset Selection → Data Loading → Validation → Normalization
```

### **4. Criação do Executor**
```python
# Instanciação do sistema de execução
Factory.create_executor() → SchedulerExecutor → ExecutionScheduler
```

### **5. Submissão de Tarefas**
```python
# Para cada algoritmo selecionado
Algorithm Creation → Task Submission → Queue Management
```

### **6. Agendamento e Execução**
```python
# Sistema de agendamento inteligente
Resource Check → FIFO Queue → Thread Pool → Process Monitoring
```

### **7. Monitoramento Visual (Opcional)**
```python
# Interface curses para monitoramento
CursesApp → Real-time Updates → Resource Display → Progress Tracking
```

### **8. Coleta de Resultados**
```python
# Agregação e formatação de resultados
Result Collection → Statistics Calculation → Report Generation
```

### **9. Geração de Relatórios**
```python
# Saída final do sistema
Report Generation → File Export → CSV/JSON Output → Cleanup
```

## 🔧 Componentes Principais

### **ExecutionScheduler (Coração do Sistema)**
```python
class ExecutionScheduler:
    """
    Scheduler avançado com as seguintes características:
    - Fila FIFO absoluta (ordem rigorosa)
    - Controle automático de recursos
    - Delay configurável entre execuções
    - Monitoramento de processos filhos
    - Thread safety completo
    """
```

### **IAlgorithm (Interface Padronizada)**
```python
class IAlgorithm(Protocol):
    """
    Interface que todos os algoritmos devem implementar:
    - run() -> Result: Execução principal
    - set_progress_callback(): Callback de progresso
    - set_warning_callback(): Callback de warnings
    """
```

### **TaskResult (Estrutura Padronizada)**
```python
class TaskResult:
    """
    Resultado padronizado de execução:
    - center: String centro encontrada
    - distance: Distância calculada
    - metadata: Informações adicionais
    - execution_time: Tempo de execução
    - memory_used: Memória utilizada
    """
```

## 🚀 Padrões de Design Utilizados

### **Factory Pattern**
- **ExecutorFactory**: Criação de executores apropriados
- **ConsoleFactory**: Criação de interfaces de console
- **DatasetFactory**: Criação de datasets

### **Observer Pattern**
- **Progress Callbacks**: Notificação de progresso
- **Console Listeners**: Atualizações de UI
- **Resource Monitors**: Monitoramento de recursos

### **Strategy Pattern**
- **Algorithm Strategies**: Diferentes algoritmos CSP
- **Export Strategies**: Diferentes formatos de saída
- **Validation Strategies**: Diferentes tipos de validação

### **Command Pattern**
- **Task Submission**: Comandos para execução
- **Batch Processing**: Processamento em lote
- **Undo/Redo**: Operações reversíveis (futuro)

## 📊 Métricas e Monitoramento

### **Métricas Coletadas**
- **Performance**: Tempo de execução, throughput
- **Recursos**: CPU, memória, I/O
- **Qualidade**: Distância, convergência
- **Robustez**: Taxa de sucesso, estabilidade

### **Logging Estruturado**
```python
# Exemplo de log estruturado
{
    "timestamp": "2025-01-08T12:00:00Z",
    "level": "INFO",
    "component": "scheduler",
    "task_id": "uuid-1234",
    "algorithm": "BLF-GA",
    "event": "task_start",
    "metadata": {
        "dataset_size": 100,
        "parameters": {...}
    }
}
```

## 🔒 Thread Safety e Concorrência

### **Componentes Thread-Safe**
- **ExecutionScheduler**: Uso de locks e queues thread-safe
- **ConsoleManager**: Gerenciamento seguro de saída
- **ResourceMonitor**: Monitoramento concorrente
- **TaskResult**: Estruturas imutáveis

### **Sincronização**
- **Threading.Lock**: Proteção de seções críticas
- **Queue.Queue**: Comunicação thread-safe
- **Future**: Resultados assíncronos
- **ThreadPoolExecutor**: Pool de threads gerenciado

## 🧪 Testabilidade

### **Interfaces Mock-Friendly**
- **Protocolos Python**: Fácil criação de mocks
- **Dependency Injection**: Injeção de dependências
- **Factory Pattern**: Criação controlada para testes

### **Estrutura de Testes**
```python
# Exemplo de estrutura de teste
tests/
├── unit/                 # Testes unitários
│   ├── test_scheduler.py
│   ├── test_algorithms.py
│   └── test_interfaces.py
├── integration/          # Testes de integração
│   ├── test_full_flow.py
│   └── test_batch_processing.py
└── fixtures/            # Dados de teste
    ├── sample_datasets/
    └── expected_results/
```

## 📈 Performance e Otimizações

### **Otimizações Implementadas**
- **Lazy Loading**: Carregamento sob demanda
- **Resource Pooling**: Reutilização de recursos
- **Memory Management**: Limpeza automática
- **Parallel Processing**: Execução paralela quando possível

### **Profiling e Debugging**
- **cProfile**: Profiling de performance
- **Memory Profiler**: Análise de uso de memória
- **Logging Detalhado**: Debugging avançado
- **Resource Monitoring**: Monitoramento em tempo real

## 🔮 Extensibilidade e Futuro

### **Pontos de Extensão**
- **Novos Algoritmos**: Sistema de registro automático
- **Novos Datasets**: Interface padronizada
- **Novos Formatos**: Exportadores plugáveis
- **Novas UIs**: Interfaces intercambiáveis

### **Roadmap Técnico**
- **Distributed Computing**: Execução em cluster
- **Web Interface**: Interface web complementar
- **Database Integration**: Persistência de resultados
- **Machine Learning**: Otimização automática de parâmetros

---

Esta arquitetura garante que o CSP-BLFGA seja não apenas funcional, mas também **mantível**, **extensível** e **robusta** para uso em pesquisa científica e desenvolvimento de novos algoritmos para o Closest String Problem.

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
