# Documentação Técnica - CSP-BLFGA

## 📋 Visão Geral Técnica

O CSP-BLFGA é uma plataforma experimental robusta para o Closest String Problem, implementada em Python com foco em **performance**, **extensibilidade** e **robustez**. Esta documentação fornece detalhes técnicos para desenvolvedores e pesquisadores.

## 🏗️ Arquitetura do Sistema

### Princípios Fundamentais

1. **Separação de Responsabilidades**: Cada componente tem uma função específica
2. **Inversão de Dependências**: Uso de interfaces e protocolos
3. **Composição sobre Herança**: Preferência por composição
4. **Fail-Fast**: Detecção precoce de erros
5. **Thread Safety**: Operações seguras em ambiente concorrente

### Camadas da Arquitetura

```
┌─────────────────────────────────────────────────────────────────┐
│                    PRESENTATION LAYER                           │
│  ┌─────────────────┐  ┌─────────────────┐  ┌─────────────────┐  │
│  │   CLI Interface │  │  Curses Monitor │  │   Batch Runner  │  │
│  └─────────────────┘  └─────────────────┘  └─────────────────┘  │
├─────────────────────────────────────────────────────────────────┤
│                     APPLICATION LAYER                          │
│  ┌─────────────────┐  ┌─────────────────┐  ┌─────────────────┐  │
│  │   App Router    │  │  Menu System    │  │ Report Generator│  │
│  └─────────────────┘  └─────────────────┘  └─────────────────┘  │
├─────────────────────────────────────────────────────────────────┤
│                      DOMAIN LAYER                              │
│  ┌─────────────────┐  ┌─────────────────┐  ┌─────────────────┐  │
│  │   Scheduler     │  │   Algorithms    │  │   Datasets      │  │
│  │   Executor      │  │   Registry      │  │   Loaders       │  │
│  └─────────────────┘  └─────────────────┘  └─────────────────┘  │
├─────────────────────────────────────────────────────────────────┤
│                  INFRASTRUCTURE LAYER                          │
│  ┌─────────────────┐  ┌─────────────────┐  ┌─────────────────┐  │
│  │ Resource Monitor│  │  File System    │  │    Logging      │  │
│  │ Process Watcher │  │  Serialization  │  │   Metrics       │  │
│  └─────────────────┘  └─────────────────┘  └─────────────────┘  │
└─────────────────────────────────────────────────────────────────┘
```

## 🔧 Componentes Centrais

### 1. ExecutionScheduler

**Localização**: `src/core/scheduler/scheduler.py`

**Responsabilidades**:
- Gerenciar fila FIFO de tarefas
- Controlar recursos do sistema
- Executar algoritmos com timeout
- Monitorar processos filhos

**Características Técnicas**:
```python
class ExecutionScheduler:
    def __init__(self, 
                 start_delay: float = 2.0,
                 cpu_threshold: float = 85.0,
                 mem_threshold: float = 300.0,
                 min_active_tasks: int = 1,
                 timeout: float = 300.0):
        # Configurações de recursos
        self.resource_checker = ResourceChecker(cpu_threshold, mem_threshold)
        self.process_watcher = ProcessWatcher()
        self.resource_monitor = ResourceMonitor()
        
        # Estruturas de dados thread-safe
        self.task_queue = queue.Queue()
        self.scheduled_tasks: Dict[str, ScheduledTask] = {}
        self.active_tasks: Set[str] = set()
        
        # Thread pool para execução
        self.executor = ThreadPoolExecutor(max_workers=os.cpu_count())
```

**Padrões de Design**:
- **Observer**: Para monitoramento de recursos
- **Command**: Para encapsular tarefas
- **State**: Para controlar estados de tarefa

### 2. Algorithm Registry

**Localização**: `algorithms/base.py`

**Responsabilidades**:
- Registrar algoritmos automaticamente
- Validar interfaces de algoritmos
- Fornecer factory para criação

**Implementação**:
```python
# Registry global thread-safe
global_registry: Dict[str, Type] = {}

def register_algorithm(cls: Type) -> Type:
    """Decorador para registro automático"""
    name = getattr(cls, "name", cls.__name__)
    global_registry[name] = cls
    return cls

@register_algorithm
class MeuAlgoritmo(CSPAlgorithm):
    name = "MeuAlgoritmo"
    # Implementação...
```

### 3. Interface System

**Localização**: `src/core/interfaces/`

**Protocolos Principais**:
```python
class IAlgorithm(Protocol):
    """Interface para algoritmos CSP"""
    def run(self) -> Result: ...
    def set_progress_callback(self, callback: Callable[[str], None]) -> None: ...
    def set_warning_callback(self, callback: Callable[[str], None]) -> None: ...

class IExecutor(Protocol):
    """Interface para executores"""
    def submit(self, algorithm: IAlgorithm, **kwargs) -> str: ...
    def get_result(self, task_id: str) -> TaskResult: ...
    def get_status(self, task_id: str) -> TaskStatus: ...
```

## 🚀 Fluxo de Dados Detalhado

### Ciclo de Vida de uma Tarefa

```
┌─────────────────┐     ┌─────────────────┐     ┌─────────────────┐
│   Task Created  │────▶│   Task Queued   │────▶│  Task Scheduled │
└─────────────────┘     └─────────────────┘     └─────────────────┘
                                                          │
┌─────────────────┐     ┌─────────────────┐     ┌─────────────────┐
│ Task Completed  │◀────│  Task Running   │◀────│ Resource Check  │
└─────────────────┘     └─────────────────┘     └─────────────────┘
```

### Estados de Tarefa

```python
class TaskState(Enum):
    QUEUED = "queued"        # Na fila esperando recursos
    RUNNING = "running"      # Executando ativamente
    COMPLETED = "completed"  # Concluída com sucesso
    FAILED = "failed"        # Falhou por erro
    CANCELLED = "cancelled"  # Cancelada pelo usuário
```

### Estrutura de Dados TaskResult

```python
@dataclass
class TaskResult:
    task_id: str
    algorithm_name: str
    dataset_info: Dict[str, Any]
    execution_time: float
    memory_used: float
    result: Dict[str, Any]  # center, distance, metadata
    status: TaskStatus
    error: Optional[str] = None
    warnings: List[str] = field(default_factory=list)
    
    # Métricas adicionais
    cpu_usage: Optional[float] = None
    start_time: Optional[float] = None
    end_time: Optional[float] = None
```

## 📊 Sistema de Monitoramento

### Resource Monitor

**Implementação**: `src/core/scheduler/resource_monitor.py`

```python
class ResourceMonitor:
    def __init__(self):
        self.cpu_threshold = 85.0
        self.memory_threshold = 300.0  # MB
        self.monitoring_interval = 1.0  # segundos
        
    def check_resources(self) -> bool:
        """Verifica se recursos estão disponíveis"""
        cpu_usage = psutil.cpu_percent(interval=0.1)
        memory_usage = psutil.virtual_memory().used / (1024 * 1024)
        
        return (cpu_usage < self.cpu_threshold and 
                memory_usage < self.memory_threshold)
```

### Process Watcher

```python
class ProcessWatcher:
    def __init__(self):
        self.watched_processes: Dict[str, psutil.Process] = {}
        self.timeouts: Dict[str, float] = {}
        
    def watch_process(self, task_id: str, process: psutil.Process, timeout: float):
        """Adiciona processo ao monitoramento"""
        self.watched_processes[task_id] = process
        self.timeouts[task_id] = time.time() + timeout
        
    def check_timeouts(self) -> List[str]:
        """Verifica processos que excederam timeout"""
        current_time = time.time()
        timed_out = []
        
        for task_id, deadline in self.timeouts.items():
            if current_time > deadline:
                timed_out.append(task_id)
                
        return timed_out
```

## 🔌 Sistema de Plugins

### Registrando Novo Algoritmo

```python
from algorithms.base import CSPAlgorithm, register_algorithm

@register_algorithm
class NovoAlgoritmo(CSPAlgorithm):
    name = "NovoAlgoritmo"
    default_params = {
        'param1': 10,
        'param2': 0.5
    }
    is_deterministic = False
    supports_internal_parallel = True
    
    def __init__(self, strings: List[str], alphabet: str, **params):
        super().__init__(strings, alphabet, **params)
        # Inicialização específica
        
    def run(self) -> Tuple[str, int, Dict[str, Any]]:
        # Implementação do algoritmo
        self._report_progress("Iniciando algoritmo...")
        
        # Lógica principal
        center = self.encontrar_centro()
        distance = self.calcular_distancia(center)
        
        # Metadados
        metadata = {
            'iterations': self.iterations,
            'convergence_time': self.convergence_time,
            'custom_metric': self.custom_metric
        }
        
        return center, distance, metadata
```

### Criando Novo Dataset Loader

```python
from src.datasets.base import DatasetLoader

class NovoDatasetLoader(DatasetLoader):
    def __init__(self, **params):
        super().__init__(**params)
        
    def load(self) -> Tuple[List[str], Dict[str, Any]]:
        # Implementação do carregamento
        sequences = self.carregar_sequencias()
        metadata = self.extrair_metadados()
        
        return sequences, metadata
        
    def validate(self, sequences: List[str]) -> bool:
        # Validação específica
        return all(self.validar_sequencia(seq) for seq in sequences)
```

## 🧪 Testing Framework

### Estrutura de Testes

```python
# tests/unit/test_scheduler.py
import pytest
from unittest.mock import Mock, patch
from src.core.scheduler.scheduler import ExecutionScheduler

class TestExecutionScheduler:
    def setup_method(self):
        self.scheduler = ExecutionScheduler(
            start_delay=0.1,
            cpu_threshold=90.0,
            mem_threshold=500.0
        )
        
    def test_task_submission(self):
        # Teste de submissão de tarefa
        mock_algorithm = Mock()
        task_id = self.scheduler.submit(mock_algorithm.run)
        
        assert task_id is not None
        assert task_id in self.scheduler.scheduled_tasks
        
    @patch('src.core.scheduler.resource_monitor.psutil.cpu_percent')
    def test_resource_check(self, mock_cpu):
        # Teste de verificação de recursos
        mock_cpu.return_value = 50.0
        
        assert self.scheduler.resource_checker.check_resources() == True
        
        mock_cpu.return_value = 95.0
        assert self.scheduler.resource_checker.check_resources() == False
```

### Fixtures para Testes

```python
# tests/fixtures/datasets.py
import pytest

@pytest.fixture
def sample_sequences():
    return [
        "ACGTACGTACGT",
        "ACGTACGTACGA",
        "ACGTACGTACGC",
        "ACGTACGTACGG"
    ]

@pytest.fixture
def sample_algorithm_params():
    return {
        'population_size': 50,
        'generations': 100,
        'mutation_rate': 0.01
    }
```

## 📈 Performance e Otimizações

### Profiling de Performance

```python
import cProfile
import pstats
from src.ui.cli.app import main

def profile_main():
    """Profiling da aplicação principal"""
    profiler = cProfile.Profile()
    profiler.enable()
    
    # Executar aplicação
    main()
    
    profiler.disable()
    stats = pstats.Stats(profiler)
    stats.sort_stats('cumulative')
    stats.print_stats(20)
```

### Memory Profiling

```python
from memory_profiler import profile

@profile
def memory_intensive_function():
    # Função que usa muita memória
    large_data = [i for i in range(1000000)]
    return process_data(large_data)
```

### Otimizações Implementadas

1. **Lazy Loading**: Carregamento de dados sob demanda
2. **Object Pooling**: Reutilização de objetos custosos
3. **Memory Cleanup**: Limpeza automática de memória
4. **Resource Throttling**: Controle de uso de recursos

## 🔐 Segurança e Robustez

### Tratamento de Erros

```python
class CSPError(Exception):
    """Exceção base para erros CSP"""
    pass

class AlgorithmError(CSPError):
    """Erro específico de algoritmo"""
    pass

class ResourceError(CSPError):
    """Erro de recursos insuficientes"""
    pass

class TimeoutError(CSPError):
    """Erro de timeout"""
    pass
```

### Validação de Entrada

```python
def validate_sequences(sequences: List[str]) -> None:
    """Valida lista de sequências"""
    if not sequences:
        raise ValueError("Lista de sequências vazia")
        
    if not all(isinstance(seq, str) for seq in sequences):
        raise TypeError("Todas as sequências devem ser strings")
        
    # Verificar se todas têm mesmo comprimento
    length = len(sequences[0])
    if not all(len(seq) == length for seq in sequences):
        raise ValueError("Todas as sequências devem ter o mesmo comprimento")
```

### Logging Structured

```python
import logging
import json

class StructuredLogger:
    def __init__(self, name: str):
        self.logger = logging.getLogger(name)
        
    def log_event(self, event_type: str, **kwargs):
        """Log estruturado em formato JSON"""
        log_data = {
            'timestamp': time.time(),
            'event_type': event_type,
            'data': kwargs
        }
        self.logger.info(json.dumps(log_data))
```

## 🚀 Deployment e Configuração

### Configuração de Ambiente

```python
# config/production.py
import os

class ProductionConfig:
    # Scheduler
    SCHEDULER_START_DELAY = 2.0
    SCHEDULER_CPU_THRESHOLD = 80.0
    SCHEDULER_MEMORY_THRESHOLD = 500.0
    
    # Logging
    LOG_LEVEL = os.getenv('LOG_LEVEL', 'INFO')
    LOG_FILE = os.getenv('LOG_FILE', 'outputs/logs/production.log')
    
    # Resources
    MAX_WORKERS = int(os.getenv('MAX_WORKERS', '4'))
    ALGORITHM_TIMEOUT = int(os.getenv('ALGORITHM_TIMEOUT', '300'))
```

### Docker Configuration

```dockerfile
# Dockerfile
FROM python:3.11-slim

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

CMD ["python", "main.py"]
```

### CI/CD Pipeline

```yaml
# .github/workflows/ci.yml
name: CI/CD Pipeline

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
    - uses: actions/checkout@v2
    - name: Set up Python
      uses: actions/setup-python@v2
      with:
        python-version: 3.11
    - name: Install dependencies
      run: |
        pip install -r requirements.txt
        pip install pytest pytest-cov
    - name: Run tests
      run: pytest tests/ --cov=src --cov-report=xml
    - name: Upload coverage
      uses: codecov/codecov-action@v1
```

## 📝 Padrões de Código

### Naming Conventions

- **Classes**: PascalCase (`ExecutionScheduler`)
- **Funções**: snake_case (`submit_task`)
- **Constantes**: UPPER_SNAKE_CASE (`MAX_WORKERS`)
- **Variáveis**: snake_case (`task_id`)

### Code Style

```python
# Usar type hints sempre
def process_sequences(sequences: List[str], alphabet: str) -> Tuple[str, int]:
    """Processa sequências e retorna resultado."""
    pass

# Documentar classes e métodos
class AlgorithmExecutor:
    """
    Executor para algoritmos CSP.
    
    Attributes:
        scheduler: Scheduler para execução
        timeout: Timeout padrão em segundos
    """
    
    def __init__(self, scheduler: ExecutionScheduler, timeout: float = 300.0):
        self.scheduler = scheduler
        self.timeout = timeout
```

### Error Handling

```python
def safe_execute(func: Callable) -> Tuple[Any, Optional[str]]:
    """Executa função com tratamento de erro seguro."""
    try:
        result = func()
        return result, None
    except Exception as e:
        logger.error(f"Erro na execução: {e}")
        return None, str(e)
```

## 🔧 Debugging e Troubleshooting

### Common Issues

1. **Memory Leaks**: Usar `tracemalloc` para rastrear
2. **Deadlocks**: Verificar ordem de locks
3. **Resource Exhaustion**: Monitorar limites
4. **Performance Issues**: Profiling regular

### Debug Tools

```python
# Debug decorator
def debug_execution(func):
    def wrapper(*args, **kwargs):
        logger.debug(f"Executing {func.__name__} with args={args}, kwargs={kwargs}")
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        logger.debug(f"Execution took {end_time - start_time:.2f} seconds")
        return result
    return wrapper
```

---

Esta documentação técnica fornece uma base sólida para desenvolvedores trabalharem com o CSP-BLFGA, garantindo manutenibilidade e extensibilidade do código.
