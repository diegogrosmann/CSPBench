# CSPBench: Framework Experimental para Teste de Algoritmos CSP

**CSPBench** é um framework experimental robusto e extensível para teste, comparação e análise de algoritmos do **Closest String Problem (CSP)**. O framework oferece uma plataforma unificada para desenvolvimento, execução e avaliação de algoritmos CSP com recursos avançados de paralelização, monitoramento e relatórios.

## 🎯 Visão Geral

O **CSPBench** é um framework científico avançado para experimentação com algoritmos do Closest String Problem, oferecendo:

### 🧬 **Biblioteca de Algoritmos CSP**
- **Sistema de Registro Automático**: Algoritmos são detectados automaticamente via decoradores
- **Interface Padronizada**: Todos os algoritmos seguem contratos bem definidos
- **Algoritmos Incluídos**: Baseline, BLF-GA, CSC, DP-CSP, H³-CSP
- **Extensibilidade**: Fácil adição de novos algoritmos

### 📊 **Gestão Avançada de Datasets**
- **Geração Sintética**: Datasets parametrizáveis com controle de ruído
- **Carregamento de Arquivos**: Suporte a formatos FASTA e texto
- **Download Automático**: Integração com NCBI/Entrez para dados reais
- **Processamento em Lote**: Múltiplos datasets com configuração YAML

### 🚀 **Sistema de Execução Inteligente**
- **Scheduler Avançado**: Controle de recursos com fila FIFO rigorosa
- **Execução Paralela**: Balanceamento dinâmico de workers
- **Monitoramento Visual**: Interface curses em tempo real
- **Controle de Timeout**: Prevenção de execuções infinitas
- **Logging Estruturado**: Rastreamento completo de operações

### 📈 **Análise e Otimização**
- **Relatórios Automáticos**: Geração de relatórios JSON/CSV/HTML detalhados
- **Análise Comparativa**: Comparação estatística entre algoritmos
- **Otimização de Hiperparâmetros**: Integração com Optuna para tuning automático
- **Análise de Sensibilidade**: Estudo do impacto de parâmetros (SALib)

## 🔄 Fluxos do Framework

### 1. **Desenvolvimento de Algoritmos**
```
Implementação → Registro Automático → Integração → Testes
```

### 2. **Experimentação Científica**
```
Datasets → Configuração → Execução Paralela → Análise Estatística
```

### 3. **Otimização de Performance**
```
Algoritmo → Espaço de Parâmetros → Otuna/SALib → Melhores Configurações
```

### 4. **Análise Comparativa**
```
Múltiplos Algoritmos → Execução Controlada → Métricas → Relatórios
```

## 🏗️ Arquitetura do CSPBench

### **Visão Geral da Arquitetura**

```
┌─────────────────────────────────────────────────────────────────┐
│                           CSPBench                              │
│                Framework de Teste CSP                           │
├─────────────────────────────────────────────────────────────────┤
│  🖥️  INTERFACE DE USUÁRIO                                       │
│  ├── CLI Interativa          # Menus e wizards                 │
│  ├── CLI Silenciosa         # Automação e scripts              │
│  ├── Execução em Lote       # Configurações YAML              │
│  └── Monitoramento Curses   # Interface visual em tempo real  │
├─────────────────────────────────────────────────────────────────┤
│  🧮 NÚCLEO DO FRAMEWORK                                         │
│  ├── 🎯 Sistema de Execução                                     │
│  │   ├── ExecutionScheduler   # Fila FIFO + Controle recursos │
│  │   ├── ResourceMonitor     # CPU/Memória em tempo real      │
│  │   ├── TaskManager         # Gerenciamento de tarefas       │
│  │   └── ProcessWatcher      # Controle de processos filhos   │
│  ├── 🔌 Interfaces Padronizadas                                │
│  │   ├── IAlgorithm          # Contrato para algoritmos       │
│  │   ├── IExecutor           # Contrato para executores       │
│  │   ├── IDataset            # Contrato para datasets         │
│  │   └── IConsole            # Contrato para interfaces       │
│  ├── 📊 Gestão de Dados                                        │
│  │   ├── TaskResult          # Estruturas de resultados       │
│  │   ├── TaskHandle          # Controle de tarefas            │
│  │   └── MetadataCollector   # Coleta de métricas             │
│  └── 📈 Sistema de Relatórios                                  │
│      ├── ResultsFormatter    # Formatação multi-formato       │
│      ├── StatisticalAnalyzer # Análise estatística            │
│      └── ReportGenerator     # Geração de relatórios          │
├─────────────────────────────────────────────────────────────────┤
│  🧬 BIBLIOTECA DE ALGORITMOS                                    │
│  ├── Sistema de Registro     # Descoberta automática          │
│  ├── Baseline               # Consenso ganancioso             │
│  ├── BLF-GA                 # Algoritmo genético híbrido      │
│  ├── CSC                    # Consensus String Clustering     │
│  ├── DP-CSP                 # Programação dinâmica exata      │
│  ├── H³-CSP                 # Busca hierárquica híbrida       │
│  └── [Extensível]           # Interface para novos algoritmos │
├─────────────────────────────────────────────────────────────────┤
│  📁 GESTÃO DE DATASETS                                          │
│  ├── Geração Sintética      # Datasets parametrizáveis        │
│  ├── Carregamento de Arquivos # FASTA, texto, formatos customizados │
│  ├── Download NCBI/Entrez   # Dados biológicos reais          │
│  ├── Processamento em Lote  # Múltiplos datasets              │
│  └── Cache Inteligente      # Otimização de acesso            │
├─────────────────────────────────────────────────────────────────┤
│  🔧 ANÁLISE E OTIMIZAÇÃO                                        │
│  ├── Optuna Integration     # Otimização de hiperparâmetros    │
│  ├── SALib Integration      # Análise de sensibilidade        │
│  ├── Statistical Analysis   # Análise estatística avançada    │
│  ├── Visualization Tools    # Gráficos e visualizações        │
│  └── Parallel Processing    # Execução paralela otimizada     │
└─────────────────────────────────────────────────────────────────┘
```

### **Fluxo de Dados**

```
📊 Dataset → 🎯 Scheduler → 🧬 Algoritmos → 📈 Análise → 💾 Relatórios
    ↓             ↓             ↓             ↓             ↓
Validação    Fila FIFO    Execução     Coleta de      Análise
Normalização  Recursos    Paralela     Métricas      Estatística
Cache        Timeouts    Callbacks    Agregação     Visualização
```

### **Estrutura de Diretórios**

```
CSPBench/
├── algorithms/                    # 🧬 Biblioteca de Algoritmos CSP
│   ├── base.py                   # Interface base e registro automático
│   ├── baseline/                 # Algoritmo de consenso ganancioso
│   ├── blf_ga/                   # Blockwise Learning Fusion + GA
│   ├── csc/                      # Consensus String Clustering
│   ├── dp_csp/                   # Programação Dinâmica (exato)
│   ├── h3_csp/                   # Hybrid Hierarchical Hamming Search
│   └── README.md                 # Guia para desenvolvimento de algoritmos
├── src/                          # � Núcleo do Framework
│   ├── core/                     # Sistema central
│   │   ├── interfaces/           # Contratos e protocolos
│   │   ├── scheduler/            # Sistema de execução
│   │   ├── io/                   # Entrada/saída de dados
│   │   └── report/               # Geração de relatórios
│   ├── datasets/                 # Gestão de datasets
│   ├── optimization/             # Otimização e análise
│   ├── ui/                       # Interfaces de usuário
│   └── utils/                    # Utilitários gerais
├── batch_configs/                # ⚙️ Configurações de Lote
├── tests/                        # 🧪 Testes automatizados
├── docs/                         # � Documentação
├── outputs/                      # � Resultados e logs
│   ├── reports/                  # Relatórios gerados
│   └── logs/                     # Logs de execução
├── saved_datasets/               # 💾 Datasets salvos
├── main.py                       # 🎯 Ponto de entrada
└── requirements.txt              # 📋 Dependências
```
```

## 🚀 Instalação e Configuração

### **Pré-requisitos**

- **Python 3.8+** (recomendado: 3.10 ou 3.11)
- **Sistema Operacional**: Linux (recomendado) / macOS / Windows
- **Terminal**: Suporte a cores (para interface curses)
- **Memória**: Mínimo 4GB RAM (8GB+ recomendado)
- **CPU**: Multi-core recomendado para paralelização

### **Instalação**

#### **Método 1: Ambiente Virtual (Recomendado)**

```bash
# Clonar o repositório
git clone https://github.com/diegogrosmann/CSPBench.git
cd CSPBench

# Criar ambiente virtual
python -m venv .venv

# Ativar ambiente virtual
source .venv/bin/activate  # Linux/macOS
# ou
.venv\Scripts\activate     # Windows

# Instalar dependências
pip install --upgrade pip
pip install -r requirements.txt

# Instalar dependências de desenvolvimento (opcional)
pip install -e .[dev]
```

#### **Método 2: Instalação Direta**

```bash
git clone https://github.com/diegogrosmann/CSPBench.git
cd CSPBench
pip install -r requirements.txt
```

### **Verificação da Instalação**

```bash
# Executar testes básicos
python -m pytest tests/ -v

# Verificar algoritmos disponíveis
python -c "from algorithms.base import global_registry; print(list(global_registry.keys()))"

# Teste de funcionalidade básica
python main.py --help
```

## 🎮 Como Usar o CSPBench

### **Execução Básica**

```bash
# Executar interface interativa
python main.py

# Executar com monitoramento visual
python main.py --visual

# Executar com logging detalhado
python main.py --debug
```

### **Execução Silenciosa (Automação)**

```bash
# Execução automatizada para benchmarks
python main.py --silent --dataset synthetic --algorithms Baseline BLF-GA CSC --num-execs 10

# Execução com configuração específica
python main.py --silent --dataset file --algorithms H3-CSP --timeout 600

# Execução com paralelismo configurado
python main.py --silent --dataset synthetic --algorithms "BLF-GA" --workers 8
```

### **Execução em Lote**

```bash
# Execução com arquivo de configuração YAML
python main.py --batch batch_configs/benchmark_completo.yaml

# Execução em lote silenciosa
python main.py --batch batch_configs/otimizacao_exemplo.yaml --silent

# Execução de análise de sensibilidade
python main.py --batch batch_configs/sensibilidade_unificado.yaml
```

### **Parâmetros da CLI**

| Parâmetro | Descrição | Exemplo |
|-----------|-----------|---------|
| `--silent` | Modo silencioso (sem interação) | `--silent` |
| `--batch FILE` | Arquivo de configuração YAML | `--batch config.yaml` |
| `--dataset TYPE` | Tipo de dataset (`synthetic`, `file`, `entrez`) | `--dataset synthetic` |
| `--algorithms ALG [ALG ...]` | Lista de algoritmos | `--algorithms Baseline BLF-GA` |
| `--num-execs N` | Número de execuções por algoritmo | `--num-execs 10` |
| `--timeout N` | Timeout por execução (segundos) | `--timeout 300` |
| `--workers N` | Número de workers paralelos | `--workers 4` |
| `--visual` | Interface visual com curses | `--visual` |
| `--debug` | Logging detalhado | `--debug` |

### **Interface Interativa**

A interface interativa oferece menus guiados para:

1. **Seleção de Dataset**
   - Geração sintética com parâmetros customizáveis
   - Carregamento de arquivos FASTA/texto
   - Download de dados do NCBI/Entrez

2. **Configuração de Algoritmos**
   - Seleção múltipla de algoritmos
   - Configuração de parâmetros específicos
   - Otimização de hiperparâmetros

3. **Configuração de Execução**
   - Número de execuções por algoritmo
   - Configuração de timeouts
   - Opções de paralelização

4. **Monitoramento e Resultados**
   - Progresso em tempo real
   - Relatórios automáticos
   - Análise estatística

## 🧬 Desenvolvimento de Algoritmos

O CSPBench oferece um sistema extensível para desenvolvimento de novos algoritmos CSP.

### **Adicionando um Novo Algoritmo**

#### **1. Estrutura de Diretórios**

```
algorithms/
  meu_algoritmo/
    __init__.py           # Módulo Python
    algorithm.py          # Wrapper do framework
    config.py            # Configurações padrão
    implementation.py    # Implementação core
    README.md           # Documentação específica
```

#### **2. Implementação da Interface**

```python
# algorithms/meu_algoritmo/algorithm.py
from algorithms.base import CSPAlgorithm, register_algorithm
from .config import MEU_ALGORITMO_DEFAULTS
from .implementation import MeuAlgoritmoCore

@register_algorithm
class MeuAlgoritmo(CSPAlgorithm):
    """Meu algoritmo personalizado para CSP."""
    
    name = "MeuAlgoritmo"
    default_params = MEU_ALGORITMO_DEFAULTS
    is_deterministic = False  # ou True se determinístico
    supports_internal_parallel = True  # se suporta paralelismo
    
    def __init__(self, strings: list[str], alphabet: str, **params):
        super().__init__(strings, alphabet, **params)
        self.core = MeuAlgoritmoCore(strings, alphabet, **self.params)
    
    def run(self) -> tuple[str, int, dict]:
        """Executa o algoritmo."""
        self._report_progress("Iniciando algoritmo...")
        
        # Executar implementação core
        center = self.core.solve()
        distance = self.core.calculate_distance(center)
        metadata = self.core.get_metadata()
        
        return center, distance, metadata
```

#### **3. Configuração de Parâmetros**

```python
# algorithms/meu_algoritmo/config.py
MEU_ALGORITMO_DEFAULTS = {
    "max_iterations": 1000,
    "convergence_threshold": 1e-6,
    "population_size": 50,
    "mutation_rate": 0.1,
    # ... outros parâmetros
}
```

#### **4. Implementação Core**

```python
# algorithms/meu_algoritmo/implementation.py
class MeuAlgoritmoCore:
    """Implementação core do algoritmo."""
    
    def __init__(self, strings, alphabet, **params):
        self.strings = strings
        self.alphabet = alphabet
        self.params = params
    
    def solve(self) -> str:
        """Resolve o CSP e retorna a string central."""
        # Implementação do algoritmo
        pass
    
    def calculate_distance(self, center: str) -> int:
        """Calcula distância máxima."""
        pass
    
    def get_metadata(self) -> dict:
        """Retorna metadados da execução."""
        pass
```

#### **5. Registro Automático**

O sistema detecta automaticamente novos algoritmos através do decorador `@register_algorithm`. Não é necessário modificar nenhum arquivo do framework.

### **Características Avançadas**

#### **Callbacks de Progresso**
```python
def run(self):
    self._report_progress("Fase 1: Inicialização")
    # ... código ...
    self._report_progress("Fase 2: Otimização")
    # ... código ...
    self._report_progress("Fase 3: Refinamento")
```

#### **Suporte a Paralelismo**
```python
supports_internal_parallel = True  # Algoritmo usa paralelismo interno

def run(self):
    # Acessar número de workers internos
    internal_workers = os.environ.get('INTERNAL_WORKERS', '1')
    # Configurar paralelismo interno
```

#### **Tratamento de Warnings**
```python
def run(self):
    if self.params['population_size'] < 10:
        self._report_warning("População muito pequena, performance pode ser afetada")
```

### **Testes e Validação**

```python
# tests/test_meu_algoritmo.py
import pytest
from algorithms.meu_algoritmo.algorithm import MeuAlgoritmo

def test_meu_algoritmo_basic():
    strings = ["ACGT", "AGCT", "ATCT"]
    alg = MeuAlgoritmo(strings, "ACGT")
    center, distance, metadata = alg.run()
    
    assert isinstance(center, str)
    assert isinstance(distance, int)
    assert isinstance(metadata, dict)
```

### **Documentação do Algoritmo**

Cada algoritmo deve incluir um README.md detalhado com:

- **Descrição**: O que o algoritmo faz e como funciona
- **Heurísticas**: Estratégias e técnicas utilizadas
- **Parâmetros**: Descrição de todos os parâmetros configuráveis
- **Casos de Uso**: Quando usar este algoritmo
- **Limitações**: Restrições e cenários não recomendados
- **Exemplos**: Código de exemplo e casos de uso

Consulte `algorithms/*/README.md` para exemplos de documentação.

## � Análise e Otimização

### **Otimização de Hiperparâmetros**

O CSPBench integra o **Optuna** para otimização automática de hiperparâmetros:

```bash
# Otimização via CLI
python main.py --optimize --algorithm BLF-GA --dataset synthetic --trials 100

# Otimização em lote via YAML
python main.py --batch batch_configs/otimizacao_exemplo.yaml
```

**Configuração YAML para Otimização:**
```yaml
task:
  type: "optimization"
  optimization:
    studies:
      - nome: "Otimização BLF-GA"
        datasets: [dataset_1]
        n_trials: 100
        timeout_per_trial: 60
        param_space:
          "BLF-GA":
            population_size: ["int", 50, 300]
            max_generations: ["int", 100, 500]
            mutation_rate: ["uniform", 0.01, 0.3]
```

### **Análise de Sensibilidade**

Integração com **SALib** para análise de sensibilidade de parâmetros:

```bash
# Análise de sensibilidade via CLI
python main.py --sensitivity --algorithm BLF-GA --param population_size --range 50,200

# Análise em lote via YAML
python main.py --batch batch_configs/sensibilidade_unificado.yaml
```

**Configuração YAML para Sensibilidade:**
```yaml
task:
  type: "sensitivity"
  sensitivity:
    analyses:
      - nome: "Sensibilidade BLF-GA"
        datasets: [dataset_1]
        n_samples: 1000
        method: "morris"  # morris, sobol, fast
        param_space:
          "BLF-GA": ["population_size", "max_generations", "mutation_rate"]
```

### **Relatórios e Visualizações**

#### **Tipos de Relatórios Gerados**

1. **Relatórios JSON**: Dados estruturados completos
2. **Relatórios CSV**: Dados tabulares para análise externa
3. **Relatórios HTML**: Visualizações interativas
4. **Logs Detalhados**: Rastreamento completo de execuções

#### **Estrutura de Resultados**

```json
{
  "execution_info": {
    "timestamp": "2025-01-09T12:00:00Z",
    "framework_version": "1.0.0",
    "total_execution_time": 125.43
  },
  "dataset_info": {
    "type": "synthetic",
    "n_strings": 50,
    "string_length": 100,
    "alphabet": "ACGT",
    "noise_level": 0.15
  },
  "algorithm_results": {
    "BLF-GA": {
      "executions": [
        {
          "execution_id": 1,
          "center": "ACGTACGTACGT...",
          "distance": 15,
          "execution_time": 45.2,
          "memory_used": 128.5,
          "metadata": {
            "generations": 245,
            "convergence_generation": 201,
            "final_population_diversity": 0.85
          }
        }
      ],
      "statistics": {
        "mean_distance": 14.8,
        "std_distance": 1.2,
        "best_distance": 13,
        "worst_distance": 17,
        "success_rate": 0.95,
        "mean_execution_time": 44.1
      }
    }
  }
}
```

### **Métricas Coletadas**

- **Performance**: Tempo de execução, uso de memória, throughput
- **Qualidade**: Distância encontrada, taxa de convergência
- **Robustez**: Taxa de sucesso, estabilidade entre execuções
- **Recursos**: Uso de CPU, memória, I/O
- **Algoritmo-específicas**: Gerações, iterações, diversidade da população

## � Recursos Avançados

### **Sistema de Execução Inteligente**

#### **Scheduler com Fila FIFO**
- **Ordem Rigorosa**: Execução estritamente sequencial das tarefas
- **Controle de Recursos**: Monitoramento automático de CPU e memória
- **Balanceamento Dinâmico**: Ajuste automático do número de workers
- **Timeout Configurável**: Prevenção de execuções infinitas
- **Delay Inteligente**: Espaçamento entre execuções para estabilidade

#### **Execução Paralela Avançada**
- **Detecção Automática**: Identifica algoritmos com paralelismo interno
- **Configuração Inteligente**: Ajuste automático de workers internos/externos
- **Prevenção de Oversubscription**: Controle de recursos baseado em núcleos disponíveis
- **Execução Heterogênea**: Suporte a algoritmos determinísticos e estocásticos

### **Monitoramento e Logging**

#### **Interface Curses**
- **Monitoramento em Tempo Real**: Progresso visual das execuções
- **Métricas Live**: CPU, memória, tempo decorrido
- **Status de Tarefas**: Fila, execução, completadas
- **Alertas Visuais**: Warnings e erros destacados

#### **Logging Estruturado**
- **Rastreamento Completo**: Logs detalhados de todas as operações
- **Níveis Configuráveis**: DEBUG, INFO, WARNING, ERROR
- **Rotação Automática**: Gestão inteligente de arquivos de log
- **Análise Post-Mortem**: Investigação de falhas e problemas

### **Gestão de Datasets**

#### **Datasets Sintéticos**
```python
# Geração programática
from src.datasets.dataset_synthetic import generate_dataset

sequences, params = generate_dataset(
    n=100,           # Número de sequências
    L=200,           # Comprimento das sequências
    alphabet='ACGT', # Alfabeto
    noise=0.15,      # Nível de ruído
    seed=42          # Reprodutibilidade
)
```

#### **Datasets Reais**
```python
# Download automático do NCBI
from src.datasets.dataset_entrez import fetch_dataset

sequences, info = fetch_dataset(
    query="COVID-19 spike protein",
    database="protein",
    retmax=100
)
```

#### **Cache Inteligente**
- **Persistência Automática**: Datasets salvos automaticamente
- **Recuperação Rápida**: Reload instantâneo de datasets processados
- **Validação de Integridade**: Verificação de consistência dos dados
- **Gestão de Espaço**: Limpeza automática de cache antigo

### **Análise Estatística Avançada**

#### **Testes Estatísticos**
- **Normalidade**: Shapiro-Wilk, Anderson-Darling
- **Comparação**: t-test, Mann-Whitney U, Kruskal-Wallis
- **Correlação**: Pearson, Spearman
- **Regressão**: Linear, polinomial

#### **Visualizações**
- **Box Plots**: Distribuição de resultados
- **Scatter Plots**: Correlação entre métricas
- **Heatmaps**: Matriz de comparação entre algoritmos
- **Convergência**: Gráficos de progresso temporal

### **Configuração Flexível**

#### **Configuração em Lote (YAML)**
```yaml
# Configuração unificada
batch_info:
  nome: "Benchmark Comparativo"
  descricao: "Comparação de algoritmos CSP"
  timeout_global: 3600

datasets:
  - id: synthetic_small
    tipo: synthetic
    parametros:
      n: 20
      L: 50
      noise: 0.1

task:
  type: "execution"
  execution:
    executions:
      - nome: "Benchmark Principal"
        dataset: synthetic_small
        runs_per_algorithm_per_base: 10
        timeout: 300

algorithms: ["Baseline", "BLF-GA", "CSC", "H3-CSP"]
```

#### **Configuração Programática**
```python
# Configuração via Python
from src.core.config import CSPBenchConfig

config = CSPBenchConfig(
    algorithms=["BLF-GA", "CSC"],
    dataset_config={
        "type": "synthetic",
        "n": 50,
        "L": 100,
        "noise": 0.15
    },
    execution_config={
        "num_executions": 10,
        "timeout": 300,
        "parallel_workers": 4
    }
)
```

## 📚 Algoritmos Incluídos

O CSPBench inclui uma biblioteca de algoritmos CSP implementados e validados. Cada algoritmo possui documentação detalhada em seu diretório específico.

### **Algoritmos Disponíveis**

| Algoritmo | Tipo | Determinístico | Paralelismo | Uso Ideal |
|-----------|------|----------------|-------------|-----------|
| **Baseline** | Consenso Ganancioso | ✅ | ❌ | Baseline rápido |
| **BLF-GA** | Algoritmo Genético Híbrido | ❌ | ✅ | Instâncias médias/grandes |
| **CSC** | Clustering + Consenso | ❌ | ✅ | Dados estruturados |
| **DP-CSP** | Programação Dinâmica | ✅ | ❌ | Instâncias pequenas (ótimo) |
| **H³-CSP** | Busca Hierárquica Híbrida | ✅ | ❌ | Instâncias médias |

### **Detalhes dos Algoritmos**

Cada algoritmo possui documentação completa em seu diretório:

- **`algorithms/baseline/README.md`**: Algoritmo de consenso ganancioso determinístico
- **`algorithms/blf_ga/README.md`**: Blockwise Learning Fusion + Genetic Algorithm
- **`algorithms/csc/README.md`**: Consensus String Clustering
- **`algorithms/dp_csp/README.md`**: Programação Dinâmica Exata
- **`algorithms/h3_csp/README.md`**: Hybrid Hierarchical Hamming Search

### **Seleção de Algoritmos**

#### **Para Benchmarks Rápidos**
```bash
python main.py --algorithms Baseline
```

#### **Para Análise Completa**
```bash
python main.py --algorithms Baseline BLF-GA CSC H3-CSP
```

#### **Para Instâncias Pequenas com Solução Ótima**
```bash
python main.py --algorithms DP-CSP Baseline
```

#### **Para Instâncias Grandes**
```bash
python main.py --algorithms BLF-GA CSC --workers 8
```

## � Debugging e Troubleshooting

### **Logs e Diagnósticos**

#### **Ativando Logging Detalhado**
```bash
# Logging completo
python main.py --debug

# Verificar logs em tempo real
tail -f outputs/logs/$(ls -t outputs/logs/ | head -1)

# Análise de logs específicos
grep "ERROR\|WARNING" outputs/logs/20250709_120000_abcd1234.log
```

#### **Monitoramento de Recursos**
```bash
# Monitoramento em tempo real
python main.py --visual --resource-monitor

# Configuração de limites de recursos
python main.py --cpu-limit 80 --memory-limit 4096
```

### **Problemas Comuns e Soluções**

| Problema | Causa Provável | Solução |
|----------|----------------|---------|
| **Timeout Constante** | Algoritmo muito lento ou dataset grande | Reduzir dataset ou aumentar `--timeout` |
| **Memória Insuficiente** | Dataset muito grande ou leak de memória | Usar dataset menor ou `--memory-limit` |
| **Resultados Inconsistentes** | Algoritmo não-determinístico sem seed | Configurar seed fixo nos parâmetros |
| **Interface Travada** | Problema com terminal/curses | Usar `--silent` ou verificar terminal |
| **Algoritmo Não Encontrado** | Erro no registro ou importação | Verificar implementação e `@register_algorithm` |
| **Paralelismo Ineficiente** | Oversubscription de recursos | Ajustar `--workers` baseado em núcleos |

### **Validação da Instalação**

```bash
# Verificar algoritmos registrados
python -c "from algorithms.base import global_registry; print(list(global_registry.keys()))"

# Teste básico de funcionalidade
python -c "
from algorithms.baseline.algorithm import BaselineAlg
alg = BaselineAlg(['ACGT', 'AGCT'], 'ACGT')
print(alg.run())
"

# Verificar dependências
python -c "
import numpy, biopython, optuna, matplotlib
print('Todas as dependências OK')
"
```

### **Performance e Profiling**

#### **Benchmark de Performance**
```bash
# Benchmark automático
python benchmark/benchmark_parallel.py --verbose

# Profiling detalhado
python -m cProfile -o profile.stats main.py --silent --dataset synthetic --algorithms BLF-GA

# Análise de profiling
python -c "
import pstats
p = pstats.Stats('profile.stats')
p.sort_stats('cumulative').print_stats(20)
"
```

#### **Otimização de Performance**
- **Paralelismo**: Usar `--workers` próximo ao número de núcleos
- **Cache**: Reutilizar datasets salvos em `saved_datasets/`
- **Timeout**: Configurar timeouts realistas para evitar espera desnecessária
- **Logging**: Usar nível `INFO` em produção, `DEBUG` apenas para desenvolvimento

## 🤝 Contribuindo

### **Como Contribuir**

O CSPBench é um projeto open-source que aceita contribuições da comunidade científica.

#### **1. Preparação do Ambiente**
```bash
# Fork e clone do repositório
git fork https://github.com/diegogrosmann/CSPBench.git
git clone https://github.com/seu-usuario/CSPBench.git
cd CSPBench

# Configurar ambiente de desenvolvimento
python -m venv .venv
source .venv/bin/activate
pip install -e .[dev]

# Instalar hooks de pre-commit
pre-commit install
```

#### **2. Desenvolvimento**
```bash
# Criar branch para sua feature
git checkout -b feature/nova-funcionalidade

# Desenvolver e testar
python -m pytest tests/ -v
python -m ruff check .
python -m black --check .

# Commit e push
git commit -m 'feat: adiciona nova funcionalidade'
git push origin feature/nova-funcionalidade
```

#### **3. Pull Request**
- Abrir PR com descrição detalhada
- Incluir testes para novas funcionalidades
- Documentar mudanças no CHANGELOG.md
- Seguir guidelines de código

### **Tipos de Contribuições**

#### **🧬 Novos Algoritmos**
- Implementar novos algoritmos CSP
- Seguir interface `CSPAlgorithm`
- Incluir documentação completa
- Adicionar testes unitários

#### **📊 Datasets**
- Novos tipos de datasets
- Melhorias no download automático
- Formatos de arquivo adicionais
- Validação de dados

#### **🔧 Melhorias no Framework**
- Otimizações de performance
- Novos recursos de monitoramento
- Melhores visualizações
- Paralelização avançada

#### **📚 Documentação**
- Tutoriais e exemplos
- Documentação de APIs
- Guides de boas práticas
- Tradução de documentos

### **Diretrizes de Código**

#### **Estilo de Código**
- **Python**: PEP 8 (automatizado com `black`)
- **Docstrings**: Estilo Google
- **Type Hints**: Obrigatório para APIs públicas
- **Imports**: Organizados com `isort`

#### **Testes**
- **Cobertura**: Mínimo 80% para novo código
- **Framework**: pytest
- **Mocks**: Para dependências externas
- **CI/CD**: Todos os testes devem passar

#### **Documentação**
- **README**: Para cada módulo/algoritmo
- **Docstrings**: Para todas as funções públicas
- **Exemplos**: Código funcional nos docs
- **CHANGELOG**: Documentar todas as mudanças

### **Estrutura de Commits**

#### **Conventional Commits**
```bash
feat: adiciona novo algoritmo XYZ
fix: corrige bug no scheduler
docs: atualiza documentação da API
test: adiciona testes para CSC
refactor: melhora estrutura do módulo de datasets
perf: otimiza execução paralela
style: aplica formatação black
ci: atualiza pipeline de CI
```

#### **Checklist de PR**
- [ ] Código segue PEP 8
- [ ] Testes passam e cobertura mantida
- [ ] Documentação atualizada
- [ ] CHANGELOG.md atualizado
- [ ] Commit messages seguem convenção
- [ ] Performance não degradada
- [ ] Funcionalidade testada manualmente

## � Suporte e Comunidade

### **Contato**

- **Email**: diegogrosmann@gmail.com
- **GitHub Issues**: [Reportar problemas](https://github.com/diegogrosmann/CSPBench/issues)
- **Discussions**: [Fórum da comunidade](https://github.com/diegogrosmann/CSPBench/discussions)

### **Documentação Adicional**

- **API Reference**: `docs/api/`
- **Tutoriais**: `docs/tutorials/`
- **Developer Guide**: `docs/DEVELOPER_GUIDE.md`
- **Architecture Guide**: `docs/ARCHITECTURE.md`

### **Roadmap**

#### **Próximas Versões**

**v1.1.0 - Visualizações Avançadas**
- [ ] Dashboard web interativo
- [ ] Gráficos de convergência em tempo real
- [ ] Visualização de paisagem de fitness
- [ ] Comparação visual entre algoritmos

**v1.2.0 - Algoritmos Avançados**
- [ ] Algoritmos evolutivos multi-objetivo
- [ ] Algoritmos de enxame (PSO, ACO)
- [ ] Machine Learning para CSP
- [ ] Algoritmos aproximados

**v1.3.0 - Infraestrutura Distribuída**
- [ ] Execução em cluster
- [ ] API REST para execução remota
- [ ] Integração com Kubernetes
- [ ] Cache distribuído

### **Histórico de Versões**

#### **v1.0.0 - Versão Inicial**
- ✅ Framework base implementado
- ✅ 5 algoritmos CSP incluídos
- ✅ Sistema de execução paralela
- ✅ Integração Optuna/SALib
- ✅ Interface curses
- ✅ Relatórios automatizados

## 📄 Licença

Este projeto está licenciado sob a **MIT License** - veja o arquivo [LICENSE](LICENSE) para detalhes.

```
MIT License

Copyright (c) 2025 Diego Grosmann

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## 🙏 Agradecimentos

- **Comunidade Python**: Pelas excelentes bibliotecas científicas
- **Pesquisadores em Bioinformática**: Pela inspiração e validação científica
- **Contribuidores Open Source**: Por melhorar continuamente o projeto
- **Universidades e Institutos**: Pelo suporte à pesquisa científica

---

**CSPBench** - Framework experimental robusto para teste e análise de algoritmos do Closest String Problem.

*Desenvolvido com ❤️ para a comunidade científica*


