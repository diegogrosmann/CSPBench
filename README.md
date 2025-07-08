# CSP-BLFGA: Plataforma Experimental para o Closest String Problem

Este projeto implementa uma **plataforma experimental completa e robusta** para resolver o Closest String Problem (CSP) usando diferentes algoritmos, com foco especial no **BLF-GA (Blockwise Learning Fusion + Genetic Algorithm)**.

## 🎯 Visão Geral

A aplicação é uma ferramenta científica avançada que permite:

### 📊 **Gestão de Datasets**
- **Geração sintética**: Criação de datasets parametrizáveis com ruído controlado
- **Carregamento de arquivos**: Suporte a formatos FASTA e texto
- **Download automático**: Integração com NCBI/Entrez para dados reais
- **Execução em lote**: Processamento de múltiplos datasets simultaneamente

### 🧬 **Algoritmos Implementados**
- **Baseline**: Algoritmo de consenso ganancioso (determinístico)
- **BLF-GA**: Algoritmo genético híbrido com aprendizado por blocos
- **CSC**: Constraint Satisfaction with Clustering
- **DP-CSP**: Programação Dinâmica para CSP
- **H3-CSP**: Heurística H3 para CSP
- **Sistema Extensível**: Fácil adição de novos algoritmos

### 🚀 **Execução e Monitoramento**
- **Scheduler Avançado**: Controle inteligente de recursos (CPU/memória)
- **Execução Paralela**: Múltiplas execuções simultâneas com balanceamento
- **Monitoramento Visual**: Interface curses com progresso em tempo real
- **Controle de Timeout**: Prevenção de execuções infinitas
- **Logging Detalhado**: Rastreamento completo de execuções

### 📈 **Análise e Relatórios**
- **Relatórios Automáticos**: Geração de relatórios detalhados em JSON/CSV
- **Análise Comparativa**: Comparação automática entre algoritmos
- **Visualizações**: Gráficos de performance e convergência
- **Exportação Flexível**: Múltiplos formatos de saída

## 🔄 Fluxo Principal da Aplicação

### 1. **Preparação dos Dados**
```
Dataset → Validação → Normalização → Estruturas de Dados
```

### 2. **Configuração de Execução**
```
Algoritmos → Parâmetros → Recursos → Scheduler
```

### 3. **Execução Controlada**
```
Fila FIFO → Controle de Recursos → Execução Paralela → Monitoramento
```

### 4. **Análise e Relatórios**
```
Coleta de Resultados → Análise Estatística → Relatórios → Visualizações
```

## 🏗️ Arquitetura do Sistema

### **Componentes Principais**

```
┌─────────────────────────────────────────────────────────────────┐
│                        CSP-BLFGA Platform                       │
├─────────────────────────────────────────────────────────────────┤
│  🖥️  USER INTERFACE (UI)                                        │
│  ├── CLI (Command Line Interface)                               │
│  │   ├── app.py           # Orquestrador principal             │
│  │   ├── menu.py          # Menus interativos                  │
│  │   └── console_manager.py # Gerenciamento thread-safe        │
│  └── curses_integration.py # Monitoramento visual em tempo real│
├─────────────────────────────────────────────────────────────────┤
│  ⚙️  CORE SYSTEM                                                │
│  ├── 🎯 Scheduler (Agendamento Inteligente)                     │
│  │   ├── ExecutionScheduler  # Fila FIFO + Controle recursos   │
│  │   ├── ResourceMonitor    # Monitoramento CPU/Memória        │
│  │   └── ProcessWatcher     # Controle de processos filhos     │
│  ├── 🔌 Interfaces (Contratos do Sistema)                       │
│  │   ├── IAlgorithm        # Interface para algoritmos         │
│  │   ├── IExecutor         # Interface para executores         │
│  │   └── IConsole          # Interface para console            │
│  ├── 📊 Data Management                                         │
│  │   ├── TaskResult        # Estrutura de resultados           │
│  │   └── TaskHandle        # Controle de tarefas               │
│  └── 📈 Reporting                                               │
│      ├── ResultsFormatter  # Formatação de resultados          │
│      └── CSPExporter       # Exportação para CSV/JSON          │
├─────────────────────────────────────────────────────────────────┤
│  🧬 ALGORITHMS (Algoritmos CSP)                                 │
│  ├── Baseline         # Consenso ganancioso determinístico     │
│  ├── BLF-GA          # Algoritmo genético híbrido              │
│  ├── CSC             # Constraint Satisfaction Clustering      │
│  ├── DP-CSP          # Programação Dinâmica                    │
│  ├── H3-CSP          # Heurística H3                           │
│  └── [Extensível]    # Sistema de registro automático          │
├─────────────────────────────────────────────────────────────────┤
│  📁 DATASETS (Gestão de Dados)                                  │
│  ├── Synthetic       # Geração sintética parametrizável        │
│  ├── File Loader     # Carregamento de arquivos FASTA          │
│  ├── NCBI/Entrez     # Download automático de dados reais      │
│  └── Batch Processing # Processamento em lote                  │
├─────────────────────────────────────────────────────────────────┤
│  🔧 UTILITIES (Utilitários)                                     │
│  ├── Logging         # Sistema de logging padronizado          │
│  ├── Config          # Configurações centralizadas             │
│  ├── Distance        # Funções de distância (Hamming, etc.)    │
│  └── Resource Monitor # Monitoramento de recursos do sistema   │
└─────────────────────────────────────────────────────────────────┘
```

### **Fluxo de Dados**

```
📊 Dataset → 🎯 Scheduler → 🧬 Algoritmos → 📈 Resultados → 💾 Relatórios
    ↓             ↓             ↓             ↓             ↓
Validação    Fila FIFO    Execução     Coleta de     Análise
Normalização  Recursos    Paralela     Métricas     Estatística
```

algorithms/                   # 🧬 Implementações dos algoritmos
├── baseline/                # Algoritmo de consenso ganancioso
├── blf_ga/                  # BLF-GA: Blockwise Learning Fusion + GA
├── csc/                     # CSC: Consensus String Clustering
├── dp_csp/                  # DP-CSP: Programação dinâmica exata
├── h3_csp/                  # H³-CSP: Híbrido hierárquico
└── README.md               # Guia para adicionar novos algoritmos

datasets/                     # 📊 Gerenciamento de datasets
├── dataset_file.py          # Leitura de arquivos
├── dataset_entrez.py        # Download do NCBI
├── dataset_synthetic.py     # Geração sintética
└── dataset_utils.py         # Utilitários

tests/                        # 🧪 Testes automatizados
main.py                      # 🚀 Ponto de entrada principal
outputs/                     # 📈 Saídas organizadas
  ├── reports/               # Relatórios gerados
  └── logs/                  # Logs do sistema
logs/                        # 📝 Logs de execução
saved_datasets/              # 💾 Datasets salvos
batch_configs/               # ⚙️ Configurações de lote
```

## Como Executar

### 1. Instale o Python 3.10+ (recomendado: 3.10 ou 3.11)

No Ubuntu/Debian:
```bash
sudo apt update
sudo apt install python3.10 python3.10-venv python3.10-dev
```
No Windows:
Baixe em https://www.python.org/downloads/

### 2. Crie um ambiente virtual

No terminal (Linux/macOS):
```bash
python3.10 -m venv .venv
source .venv/bin/activate
```
No Windows (cmd):
```cmd
python -m venv .venv
.venv\Scripts\activate
```

### 3. Instale as dependências

```bash
pip install --upgrade pip
pip install -r requirements.txt
```

### 4. Execute a aplicação

```bash
python main.py
```

- Siga o menu interativo para escolher dataset, algoritmos e parâmetros.
- Para execução em lote, escolha a opção correspondente e selecione um arquivo YAML/XML em `batch_configs/`.

## Adicionando Novos Algoritmos

1. Crie uma nova pasta em `algorithms/` seguindo o padrão:
    ```
    algorithms/
      meu_algoritmo/
        __init__.py
        algorithm.py
        config.py
        implementation.py
    ```
2. Implemente a interface `Algorithm` e use o decorador `@register_algorithm`.
3. O algoritmo aparecerá automaticamente no menu, sem necessidade de alterar o main.

Veja `algorithms/README.md` para exemplos detalhados.

## Requisitos

- Python 3.10+
- Biopython
- scikit-learn
- tabulate

Instale dependências com:
```bash
pip install -r requirements.txt
```

## Relatórios e Resultados

- Relatórios detalhados são salvos em `outputs/reports/` após cada execução.
- Resumos rápidos são exibidos no console.
- Execuções em lote geram relatórios consolidados.

## Interface Gráfica (Futuro)

O projeto está preparado para interface gráfica que será implementada futuramente:

```python
from src.ui.widgets import run_gui
# run_gui()  # Será implementado em versões futuras
```

## Suporte e Documentação

- Cada algoritmo possui README próprio explicando heurísticas, funcionamento e parâmetros.
- Consulte `REESTRUTURACAO.md` para detalhes sobre a arquitetura e modularização.

---

### Observações sobre o main.py

O arquivo `main.py` está totalmente documentado com docstrings no estilo Google, detalhando o fluxo, parâmetros e retornos de cada função. Consulte o código para detalhes de uso programático e integração.

## Uso

### Execução Básica

```bash
python main.py
```

### Parâmetros da CLI

```bash
python main.py --help
```

Principais opções:
- `--silent`: Modo silencioso (sem interação)
- `--dataset {synthetic,file,entrez,batch}`: Tipo de dataset
- `--algorithms ALGS [ALGS ...]`: Algoritmos a executar
- `--num-execs N`: Número de execuções por algoritmo
- `--timeout N`: Timeout por execução (segundos)
- `--workers N` ou `-w N`: Número de workers paralelos (padrão: 4)

### Configuração de Paralelismo

O sistema detecta automaticamente algoritmos que suportam paralelismo interno:
- **Algoritmos sem paralelismo interno**: Usa múltiplos workers externos
- **Algoritmos com paralelismo interno**: Usa 1 worker externo e configura workers internos

Exemplo com múltiplos algoritmos:
```bash
python main.py --dataset synthetic --algorithms Baseline CSC BLF-GA --workers 4
```

A variável de ambiente `INTERNAL_WORKERS` é configurada automaticamente baseada no número de CPUs disponíveis.

## 🚀 Características Avançadas

### **Scheduler Inteligente**
- **Fila FIFO Absoluta**: Ordem de execução rigorosamente respeitada
- **Controle de Recursos**: Monitoramento automático de CPU e memória
- **Balanceamento Dinâmico**: Ajuste automático do número de workers
- **Timeout Configurável**: Prevenção de execuções infinitas
- **Delay Inteligente**: Espaçamento entre execuções para estabilidade

### **Sistema de Monitoramento**
- **Interface Curses**: Monitoramento visual em tempo real
- **Logging Estruturado**: Rastreamento detalhado de todas as operações
- **Métricas de Performance**: Coleta automática de estatísticas
- **Detecção de Anomalias**: Identificação de problemas durante execução

### **Extensibilidade**
- **Sistema de Registro**: Algoritmos registrados automaticamente
- **Interfaces Padronizadas**: Contratos bem definidos para todos os componentes
- **Arquitetura Modular**: Componentes independentes e reutilizáveis
- **Configuração Flexível**: Parâmetros ajustáveis em tempo de execução

## 📋 Pré-requisitos

### **Sistema Operacional**
- Linux (recomendado) / macOS / Windows
- Python 3.8+ (testado com Python 3.12)
- Terminal com suporte a cores (para interface curses)

### **Dependências Principais**
- **NumPy**: Computação numérica eficiente
- **Biopython**: Processamento de sequências biológicas
- **Optuna**: Otimização de hiperparâmetros
- **Matplotlib**: Visualização de dados
- **Rich**: Interface de terminal rica

## 🛠️ Instalação

### **Método 1: Ambiente Virtual (Recomendado)**
```bash
# Clonar o repositório
git clone https://github.com/diegogrosmann/CSP.git
cd CSP

# Criar ambiente virtual
python -m venv .venv

# Ativar ambiente virtual
source .venv/bin/activate  # Linux/macOS
# ou
.venv\Scripts\activate     # Windows

# Instalar dependências
pip install -r requirements.txt

# Instalar dependências de desenvolvimento (opcional)
pip install -e .[dev]
```

### **Método 2: Instalação Direta**
```bash
git clone https://github.com/seu-usuario/csp-blfga.git
cd csp-blfga
pip install -r requirements.txt
```

## 🎮 Como Usar

### **Execução Interativa**
```bash
# Executar com interface completa
python main.py

# Executar com monitoramento visual
python main.py --visual

# Executar com logging detalhado
python main.py --debug
```

### **Execução Automatizada**
```bash
# Execução silenciosa para testes
python main.py --silent --dataset synthetic --algorithms Baseline BLF-GA --num-execs 5

# Execução em lote com arquivo de configuração
python main.py --batch batch_configs/otimizacao_completa.yaml

# Execução com timeout personalizado
python main.py --timeout 600 --dataset file --algorithms BLF-GA
```

### **Otimização de Hiperparâmetros**
```bash
# Otimização automática com Optuna
python main.py --optimize --algorithm BLF-GA --dataset synthetic --trials 100

# Análise de sensibilidade
python main.py --sensitivity --algorithm BLF-GA --param population_size --range 50,200
```

## 📊 Exemplos de Uso

### **Exemplo 1: Comparação Simples**
```python
from main import main
import sys

# Configurar argumentos
sys.argv = ['main.py', '--dataset', 'synthetic', 
           '--algorithms', 'Baseline', 'BLF-GA', 
           '--num-execs', '10']

# Executar
main()
```

### **Exemplo 2: Dataset Customizado**
```python
from src.datasets.dataset_synthetic import generate_dataset
from algorithms.blf_ga.algorithm import BLFGAAlgorithm

# Gerar dataset
sequences, params = generate_dataset(
    n=50, L=200, alphabet='ACGT', noise=0.15
)

# Executar algoritmo
algorithm = BLFGAAlgorithm(sequences, 'ACGT', 
                          population_size=100, 
                          generations=500)
center, distance, metadata = algorithm.run()
```

### **Exemplo 3: Execução em Lote**
```yaml
# batch_config.yaml
datasets:
  - type: synthetic
    params:
      n: 30
      L: 100
      noise: 0.1
  - type: file
    params:
      filepath: "data/sequences.fasta"

algorithms:
  - name: Baseline
  - name: BLF-GA
    params:
      population_size: 100
      generations: 300

execution:
  num_execs: 10
  timeout: 300
  visual: true
```

## 🔧 Adicionando Novos Algoritmos

Para adicionar um novo algoritmo ao sistema:

### **1. Criar a Classe do Algoritmo**
```python
from algorithms.base import CSPAlgorithm, register_algorithm

@register_algorithm
class MeuAlgoritmo(CSPAlgorithm):
    name = "MeuAlgoritmo"
    default_params = {
        'parametro1': 10,
        'parametro2': 0.5
    }
    is_deterministic = False
    supports_internal_parallel = True
    
    def __init__(self, strings, alphabet, **params):
        super().__init__(strings, alphabet, **params)
        # Inicialização específica
    
    def run(self):
        # Implementação do algoritmo
        center = self.meu_algoritmo_core()
        distance = self.calculate_distance(center)
        metadata = self.get_metadata()
        
        return center, distance, metadata
```

### **2. Registrar o Algoritmo**
```python
# O decorador @register_algorithm já registra automaticamente
# Não é necessário código adicional!
```

### **3. Usar o Novo Algoritmo**
```bash
python main.py --algorithms MeuAlgoritmo
```

## 📈 Análise de Resultados

### **Relatórios Automáticos**
O sistema gera automaticamente:
- **Relatórios JSON**: Dados estruturados completos
- **Relatórios CSV**: Dados tabulares para análise
- **Logs Detalhados**: Rastreamento completo de execuções
- **Gráficos**: Visualizações de performance (quando disponível)

### **Métricas Coletadas**
- **Performance**: Tempo de execução, memória utilizada
- **Qualidade**: Distância encontrada, convergência
- **Robustez**: Taxa de sucesso, estabilidade
- **Recursos**: Uso de CPU, memória, I/O

### **Estrutura de Dados dos Resultados**
```json
{
  "algorithm": "BLF-GA",
  "dataset": "synthetic_n30_L100",
  "executions": [
    {
      "execution_id": 1,
      "center": "ACGTACGTACGT...",
      "distance": 15,
      "execution_time": 45.2,
      "memory_used": 128.5,
      "metadata": {
        "iterations": 245,
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
    "success_rate": 0.95
  }
}
```

## 🔍 Debugging e Troubleshooting

### **Logs Detalhados**
```bash
# Ativar logging detalhado
python main.py --debug

# Verificar logs
tail -f outputs/logs/20250708_120000_abcd1234.log
```

### **Monitoramento de Recursos**
```bash
# Verificar uso de recursos em tempo real
python main.py --visual --resource-monitor

# Ajustar limites de recursos
python main.py --cpu-limit 80 --memory-limit 2048
```

### **Problemas Comuns**

| Problema | Causa | Solução |
|----------|-------|---------|
| Timeout constante | Algoritmo muito lento | Reduzir tamanho do dataset ou aumentar timeout |
| Memória insuficiente | Dataset muito grande | Usar dataset menor ou aumentar limite de memória |
| Resultados inconsistentes | Algoritmo não-determinístico | Usar seed fixo ou algoritmo determinístico |
| Interface travada | Problema com curses | Usar modo `--silent` ou verificar terminal |

## 🤝 Contribuindo

### **Como Contribuir**
1. **Fork** do repositório
2. **Criar branch** para sua feature: `git checkout -b feature/nova-feature`
3. **Commit** suas mudanças: `git commit -m 'Adiciona nova feature'`
4. **Push** para a branch: `git push origin feature/nova-feature`
5. **Abrir Pull Request**

### **Diretrizes de Contribuição**
- **Código**: Seguir PEP 8 (usar `black` para formatação)
- **Testes**: Adicionar testes para novas funcionalidades
- **Documentação**: Documentar todas as funções e classes
- **Commits**: Mensagens claras e descritivas

### **Estrutura de Desenvolvimento**
```bash
# Executar testes
python -m pytest tests/ -v

# Verificar qualidade do código
python -m ruff check .
python -m black --check .

# Executar análise de cobertura
python -m pytest tests/ --cov=src --cov-report=html
```

## 📄 Licença

Este projeto está licenciado sob a [MIT License](LICENSE).

## 🙏 Agradecimentos

- **Comunidade Python**: Pelas excelentes bibliotecas
- **Pesquisadores em Bioinformática**: Pela inspiração e dados
- **Contribuidores**: Por melhorar continuamente o projeto

## 📞 Contato

Para dúvidas, sugestões ou colaborações:
- **Email**: diegogrosmann@gmail.com
- **GitHub**: [Abrir uma issue](https://github.com/seu-usuario/csp-blfga/issues)

---

**CSP-BLFGA** - Uma plataforma robusta e extensível para experimentação com algoritmos do Closest String Problem.
