# Sistema de Batch Unificado CSP-BLFGA

## Visão Geral

O sistema de batch unificado permite executar três tipos diferentes de tarefas usando uma estrutura padronizada de configuração YAML:

1. **Execução (execution)**: Executa algoritmos em datasets
2. **Otimização (optimization)**: Otimiza hiperparâmetros usando Optuna
3. **Análise de Sensibilidade (sensitivity)**: Analisa sensibilidade de parâmetros

## Uso

### Via linha de comando

```bash
# Executar batch diretamente
python main.py --batch batch_configs/exemplo_unificado.yaml --silent

# Executar em modo interativo
python main.py --batch batch_configs/exemplo_unificado.yaml

# Executar via interface
python main.py
# Escolher opção 4 (Execução em lote unificada)
```

### Estrutura do arquivo de configuração

```yaml
# =====================================================================
# ESTRUTURA PADRONIZADA CSP-BLFGA (v1.3)
# =====================================================================

batch_info:
  nome: "Nome do batch"
  descricao: "Descrição do objetivo"
  autor: "Autor"
  versao: "1.3"
  timeout_global: 1800          # Timeout global em segundos

# ---------------------------------------------------------------------
# DATASETS - Definição dos datasets utilizados
# ---------------------------------------------------------------------
datasets:
  - id: dataset_1
    nome: "Nome do Dataset"
    tipo: "synthetic"           # synthetic | file | entrez
    parametros:
      # Para synthetic
      n: 8
      L: 30
      alphabet: "ACGT"
      noise: 0.10
      fully_random: false
      seed: 12345
      
      # Para file
      # filename: "caminho/para/arquivo.fasta"
      
      # Para entrez
      # query: "termo_de_busca"
      # db: "nucleotide"
      # retmax: 100

# ---------------------------------------------------------------------
# TAREFA - Escolha UM tipo por arquivo
# ---------------------------------------------------------------------
task:
  type: "execution"             # execution | sensitivity | optimization

  # ===================== EXECUTION ====================================
  execution:
    executions:
      - nome: "Nome da Execução"
        dataset: dataset_1      # ID do dataset
        runs_per_algorithm_per_base: 3
        num_bases: 1
        timeout: 60

  # ===================== SENSITIVITY ==================================
  sensitivity:
    analyses:
      - nome: "Nome da Análise"
        datasets: [dataset_1]   # Lista de IDs de datasets
        n_samples: 100
        timeout_per_sample: 30
        method: "morris"        # morris | sobol | fast
        param_space:
          "BLF-GA": ["pop_size", "max_gens", "cross_prob"]

  # ===================== OPTIMIZATION ================================
  optimization:
    studies:
      - nome: "Nome do Estudo"
        datasets: [dataset_1]
        n_trials: 50
        timeout_per_trial: 60
        direction: "minimize"
        sampler: "TPE"          # TPE | Random | Grid
        param_space:
          "BLF-GA":
            pop_size: ["loguniform", 50, 400]
            max_gens: ["int", 100, 500]

# ---------------------------------------------------------------------
# ALGORITMOS E CONFIGURAÇÕES
# ---------------------------------------------------------------------
algorithms: ["Baseline", "BLF-GA", "CSC"]

algorithm_params:
  "BLF-GA":
    tournament_k: 3

# ---------------------------------------------------------------------
# SAÍDA
# ---------------------------------------------------------------------
output:
  save_results: true
  save_detailed_results: true
  save_plots: true
  plot_format: "png"
  results_dir: "outputs/batch_results"

# ---------------------------------------------------------------------
# OPÇÕES GLOBAIS
# ---------------------------------------------------------------------
advanced:
  use_curses: true
  parallel:
    n_jobs: -1
  log_level: "INFO"
```

## Arquivos de Exemplo

### Execução de Algoritmos
- `batch_configs/unificado_processamento.yaml`
- Executa algoritmos em múltiplos datasets

### Otimização de Hiperparâmetros
- `batch_configs/unificado_otimizacao.yaml`
- Otimiza parâmetros usando Optuna

### Análise de Sensibilidade
- `batch_configs/unificado_sensibilidade.yaml`
- Analisa sensibilidade de parâmetros

## Componentes

### BatchConfigExtractor
- Carrega e valida configurações YAML
- Extrai informações para cada tipo de tarefa
- Gera datasets baseado nas configurações

### UnifiedBatchProcessor
- Processa todas as tarefas de forma unificada
- Suporte a execução com/sem interface curses
- Gerencia timeouts e paralelismo

### Menu Unificado
- Interface única para seleção de arquivos
- Suporte a arquivos personalizados
- Integração com sistema legado

## Vantagens

1. **Padronização**: Estrutura única para todos os tipos de tarefa
2. **Flexibilidade**: Suporte a múltiplos datasets e algoritmos
3. **Facilidade**: Execução via linha de comando ou interface
4. **Extensibilidade**: Fácil adição de novos tipos de tarefa
5. **Compatibilidade**: Mantém compatibilidade com sistema legado

## Migração

Para migrar configurações antigas:

1. Copie a estrutura do arquivo exemplo
2. Adapte as configurações específicas
3. Defina o tipo de tarefa apropriado
4. Ajuste os parâmetros conforme necessário

O sistema legado (opções 5-6 no menu) ainda funciona para compatibilidade.
