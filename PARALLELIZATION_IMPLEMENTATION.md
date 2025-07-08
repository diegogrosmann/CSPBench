# Paralelização do Sistema de Otimização CSP-BLFGA

## Resumo das Implementações

Este documento descreve as melhorias implementadas para paralelizar o sistema de otimização e análise de sensibilidade do CSP-BLFGA.

### 1. Expandir arquivo de configuração YAML ✅

**Arquivos modificados:**
- `batch_configs/otimizacao_exemplo.yaml`
- `batch_configs/otimizacao_completa.yaml`
- `batch_configs/sensibilidade_exemplo.yaml`
- `batch_configs/sensibilidade_teste_rapido.yaml`

**Adicionado:**
```yaml
# Em optimization_config:
parallel:
  n_jobs: 4              # Número de jobs paralelos para Optuna
  storage: null          # Storage para Optuna (ex: "sqlite:///optuna.db")
  study_name: null       # Nome do estudo

# Em sensitivity_config:
parallel:
  n_jobs: 4              # Número de jobs paralelos para SALib
```

### 2. Paralelizar o módulo OptunaOptimizer ✅

**Arquivos modificados:**
- `src/optimization/optuna_optimizer.py`
- `src/optimization/batch_optimizer.py`

**Principais mudanças:**
- Adicionado campos `n_jobs` e `internal_workers` em `OptimizationConfig`
- Implementado suporte a storage seguro para multiprocessamento
- Configuração automática da variável `INTERNAL_WORKERS` para algoritmos
- Chamada `study.optimize()` com parâmetro `n_jobs`
- Integração com `worker_calculator` para cálculo dinâmico de workers

### 3. Paralelizar o módulo SensitivityAnalyzer ✅

**Arquivos modificados:**
- `src/optimization/sensitivity_analyzer.py`
- `src/optimization/batch_sensitivity.py`

**Principais mudanças:**
- Adicionado campo `n_jobs` em `SensitivityConfig`
- Implementado `ProcessPoolExecutor` para avaliação paralela de amostras
- Separação entre execução serial e paralela
- Preservação da ordem dos resultados
- Manutenção da barra de progresso mesmo no modo paralelo

### 4. Centralizar cálculo de workers ✅

**Arquivo criado:**
- `src/utils/worker_calculator.py`

**Funcionalidades:**
- `get_cpu_count()`: Detecção automática de CPUs disponíveis
- `calculate_optuna_workers()`: Cálculo de workers para Optuna
- `calculate_salib_workers()`: Cálculo de workers para SALib
- `calculate_internal_workers()`: Cálculo de workers internos para algoritmos
- `get_worker_config()`: Configuração completa baseada em contexto

**Heurísticas implementadas:**
- Optuna: 75% dos CPUs disponíveis
- SALib: Limitado por número de amostras (min 10 amostras/worker)
- Workers internos: Ajustados baseado em paralelismo externo

### 5. Criar script de benchmark rápido ✅

**Arquivo criado:**
- `benchmark/benchmark_parallel.py`

**Funcionalidades:**
- Benchmark comparativo entre execução serial e paralela
- Medição de speedup (meta: ≥ 2x)
- Testes para Optuna e SALib
- Datasets de teste configuráveis (small/medium)
- Relatório detalhado de performance

**Uso:**
```bash
# Executar do diretório raiz do projeto
python benchmark/benchmark_parallel.py --verbose
python benchmark/benchmark_parallel.py --skip-optuna  # Apenas SALib
python benchmark/benchmark_parallel.py --skip-salib   # Apenas Optuna
python benchmark/benchmark_parallel.py --dataset-size medium
```

### 6. Limpar código legado ✅

**Ações realizadas:**
- Atualização das docstrings para mencionar paralelização
- Integração consistente do `yaml_config` em todas as funções
- Remoção de comentários sobre paralelização futura
- Padronização dos imports e estruturas

## Arquitetura da Paralelização

### Fluxo Optuna
```
YAML Config → WorkerCalculator → OptimizationConfig → OptunaOptimizer
     ↓              ↓                    ↓                  ↓
  n_jobs      optuna_workers      configuração      study.optimize()
                internal_workers   INTERNAL_WORKERS      (n_jobs)
```

### Fluxo SALib
```
YAML Config → WorkerCalculator → SensitivityConfig → SensitivityAnalyzer
     ↓              ↓                    ↓                    ↓
  n_jobs      salib_workers         configuração      ProcessPoolExecutor
                                                          (max_workers)
```

### Controle de Recursos
```
Total CPUs
    ↓
WorkerCalculator
    ↓
┌─────────────────┬─────────────────┐
│   External      │   Internal      │
│   Workers       │   Workers       │
├─────────────────┼─────────────────┤
│ Optuna: 75%     │ BLF-GA: dinâmic │
│ SALib: 100%     │ CSC: dinâmico   │
│                 │ Outros: 1       │
└─────────────────┴─────────────────┘
```

## Benefícios Esperados

1. **Performance**: Speedup ≥ 2x em sistemas multi-core
2. **Escalabilidade**: Uso eficiente de recursos em servidores
3. **Flexibilidade**: Configuração através de YAML
4. **Robustez**: Prevenção de oversubscription
5. **Monitoramento**: Logging detalhado de execução

## Configuração Recomendada

### Para desenvolvimento (4 CPUs):
```yaml
optimization_config:
  parallel:
    n_jobs: 3
    storage: "sqlite:///outputs/optuna_dev.db"

sensitivity_config:
  parallel:
    n_jobs: 4
```

### Para produção (16+ CPUs):
```yaml
optimization_config:
  parallel:
    n_jobs: 12
    storage: "sqlite:///outputs/optuna_prod.db"

sensitivity_config:
  parallel:
    n_jobs: 16
```

## Testes de Validação

Execute o benchmark para validar a implementação:

```bash
# Teste rápido
python benchmark/benchmark_parallel.py

# Teste completo
python benchmark/benchmark_parallel.py --verbose --dataset-size medium

# Teste apenas otimização
python benchmark/benchmark_parallel.py --skip-salib --verbose
```

## Notas Técnicas

### Limitações
- ProcessPoolExecutor tem overhead de inicialização
- SQLite storage tem limitações de concorrência
- Algoritmos determinísticos não se beneficiam de paralelismo interno

### Recomendações
- Use SQLite storage para estudos Optuna persistentes
- Configure `n_jobs: -1` para usar todos os CPUs
- Monitore uso de memória em datasets grandes
- Ajuste `timeout_per_trial` para execuções paralelas

### Troubleshooting
- Se speedup < 2x: Verifique se há gargalos de I/O
- Se erro de storage: Certifique-se que diretório existe
- Se timeout: Aumente `timeout_per_trial` no YAML
