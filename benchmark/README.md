# Benchmarks CSP-BLFGA

Esta pasta contém scripts de benchmark para testar performance e paralelização do sistema CSP-BLFGA.

## Scripts Disponíveis

### benchmark_parallel.py
Benchmark rápido para testar paralelização do sistema de otimização e análise de sensibilidade.

**Uso:**
```bash
# Executar do diretório raiz do projeto
cd /home/diego_grosmann/csp-blfga

# Benchmark completo
python benchmark/benchmark_parallel.py --verbose

# Apenas otimização Optuna
python benchmark/benchmark_parallel.py --skip-salib --verbose

# Apenas análise de sensibilidade SALib
python benchmark/benchmark_parallel.py --skip-optuna --verbose

# Dataset maior para teste mais realista
python benchmark/benchmark_parallel.py --dataset-size medium --verbose
```

**Opções:**
- `--verbose, -v`: Saída detalhada
- `--skip-optuna`: Pular benchmark Optuna
- `--skip-salib`: Pular benchmark SALib
- `--dataset-size {small,medium}`: Tamanho do dataset de teste

**Meta de Performance:**
- Speedup ≥ 2x com paralelização
- Medição de tempo serial vs paralelo
- Relatório detalhado de performance

## Resultados Esperados

O benchmark deve mostrar:
- Tempo de execução serial vs paralelo
- Speedup calculado
- Informações do sistema (CPUs, Python)
- Status se meta foi alcançada

### Exemplo de Saída
```
🚀 Iniciando benchmark de paralelização
============================================================
🖥️  Sistema: 8 CPUs, Python 3.12.0
📊 Criando dataset de teste (small)...
   Dataset: 8 sequências, tamanho 25

📊 Benchmark Optuna (15 trials)
  🔄 Executando modo serial...
    ✅ Serial: 45.23s, melhor=12.00
  🚀 Executando modo paralelo (4 jobs)...
    ✅ Paralelo: 18.45s, melhor=11.50

============================================================
📈 RESULTADOS DO BENCHMARK
============================================================
🖥️  Sistema:
   CPUs: 8
   Python: 3.12.0

🔧 Optuna (Otimização):
   Serial:   45.23s
   Paralelo: 18.45s
   Speedup:  2.45x
   ✅ Meta alcançada (≥2x)

📊 Speedup médio: 2.45x
🎉 Paralelização bem-sucedida!
============================================================
```

## Requisitos

- Python 3.8+
- Todas as dependências do projeto CSP-BLFGA
- Sistema multi-core para testar paralelização
- Pelo menos 4 CPUs para resultados significativos

## Troubleshooting

### Erro de Import
Se houver erros de import, certifique-se de executar do diretório raiz:
```bash
cd /home/diego_grosmann/csp-blfga
python benchmark/benchmark_parallel.py
```

### Speedup Baixo
Se o speedup for menor que 2x:
- Verifique se há gargalos de I/O
- Aumente o tamanho do dataset (`--dataset-size medium`)
- Verifique se outros processos estão consumindo CPU

### Timeout
Se houver timeouts:
- Reduza o número de trials/amostras
- Aumente o timeout no código
- Use dataset menor (`--dataset-size small`)

## Desenvolvimento

Para adicionar novos benchmarks:
1. Crie um novo arquivo na pasta benchmark/
2. Siga o padrão de estrutura do benchmark_parallel.py
3. Documente no README.md
4. Teste com diferentes tamanhos de dataset
