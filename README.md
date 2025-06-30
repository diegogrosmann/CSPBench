# Closest String Problem (CSP) - Plataforma Experimental

Este projeto implementa uma arquitetura modular e extensível para experimentação com algoritmos de solução do Closest String Problem (CSP), incluindo baseline, heurísticas avançadas e métodos exatos.

## Visão Geral

A aplicação permite:
- Geração, leitura e download de datasets (sintéticos, arquivos, NCBI)
- Execução interativa ou em lote de múltiplos algoritmos
- Relatórios detalhados e comparativos automáticos
- Adição de novos algoritmos sem modificar o código principal

## Estrutura do Projeto

```
algorithms/
    baseline/         # Algoritmo de consenso ganancioso (baseline)
    blf_ga/           # BLF-GA: Blockwise Learning Fusion + GA
    csc/              # CSC: Consensus String Clustering
    dp_csp/           # DP-CSP: Programação dinâmica exata
    h3_csp/           # H³-CSP: Híbrido hierárquico
    README.md         # Guia para adicionar novos algoritmos

datasets/
    dataset_file.py      # Leitura de datasets de arquivo
    dataset_entrez.py    # Download de datasets do NCBI
    dataset_synthetic.py # Geração de datasets sintéticos
    dataset_utils.py     # Utilitários para datasets

utils/
    config.py        # Parâmetros e configurações globais
    distance.py      # Funções de distância de Hamming
    logging_utils.py # Logging customizado
    resource_monitor.py # Monitoramento de recursos

src/
    menu.py              # Menus interativos
    runner.py            # Execução e controle de algoritmos
    report_utils.py      # Relatórios e resumos rápidos
    results_formatter.py # Formatação de relatórios detalhados
    console_manager.py   # Saída thread-safe
    batch_executor.py    # Execução em lote

main.py             # Interface principal (ponto de entrada)
results/            # Relatórios gerados
logs/               # Logs de execução
saved_datasets/     # Datasets salvos
batch_configs/      # Configurações de execução em lote
```

## Como Executar

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

- Relatórios detalhados são salvos em `results/` após cada execução.
- Resumos rápidos são exibidos no console.
- Execuções em lote geram relatórios consolidados.

## Suporte e Documentação

- Cada algoritmo possui README próprio explicando heurísticas, funcionamento e parâmetros.
- Consulte `REESTRUTURACAO.md` para detalhes sobre a arquitetura e modularização.
