"""
README do projeto Closest String Problem (CSP)

Este projeto implementa uma arquitetura modular para experimentação com algoritmos de solução do Closest String Problem (CSP), incluindo baseline, heurísticas e métodos exatos.

## Estrutura do Projeto

```
algorithms/
    baseline/         # Algoritmo de consenso guloso (baseline)
    blf_ga/           # Algoritmo BLF-GA (Blockwise Learning Fusion + GA)
    csc/              # Algoritmo Cluster & String Consensus
    dp_csp/           # Algoritmo exato por programação dinâmica
    h3_csp/           # Algoritmo híbrido hierárquico

datasets/
    dataset_file.py      # Leitura de datasets de arquivo
    dataset_entrez.py    # Download de datasets do NCBI
    dataset_synthetic.py # Geração de datasets sintéticos
    dataset_utils.py     # Utilitários para datasets

utils/
    config.py        # Parâmetros e configurações globais
    distance.py      # Funções de distância
    logging_utils.py # Logging customizado

src/
    menu.py              # Menus interativos
    runner.py            # Execução e controle de algoritmos
    report_utils.py      # Relatórios e resumos
    results_formatter.py # Formatação de resultados
    console_manager.py   # Saída thread-safe

main.py             # Interface principal (ponto de entrada)
results/            # Resultados das execuções
logs/               # Logs de execução
saved_datasets/     # Datasets salvos
```

## Como executar

```bash
python main.py
```

## Descrição
- **algorithms/**: Algoritmos de solução do CSP, cada um modularizado em sua pasta.
- **datasets/**: Leitura, geração e manipulação de datasets.
- **utils/**: Configurações e utilitários gerais.
- **src/**: Componentes auxiliares, menus, runner, relatórios e formatação.
- **main.py**: Interface interativa para escolha de dataset e execução dos algoritmos.

## Requisitos
- Python 3.10+
- Biopython
- scikit-learn
- tabulate

Instale dependências com:
```bash
pip install -r requirements.txt
```

## Adicionando Novos Algoritmos

Basta criar uma nova pasta em `algorithms/` seguindo o padrão:

```
algorithms/
  meu_algoritmo/
    __init__.py
    algorithm.py
    config.py
    implementation.py
```

Implemente a interface `Algorithm` e use o decorador `@register_algorithm` para registro automático. O algoritmo aparecerá no menu sem necessidade de alterar o main.

Veja `algorithms/README.md` para exemplos detalhados.
```
