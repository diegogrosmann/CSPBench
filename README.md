# Closest String Problem – BLF-GA

## Estrutura do Projeto

```
algorithms/
    baseline.py      # Algoritmo de consenso guloso (baseline)
    blf_ga.py        # Algoritmo BLF-GA (Blockwise Learning Fusion + GA)
    csc.py           # Algoritmo Cluster & String Consensus

datasets/
    dataset_file.py      # Leitura de datasets de arquivo
    dataset_entrez.py    # Download de datasets do NCBI
    dataset_synthetic.py # Geração de datasets sintéticos
    dataset_utils.py     # Utilitários para datasets

utils/
    config.py        # Parâmetros e configurações globais

main.py             # Interface principal (ponto de entrada)
results.txt         # Resultados das execuções
```

## Como executar

```bash
python main.py
```

## Descrição
- **algorithms/**: Algoritmos de solução do CSP.
- **datasets/**: Leitura, geração e manipulação de datasets.
- **utils/**: Configurações e utilitários gerais.
- **main.py**: Interface interativa para escolha de dataset e execução dos algoritmos.

## Requisitos
- Python 3.10+
- Biopython
- scikit-learn

Instale dependências com:
```bash
pip install -r requirements.txt
```
