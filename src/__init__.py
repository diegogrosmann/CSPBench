"""
CSPBench: Framework Experimental para Algoritmos do Closest String Problem

CSPBench é um framework robusto e extensível para teste, comparação e análise
de algoritmos do Closest String Problem (CSP). Oferece uma plataforma unificada
para desenvolvimento, execução e avaliação de algoritmos CSP com recursos
avançados de paralelização, monitoramento e relatórios.

CARACTERÍSTICAS PRINCIPAIS:
==========================
- Biblioteca de algoritmos CSP (Baseline, BLF-GA, CSC, DP-CSP, H³-CSP)
- Sistema de execução paralela com scheduler avançado
- Gestão inteligente de datasets (sintéticos, arquivos, download automático)
- Análise e otimização (Optuna, SALib, visualizações)
- Interface curses para monitoramento em tempo real
- Relatórios automatizados e análise estatística
- Sistema extensível para novos algoritmos

COMPONENTES PRINCIPAIS:
======================
- algorithms/: Biblioteca de algoritmos CSP
- src/core/: Núcleo do framework (scheduler, interfaces, relatórios)
- src/datasets/: Gestão de datasets
- src/optimization/: Otimização e análise
- src/ui/: Interfaces de usuário
- src/utils/: Utilitários gerais

EXEMPLO DE USO:
==============
```python
from algorithms.blf_ga import BLFGAAlgorithm
from src.datasets.dataset_synthetic import generate_dataset

# Gerar dataset sintético
strings, _ = generate_dataset(n=50, L=100, alphabet="ACGT")

# Executar algoritmo
algorithm = BLFGAAlgorithm(strings, "ACGT")
center, distance, metadata = algorithm.run()

print(f"Centro encontrado: {center}")
print(f"Distância máxima: {distance}")
```

ARQUITETURA:
===========
CSPBench segue uma arquitetura modular com interfaces bem definidas:
- Algoritmos implementam a interface CSPAlgorithm
- Datasets seguem protocolos padronizados
- Executores gerenciam paralelização e recursos
- Relatórios são gerados automaticamente

EXTENSIBILIDADE:
===============
Novos algoritmos podem ser adicionados facilmente:
1. Implementar a interface CSPAlgorithm
2. Usar o decorador @register_algorithm
3. O sistema detecta automaticamente o novo algoritmo

Para mais informações, consulte a documentação completa.
"""

__version__ = "1.0.0"
__author__ = "Diego Grosmann"
__project__ = "CSPBench"
__description__ = "Framework experimental para algoritmos do Closest String Problem"
