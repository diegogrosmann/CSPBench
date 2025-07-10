"""
Pacote de algoritmos CSP.

Este módulo inicializa o pacote algorithms e implementa o sistema de
auto-descoberta de algoritmos. Todos os algoritmos implementados são
automaticamente registrados no registry global.

ALGORITMOS IMPLEMENTADOS:

1. **BLF-GA (Blockwise Learning Fusion + Genetic Algorithm)**:
   - Metaheurística híbrida avançada
   - Combina aprendizado por blocos com evolução genética
   - Mecanismos adaptativos para diversidade e convergência
   - Paralelização eficiente e configurabilidade avançada

2. **Baseline (Algoritmos Simples)**:
   - Greedy Consensus: Consenso por votação majoritária
   - Max Distance: Busca exaustiva limitada
   - Baselines para comparação de performance

3. **CSC (Consensus String Clustering)**:
   - Estratégia de divisão e conquista
   - Clustering hierárquico de strings
   - Otimização local com recombinação global

4. **DP-CSP (Dynamic Programming CSP)**:
   - Programação dinâmica com poda
   - Garantia de otimalidade dentro de limites de recurso
   - Modelagem de estados eficiente

5. **H³-CSP (Hybrid Hierarchical Hamming Search)**:
   - Abordagem hierárquica em três camadas
   - Seleção adaptativa de técnicas por bloco
   - Balanceamento automático entre qualidade e eficiência

AUTO-DESCOBERTA:
O sistema automaticamente descobre e registra todos os algoritmos
implementados nos subpacotes, permitindo uso dinâmico através do
registry global.

EXEMPLO DE USO:
```python
from algorithms import global_registry

# Listar algoritmos disponíveis
print("Algoritmos disponíveis:")
for name, cls in global_registry.items():
    print(f"  {name}: {cls.__doc__.split('.')[0]}")

# Instanciar algoritmo dinamicamente
algorithm_class = global_registry["BLF-GA"]
algorithm = algorithm_class(strings, alphabet, **params)
solution, distance, metadata = algorithm.run()
```

Exports:
    global_registry: Dicionário com todos os algoritmos registrados
    register_algorithm: Decorador para registro de novos algoritmos
"""

# algorithms package

import importlib
import pkgutil

from .base import global_registry, register_algorithm

# Auto-import subpackages to register algorithms
for _, name, ispkg in pkgutil.iter_modules(__path__):
    if ispkg:
        importlib.import_module(f"{__name__}.{name}")
