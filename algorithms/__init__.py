"""
Pacote de algoritmos CSP.

Este módulo inicializa o pacote algorithms e implementa o sistema de
auto-descoberta de algoritmos. Todos os algoritmos implementados são
automaticamente registrados no registry global.

ALGORITMOS IMPLEMENTADOS:

1. **Baseline (Algoritmos Simples)**:
   - Greedy Consensus: Consenso por votação majoritária
   - Baselines para comparação de performance

2. **BLF-GA (Blockwise Learning Fusion + Genetic Algorithm)**:
   - Metaheurística híbrida avançada
   - Combina aprendizado por blocos com evolução genética
   - Mecanismos adaptativos para diversidade e convergência
   - Paralelização eficiente e configurabilidade avançada

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
from cspbench.domain.algorithms import global_registry

# Listar algoritmos disponíveis
print("Algoritmos disponíveis:")
for name, cls in global_registry.items():
    print(f"  {name}: {cls.__doc__.split('.')[0]}")

# Usar um algoritmo específico
algorithm_class = global_registry["Baseline"]
algorithm = algorithm_class(strings=["ATCG", "ATCC"], alphabet="ATCG")
result_string, max_distance, metadata = algorithm.run()
```

ESTRUTURA:
Cada algoritmo deve estar em seu próprio subpacote com:
- __init__.py: Exposição da classe principal
- algorithm.py: Wrapper com decorador @register_algorithm
- implementation.py: Lógica específica do algoritmo
- config.py: Configurações e parâmetros padrão
- README.md: Documentação detalhada
"""

import importlib
import pkgutil
from pathlib import Path

# Importa o registry do domínio
from src.domain.algorithms import global_registry, register_algorithm


# Auto-descoberta e importação de algoritmos
def _discover_algorithms():
    """Descobre e importa automaticamente todos os algoritmos."""
    algorithms_path = Path(__file__).parent

    for importer, modname, ispkg in pkgutil.iter_modules([str(algorithms_path)]):
        if ispkg and not modname.startswith("_"):
            try:
                # Importa o subpacote para ativar o registro automático
                importlib.import_module(f"algorithms.{modname}")
            except ImportError as e:
                # Algoritmo pode ter dependências opcionais
                print(f"Aviso: Algoritmo '{modname}' não pôde ser carregado: {e}")


# Executa auto-descoberta na importação
_discover_algorithms()

__all__ = ["global_registry", "register_algorithm"]
