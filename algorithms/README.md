# Como Adicionar Novos Algoritmos

A aplicação possui uma interface padronizada que permite adicionar novos algoritmos sem modificar o código principal (`main.py`).

## Estrutura de um Algoritmo

Cada algoritmo deve estar em sua própria pasta dentro de `algorithms/` com a seguinte estrutura:

```
algorithms/
├── meu_algoritmo/
│   ├── __init__.py          # Expõe o algoritmo
│   ├── algorithm.py         # Wrapper que implementa a interface Algorithm
│   ├── config.py            # Configurações específicas do algoritmo
│   └── implementation.py    # Implementação real do algoritmo
```

## Passo a Passo

### 1. Criar a pasta do algoritmo
```bash
mkdir algorithms/meu_algoritmo
```

### 2. Criar o arquivo de configuração (`config.py`)
```python
# Configurações do Meu Algoritmo
MEU_ALGORITMO_DEFAULTS = {
    'param1': 'valor_padrao',
    'param2': 42,
    'max_time': 300.0,
}
```

### 3. Criar a implementação (`implementation.py`)
```python
# Sua implementação real do algoritmo
def meu_algoritmo_funcao(strings, alphabet, **params):
    # Sua lógica aqui
    center = "ACGT"  # exemplo
    return center
```

### 4. Criar o wrapper (`algorithm.py`)
```python
from algorithms.base import Algorithm, register_algorithm
from .config import MEU_ALGORITMO_DEFAULTS
from .implementation import meu_algoritmo_funcao
from utils.distance import max_hamming

@register_algorithm
class MeuAlgoritmo(Algorithm):
    """Descrição do meu algoritmo"""
    name = "Meu Algoritmo"
    default_params = MEU_ALGORITMO_DEFAULTS

    def __init__(self, strings: list[str], alphabet: str, **params):
        self.strings = strings
        self.alphabet = alphabet
        self.params = {**self.default_params, **params}

    def run(self) -> tuple[str, int]:
        """Execute o algoritmo"""
        center = meu_algoritmo_funcao(
            self.strings, 
            self.alphabet, 
            **self.params
        )
        dist = max_hamming(center, self.strings)
        return center, dist
```

### 5. Criar o `__init__.py`
```python
from .algorithm import MeuAlgoritmo
```

## Importante

- **O decorator `@register_algorithm`** é obrigatório para que o algoritmo apareça no menu
- **O nome do algoritmo** será exibido no menu de seleção
- **A interface `run()`** deve retornar uma tupla `(center_string, max_distance)`
- **Configurações específicas** devem ficar no `config.py` do algoritmo, não no `utils/config.py`

## Exemplo Completo

Veja os algoritmos existentes em:
- `algorithms/baseline/` - Exemplo mais simples
- `algorithms/blf_ga/` - Exemplo mais complexo
- `algorithms/csc/` - Exemplo com lógica de fallback

Após criar seu algoritmo seguindo essa estrutura, ele aparecerá automaticamente no menu principal sem necessidade de modificar o `main.py`.
