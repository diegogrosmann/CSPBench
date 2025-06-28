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
def meu_algoritmo_funcao(strings, alphabet, progress_callback=None, **params):
    """
    Exemplo de função de implementação de algoritmo customizado.

    Args:
        strings (list[str]): Lista de strings de entrada.
        alphabet (str): Alfabeto utilizado.
        progress_callback (callable, opcional): Callback de progresso.
        **params: Parâmetros customizados.

    Returns:
        str: String centro encontrada.
    """
    if progress_callback:
        progress_callback("Iniciando processamento...")
    # Sua lógica aqui
    center = "ACGT"  # exemplo
    if progress_callback:
        progress_callback("Finalizando...")
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
        self.progress_callback = None

    def set_progress_callback(self, callback):
        """Define callback para reportar progresso"""
        self.progress_callback = callback

    def run(self) -> tuple[str, int, dict]:
        """Execute o algoritmo e retorna (center, distance, info)"""
        try:
            center = meu_algoritmo_funcao(
                self.strings, 
                self.alphabet, 
                progress_callback=self.progress_callback,
                **self.params
            )
            dist = max_hamming(center, self.strings)
            
            # Informações adicionais para auditoria
            info = {
                'iteracoes': 1,
                'melhor_string': center
            }
            
            return center, dist, info
            
        except Exception as e:
            # Retorna erro estruturado
            return None, float('inf'), {'erro': str(e)}
```

### 5. Expor no `__init__.py`
```python
from .algorithm import MeuAlgoritmo
```

## Dicas
- Use o decorador `@register_algorithm` para registro automático.
- O algoritmo aparecerá no menu sem necessidade de alterar o main.
- Veja exemplos em `algorithms/baseline/`, `algorithms/blf_ga/`, etc.
