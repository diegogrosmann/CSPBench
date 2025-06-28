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
    # Use progress_callback para reportar progresso
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

### 5. Criar o `__init__.py`
```python
from .algorithm import MeuAlgoritmo
```

## Interface Obrigatória

### Método `run()`
Deve retornar uma tupla com 3 elementos:
- **center** (str): String centro encontrada (ou None em caso de erro)
- **distance** (int): Distância máxima (ou float('inf') em caso de erro)  
- **info** (dict): Informações adicionais incluindo:
  - `iteracoes`: Número de iterações realizadas
  - `melhor_string`: String resultado (mesmo que center)
  - `erro`: Mensagem de erro (se aplicável)

### Callback de Progresso
- **Método `set_progress_callback(callback)`**: Define função para reportar progresso
- **Uso**: Chame `callback(mensagem)` para enviar updates de progresso
- **Importante**: NUNCA use `print()` diretamente nos algoritmos

### Tratamento de Erros
- **Capture todas as exceções** no método `run()`
- **Retorne erros estruturados** via campo `erro` no dict info
- **Não deixe exceções escaparem** do algoritmo

## Importante

- **O decorator `@register_algorithm`** é obrigatório para que o algoritmo apareça no menu
- **O nome do algoritmo** será exibido no menu de seleção
- **Callback de progresso** deve ser usado em vez de prints diretos
- **Tratamento de erro** deve ser feito internamente e retornado estruturado
- **Configurações específicas** devem ficar no `config.py` do algoritmo, não no `utils/config.py`

## Exemplo Completo

Veja os algoritmos existentes em:
- `algorithms/baseline/` - Exemplo mais simples
- `algorithms/blf_ga/` - Exemplo mais complexo com progresso
- `algorithms/csc/` - Exemplo com lógica de fallback

Após criar seu algoritmo seguindo essa estrutura, ele aparecerá automaticamente no menu principal sem necessidade de modificar o `main.py`.
