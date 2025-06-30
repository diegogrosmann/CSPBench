# Como Adicionar Novos Algoritmos ao CSP

A plataforma CSP possui uma interface padronizada que permite adicionar novos algoritmos sem modificar o código principal (`main.py`). O registro é automático via decorador.

## Estrutura Recomendada

Cada algoritmo deve estar em sua própria pasta dentro de `algorithms/`:

```
algorithms/
├── meu_algoritmo/
│   ├── __init__.py          # Expõe o algoritmo
│   ├── algorithm.py         # Wrapper que implementa a interface Algorithm
│   ├── config.py            # Configurações específicas do algoritmo
│   └── implementation.py    # Implementação real do algoritmo
```

## Passo a Passo

1. **Criar a pasta do algoritmo**
    ```bash
    mkdir algorithms/meu_algoritmo
    ```

2. **Criar o arquivo de configuração (`config.py`)**
    ```python
    # Exemplo
    MEU_ALGORITMO_DEFAULTS = {
        'param1': 'valor_padrao',
        'param2': 42,
        'max_time': 300.0,
    }
    ```

3. **Implementar a lógica principal (`implementation.py`)**
    ```python
    def meu_algoritmo_funcao(strings, alphabet, progress_callback=None, **params):
        # Sua lógica aqui
        if progress_callback:
            progress_callback("Iniciando processamento...")
        center = "ACGT"  # Exemplo
        if progress_callback:
            progress_callback("Finalizando...")
        return center
    ```

4. **Criar o wrapper (`algorithm.py`)**
    ```python
    from algorithms.base import Algorithm, register_algorithm
    from .config import MEU_ALGORITMO_DEFAULTS
    from .implementation import meu_algoritmo_funcao
    from utils.distance import max_hamming

    @register_algorithm
    class MeuAlgoritmo(Algorithm):
        name = "Meu Algoritmo"
        default_params = MEU_ALGORITMO_DEFAULTS

        def __init__(self, strings, alphabet, **params):
            self.strings = strings
            self.alphabet = alphabet
            self.params = {**self.default_params, **params}
            self.progress_callback = None

        def set_progress_callback(self, callback):
            self.progress_callback = callback

        def run(self):
            try:
                center = meu_algoritmo_funcao(
                    self.strings, self.alphabet,
                    progress_callback=self.progress_callback,
                    **self.params
                )
                dist = max_hamming(center, self.strings)
                info = {'iteracoes': 1, 'melhor_string': center}
                return center, dist, info
            except Exception as e:
                return None, float('inf'), {'erro': str(e)}
    ```

5. **Expor no `__init__.py`**
    ```python
    from .algorithm import MeuAlgoritmo
    ```

## Dicas

- Use o decorador `@register_algorithm` para registro automático.
- O algoritmo aparecerá no menu sem necessidade de alterar o main.
- Consulte exemplos em `algorithms/baseline/`, `algorithms/blf_ga/`, etc.

## Interface Esperada

- O método `run()` deve retornar `(center, distancia, info)`:
    - `center`: string encontrada
    - `distancia`: distância máxima
    - `info`: dicionário com detalhes (iteracoes, melhor_string, erro, etc.)

- Para algoritmos não determinísticos, utilize o parâmetro de execuções múltiplas.

## Documentação

Cada algoritmo deve possuir um `README.md` explicando:
- Estratégia e heurísticas principais
- Parâmetros configuráveis
- Fluxo de execução
- Limitações e cenários ideais

Veja os READMEs dos algoritmos existentes para inspiração.
