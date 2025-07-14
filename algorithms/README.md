# Guia para Adição de Novos Algoritmos ao CSPBench

O framework CSPBench possui uma interface padronizada e documentação detalhada (Google style docstrings) em todos os arquivos, facilitando a integração de novos algoritmos sem modificar o código principal (`main.py`). O registro é automático via decorador.

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

## Passo a Passo para Integração

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
        """
        Executa o algoritmo MeuAlgoritmo.
        Args:
            strings (list[str]): Lista de strings de entrada.
            alphabet (str): Alfabeto utilizado.
            progress_callback (callable, opcional): Callback de progresso.
            **params: Parâmetros do algoritmo.
        Returns:
            str: String central encontrada.
        """
        if progress_callback:
            progress_callback("Iniciando processamento...")
        center = "ACGT"  # Exemplo
        if progress_callback:
            progress_callback("Finalizando...")
        return center
    ```

4. **Criar o wrapper (`algorithm.py`)**
    ```python
    from cspbench.domain.algorithms import CSPAlgorithm, register_algorithm
    from .config import MEU_ALGORITMO_DEFAULTS
    from .implementation import meu_algoritmo_funcao
    from utils.distance import max_hamming

    @register_algorithm
    class MeuAlgoritmo(CSPAlgorithm):
        """
        Wrapper para integração do MeuAlgoritmo ao framework CSPBench.
        """
        name = "Meu Algoritmo"
        default_params = MEU_ALGORITMO_DEFAULTS

        def __init__(self, strings, alphabet, **params):
            self.strings = strings
            self.alphabet = alphabet
            self.params = {**self.default_params, **params}

        def run(self):
            center = meu_algoritmo_funcao(self.strings, self.alphabet, **self.params)
            dist = max_hamming(center, self.strings)
            return center, dist
    ```

5. **Expor o algoritmo no `__init__.py`**
    ```python
    from .algorithm import MeuAlgoritmo
    ```

6. **(Opcional) Adicionar README.md explicando a heurística, parâmetros e uso.**

## Observações

- Todos os algoritmos devem implementar a interface `Algorithm` e usar o decorador `@register_algorithm`.
- O algoritmo aparecerá automaticamente no menu do sistema.
- Consulte os READMEs dos algoritmos existentes para exemplos detalhados e padrões de documentação.
- Utilize docstrings no estilo Google para facilitar a geração de documentação automática.
