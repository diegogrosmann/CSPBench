"""
Base e registry global para algoritmos CSP.

Este módulo define a infraestrutura fundamental para todos os algoritmos CSP
do projeto, estabelecendo interfaces padronizadas, sistema de registro e
funcionalidades compartilhadas.

ARQUITETURA DO SISTEMA:

┌─────────────────────────────────────────────────────────────────────────────────┐
│                        INFRAESTRUTURA DE ALGORITMOS CSP                         │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 1. REGISTRY GLOBAL                                                             │
│   ├── Registro automático de algoritmos via decorador                          │
│   ├── Descoberta dinâmica de algoritmos disponíveis                            │
│   ├── Lookup por nome para instanciação                                        │
│   └── Validação de conformidade com interfaces                                 │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 2. INTERFACE PADRÃO (CSPAlgorithm)                                            │
│   ├── Inicialização consistente (strings, alfabeto, parâmetros)               │
│   ├── Execução padronizada com retorno estruturado                             │
│   ├── Sistema de callbacks para progresso e warnings                           │
│   ├── Suporte a timeout e paralelismo interno                                  │
│   ├── Metadados estruturados para análise                                      │
│   └── Configuração dinâmica de parâmetros                                      │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 3. COMPATIBILIDADE LEGACY                                                      │
│   ├── Interface Algorithm para código legado                                   │
│   ├── Migração gradual para nova interface                                     │
│   └── Preservação de funcionalidades existentes                                │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 4. FUNCIONALIDADES AVANÇADAS                                                   │
│   ├── Callbacks de progresso para UIs interativas                              │
│   ├── Sistema de warnings para alertas não-críticos                            │
│   ├── Metadados ricos para análise de performance                              │
│   ├── Configuração dinâmica de parâmetros                                      │
│   └── Compatibilidade com ProcessPoolExecutor                                  │
└─────────────────────────────────────────────────────────────────────────────────┘

VANTAGENS DA ARQUITETURA:

• **EXTENSIBILIDADE**: Facilita adição de novos algoritmos
• **CONSISTÊNCIA**: Interface padronizada para todos os algoritmos
• **MONITORAMENTO**: Sistema de callbacks para acompanhar execução
• **FLEXIBILIDADE**: Configuração dinâmica de parâmetros
• **INTEROPERABILIDADE**: Compatibilidade com ferramentas de paralelismo
• **MANUTENIBILIDADE**: Código legado suportado durante migração

DESIGN PATTERNS IMPLEMENTADOS:

• **Registry Pattern**: Para descoberta e instanciação de algoritmos
• **Template Method**: Interface abstrata com comportamento compartilhado
• **Observer Pattern**: Sistema de callbacks para notificações
• **Strategy Pattern**: Algoritmos intercambiáveis com interface comum

EXEMPLO DE USO:

```python
# Registro automático de algoritmo
@register_algorithm
class MeuAlgoritmo(CSPAlgorithm):
    name = "meu_algoritmo"
    default_params = {"param1": 10, "param2": 0.5}
    is_deterministic = True

    def __init__(self, strings, alphabet, **params):
        super().__init__(strings, alphabet, **params)
        # Inicialização específica

    def run(self):
        self._report_progress("Iniciando algoritmo...")
        # Implementação do algoritmo
        return solution, distance, metadata

# Descoberta e uso
algorithm_class = global_registry["meu_algoritmo"]
algorithm = algorithm_class(strings, alphabet, param1=20)
algorithm.set_progress_callback(lambda msg: print(f"Progresso: {msg}"))
solution, distance, metadata = algorithm.run()
```

Classes:
    CSPAlgorithm: Interface abstrata moderna para algoritmos CSP.
    Algorithm: Interface abstrata legacy (mantida para compatibilidade).

Funções:
    register_algorithm(cls): Decorador para registrar algoritmos no registry global.

Atributos:
    global_registry (dict): Registry global de algoritmos disponíveis.

Author: Framework de algoritmos CSP
Version: Interface moderna com compatibilidade legacy
"""

from abc import ABC, abstractmethod
from collections.abc import Callable
from typing import Any

# Registry for all algorithms
global_registry: dict[str, type] = {}


def register_algorithm(cls: type) -> type:
    """
    Decorador para registrar uma classe de algoritmo no registry global.

    Esta função implementa o padrão Registry, permitindo descoberta automática
    de algoritmos disponíveis no sistema. Quando uma classe é decorada com
    @register_algorithm, ela é automaticamente adicionada ao registry global
    e pode ser instanciada dinamicamente por nome.

    FUNCIONAMENTO:
    1. Extrai o nome do algoritmo da classe (atributo 'name' ou nome da classe)
    2. Registra a classe no dicionário global_registry
    3. Retorna a classe original (comportamento de decorador)

    VANTAGENS:
    - Descoberta automática de algoritmos
    - Instanciação dinâmica por nome
    - Desacoplamento entre uso e implementação
    - Facilita extensibilidade do sistema

    EXEMPLO DE USO:
    ```python
    @register_algorithm
    class MeuAlgoritmo(CSPAlgorithm):
        name = "meu_algoritmo"
        # ...implementação...

    # Uso posterior
    AlgClass = global_registry["meu_algoritmo"]
    algoritmo = AlgClass(strings, alphabet, **params)
    ```

    Args:
        cls (type): Classe do algoritmo a ser registrada.
                   Deve implementar a interface CSPAlgorithm.

    Returns:
        type: A própria classe, permitindo uso como decorador.

    Note:
        Se a classe não possui atributo 'name', usa cls.__name__ como nome.
        Classes com mesmo nome sobrescrevem registros anteriores.
    """
    # Extrai nome da classe: prioriza atributo 'name', fallback para __name__
    algorithm_name = getattr(cls, "name", cls.__name__)

    # Registra classe no registry global
    global_registry[algorithm_name] = cls

    # Retorna classe original para uso como decorador
    return cls


class CSPAlgorithm(ABC):
    """
    Classe base abstrata moderna para todos os algoritmos CSP.

    Esta classe define a interface padrão que todos os algoritmos CSP devem
    implementar, fornecendo funcionalidades comuns e garantindo consistência
    na integração com o sistema.

    FUNCIONALIDADES PRINCIPAIS:

    1. **INICIALIZAÇÃO PADRONIZADA**:
       - Recebe strings, alfabeto e parâmetros de forma consistente
       - Merge de parâmetros padrão com parâmetros específicos
       - Configuração de callbacks de progresso e warnings

    2. **SISTEMA DE CALLBACKS**:
       - Callback de progresso para UIs interativas
       - Callback de warnings para alertas não-críticos
       - Métodos auxiliares para reportar eventos

    3. **CONFIGURAÇÃO DINÂMICA**:
       - Atualização de parâmetros durante execução
       - Propagação automática para instâncias internas
       - Flexibilidade para ajuste fino

    4. **METADADOS ESTRUTURADOS**:
       - Informações sobre o algoritmo e sua configuração
       - Dados de entrada (tamanho, comprimento, alfabeto)
       - Características do algoritmo (determinismo, paralelismo)

    5. **COMPATIBILIDADE COM PARALELISMO**:
       - Suporte a ProcessPoolExecutor
       - Indicação de capacidade de paralelismo interno
       - Interface thread-safe para callbacks

    ATRIBUTOS OBRIGATÓRIOS:
    - name (str): Nome único do algoritmo
    - default_params (dict): Parâmetros padrão
    - is_deterministic (bool): Se produz resultados reproduzíveis
    - supports_internal_parallel (bool): Se suporta paralelismo interno

    MÉTODOS OBRIGATÓRIOS:
    - __init__: Inicialização com strings, alfabeto e parâmetros
    - run: Execução principal retornando (solução, distância, metadata)

    EXEMPLO DE IMPLEMENTAÇÃO:
    ```python
    @register_algorithm
    class MeuAlgoritmo(CSPAlgorithm):
        name = "meu_algoritmo"
        default_params = {"iterations": 100, "threshold": 0.1}
        is_deterministic = True
        supports_internal_parallel = True

        def __init__(self, strings, alphabet, **params):
            super().__init__(strings, alphabet, **params)
            self.iterations = self.params["iterations"]
            self.threshold = self.params["threshold"]

        def run(self):
            self._report_progress("Iniciando...")
            # Implementação do algoritmo
            solution = "ACGT"
            distance = 2
            metadata = {"iterations_used": 50}
            return solution, distance, metadata
    ```

    Atributos obrigatórios:
        name (str): Nome único do algoritmo para registro.
        default_params (dict): Parâmetros padrão do algoritmo.
        is_deterministic (bool): Se é determinístico (padrão: False).
        supports_internal_parallel (bool): Se suporta paralelismo interno (padrão: False).

    Atributos de instância:
        strings (list[str]): Lista de strings do dataset.
        alphabet (str): Alfabeto utilizado.
        params (dict): Parâmetros atuais do algoritmo.
        progress_callback: Callback para reportar progresso.
        warning_callback: Callback para reportar warnings.
    """

    # Atributos de classe obrigatórios (devem ser definidos pelas subclasses)
    name: str
    default_params: dict
    is_deterministic: bool = False
    supports_internal_parallel: bool = False

    @abstractmethod
    def __init__(self, strings: list[str], alphabet: str, **params):
        """
        Inicializa o algoritmo com as strings e alfabeto.

        Configura o estado inicial do algoritmo, incluindo merge de parâmetros
        padrão com parâmetros específicos e inicialização de callbacks.

        PROCESSO DE INICIALIZAÇÃO:
        1. Armazena strings e alfabeto de entrada
        2. Merge parâmetros padrão com parâmetros fornecidos
        3. Inicializa callbacks como None
        4. Subclasses podem fazer inicialização adicional

        Args:
            strings (list[str]): Lista de strings do dataset.
                                Todas devem ter mesmo comprimento.
            alphabet (str): Alfabeto utilizado (ex: "ACGT" para DNA).
            **params: Parâmetros específicos do algoritmo.
                     Sobrescrevem valores em default_params.

        Example:
            >>> alg = MeuAlgoritmo(["ACGT", "AGCT"], "ACGT", iterations=50)
            >>> alg.params["iterations"]
            50
        """
        self.strings = strings
        self.alphabet = alphabet

        # Merge parâmetros padrão com parâmetros específicos
        self.params = {**self.default_params, **params}

        # Inicializa callbacks como None
        self.progress_callback: Callable[[str], None] | None = None
        self.warning_callback: Callable[[str], None] | None = None

    def set_progress_callback(self, callback: Callable[[str], None]) -> None:
        """
        Define callback para relatar progresso do algoritmo.

        O callback será chamado periodicamente durante a execução do algoritmo
        com mensagens descritivas sobre o progresso atual. Permite implementar
        barras de progresso, logs detalhados ou atualizações de UI.

        Args:
            callback (Callable[[str], None]): Função que recebe mensagens
                                             de progresso como string.

        Example:
            >>> def progress_handler(msg):
            ...     print(f"[{datetime.now()}] {msg}")
            >>> algorithm.set_progress_callback(progress_handler)
        """
        self.progress_callback = callback

    def set_warning_callback(self, callback: Callable[[str], None]) -> None:
        """
        Define callback para relatar warnings do algoritmo.

        O callback será chamado quando o algoritmo encontrar situações
        não-críticas que merecem atenção, como timeouts, convergência
        prematura ou parâmetros subótimos.

        Args:
            callback (Callable[[str], None]): Função que recebe mensagens
                                             de warning como string.

        Example:
            >>> def warning_handler(msg):
            ...     logging.warning(f"Algorithm warning: {msg}")
            >>> algorithm.set_warning_callback(warning_handler)
        """
        self.warning_callback = callback

    def _report_progress(self, message: str) -> None:
        """
        Relata progresso se callback estiver definido.

        Método auxiliar para reportar progresso de forma segura.
        Não faz nada se callback não estiver configurado.

        Args:
            message (str): Mensagem de progresso a ser reportada.
        """
        if self.progress_callback:
            self.progress_callback(message)

    def _report_warning(self, message: str) -> None:
        """
        Relata warning se callback estiver definido.

        Método auxiliar para reportar warnings de forma segura.
        Não faz nada se callback não estiver configurado.

        Args:
            message (str): Mensagem de warning a ser reportada.
        """
        if self.warning_callback:
            self.warning_callback(message)

    @abstractmethod
    def run(self) -> tuple[str, int, dict[str, Any]]:
        """
        Executa o algoritmo e retorna resultado estruturado.

        Método principal que implementa a lógica do algoritmo. Deve ser
        implementado por todas as subclasses.

        CONTRATO DE RETORNO:
        - Primeiro elemento: String central encontrada
        - Segundo elemento: Distância máxima (valor objetivo)
        - Terceiro elemento: Metadados da execução

        METADADOS SUGERIDOS:
        - iterations: Número de iterações executadas
        - time: Tempo de execução em segundos
        - converged: Se convergiu ou foi interrompido
        - best_fitness_history: Histórico de melhores fitness

        Returns:
            tuple[str, int, dict[str, Any]]: Tupla contendo:
                - str: String central encontrada
                - int: Distância máxima (valor objetivo)
                - dict: Metadados da execução

        Example:
            >>> solution, distance, metadata = algorithm.run()
            >>> print(f"Solução: {solution}")
            >>> print(f"Distância: {distance}")
            >>> print(f"Iterações: {metadata['iterations']}")
        """

    def set_params(self, **params) -> None:
        """
        Define novos parâmetros para o algoritmo.

        Permite configuração dinâmica de parâmetros após inicialização.
        Automaticamente tenta propagar mudanças para instâncias internas
        conhecidas dos algoritmos.

        PROPAGAÇÃO AUTOMÁTICA:
        Tenta atualizar instâncias internas com nomes conhecidos:
        - blf_ga_instance: Para algoritmos BLF-GA
        - csc_instance: Para algoritmos CSC
        - h3_instance: Para algoritmos H³-CSP
        - dp_instance: Para algoritmos DP-CSP

        Args:
            **params: Parâmetros a serem atualizados.

        Example:
            >>> algorithm.set_params(iterations=200, threshold=0.05)
            >>> algorithm.params["iterations"]
            200
        """
        # Atualiza parâmetros principais
        self.params.update(params)

        # Tenta atualizar instâncias internas se existirem
        internal_instances = [
            "blf_ga_instance",
            "csc_instance",
            "h3_instance",
            "dp_instance",
        ]

        for instance_name in internal_instances:
            if hasattr(self, instance_name):
                instance = getattr(self, instance_name)
                # Atualiza apenas atributos que existem na instância
                for key, value in params.items():
                    if hasattr(instance, key):
                        setattr(instance, key, value)

    def get_metadata(self) -> dict[str, Any]:
        """
        Retorna metadados do algoritmo.

        Fornece informações estruturadas sobre o algoritmo e sua configuração
        atual. Útil para análise, debugging e relatórios.

        METADADOS INCLUÍDOS:
        - name: Nome do algoritmo
        - params: Parâmetros atuais
        - is_deterministic: Se é determinístico
        - supports_internal_parallel: Se suporta paralelismo interno
        - input_size: Número de strings de entrada
        - string_length: Comprimento das strings
        - alphabet_size: Tamanho do alfabeto

        Returns:
            dict[str, Any]: Dicionário com informações sobre o algoritmo.

        Example:
            >>> metadata = algorithm.get_metadata()
            >>> print(f"Algoritmo: {metadata['name']}")
            >>> print(f"Entrada: {metadata['input_size']} strings")
            >>> print(f"Determinístico: {metadata['is_deterministic']}")
        """
        return {
            "name": self.name,
            "params": self.params.copy(),  # Cópia para evitar modificações
            "is_deterministic": self.is_deterministic,
            "supports_internal_parallel": self.supports_internal_parallel,
            "input_size": len(self.strings),
            "string_length": len(self.strings[0]) if self.strings else 0,
            "alphabet_size": len(self.alphabet),
        }
