"""
Adaptador para algoritmos legados para a interface IAlgorithm.

Converte algoritmos existentes para a nova interface padronizada,
mantendo compatibilidade com o código existente.
"""

import logging
from typing import Any, Callable, Optional

from ..algorithm import IAlgorithm, Result

logger = logging.getLogger(__name__)


class AlgorithmAdapter:
    """
    Adaptador que converte algoritmos legados para IAlgorithm.

    Permite que algoritmos existentes funcionem com a nova interface
    padronizada sem modificações no código original.
    """

    def __init__(self, algorithm_instance: Any):
        """
        Inicializa o adaptador.

        Args:
            algorithm_instance: Instância do algoritmo legado
        """
        self._algorithm = algorithm_instance
        self._progress_callback: Optional[Callable[[str], None]] = None
        self._warning_callback: Optional[Callable[[str], None]] = None

        # Verificar se o algoritmo tem o atributo is_deterministic
        self.is_deterministic = getattr(algorithm_instance, "is_deterministic", False)

        logger.debug(
            f"AlgorithmAdapter criado para {type(algorithm_instance).__name__}"
        )

    def run(self) -> Result:
        """
        Executa o algoritmo legado e adapta o resultado.

        Returns:
            Result: Resultado padronizado
        """
        try:
            # Configurar callbacks se o algoritmo suporta
            if self._progress_callback and hasattr(
                self._algorithm, "set_progress_callback"
            ):
                self._algorithm.set_progress_callback(self._progress_callback)

            if self._warning_callback and hasattr(
                self._algorithm, "set_warning_callback"
            ):
                self._algorithm.set_warning_callback(self._warning_callback)

            # Executar algoritmo
            result = self._algorithm.run()

            # Adaptar resultado baseado no formato retornado
            if isinstance(result, dict) and "center" in result:
                # Já está no formato correto
                return result
            elif isinstance(result, tuple):
                if len(result) == 3:
                    # Formato: (center, distance, metadata)
                    center, distance, metadata = result
                    return {
                        "center": center,
                        "distance": distance,
                        "metadata": metadata or {},
                    }
                elif len(result) == 2:
                    # Formato legado: (center, distance)
                    center, distance = result
                    metadata = {}

                    # Tentar extrair informações adicionais do algoritmo
                    if hasattr(self._algorithm, "geracao"):
                        metadata["iteracoes"] = self._algorithm.geracao
                    elif hasattr(self._algorithm, "iterations"):
                        metadata["iteracoes"] = self._algorithm.iterations
                    elif hasattr(self._algorithm, "num_iteracoes"):
                        metadata["iteracoes"] = self._algorithm.num_iteracoes

                    return {
                        "center": center,
                        "distance": distance,
                        "metadata": metadata,
                    }
                else:
                    raise ValueError(
                        f"Formato de resultado inválido: {len(result)} elementos"
                    )
            else:
                # Formato desconhecido - tentar extrair informações
                center = str(result) if result is not None else ""
                distance = float("inf")
                metadata = {"erro": "Formato de resultado desconhecido"}

                return {"center": center, "distance": distance, "metadata": metadata}

        except Exception as e:
            logger.error("Erro na execução do algoritmo adaptado: %s", e)
            return {
                "center": "",
                "distance": float("inf"),
                "metadata": {"erro": str(e)},
            }

    def set_progress_callback(self, callback: Callable[[str], None]) -> None:
        """
        Define callback para progresso.

        Args:
            callback: Função de callback para progresso
        """
        self._progress_callback = callback

        # Se o algoritmo suporta, configurar imediatamente
        if hasattr(self._algorithm, "set_progress_callback"):
            self._algorithm.set_progress_callback(callback)

    def set_warning_callback(self, callback: Callable[[str], None]) -> None:
        """
        Define callback para avisos.

        Args:
            callback: Função de callback para avisos
        """
        self._warning_callback = callback

        # Se o algoritmo suporta, configurar imediatamente
        if hasattr(self._algorithm, "set_warning_callback"):
            self._algorithm.set_warning_callback(callback)

    def get_original_algorithm(self) -> Any:
        """
        Retorna a instância original do algoritmo.

        Returns:
            Any: Instância do algoritmo original
        """
        return self._algorithm

    def __getattr__(self, name: str) -> Any:
        """
        Delega atributos não encontrados para o algoritmo original.

        Args:
            name: Nome do atributo

        Returns:
            Any: Valor do atributo do algoritmo original
        """
        return getattr(self._algorithm, name)


def adapt_algorithm(algorithm_instance: Any) -> IAlgorithm:
    """
    Adapta um algoritmo legado para a interface IAlgorithm.

    Args:
        algorithm_instance: Instância do algoritmo legado

    Returns:
        IAlgorithm: Algoritmo adaptado
    """
    # Verificar se já implementa a interface
    if (
        hasattr(algorithm_instance, "run")
        and hasattr(algorithm_instance, "set_progress_callback")
        and hasattr(algorithm_instance, "set_warning_callback")
    ):

        # Verificar se run() retorna Result
        try:
            # Não executar, apenas verificar assinatura
            import inspect

            sig = inspect.signature(algorithm_instance.run)
            # Se chegou até aqui, provavelmente já está compatível
            return algorithm_instance
        except:
            pass

    # Criar adaptador
    return AlgorithmAdapter(algorithm_instance)


def is_algorithm_adapted(algorithm_instance: Any) -> bool:
    """
    Verifica se um algoritmo já foi adaptado.

    Args:
        algorithm_instance: Instância do algoritmo

    Returns:
        bool: True se é um algoritmo adaptado
    """
    return isinstance(algorithm_instance, AlgorithmAdapter)
