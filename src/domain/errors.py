"""
Exceções do Domínio CSPBench

Define exceções customizadas para erros de domínio e aplicação.
"""


class CSPBenchError(Exception):
    """Exceção base para todas as exceções do CSPBench."""

    pass


class DomainError(CSPBenchError):
    """Exceção base para erros de domínio."""

    pass


class ApplicationError(CSPBenchError):
    """Exceção base para erros de aplicação."""

    pass


# Exceções de Dataset
class DatasetError(DomainError):
    """Erro relacionado a datasets."""

    pass


class DatasetNotFoundError(DatasetError):
    """Dataset não encontrado."""

    pass


class DatasetValidationError(DatasetError):
    """Erro de validação de dataset."""

    pass


class DatasetEmptyError(DatasetError):
    """Dataset vazio."""

    pass


# Exceções de Algoritmo
class AlgorithmError(DomainError):
    """Erro relacionado a algoritmos."""

    pass


class AlgorithmNotFoundError(AlgorithmError):
    """Algoritmo não encontrado."""

    pass


class AlgorithmExecutionError(AlgorithmError):
    """Erro na execução de algoritmo."""

    pass


class AlgorithmRegistrationError(AlgorithmError):
    """Erro no registro de algoritmo."""

    pass


class AlgorithmTimeoutError(AlgorithmError):
    """Timeout na execução de algoritmo."""

    pass


class AlgorithmParameterError(AlgorithmError):
    """Erro nos parâmetros do algoritmo."""

    pass


# Exceções de Configuração
class ConfigurationError(ApplicationError):
    """Erro de configuração."""

    pass


class BatchConfigurationError(ConfigurationError):
    """Erro na configuração de batch."""

    pass


class OptimizationConfigurationError(ConfigurationError):
    """Erro na configuração de otimização."""

    pass


class SensitivityConfigurationError(ConfigurationError):
    """Erro na configuração de análise de sensibilidade."""

    pass


# Exceções de Execução
class ExecutionError(ApplicationError):
    """Erro de execução."""

    pass


class BatchExecutionError(ExecutionError):
    """Erro na execução de batch."""

    pass


class OptimizationExecutionError(ExecutionError):
    """Erro na execução de otimização."""

    pass


class SensitivityExecutionError(ExecutionError):
    """Erro na execução de análise de sensibilidade."""

    pass


# Exceções de Export
class ExportError(ApplicationError):
    """Erro de exportação."""

    pass


class UnsupportedFormatError(ExportError):
    """Formato de exportação não suportado."""

    pass


class ExportDestinationError(ExportError):
    """Erro no destino de exportação."""

    pass


# Exceções de Repositório
class RepositoryError(ApplicationError):
    """Erro de repositório."""

    pass


class RepositoryConnectionError(RepositoryError):
    """Erro de conexão com repositório."""

    pass


class RepositoryPermissionError(RepositoryError):
    """Erro de permissão no repositório."""

    pass
