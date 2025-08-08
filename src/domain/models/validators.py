"""
Domain Model Validators

Utility functions for validating domain models.
"""

from typing import Any, Dict, List
from .enums import DatasetType, BatchType


def validate_dataset_parameters(dataset_type: DatasetType, parameters: Dict[str, Any]) -> None:
    """
    Validate dataset parameters based on type.
    
    Args:
        dataset_type: Type of dataset to validate
        parameters: Dictionary of parameters to validate
        
    Raises:
        ValueError: If required parameters are missing or invalid
    """
    if dataset_type == DatasetType.SYNTHETIC:
        required = {"n_strings", "string_length", "alphabet"}
        missing = required - set(parameters.keys())
        if missing:
            raise ValueError(f"Synthetic dataset missing: {missing}")
        
        # Validar tipos e valores
        if not isinstance(parameters.get("n_strings"), int) or parameters["n_strings"] <= 0:
            raise ValueError("n_strings must be positive integer")
        if not isinstance(parameters.get("string_length"), int) or parameters["string_length"] <= 0:
            raise ValueError("string_length must be positive integer")
        if not isinstance(parameters.get("alphabet"), str) or len(parameters["alphabet"]) == 0:
            raise ValueError("alphabet must be non-empty string")
    
    elif dataset_type == DatasetType.FILE:
        if "filepath" not in parameters:
            raise ValueError("File dataset requires 'filepath'")
        if not isinstance(parameters["filepath"], str):
            raise ValueError("filepath must be string")
    
    elif dataset_type == DatasetType.ENTREZ:
        required = {"database", "query"}
        missing = required - set(parameters.keys())
        if missing:
            raise ValueError(f"Entrez dataset missing: {missing}")


def validate_repetitions(repetitions: int) -> None:
    """
    Validate repetitions parameter.
    
    Args:
        repetitions: Number of repetitions to validate
        
    Raises:
        ValueError: If repetitions is not a positive integer
    """
    if not isinstance(repetitions, int) or repetitions <= 0:
        raise ValueError("repetitions must be positive integer")


def validate_timeout(timeout: int, field_name: str = "timeout") -> None:
    """
    Validate timeout parameters.
    
    Args:
        timeout: Timeout value in seconds
        field_name: Name of the field being validated (for error messages)
        
    Raises:
        ValueError: If timeout is not a positive integer
    """
    if not isinstance(timeout, int) or timeout <= 0:
        raise ValueError(f"{field_name} must be positive integer")


def validate_memory_limit(memory_gb: float, field_name: str = "memory_limit") -> None:
    """
    Validate memory limit parameters.
    
    Args:
        memory_gb: Memory limit in GB
        field_name: Name of the field being validated (for error messages)
        
    Raises:
        ValueError: If memory limit is not a positive number
    """
    if not isinstance(memory_gb, (int, float)) or memory_gb <= 0:
        raise ValueError(f"{field_name} must be positive number")


def validate_cores(cores: int, field_name: str = "cores") -> None:
    """
    Validate CPU cores parameter.
    
    Args:
        cores: Number of CPU cores
        field_name: Name of the field being validated (for error messages)
        
    Raises:
        ValueError: If cores is not a positive integer
    """
    if not isinstance(cores, int) or cores <= 0:
        raise ValueError(f"{field_name} must be positive integer")


def validate_n_samples(n_samples: int, minimum: int = 1) -> None:
    """
    Validate number of samples for sensitivity analysis.
    
    Args:
        n_samples: Number of samples
        minimum: Minimum required samples
        
    Raises:
        ValueError: If n_samples is not sufficient
    """
    if not isinstance(n_samples, int) or n_samples < minimum:
        raise ValueError(f"n_samples must be integer >= {minimum}")


def validate_export_formats(formats: List[str]) -> None:
    """
    Validate export formats list.
    
    Args:
        formats: List of export format strings
        
    Raises:
        ValueError: If formats list is empty or contains invalid formats
    """
    if not formats:
        raise ValueError("export formats list cannot be empty")
    
    valid_formats = {"csv", "json", "parquet", "pickle"}
    invalid = set(formats) - valid_formats
    if invalid:
        raise ValueError(f"Invalid export formats: {invalid}. Valid formats: {valid_formats}")
