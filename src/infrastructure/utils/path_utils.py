"""
Path utilities for consistent path handling across CSPBench.

This module provides utilities for handling paths from environment variables,
automatically detecting relative vs absolute paths and ensuring consistent
behavior across CLI and web interfaces.
"""

import os
from pathlib import Path
from typing import Union


def get_env_path(env_var: str, default: str, create_if_missing: bool = True) -> Path:
    """
    Get a path from environment variable with automatic relative/absolute detection.
    
    Args:
        env_var: Environment variable name
        default: Default path if environment variable is not set
        create_if_missing: Whether to create the directory if it doesn't exist
        
    Returns:
        Path object (resolved to absolute path)
        
    Examples:
        >>> get_env_path("DATASET_DIRECTORY", "./datasets")
        PosixPath('/workspaces/CSPBench/data/datasets')
        
        >>> get_env_path("OUTPUT_BASE_DIRECTORY", "./outputs") 
        PosixPath('/workspaces/CSPBench/data/outputs')
    """
    path_str = os.getenv(env_var, default)
    path = Path(path_str)
    
    # Convert to absolute path
    if not path.is_absolute():
        # Relative paths are resolved relative to current working directory
        path = path.resolve()
    
    # Create directory if requested and doesn't exist
    if create_if_missing:
        path.mkdir(parents=True, exist_ok=True)
    
    return path


def get_dataset_directory() -> Path:
    """Get the dataset directory from environment variable."""
    return get_env_path("DATASET_DIRECTORY", "./datasets")


def get_batch_directory() -> Path:
    """Get the batch directory from environment variable."""
    return get_env_path("BATCH_DIRECTORY", "./batches")


def get_output_base_directory() -> Path:
    """Get the output base directory from environment variable."""
    return get_env_path("OUTPUT_BASE_DIRECTORY", "./outputs")


def get_work_db_path() -> Path:
    """Get the work database path from environment variable."""
    db_path_str = os.getenv("WORK_DB_PATH", "./data/work_manager.db")
    db_path = Path(db_path_str)
    
    # Convert to absolute path
    if not db_path.is_absolute():
        db_path = db_path.resolve()
    
    # Create parent directory if it doesn't exist
    db_path.parent.mkdir(parents=True, exist_ok=True)
    
    return db_path


def is_relative_path(path: Union[str, Path]) -> bool:
    """
    Check if a path is relative.
    
    Args:
        path: Path string or Path object
        
    Returns:
        True if path is relative, False if absolute
    """
    return not Path(path).is_absolute()


def normalize_path(path: Union[str, Path], base_dir: Union[str, Path, None] = None) -> Path:
    """
    Normalize a path, converting relative paths to absolute.
    
    Args:
        path: Path to normalize
        base_dir: Base directory for relative paths (default: current working directory)
        
    Returns:
        Normalized absolute Path object
    """
    path_obj = Path(path)
    
    if path_obj.is_absolute():
        return path_obj
    
    if base_dir:
        return (Path(base_dir) / path_obj).resolve()
    
    return path_obj.resolve()
