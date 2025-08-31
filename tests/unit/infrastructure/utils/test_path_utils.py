"""
Tests for path utilities to ensure consistent path handling.
"""

import os
import tempfile
from pathlib import Path
from unittest import mock

import pytest

from src.infrastructure.utils.path_utils import (
    get_env_path,
    get_dataset_directory,
    get_batch_directory,
    get_output_base_directory,
    get_work_db_path,
    is_relative_path,
    normalize_path,
)


class TestPathUtils:
    """Test path utility functions."""

    def test_get_env_path_with_absolute_path(self):
        """Test get_env_path with absolute path."""
        with tempfile.TemporaryDirectory() as temp_dir:
            test_path = Path(temp_dir) / "test_dir"
            
            with mock.patch.dict(os.environ, {"TEST_PATH": str(test_path)}):
                result = get_env_path("TEST_PATH", "./default")
                
                assert result.is_absolute()
                assert result.exists()  # Should be created
                assert result == test_path.resolve()

    def test_get_env_path_with_relative_path(self):
        """Test get_env_path with relative path."""
        with mock.patch.dict(os.environ, {"TEST_PATH": "./relative/path"}):
            result = get_env_path("TEST_PATH", "./default")
            
            assert result.is_absolute()
            assert result.exists()  # Should be created
            # Should be resolved relative to cwd
            expected = (Path.cwd() / "relative/path").resolve()
            assert result == expected

    def test_get_env_path_with_default(self):
        """Test get_env_path falls back to default when env var not set."""
        # Ensure TEST_PATH is not set
        env_vars = {k: v for k, v in os.environ.items() if k != "TEST_PATH"}
        
        with mock.patch.dict(os.environ, env_vars, clear=True):
            with tempfile.TemporaryDirectory() as temp_dir:
                default_path = Path(temp_dir) / "default_dir"
                
                result = get_env_path("TEST_PATH", str(default_path))
                
                assert result.is_absolute()
                assert result.exists()  # Should be created
                assert result == default_path.resolve()

    def test_get_env_path_no_create(self):
        """Test get_env_path with create_if_missing=False."""
        with tempfile.TemporaryDirectory() as temp_dir:
            test_path = Path(temp_dir) / "nonexistent"
            
            with mock.patch.dict(os.environ, {"TEST_PATH": str(test_path)}):
                result = get_env_path("TEST_PATH", "./default", create_if_missing=False)
                
                assert result.is_absolute()
                assert not result.exists()  # Should not be created
                assert result == test_path.resolve()

    def test_get_dataset_directory(self):
        """Test get_dataset_directory uses correct env var."""
        with tempfile.TemporaryDirectory() as temp_dir:
            dataset_path = Path(temp_dir) / "datasets"
            
            with mock.patch.dict(os.environ, {"DATASET_DIRECTORY": str(dataset_path)}):
                result = get_dataset_directory()
                
                assert result == dataset_path.resolve()
                assert result.exists()

    def test_get_batch_directory(self):
        """Test get_batch_directory uses correct env var."""
        with tempfile.TemporaryDirectory() as temp_dir:
            batch_path = Path(temp_dir) / "batches"
            
            with mock.patch.dict(os.environ, {"BATCH_DIRECTORY": str(batch_path)}):
                result = get_batch_directory()
                
                assert result == batch_path.resolve()
                assert result.exists()

    def test_get_output_base_directory(self):
        """Test get_output_base_directory uses correct env var."""
        with tempfile.TemporaryDirectory() as temp_dir:
            output_path = Path(temp_dir) / "outputs"
            
            with mock.patch.dict(os.environ, {"OUTPUT_BASE_DIRECTORY": str(output_path)}):
                result = get_output_base_directory()
                
                assert result == output_path.resolve()
                assert result.exists()

    def test_get_work_db_path(self):
        """Test get_work_db_path uses correct env var and creates parent dir."""
        with tempfile.TemporaryDirectory() as temp_dir:
            db_path = Path(temp_dir) / "subdir" / "work.db"
            
            with mock.patch.dict(os.environ, {"WORK_DB_PATH": str(db_path)}):
                result = get_work_db_path()
                
                assert result == db_path.resolve()
                assert result.parent.exists()  # Parent directory should be created
                # Note: db file itself should not exist yet

    def test_is_relative_path(self):
        """Test is_relative_path detection."""
        assert is_relative_path("./relative/path")
        assert is_relative_path("relative/path")
        assert is_relative_path("../parent/path")
        
        assert not is_relative_path("/absolute/path")
        assert not is_relative_path(Path("/absolute/path"))

    def test_normalize_path_absolute(self):
        """Test normalize_path with absolute path."""
        abs_path = Path("/absolute/path")
        result = normalize_path(abs_path)
        
        assert result.is_absolute()
        assert result == abs_path

    def test_normalize_path_relative_no_base(self):
        """Test normalize_path with relative path and no base."""
        rel_path = Path("relative/path")
        result = normalize_path(rel_path)
        
        assert result.is_absolute()
        assert result == (Path.cwd() / rel_path).resolve()

    def test_normalize_path_relative_with_base(self):
        """Test normalize_path with relative path and base directory."""
        with tempfile.TemporaryDirectory() as temp_dir:
            base_dir = Path(temp_dir)
            rel_path = Path("relative/path")
            
            result = normalize_path(rel_path, base_dir)
            
            assert result.is_absolute()
            assert result == (base_dir / rel_path).resolve()

    def test_path_consistency_across_functions(self):
        """Test that all path functions return consistent absolute paths."""
        # Set relative paths in environment
        test_env = {
            "DATASET_DIRECTORY": "./test_datasets",
            "BATCH_DIRECTORY": "./test_batches", 
            "OUTPUT_BASE_DIRECTORY": "./test_outputs",
            "WORK_DB_PATH": "./test_data/work.db",
        }
        
        with mock.patch.dict(os.environ, test_env):
            dataset_dir = get_dataset_directory()
            batch_dir = get_batch_directory()
            output_dir = get_output_base_directory()
            work_db = get_work_db_path()
            
            # All should be absolute paths
            assert dataset_dir.is_absolute()
            assert batch_dir.is_absolute()
            assert output_dir.is_absolute()
            assert work_db.is_absolute()
            
            # All should exist (except db file)
            assert dataset_dir.exists()
            assert batch_dir.exists()
            assert output_dir.exists()
            assert work_db.parent.exists()
            
            # Should be resolved relative to current working directory
            cwd = Path.cwd()
            assert dataset_dir == (cwd / "test_datasets").resolve()
            assert batch_dir == (cwd / "test_batches").resolve()
            assert output_dir == (cwd / "test_outputs").resolve()
            assert work_db == (cwd / "test_data/work.db").resolve()
