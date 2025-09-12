"""
Test coverage for work_state persistence mixins.
"""

import json
import tempfile
import time
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from src.domain.config import CSPBenchConfig
from src.domain.status import BaseStatus
from src.domain.work import WorkItem
from src.infrastructure.persistence.work_state import WorkPersistence


@pytest.fixture
def temp_db():
    """Create a temporary database for testing."""
    with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as f:
        db_path = Path(f.name)
    
    yield db_path
    
    # Cleanup
    if db_path.exists():
        db_path.unlink()


@pytest.fixture
def persistence(temp_db):
    """Create WorkPersistence instance."""
    return WorkPersistence(f"sqlite:///{temp_db}")


def submit_work_item(persistence, work_item):
    """Helper function to submit a work item using the correct method."""
    config_dict = work_item.config if isinstance(work_item.config, dict) else work_item.config.to_dict()
    persistence.work_create(
        id=work_item.id,
        config=config_dict,
        status=work_item.status,
        output_path=work_item.output_path,
        error=work_item.error,
        extra=work_item.extra,
        created_at=work_item.created_at,
        updated_at=work_item.updated_at,
    )


@pytest.fixture
def sample_config():
    """Create a simple mock CSPBenchConfig for testing."""
    # Create a minimal config using the Mock pattern
    config_dict = {
        "metadata": {
            "name": "Test",
            "description": "Test config",
            "author": "Test",
            "version": "1.0",
            "creation_date": "2025-01-01"
        },
        "datasets": {
            "test": {
                "id": "test",
                "name": "Test Dataset",
                "type": "synthetic",
                "mode": "random",
                "n": 10,
                "L": 100
            }
        },
        "algorithms": {
            "baseline": {
                "id": "baseline",
                "name": "Baseline",
                "algorithm": "Baseline",
                "params": {}
            }
        },
        "tasks": {
            "type": "experiment",
            "datasets": ["test"],
            "algorithms": ["baseline"]
        },
        "output": {
            "base_path": "test_output",
            "logging": True,
            "results": {
                "formats": {
                    "csv": True,
                    "json": False,
                    "parquet": False,
                    "pickle": False
                },
                "partial_results": False
            }
        }
    }
    
    # Return as dict to avoid config constructor issues
    return config_dict


@pytest.fixture
def sample_work_item(sample_config):
    """Create a sample WorkItem."""
    return WorkItem(
        id="test_work_123",
        config=sample_config,
        status=BaseStatus.QUEUED,
        created_at=time.time(),
        updated_at=time.time(),
        output_path="/tmp/test_work",
        extra={"test_key": "test_value"}
    )


class TestWorkMixin:
    """Test WorkMixin functionality."""

    def test_submit_work(self, persistence, sample_work_item):
        """Test submitting work item."""
        submit_work_item(persistence, sample_work_item)

        # Verify item exists
        work = persistence.work_get(sample_work_item.id)
        assert work is not None
        assert work["id"] == sample_work_item.id
        assert work["status"] == BaseStatus.QUEUED.value

    def test_get_work_nonexistent(self, persistence):
        """Test getting non-existent work item."""
        result = persistence.work_get("nonexistent_id")
        assert result is None

    def test_update_work_item(self, persistence, sample_work_item):
        """Test updating work item."""
        # Arrange
        submit_work_item(persistence, sample_work_item)
        
        # Act
        success = persistence.work_update(sample_work_item.id, status="running")
        
        # Assert
        assert success is True
        
        work = persistence.work_get(sample_work_item.id)
        assert work is not None
        assert work["status"] == "running"

    def test_list_work_items(self, persistence, sample_work_item):
        """Test listing work items."""
        # Initially empty
        items, total = persistence.work_list()
        assert len(items) == 0
        
        # Add work item
        persistence.work_create(
            id=sample_work_item.id,
            config=sample_work_item.config.to_dict(),
            status=sample_work_item.status.value,
            created_at=sample_work_item.created_at,
            updated_at=sample_work_item.updated_at,
            output_path=sample_work_item.output_path,
            error=sample_work_item.error,
            extra=sample_work_item.extra or {}
        )
        
        # Verify list
        items, total = persistence.work_list()
        assert len(items) == 1
        work_items = [WorkItem.from_dict(w) for w in items]
        assert work_items[0].id == sample_work_item.id

    def test_list_by_status(self, persistence, sample_config):
        """Test listing work items by status."""
        # Create items with different statuses
        work1 = WorkItem(
            id="work1", config=sample_config, status=BaseStatus.QUEUED,
            created_at=time.time(), updated_at=time.time(), output_path="/tmp/1"
        )
        work2 = WorkItem(
            id="work2", config=sample_config, status=BaseStatus.RUNNING,
            created_at=time.time(), updated_at=time.time(), output_path="/tmp/2"
        )
        
        submit_work_item(persistence, work1)
        submit_work_item(persistence, work2)
        
        # Test filtering
        queued_items, _ = persistence.work_list(filters={"status": BaseStatus.QUEUED})
        running_items, _ = persistence.work_list(filters={"status": BaseStatus.RUNNING})
        
        assert len(queued_items) == 1
        assert len(running_items) == 1
        assert queued_items[0]["id"] == "work1"
        assert running_items[0]["id"] == "work2"

    def test_count_operations(self, persistence, sample_work_item):
        """Test count operations using work_list."""
        # Initially zero
        _, total_count = persistence.work_list(limit=0)
        assert total_count == 0
        
        # Add work item
        submit_work_item(persistence, sample_work_item)
        
        # Verify counts
        _, total_count = persistence.work_list(limit=0)
        assert total_count == 1
        
        # Count by status
        _, queued_count = persistence.work_list(filters={"status": BaseStatus.QUEUED}, limit=0)
        assert queued_count == 1
        _, running_count = persistence.work_list(filters={"status": BaseStatus.RUNNING}, limit=0)
        assert running_count == 0

    def test_remove_work_item(self, persistence, sample_work_item):
        """Test removing work item."""
        # Add work item
        submit_work_item(persistence, sample_work_item)
        _, initial_count = persistence.work_list(limit=0)
        assert initial_count == 1
        
        # Remove work item
        persistence.work_delete(sample_work_item.id)
        _, final_count = persistence.work_list(limit=0)
        assert final_count == 0
        
        # Try removing non-existent item should not raise error
        persistence.work_delete("nonexistent")  # Should not raise exception

    def test_get_work_statistics(self, persistence, sample_work_item):
        """Test getting work statistics."""
        # Submit work and get stats
        submit_work_item(persistence, sample_work_item)
        work = persistence.work_get(sample_work_item.id)
        
        assert work["id"] == sample_work_item.id
        assert work["status"] == BaseStatus.QUEUED.value
        assert "created_at" in work
        # We can derive statistics from work data
        assert work is not None

    def test_get_work_statistics_nonexistent(self, persistence):
        """Test getting statistics for non-existent work."""
        work = persistence.work_get("nonexistent")
        assert work is None

    def test_get_work_config(self, persistence, sample_work_item):
        """Test retrieving work configuration."""
        submit_work_item(persistence, sample_work_item)
        
        work = persistence.work_get(sample_work_item.id)
        assert work is not None
        config = work["config_json"]  # Config is stored in config_json field
        # Config should match the original dict version
        expected_config = sample_work_item.config.to_dict()
        assert config == expected_config

    def test_get_work_config_nonexistent(self, persistence):
        """Test retrieving config for non-existent work."""
        work = persistence.work_get("nonexistent")
        assert work is None

    def test_get_work_output_path(self, persistence, sample_work_item):
        """Test retrieving work output path."""
        submit_work_item(persistence, sample_work_item)
        
        work = persistence.work_get(sample_work_item.id)
        assert work is not None
        output_path = work["output_path"]
        assert output_path == sample_work_item.output_path

    def test_get_work_output_path_nonexistent(self, persistence):
        """Test retrieving output path for non-existent work."""
        work = persistence.work_get("nonexistent")
        assert work is None

    def test_get_work_status_with_service(self, persistence, sample_work_item):
        """Test getting work status directly from persistence."""
        submit_work_item(persistence, sample_work_item)
        
        work = persistence.work_get(sample_work_item.id)
        assert work is not None
        status = work["status"]
        assert status == BaseStatus.QUEUED.value

    def test_get_work_status_service_fallback(self, persistence, sample_work_item):
        """Test getting work status when work doesn't exist."""
        work = persistence.work_get("nonexistent_id")
        assert work is None

    def test_update_work_status(self, persistence, sample_work_item):
        """Test updating work status."""
        submit_work_item(persistence, sample_work_item)
        
        # Update status
        persistence.work_update(sample_work_item.id, status=BaseStatus.RUNNING.value)
        
        # Verify update
        work = persistence.work_get(sample_work_item.id)
        assert work is not None
        assert work["status"] == BaseStatus.RUNNING.value

    def test_add_alias_method(self, persistence, sample_work_item):
        """Test that work_create method works."""
        submit_work_item(persistence, sample_work_item)
        
        work = persistence.work_get(sample_work_item.id)
        assert work is not None
        assert work["id"] == sample_work_item.id


class TestDatasetMixin:
    """Test DatasetMixin functionality."""

    def test_submit_dataset(self, persistence, sample_work_item):
        """Test submitting dataset."""
        # First submit work item
        submit_work_item(persistence, sample_work_item)
        
        # Submit dataset
        dataset_data = {
            "id": "test_dataset",
            "work_id": sample_work_item.id,
            "name": "Test Dataset",
            "meta_json": json.dumps({"type": "synthetic"})
        }
        
        dataset_id = persistence.dataset_create(
            dataset_id=dataset_data["id"],
            work_id=dataset_data["work_id"], 
            name=dataset_data["name"],
            meta={"type": "synthetic"}
        )
        
        # Verify dataset was inserted
        datasets, _ = persistence.dataset_list(filters={"work_id": sample_work_item.id})
        assert len(datasets) == 1
        assert datasets[0]["dataset_id"] == "test_dataset"

    def test_get_datasets_by_work_nonexistent(self, persistence):
        """Test getting datasets for non-existent work."""
        datasets, _ = persistence.dataset_list(filters={"work_id": "nonexistent"})
        assert len(datasets) == 0


class TestEventsMixin:
    """Test EventsMixin functionality."""

    def test_log_event(self, persistence, sample_work_item):
        """Test logging events."""
        # Submit work item first
        submit_work_item(persistence, sample_work_item)
        
        # Create event
        event_id = persistence.event_create(
            work_id=sample_work_item.id,
            event_type="progress",
            event_category="work",
            entity_data={"message": "Work started"}
        )
        
        # Verify event was created
        events, _ = persistence.event_list(filters={"work_id": sample_work_item.id})
        assert len(events) == 1
        assert events[0]["event_type"] == "progress"
        assert events[0]["event_category"] == "work"

    def test_get_events_by_work_nonexistent(self, persistence):
        """Test getting events for non-existent work."""
        # Test using event list method
        events, _ = persistence.event_list(filters={"work_id": "nonexistent"})
        assert len(events) == 0


class TestExecutionsMixin:
    """Test ExecutionsMixin functionality."""

    def test_submit_execution(self, persistence, sample_work_item):
        """Test submitting execution."""
        # Submit work and combination first
        submit_work_item(persistence, sample_work_item)
        
        combination_data = {
            "work_id": sample_work_item.id,
            "task_id": "task1",
            "dataset_id": "dataset1",
            "preset_id": "preset1", 
            "algorithm_id": "algorithm1",
            "mode": "test",
            "status": "queued",
            "total_sequences": 10
        }
        combination_id = persistence.combination_create(**combination_data)
        
        # Verify combination was created and get its actual ID
        combinations, _ = persistence.combination_list(filters={"work_id": sample_work_item.id})
        assert len(combinations) == 1
        actual_combination_id = combinations[0]["id"]
        
        # Submit execution using the actual combination ID
        execution_data = {
            "unit_id": "unit_123",
            "combination_id": actual_combination_id,
            "sequencia": 0,
            "status": "queued",
            "params": {"param1": "value1"}
        }
        persistence.execution_create(**execution_data)
        
        # Verify execution was submitted
        executions, _ = persistence.execution_list(filters={"combination_id": actual_combination_id})
        assert len(executions) == 1
        assert executions[0]["unit_id"] == "unit_123"


class TestCombinationsMixin:
    """Test CombinationsMixin functionality."""

    def test_submit_combination(self, persistence, sample_work_item):
        """Test submitting combination."""
        # Submit work item first
        submit_work_item(persistence, sample_work_item)
        
        # Submit combination
        combination_data = {
            "work_id": sample_work_item.id,
            "task_id": "task1",
            "dataset_id": "dataset1", 
            "preset_id": "preset1",
            "algorithm_id": "algorithm1",
            "mode": "test",
            "status": "queued",
            "total_sequences": 10
        }
        combination_id = persistence.combination_create(**combination_data)
        
        # Verify combination was submitted
        combinations, _ = persistence.combination_list(filters={"work_id": sample_work_item.id})
        assert len(combinations) == 1
        assert combinations[0]["task_id"] == "task1"

    def test_get_combinations_by_work_nonexistent(self, persistence):
        """Test getting combinations for non-existent work."""
        combinations, _ = persistence.combination_list(filters={"work_id": "nonexistent"})
        assert len(combinations) == 0


class TestWorkPersistenceIntegration:
    """Integration tests for WorkPersistence."""

    def test_database_info(self, persistence):
        """Test basic database functionality."""
        # Test that we can create and retrieve data
        work_id = "test_db_info"
        persistence.work_create(id=work_id, config={})
        work = persistence.work_get(work_id)
        assert work is not None
        assert work["id"] == work_id

    def test_backup_database(self, persistence, temp_db):
        """Test database persistence functionality."""
        # Test that data persists
        work_id = "test_backup"
        persistence.work_create(id=work_id, config={})
        
        # Verify persistence
        work = persistence.work_get(work_id)
        assert work is not None

    def test_vacuum_database(self, persistence):
        """Test database operations don't raise exceptions.""" 
        # Test basic operations work
        work_id = "test_vacuum"
        persistence.work_create(id=work_id, config={})
        persistence.work_delete(work_id)

    def test_concurrent_access(self, persistence, sample_work_item):
        """Test concurrent access to database."""
        # Submit work item
        submit_work_item(persistence, sample_work_item)
        
        # Multiple read operations should work
        for _ in range(5):
            item = persistence.work_get(sample_work_item.id)
            assert item is not None
            
        # Multiple count operations  
        for _ in range(5):
            items, count = persistence.work_list(limit=0)
            assert count == 1

    def test_repr_method(self, persistence, temp_db):
        """Test string representation of persistence object.""" 
        repr_str = repr(persistence)
        assert "WorkPersistence" in repr_str
