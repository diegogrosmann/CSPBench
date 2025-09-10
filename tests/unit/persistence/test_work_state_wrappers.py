"""
Test coverage for work_state persistence wrappers.
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
from src.infrastructure.persistence.work_state.wrappers.work_scoped import WorkScopedPersistence
from src.infrastructure.persistence.work_state.wrappers.combination_scoped import CombinationScopedPersistence
from src.infrastructure.persistence.work_state.wrappers.execution_scoped import ExecutionScopedPersistence


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


@pytest.fixture
def work_scoped_persistence(persistence, sample_work_item):
    """Create WorkScopedPersistence instance with submitted work."""
    persistence.submit_work(sample_work_item)
    return WorkScopedPersistence(sample_work_item.id, persistence)


@pytest.fixture
def combination_data(sample_work_item):
    """Create sample combination data."""
    return {
        "work_id": sample_work_item.id,
        "task_id": "task1",
        "dataset_id": "dataset1",
        "preset_id": "preset1",
        "algorithm_id": "algorithm1", 
        "mode": "test",
        "status": "queued",
        "total_sequences": 10
    }


@pytest.fixture
def combination_scoped_persistence(persistence, sample_work_item, combination_data):
    """Create CombinationScopedPersistence instance with submitted combination using factory method."""
    persistence.submit_work(sample_work_item)
    persistence.submit_combination(combination_data)
    
    combinations = persistence.get_combinations(work_id=sample_work_item.id)
    combination_id = combinations[0]["id"]
    
    work_store = WorkScopedPersistence(sample_work_item.id, persistence)
    return work_store.for_combination(combination_id)


@pytest.fixture
def execution_data():
    """Create sample execution data."""
    return {
        "unit_id": "unit_123",
        "sequencia": 0,
        "status": "queued",
        "params_json": json.dumps({"param1": "value1"})
    }


@pytest.fixture
def execution_scoped_persistence(persistence, sample_work_item, combination_data, execution_data):
    """Create ExecutionScopedPersistence instance with submitted execution using factory method."""
    persistence.submit_work(sample_work_item)
    persistence.submit_combination(combination_data)
    
    combinations = persistence.get_combinations(work_id=sample_work_item.id)
    combination_id = combinations[0]["id"]
    
    execution_data["combination_id"] = combination_id
    persistence.submit_execution(execution_data)
    
    executions = persistence.get_executions_by_combination(combination_id)
    execution_id = executions[0]["id"]
    
    work_store = WorkScopedPersistence(sample_work_item.id, persistence)
    return work_store.for_execution(execution_data["unit_id"])


class TestWorkScopedPersistence:
    """Test WorkScopedPersistence wrapper."""

    def test_initialization(self, persistence, sample_work_item):
        """Test proper initialization of work-scoped persistence."""
        persistence.submit_work(sample_work_item)
        
        scoped = WorkScopedPersistence(sample_work_item.id, persistence)
        assert scoped.work_id == sample_work_item.id
        assert scoped.store is persistence

    def test_initialization_nonexistent_work(self, persistence):
        """Test initialization with non-existent work ID."""
        # WorkScopedPersistence doesn't validate work existence on init
        scoped = WorkScopedPersistence("nonexistent", persistence)
        assert scoped.work_id == "nonexistent"

    def test_get_work_info(self, work_scoped_persistence, sample_work_item):
        """Test getting work information."""
        work_info = work_scoped_persistence.get_work_info()
        
        assert work_info["id"] == sample_work_item.id
        assert work_info["status"] == BaseStatus.QUEUED.value

    def test_log_event(self, work_scoped_persistence):
        """Test logging events through scoped persistence."""
        work_scoped_persistence.log_event(
            event_type="progress",
            event_category="work",
            entity_data={"message": "Test event"}
        )
        
        events = work_scoped_persistence.get_events()
        assert len(events) == 1
        assert events[0]["event_type"] == "progress"

    def test_get_events(self, work_scoped_persistence):
        """Test getting events for scoped work."""
        # Log multiple events
        work_scoped_persistence.log_event("progress", "work", {"step": 1})
        work_scoped_persistence.log_event("error", "work", {"error": "test error"})
        
        events = work_scoped_persistence.get_events()
        assert len(events) == 2

    def test_get_combinations(self, work_scoped_persistence, combination_data):
        """Test getting combinations for scoped work."""
        # Initially no combinations
        combinations = work_scoped_persistence.get_combinations()
        assert len(combinations) == 0
        
        # Add combination
        work_scoped_persistence._persistence.submit_combination(combination_data)
        
        # Verify combination is returned
        combinations = work_scoped_persistence.get_combinations()
        assert len(combinations) == 1
        assert combinations[0]["task_id"] == "task1"

    def test_update_status(self, work_scoped_persistence):
        """Test updating work status through scoped persistence."""
        work_scoped_persistence.update_status(BaseStatus.RUNNING.value)
        
        work_info = work_scoped_persistence.get_work_info()
        assert work_info["status"] == BaseStatus.RUNNING.value

    def test_get_statistics(self, work_scoped_persistence):
        """Test getting work statistics through scoped persistence."""
        stats = work_scoped_persistence.get_statistics()
        
        assert "work_id" in stats
        assert "status" in stats
        assert "runtime_seconds" in stats

    def test_get_config(self, work_scoped_persistence, sample_config):
        """Test getting work configuration through scoped persistence."""
        config = work_scoped_persistence.get_config()
        
        assert config is not None
        assert config.datasets == sample_config.datasets
        assert config.algorithms == sample_config.algorithms

    def test_get_output_path(self, work_scoped_persistence, sample_work_item):
        """Test getting output path through scoped persistence."""
        output_path = work_scoped_persistence.get_output_path()
        
        assert output_path is not None
        assert str(output_path) == sample_work_item.output_path


class TestCombinationScopedPersistence:
    """Test CombinationScopedPersistence wrapper."""

    def test_initialization(self, persistence, sample_work_item, combination_data):
        """Test proper initialization of combination-scoped persistence using factory method."""
        persistence.submit_work(sample_work_item)
        persistence.submit_combination(combination_data)
        
        combinations = persistence.get_combinations(work_id=sample_work_item.id)
        combination_id = combinations[0]["id"]
        
        work_store = WorkScopedPersistence(sample_work_item.id, persistence)
        scoped = work_store.for_combination(combination_id)
        assert scoped.combination_id == combination_id
        assert scoped.store is persistence

    def test_initialization_nonexistent_combination(self, persistence, sample_work_item):
        """Test initialization with non-existent combination ID."""
        persistence.submit_work(sample_work_item)
        work_store = WorkScopedPersistence(sample_work_item.id, persistence)
        
        with pytest.raises(ValueError, match="Combinação id=999 não encontrada"):
            work_store.for_combination(999)

    def test_get_combination_info(self, combination_scoped_persistence):
        """Test getting combination information."""
        combo_info = combination_scoped_persistence.get_combination_info()
        
        assert combo_info["task_id"] == "task1"
        assert combo_info["status"] == "queued"

    def test_update_status(self, combination_scoped_persistence):
        """Test updating combination status."""
        combination_scoped_persistence.update_status("running")
        
        combo_info = combination_scoped_persistence.get_combination_info()
        assert combo_info["status"] == "running"

    def test_get_executions(self, combination_scoped_persistence, execution_data):
        """Test getting executions for scoped combination."""
        # Initially no executions
        executions = combination_scoped_persistence.get_executions()
        assert len(executions) == 0
        
        # Add execution
        execution_data["combination_id"] = combination_scoped_persistence.combination_id
        combination_scoped_persistence._persistence.submit_execution(execution_data)
        
        # Verify execution is returned
        executions = combination_scoped_persistence.get_executions()
        assert len(executions) == 1
        assert executions[0]["unit_id"] == "unit_123"

    def test_log_event(self, combination_scoped_persistence):
        """Test logging events through scoped persistence."""
        # Get work_id for the combination
        combo_info = combination_scoped_persistence.get_combination_info()
        work_id = combo_info["work_id"]
        
        combination_scoped_persistence.log_event(
            work_id=work_id,
            event_type="progress",
            event_category="combination",
            entity_data={"combination_id": combination_scoped_persistence.combination_id}
        )
        
        events = combination_scoped_persistence._persistence.get_events(work_id=work_id)
        assert len(events) == 1
        assert events[0]["event_category"] == "combination"

    def test_submit_execution(self, combination_scoped_persistence):
        """Test submitting execution through scoped persistence."""
        execution_data = {
            "unit_id": "unit_456",
            "combination_id": combination_scoped_persistence.combination_id,
            "sequencia": 1,
            "status": "queued",
            "params_json": json.dumps({"param2": "value2"})
        }
        
        combination_scoped_persistence.submit_execution(execution_data)
        
        executions = combination_scoped_persistence.get_executions()
        assert len(executions) == 1
        assert executions[0]["unit_id"] == "unit_456"

    def test_get_progress_summary(self, combination_scoped_persistence, execution_data):
        """Test getting progress summary for combination."""
        # Add multiple executions with different statuses
        for i, status in enumerate(["queued", "running", "completed"]):
            exec_data = execution_data.copy()
            exec_data["unit_id"] = f"unit_{i}"
            exec_data["combination_id"] = combination_scoped_persistence.combination_id
            exec_data["sequencia"] = i
            exec_data["status"] = status
            combination_scoped_persistence._persistence.submit_execution(exec_data)
        
        summary = combination_scoped_persistence.get_progress_summary()
        
        assert "total_executions" in summary
        assert "status_counts" in summary
        assert summary["total_executions"] == 3


class TestExecutionScopedPersistence:
    """Test ExecutionScopedPersistence wrapper."""

    def test_initialization(self, persistence, sample_work_item, combination_data, execution_data):
        """Test proper initialization of execution-scoped persistence using factory method."""
        persistence.submit_work(sample_work_item)
        persistence.submit_combination(combination_data)
        
        combinations = persistence.get_combinations(work_id=sample_work_item.id)
        combination_id = combinations[0]["id"]
        
        execution_data["combination_id"] = combination_id
        persistence.submit_execution(execution_data)
        
        executions = persistence.get_executions_by_combination(combination_id)
        execution_id = executions[0]["id"]
        unit_id = executions[0]["unit_id"]  # Get unit_id from execution
        
        work_store = WorkScopedPersistence(sample_work_item.id, persistence)
        scoped = work_store.for_execution(unit_id)
        assert scoped.execution_id == execution_id

    def test_initialization_nonexistent_execution(self, persistence, sample_work_item):
        """Test initialization with non-existent execution ID."""
        persistence.submit_work(sample_work_item)
        work_store = WorkScopedPersistence(sample_work_item.id, persistence)
        
        with pytest.raises(ValueError, match="Unidade de execução unit_id=unit_999.*não encontrada"):
            work_store.for_execution("unit_999")

    def test_get_execution_info(self, execution_scoped_persistence):
        """Test getting execution information."""
        exec_info = execution_scoped_persistence.get_execution_info()
        
        assert exec_info["unit_id"] == "unit_123"
        assert exec_info["status"] == "queued"

    def test_update_status(self, execution_scoped_persistence):
        """Test updating execution status."""
        execution_scoped_persistence.update_status("running")
        
        exec_info = execution_scoped_persistence.get_execution_info()
        assert exec_info["status"] == "running"

    def test_update_progress(self, execution_scoped_persistence):
        """Test updating execution progress."""
        execution_scoped_persistence.update_progress(0.5, "Halfway complete")
        
        progress = execution_scoped_persistence.get_progress()
        assert len(progress) == 1
        assert progress[0]["progress"] == 0.5
        assert progress[0]["message"] == "Halfway complete"

    def test_get_progress(self, execution_scoped_persistence):
        """Test getting execution progress."""
        # Initially no progress
        progress = execution_scoped_persistence.get_progress()
        assert len(progress) == 0
        
        # Add progress entries
        execution_scoped_persistence.update_progress(0.3, "30% complete")
        execution_scoped_persistence.update_progress(0.7, "70% complete")
        
        progress = execution_scoped_persistence.get_progress()
        assert len(progress) == 2

    def test_set_result(self, execution_scoped_persistence):
        """Test setting execution result."""
        result_data = {"objective": 95.5, "details": "Optimal solution found"}
        
        execution_scoped_persistence.set_result(result_data, 95.5)
        
        exec_info = execution_scoped_persistence.get_execution_info()
        result_json = json.loads(exec_info["result_json"])
        assert result_json == result_data
        assert exec_info["objective"] == 95.5

    def test_set_params(self, execution_scoped_persistence):
        """Test setting execution parameters."""
        params = {"max_iterations": 1000, "population_size": 50}
        
        execution_scoped_persistence.set_params(params)
        
        exec_info = execution_scoped_persistence.get_execution_info()
        params_json = json.loads(exec_info["params_json"])
        assert params_json == params

    def test_get_combination_info(self, execution_scoped_persistence):
        """Test getting combination info through execution-scoped persistence."""
        combo_info = execution_scoped_persistence.get_combination_info()
        
        assert combo_info["task_id"] == "task1"
        assert "work_id" in combo_info

    def test_log_event(self, execution_scoped_persistence):
        """Test logging events through execution-scoped persistence."""
        # Get work_id through combination
        combo_info = execution_scoped_persistence.get_combination_info()
        work_id = combo_info["work_id"]
        
        execution_scoped_persistence.log_event(
            work_id=work_id,
            event_type="progress",
            event_category="unit",
            entity_data={"execution_id": execution_scoped_persistence.execution_id}
        )
        
        events = execution_scoped_persistence._persistence.get_events(work_id=work_id)
        assert len(events) == 1
        assert events[0]["event_category"] == "unit"


class TestWrapperIntegration:
    """Integration tests for persistence wrappers."""

    def test_work_to_combination_scoped_transition(self, persistence, sample_work_item, combination_data):
        """Test transitioning from work-scoped to combination-scoped persistence using factory methods."""
        # Start with work-scoped
        persistence.submit_work(sample_work_item)
        work_scoped = WorkScopedPersistence(sample_work_item.id, persistence)
        
        # Submit combination through work-scoped
        persistence.submit_combination(combination_data)
        
        # Get combination and create combination-scoped using factory method
        combinations = work_scoped.get_combinations()
        combination_id = combinations[0]["id"]
        
        combo_scoped = work_scoped.for_combination(combination_id)
        
        # Verify both scopes work
        assert work_scoped.get_work_info()["id"] == sample_work_item.id
        assert combo_scoped.get_combination_info()["task_id"] == "task1"

    def test_combination_to_execution_scoped_transition(self, combination_scoped_persistence, execution_data):
        """Test transitioning from combination-scoped to execution-scoped persistence using factory methods."""
        # Submit execution through combination-scoped
        execution_data["combination_id"] = combination_scoped_persistence.combination_id
        combination_scoped_persistence.submit_execution(execution_data)
        
        # Get execution and create execution-scoped using factory method
        executions = combination_scoped_persistence.get_executions()
        unit_id = executions[0]["unit_id"]
        
        exec_scoped = combination_scoped_persistence.for_execution(unit_id)
        
        # Verify both scopes work
        assert combination_scoped_persistence.get_combination_info()["task_id"] == "task1"
        assert exec_scoped.get_execution_info()["unit_id"] == "unit_123"

    def test_cross_scope_event_logging(self, persistence, sample_work_item, combination_data, execution_data):
        """Test event logging across different scopes using factory methods."""
        # Set up all scopes
        persistence.submit_work(sample_work_item)
        persistence.submit_combination(combination_data)
        
        combinations = persistence.get_combinations(work_id=sample_work_item.id)
        combination_id = combinations[0]["id"]
        
        execution_data["combination_id"] = combination_id
        persistence.submit_execution(execution_data)
        
        executions = persistence.get_executions_by_combination(combination_id)
        execution_id = executions[0]["id"]
        
        # Create scoped persistences using factory methods
        work_scoped = WorkScopedPersistence(sample_work_item.id, persistence)
        combo_scoped = work_scoped.for_combination(combination_id)
        exec_scoped = work_scoped.for_execution("unit_789")
        
        # Log events from different scopes
        work_scoped.log_event("progress", "work", {"message": "Work started"})
        combo_scoped.log_event(sample_work_item.id, "progress", "combination", {"combination_id": combination_id})
        exec_scoped.log_event(sample_work_item.id, "progress", "unit", {"execution_id": execution_id})
        
        # Verify all events are logged for the work
        events = work_scoped.get_events()
        assert len(events) == 3
        
        event_categories = [e["event_category"] for e in events]
        assert "work" in event_categories
        assert "combination" in event_categories
        assert "unit" in event_categories
