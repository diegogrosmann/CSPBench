"""
Simplified test coverage for work_state persistence components.
"""

import json
import tempfile
import time
from pathlib import Path

import pytest

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


class TestWorkStatePersistenceBasic:
    """Basic test coverage for WorkStatePersistence components."""

    def test_submit_work_raw_dict(self, persistence):
        """Test submitting work item using low-level method."""
        # Create a simple work dict that can bypass WorkItem validation
        work_id = "test_work_123"

        # Use work_create method instead of direct SQL
        persistence.work_create(
            id=work_id,
            config={"test": "config"},
            status=BaseStatus.QUEUED.value,
            created_at=time.time(),
            updated_at=time.time(),
            output_path="/tmp/test_work",
            error="",
            extra={"test_key": "test_value"}
        )

        # Verify it was inserted
        work_data = persistence.work_get(work_id)
        assert work_data is not None
        work = WorkItem.from_dict(work_data)
        assert work.id == work_id

    def test_get_work_nonexistent(self, persistence):
        """Test getting non-existent work."""
        result = persistence.work_get("nonexistent")
        assert result is None

    def test_list_work_items(self, persistence):
        """Test listing work items."""
        work_id = "test_work_123"
        
        # Insert test work using work_create
        persistence.work_create(
            id=work_id,
            config={"test": "config"},
            status=BaseStatus.QUEUED.value,
            created_at=time.time(),
            updated_at=time.time(),
            output_path="/tmp/test_work",
            error="",
            extra={}
        )
        
        works, total = persistence.work_list()
        assert len(works) == 1
        work_items = [WorkItem.from_dict(w) for w in works]
        assert work_items[0].id == work_id

    def test_update_work_status(self, persistence):
        """Test updating work status."""
        work_id = "test_work_123"
        
        # Insert test work
        persistence.work_create(
            id=work_id,
            config={"test": "config"},
            status=BaseStatus.QUEUED.value,
            created_at=time.time(),
            updated_at=time.time(),
            output_path="/tmp/test_work",
            error="",
            extra={}
        )
        
        persistence.work_update(work_id, status=BaseStatus.RUNNING.value)

        work_data = persistence.work_get(work_id)
        work = WorkItem.from_dict(work_data)
        assert work.status == BaseStatus.RUNNING.value

    def test_log_event(self, persistence):
        """Test logging events."""
        work_id = "test_work_123"
        
        # Insert test work first
        persistence.work_create(
            id=work_id,
            config={"test": "config"},
            status=BaseStatus.QUEUED.value,
            created_at=time.time(),
            updated_at=time.time(),
            output_path="/tmp/test_work",
            error="",
            extra={}
        )
        
        persistence.event_create(
            work_id=work_id,
            event_type="error",
            event_category="work",
            entity_data={"message": "Test event"}
        )
        
        events, total = persistence.event_list(filters={"work_id": work_id})
        assert len(events) >= 1
        # Find the test event we just created
        test_events = [e for e in events if e.get("event_type") == "error"]
        assert len(test_events) >= 1

    def test_count_operations(self, persistence):
        """Test count operations."""
        work_id = "test_work_123"
        
        # Insert test work
        persistence.work_create(
            id=work_id,
            config={"test": "config"},
            status=BaseStatus.QUEUED.value,
            created_at=time.time(),
            updated_at=time.time(),
            output_path="/tmp/test_work",
            error="",
            extra={}
        )
        
        work_list, total_count = persistence.work_list()
        assert total_count == 1

    @pytest.mark.skip("backup_database method not implemented")
    def test_backup_database(self, persistence, temp_db):
        """Test database backup functionality."""
        backup_path = temp_db.parent / f"{temp_db.name}.backup"
        
        result = persistence.backup_database(backup_path)
        assert result == backup_path  # Returns path, not boolean
        assert backup_path.exists()
        
        # Cleanup
        backup_path.unlink()

    @pytest.mark.skip("get_database_info method not implemented")
    def test_get_database_info(self, persistence):
        """Test get_database_info method."""
        info = persistence.get_database_info()
        assert isinstance(info, dict)
        assert "db_path" in info

    def test_vacuum_database(self, persistence):
        """Test database vacuum operation."""
        result = persistence.vacuum_database()
        assert result is None  # Method returns None

    def test_repr_method(self, persistence, temp_db):
        """Test string representation."""
        repr_str = repr(persistence)
        assert str(temp_db) in repr_str


class TestWorkScopedPersistenceBasic:
    """Basic test coverage for WorkScopedPersistence wrapper."""

    def test_initialization(self, persistence):
        """Test WorkScopedPersistence initialization."""
        from src.infrastructure.persistence.work_state.wrappers.work_scoped import WorkScopedPersistence
        
        work_id = "test_work_123"
        
        # Insert test work
        persistence.work_create(
            id=work_id,
            config={"test": "config"},
            status=BaseStatus.QUEUED.value,
            created_at=time.time(),
            updated_at=time.time(),
            output_path="/tmp/test_work",
            error="",
            extra={}
        )
        
        scoped = WorkScopedPersistence(work_id, persistence)
        
        assert scoped.work_id == work_id
        assert scoped.store is persistence

    def test_initialization_nonexistent_work(self, persistence):
        """Test initialization with non-existent work ID."""
        from src.infrastructure.persistence.work_state.wrappers.work_scoped import WorkScopedPersistence
        
        # WorkScopedPersistence doesn't validate work existence on init
        scoped = WorkScopedPersistence("nonexistent", persistence)
        assert scoped.work_id == "nonexistent"


class TestCombinationScopedPersistenceBasic:
    """Basic test coverage for CombinationScopedPersistence wrapper."""

    def test_initialization_nonexistent_combination(self, persistence):
        """Test initialization with non-existent combination ID using factory method."""
        from src.infrastructure.persistence.work_state.wrappers.work_scoped import WorkScopedPersistence
        
        # Need a work to use factory method
        work_id = "test_work_123"
        persistence.work_create(
            id=work_id,
            config={"test": "config"},
            status=BaseStatus.QUEUED.value,
            created_at=time.time(),
            updated_at=time.time(),
            output_path="/tmp/test_work",
            error="",
            extra={}
        )
        
        work_store = WorkScopedPersistence(work_id, persistence)
        with pytest.raises(ValueError, match="Combinação id=999 não encontrada"):
            work_store.for_combination(999)


class TestExecutionScopedPersistenceBasic:
    """Basic test coverage for ExecutionScopedPersistence wrapper."""

    def test_initialization_nonexistent_execution(self, persistence):
        """Test initialization with non-existent execution ID using factory method."""
        from src.infrastructure.persistence.work_state.wrappers.work_scoped import WorkScopedPersistence
        
        # Need a work to use factory method
        work_id = "test_work_456"
        persistence.work_create(
            id=work_id,
            config={"test": "config"},
            status=BaseStatus.QUEUED.value,
            created_at=time.time(),
            updated_at=time.time(),
            output_path="/tmp/test_work",
            error="",
            extra={}
        )
        
        work_store = WorkScopedPersistence(work_id, persistence)
        with pytest.raises(ValueError, match="Unidade de execução unit_id=unit_999.*não encontrada"):
            work_store.for_execution("unit_999")


class TestMixinCoverage:
    """Test coverage for individual mixin methods."""

    def test_dataset_mixin_basic(self, persistence):
        """Test basic DatasetMixin functionality."""
        # Test method existence
        assert hasattr(persistence, 'submit_dataset')

    def test_events_mixin_basic(self, persistence):
        """Test basic EventsMixin functionality."""
        # Test with empty work_id - should return empty list
        events = persistence.get_events(work_id="nonexistent")
        assert len(events) == 0

    def test_combinations_mixin_basic(self, persistence):
        """Test basic CombinationsMixin functionality."""
        # Test method existence
        assert hasattr(persistence, 'get_combinations')

    def test_executions_mixin_basic(self, persistence):
        """Test basic ExecutionsMixin functionality."""
        # Test method existence - check for core execution methods
        assert hasattr(persistence, 'get_executions')

    def test_work_mixin_alias_methods(self, persistence):
        """Test WorkMixin alias methods."""
        work_id = "test_work_123"
        
        # Insert test work directly
        persistence._execute(
            """
            INSERT INTO work(id, config_json, status, created_at, updated_at, output_path, error, extra_json)
            VALUES(?,?,?,?,?,?,?,?)
            """,
            (
                work_id,
                json.dumps({"test": "config"}),
                BaseStatus.QUEUED.value,
                time.time(),
                time.time(),
                "/tmp/test_work",
                "",
                json.dumps({})
            ),
        )
        
        # Verify work exists via work_get()
        work_data = persistence.work_get(work_id)
        work = WorkItem.from_dict(work_data)
        assert work.id == work_id
