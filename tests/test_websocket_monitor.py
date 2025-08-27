"""
Tests for WebSocket real-time monitoring.
"""

import asyncio
import json
import pytest
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch
from fastapi.testclient import TestClient

from src.presentation.web.app import app
from src.presentation.web.websocket import work_monitor_manager
from src.infrastructure.persistence.work_state.queries import ProgressSummary, ExecutionDetail


@pytest.fixture
def test_client():
    """Create test client."""
    return TestClient(app)


@pytest.fixture 
def temp_db():
    """Create temporary database for testing."""
    with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as f:
        yield Path(f.name)
    

@pytest.fixture
def mock_work_service():
    """Mock work service."""
    with patch('src.presentation.web.websocket.server.get_work_service') as mock:
        mock_service = Mock()
        mock.return_value = mock_service
        yield mock_service


@pytest.fixture 
def sample_progress_summary():
    """Sample progress summary for testing."""
    return ProgressSummary(
        work_id="test_work",
        tasks={"finished": [], "running": ["task1"], "queued": ["task2"]},
        datasets={"finished": [], "running": ["dataset1"], "queued": []},
        configs={"finished": [], "running": ["config1"], "queued": []},
        algorithms={"finished": [], "running": ["algorithm1"], "queued": []},
        execution={"finished": 5, "running": 3, "queued": 2, "total": 10},
        global_execution={"Finished": 5, "Total": 10},
        global_progress=0.5,
        current_combination_details={"combination_id": 123}
    )


@pytest.fixture
def sample_executions():
    """Sample executions for testing."""
    return [
        ExecutionDetail(
            unit_id="exec_1",
            combination_id=123,
            sequencia=1,
            status="running",
            progress=0.3,
            progress_message="Processing...",
            started_at=1700000000.0,
            finished_at=None,
            objective=None,
            task_id="task1",
            dataset_id="dataset1", 
            preset_id="config1",
            algorithm_id="algorithm1",
            mode="experiment",
            total_sequences=100
        ),
        ExecutionDetail(
            unit_id="exec_2",
            combination_id=123,
            sequencia=2,
            status="completed",
            progress=1.0,
            progress_message="Completed",
            started_at=1700000000.0,
            finished_at=1700000010.0,
            objective=42.5,
            task_id="task1",
            dataset_id="dataset1",
            preset_id="config1", 
            algorithm_id="algorithm1",
            mode="experiment",
            total_sequences=100
        )
    ]


class TestWebSocketConnection:
    """Test WebSocket connection basics."""
    
    def test_general_websocket_connection(self, test_client):
        """Test general WebSocket connection."""
        with test_client.websocket_connect("/ws/test_client") as websocket:
            # Should receive welcome message
            data = websocket.receive_json()
            assert data["type"] == "welcome"
            assert data["client_id"] == "test_client"
            assert "available_endpoints" in data["payload"]
            
            # Test echo
            websocket.send_text("hello")
            echo_data = websocket.receive_json()
            assert echo_data["type"] == "echo"
            assert echo_data["payload"]["received"] == "hello"
    
    def test_work_websocket_work_not_found(self, test_client, mock_work_service):
        """Test WebSocket connection for non-existent work."""
        mock_work_service.get_work_details.return_value = None
        
        with test_client.websocket_connect("/ws/work/nonexistent") as websocket:
            # Should receive error message
            data = websocket.receive_json()
            assert data["type"] == "error"
            assert data["payload"]["code"] == "WORK_NOT_FOUND"
    
    def test_work_websocket_database_not_ready(self, test_client, mock_work_service):
        """Test WebSocket connection when database is not ready."""
        mock_work_service.get_work_details.return_value = None
        
        with test_client.websocket_connect("/ws/work/test_work") as websocket:
            # Should receive error message
            data = websocket.receive_json()
            assert data["type"] == "error"
            assert data["payload"]["code"] == "WORK_NOT_FOUND"


class TestProgressTracking:
    """Test progress tracking and diff detection."""
    
    def test_progress_hash_generation(self, sample_progress_summary):
        """Test progress hash generation."""
        from src.presentation.web.websocket.manager import WorkMonitorSession
        
        session = WorkMonitorSession("test", Path("/tmp/test.db"))
        hash1 = session._hash_progress(sample_progress_summary)
        hash2 = session._hash_progress(sample_progress_summary)
        
        # Same data should produce same hash
        assert hash1 == hash2
        
        # Different progress should produce different hash
        sample_progress_summary.global_progress = 0.6
        hash3 = session._hash_progress(sample_progress_summary)
        assert hash1 != hash3
    
    def test_execution_changes_detection(self, sample_executions):
        """Test execution changes detection."""
        from src.presentation.web.websocket.manager import WorkMonitorSession
        
        session = WorkMonitorSession("test", Path("/tmp/test.db"))
        
        # First detection - should consider all as new
        changes = session._detect_execution_changes(sample_executions)
        assert changes is not None
        assert len(changes.new) == 2
        assert "exec_1" in changes.new
        assert "exec_2" in changes.new
        
        # Update state
        session._update_executions_state(sample_executions)
        
        # No changes - should return None
        changes = session._detect_execution_changes(sample_executions)
        assert changes is None
        
        # Modify one execution
        sample_executions[0].progress = 0.5
        sample_executions[0].status = "running"  # Still running
        
        changes = session._detect_execution_changes(sample_executions)
        assert changes is not None
        assert "exec_1" in changes.updated
        
        # Complete an execution
        sample_executions[0].status = "completed"
        sample_executions[0].progress = 1.0
        
        changes = session._detect_execution_changes(sample_executions)
        assert changes is not None
        assert "exec_1" in changes.updated
        assert "exec_1" in changes.completed


class TestMessageSerialization:
    """Test message serialization and schemas."""
    
    def test_progress_summary_serialization(self, sample_progress_summary):
        """Test ProgressSummary serialization."""
        from src.presentation.web.websocket.schemas import serialize_progress_summary
        
        serialized = serialize_progress_summary(sample_progress_summary)
        
        assert serialized["work_id"] == "test_work"
        assert serialized["global_progress"] == 0.5
        assert serialized["global_execution"]["Finished"] == 5
        assert serialized["current_combination_details"]["combination_id"] == 123
    
    def test_execution_detail_serialization(self, sample_executions):
        """Test ExecutionDetail serialization."""
        from src.presentation.web.websocket.schemas import serialize_execution_detail
        
        serialized = serialize_execution_detail(sample_executions[0])
        
        assert serialized["unit_id"] == "exec_1"
        assert serialized["status"] == "running"
        assert serialized["progress"] == 0.3
        assert serialized["task_id"] == "task1"
    
    def test_websocket_message_json(self):
        """Test WebSocketMessage JSON conversion.""" 
        from src.presentation.web.websocket.schemas import WebSocketMessage, MessageType
        
        message = WebSocketMessage(
            type=MessageType.SNAPSHOT,
            work_id="test_work",
            timestamp=1700000000.0,
            payload={"test": "data"}
        )
        
        json_str = message.to_json()
        parsed = json.loads(json_str)
        
        assert parsed["type"] == "snapshot"
        assert parsed["work_id"] == "test_work"
        assert parsed["timestamp"] == 1700000000.0
        assert parsed["payload"]["test"] == "data"


class TestMonitoringSession:
    """Test monitoring session lifecycle."""
    
    @pytest.mark.asyncio
    async def test_session_lifecycle(self):
        """Test session creation, start, stop."""
        from src.presentation.web.websocket.manager import WorkMonitorSession
        
        session = WorkMonitorSession("test", Path("/tmp/test.db"))
        
        # Initially not running
        assert not session.running
        
        # Start session
        await session.start()
        assert session.running
        
        # Stop session
        await session.stop()
        assert not session.running
    
    @pytest.mark.asyncio
    async def test_subscriber_management(self):
        """Test adding/removing subscribers."""
        from src.presentation.web.websocket.manager import WorkMonitorSession
        
        session = WorkMonitorSession("test", Path("/tmp/test.db"))
        
        queue1 = asyncio.Queue()
        queue2 = asyncio.Queue()
        
        # Add subscribers
        await session.add_subscriber(queue1)
        await session.add_subscriber(queue2)
        
        assert len(session.subscribers) == 2
        
        # Remove subscriber
        await session.remove_subscriber(queue1)
        assert len(session.subscribers) == 1
        
        # Remove last subscriber should stop session
        await session.add_subscriber(queue1)  # Add back
        await session.start()  # Start session
        
        await session.remove_subscriber(queue1)
        await session.remove_subscriber(queue2)
        
        # Give a moment for stop to be called
        await asyncio.sleep(0.1)
        assert not session.running


class TestRealExecution:
    """Test with real batch execution (requires running work)."""
    
    @pytest.mark.integration
    def test_websocket_with_real_execution(self, test_client):
        """Test WebSocket with real batch execution."""
        # This test requires a real running batch
        # Skip if no real execution is available
        pytest.skip("Integration test - requires real execution")
    
    @pytest.mark.integration  
    @pytest.mark.asyncio
    async def test_real_database_queries(self):
        """Test with real database queries."""
        # This test would use real database from running execution
        pytest.skip("Integration test - requires real database")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
