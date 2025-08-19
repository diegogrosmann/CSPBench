"""
Basic tests for Batch Execution API endpoints.

Tests the core functionality of batch execution including:
- Batch submission and validation
- Status monitoring
- Basic error handling
"""

import pytest
import tempfile
import yaml
from datetime import datetime
from pathlib import Path
from unittest.mock import Mock, patch

from fastapi.testclient import TestClient
from src.presentation.web.app import app

client = TestClient(app)

# Sample batch configuration for testing
SAMPLE_BATCH_CONFIG = {
    "metadata": {
        "name": "Test Batch",
        "description": "Test batch configuration",
        "author": "Test Author",
        "version": "1.0",
        "creation_date": "2025-01-01",
        "tags": ["test"]
    },
    "datasets": [
        {
            "id": "test_dataset",
            "name": "Test Dataset",
            "type": "synthetic",
            "parameters": {
                "mode": "random",
                "n": 5,
                "L": 10,
                "alphabet": "ACGT",
                "seed": 42
            }
        }
    ],
    "algorithms": [
        {
            "id": "test_config",
            "name": "Test Configuration",
            "description": "Test algorithm configuration",
            "algorithms": ["Baseline"],
            "algorithm_params": {
                "Baseline": {
                    "tie_break": "lex"
                }
            }
        }
    ],
    "task": {
        "type": "experiment",
        "parameters": {
            "title": "Test Experiment",
            "description": "Test experiment description",
            "algorithm_preset": "test_config",
            "experiment_instances": [
                {
                    "id": "test_instance",
                    "title": "Test Instance",
                    "dataset": "test_dataset",
                    "repetitions": 1,
                    "seeds": [42]
                }
            ]
        }
    }
}


def test_batch_execution_invalid_file():
    """Test batch execution with invalid file path."""
    response = client.post("/api/batch/execute", json={
        "batch_file": "nonexistent.yaml",
        "monitor_type": "log"
    })
    
    assert response.status_code == 404
    assert "not found" in response.json()["detail"]


def test_batch_execution_empty_file():
    """Test batch execution with empty file path."""
    response = client.post("/api/batch/execute", json={
        "batch_file": "",
        "monitor_type": "log"
    })
    
    assert response.status_code == 422  # Pydantic validation error
    assert "cannot be empty" in response.json()["detail"][0]["msg"]


def test_batch_execution_invalid_monitor_type():
    """Test batch execution with invalid monitor type."""
    
    response = client.post("/api/batch/execute", json={
        "batch_file": "TEMPLATE.yaml",
        "monitor_type": "invalid_type"
    })
    
    assert response.status_code == 422  # Validation error


def test_batch_execution_invalid_extension():
    """Test batch execution with invalid file extension."""
    response = client.post("/api/batch/execute", json={
        "batch_file": "somefile.txt",
        "monitor_type": "log"
    })
    
    assert response.status_code == 422  # Pydantic validation error
    assert "must have .yaml or .yml extension" in response.json()["detail"][0]["msg"]


@patch('src.presentation.web.routes.batch_execution.PipelineService.run')
@patch('src.presentation.web.routes.batch_execution.get_global_work_manager')
@patch('src.presentation.web.routes.batch_execution.submit_work')
@patch('src.presentation.web.routes.batch_execution.load_cspbench_config')
@patch('pathlib.Path.exists')
@patch('pathlib.Path.is_file')
def test_batch_execution_success(mock_is_file, mock_exists, mock_load_config, mock_submit_work, mock_get_work_manager, mock_pipeline_run):
    """Test successful batch execution submission."""
    # Mock file existence
    mock_exists.return_value = True
    mock_is_file.return_value = True
    
    # Mock configuration
    mock_config = Mock()
    mock_config.metadata.name = "Test Batch"
    mock_load_config.return_value = mock_config
    
    # Mock submit_work
    mock_submit_work.return_value = ("test_work_id", {
        "id": "test_work_id",
        "status": "queued",
        "created_at": datetime.now().timestamp(),
        "updated_at": datetime.now().timestamp(),
        "output_path": "/tmp/test_output",
        "submitted_at": datetime.now()
    })
    
    # Mock work manager
    mock_wm = Mock()
    mock_work_item = {
        "id": "test_work_id",
        "status": "queued",
        "created_at": datetime.now().timestamp(),
        "updated_at": datetime.now().timestamp(),
        "output_path": "/tmp/test_output"
    }
    mock_wm.get.return_value = mock_work_item
    mock_get_work_manager.return_value = mock_wm
    
    # Mock pipeline service
    mock_pipeline_run.return_value = "test_work_id"
    
    response = client.post("/api/batch/execute", json={
        "batch_file": "test_batch.yaml",
        "monitor_type": "log"
    })
    
    assert response.status_code == 200
    data = response.json()
    assert data["work_id"] == "test_work_id"
    assert data["status"] == "queued"
    assert "started successfully" in data["message"]
    assert data["output_path"] == "/tmp/test_output"


@patch('src.presentation.web.routes.batch_execution.get_work_details')
def test_get_batch_status_not_found(mock_get_work_details):
    """Test getting status for non-existent batch."""
    mock_get_work_details.return_value = None
    
    response = client.get("/api/batch/nonexistent_id/status")
    
    assert response.status_code == 404
    assert "not found" in response.json()["detail"]


@patch('src.presentation.web.routes.batch_execution.get_work_details')
def test_get_batch_status_success(mock_get_work_details):
    """Test getting batch status successfully."""
    mock_config = Mock()
    mock_config.metadata.name = "Test Batch"
    
    mock_work_item = {
        "id": "test_work_id",
        "status": "running",
        "created_at": datetime.now().timestamp(),
        "updated_at": datetime.now().timestamp(),
        "output_path": "/tmp/test_output",
        "error": None,
        "config": mock_config
    }
    mock_get_work_details.return_value = mock_work_item
    
    response = client.get("/api/batch/test_work_id/status")
    
    assert response.status_code == 200
    data = response.json()
    assert data["work_id"] == "test_work_id"
    assert data["status"] == "running"


@patch('src.presentation.web.routes.batch_execution.get_work_details')
def test_get_batch_results_not_complete(mock_get_work_details):
    """Test getting results for incomplete batch."""
    mock_work_item = {
        "id": "test_work_id",
        "status": "running",  # Not terminal status
        "created_at": datetime.now().timestamp(),
        "updated_at": datetime.now().timestamp()
    }
    mock_get_work_details.return_value = mock_work_item
    
    response = client.get("/api/batch/test_work_id/results")
    
    assert response.status_code == 409
    assert "not complete" in response.json()["detail"]


@patch('src.presentation.web.routes.batch_execution.get_work_details')
@patch('src.presentation.web.routes.batch_execution.control_work')
def test_cancel_batch_success(mock_control_work, mock_get_work_details):
    """Test successfully canceling a batch."""
    mock_work_item = {
        "id": "test_work_id",
        "status": "running"
    }
    mock_get_work_details.return_value = mock_work_item
    mock_control_work.return_value = True
    
    response = client.post("/api/batch/test_work_id/control", json={"action": "cancel"})
    
    assert response.status_code == 200
    data = response.json()
    assert data["work_id"] == "test_work_id"
    assert data["operation"] == "cancel"
    assert data["success"] is True


@patch('src.presentation.web.routes.batch_execution.list_all_work')
def test_list_batch_executions(mock_list_all_work):
    """Test listing batch executions."""
    mock_items = [
        {
            "id": "work_1",
            "status": "finished",
            "created_at": datetime.now().timestamp(),
            "updated_at": datetime.now().timestamp(),
            "output_path": "/tmp/work_1",
            "error": None,
            "config": None
        },
        {
            "id": "work_2", 
            "status": "running",
            "created_at": datetime.now().timestamp(),
            "updated_at": datetime.now().timestamp(),
            "output_path": "/tmp/work_2",
            "error": None,
            "config": None
        }
    ]
    mock_list_all_work.return_value = mock_items
    
    response = client.get("/api/batch/list")
    
    assert response.status_code == 200
    data = response.json()
    assert data["total"] == 2
    assert len(data["executions"]) == 2
    assert data["executions"][0]["work_id"] == "work_1"
    assert data["executions"][1]["work_id"] == "work_2"


@patch('src.presentation.web.routes.batch_execution.list_all_work')
def test_list_batch_executions_filtered(mock_list_all_work):
    """Test listing batch executions with status filter."""
    mock_items = [
        {
            "id": "work_1",
            "status": "finished",
            "created_at": datetime.now().timestamp(),
            "updated_at": datetime.now().timestamp(),
            "output_path": "/tmp/work_1",
            "error": None,
            "config": None
        },
        {
            "id": "work_2", 
            "status": "running",
            "created_at": datetime.now().timestamp(),
            "updated_at": datetime.now().timestamp(),
            "output_path": "/tmp/work_2",
            "error": None,
            "config": None
        }
    ]
    mock_list_all_work.return_value = mock_items
    
    response = client.get("/api/batch/list?status=running")
    
    assert response.status_code == 200
    data = response.json()
    assert data["total"] == 2  # Total without filter
    assert data["filtered"] == 1  # Total with filter
    assert len(data["executions"]) == 1
    assert data["executions"][0]["status"] == "running"
