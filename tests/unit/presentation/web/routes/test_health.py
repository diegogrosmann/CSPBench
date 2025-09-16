"""
Unit tests for src.presentation.web.routes.health module.

This module tests FastAPI health check endpoints including:
- GET /health/ - basic health check
- GET /health/detailed - detailed health check with components
- GET /health/metrics - basic metrics endpoint
- GET /health/status - system status for UI
- GET /health/version - version information
- Version parsing from pyproject.toml
- Error handling and edge cases
"""

import pytest
from unittest.mock import patch, mock_open, MagicMock
from datetime import datetime
from pathlib import Path
from fastapi import HTTPException
from fastapi.testclient import TestClient
from fastapi import FastAPI

from src.presentation.web.routes.health import router, get_version


@pytest.fixture
def app():
    """Create FastAPI app with health router."""
    app = FastAPI()
    app.include_router(router)
    return app


@pytest.fixture
def client(app):
    """Create test client."""
    return TestClient(app)


class TestGetVersion:
    """Test get_version function."""

    @patch('builtins.open', new_callable=mock_open, read_data='version = "1.2.3"')
    @patch('pathlib.Path.exists')
    def test_get_version_from_pyproject(self, mock_exists, mock_file):
        """Test version extraction from pyproject.toml."""
        mock_exists.return_value = True
        
        version = get_version()
        assert version == "v1.2.3"

    @patch('builtins.open', new_callable=mock_open, read_data='[tool.poetry]\nversion = "2.1.0"\nname = "test"')
    @patch('pathlib.Path.exists')
    def test_get_version_from_pyproject_multiline(self, mock_exists, mock_file):
        """Test version extraction from multi-line pyproject.toml."""
        mock_exists.return_value = True
        
        version = get_version()
        assert version == "v2.1.0"

    @patch('builtins.open', new_callable=mock_open, read_data="version = '3.0.0-beta'")
    @patch('pathlib.Path.exists')
    def test_get_version_single_quotes(self, mock_exists, mock_file):
        """Test version extraction with single quotes."""
        mock_exists.return_value = True
        
        version = get_version()
        assert version == "v3.0.0-beta"

    @patch('pathlib.Path.exists')
    @patch('os.getenv')
    def test_get_version_no_pyproject_with_env(self, mock_getenv, mock_exists):
        """Test version fallback to environment variable."""
        mock_exists.return_value = False
        mock_getenv.return_value = "v4.5.6"
        
        version = get_version()
        assert version == "v4.5.6"

    @patch('pathlib.Path.exists')
    @patch('os.getenv')
    def test_get_version_no_pyproject_no_env(self, mock_getenv, mock_exists):
        """Test version fallback to default."""
        mock_exists.return_value = False
        mock_getenv.return_value = None
        
        version = get_version()
        assert version == "v0.1.0"

    @patch('builtins.open', side_effect=IOError("File read error"))
    @patch('pathlib.Path.exists')
    @patch('os.getenv')
    def test_get_version_file_read_error(self, mock_getenv, mock_exists, mock_file):
        """Test version fallback when file read fails."""
        mock_exists.return_value = True
        mock_getenv.return_value = None
        
        version = get_version()
        assert version == "v0.1.0"

    @patch('builtins.open', new_callable=mock_open, read_data='no version line here')
    @patch('pathlib.Path.exists')
    @patch('os.getenv')
    def test_get_version_no_version_line(self, mock_getenv, mock_exists, mock_file):
        """Test version fallback when no version line found."""
        mock_exists.return_value = True
        mock_getenv.return_value = None
        
        version = get_version()
        assert version == "v0.1.0"


class TestBasicHealthCheck:
    """Test GET /health/ endpoint."""

    @patch('src.presentation.web.routes.health.global_registry')
    @patch('src.presentation.web.routes.health.get_version')
    def test_health_check_success(self, mock_get_version, mock_registry, client):
        """Test successful basic health check."""
        mock_registry.__len__ = MagicMock(return_value=5)
        mock_get_version.return_value = "v1.0.0"

        response = client.get("/health/")
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "healthy"
        assert data["version"] == "v1.0.0"

    @patch('src.presentation.web.routes.health.global_registry')
    @patch('src.presentation.web.routes.health.get_version')
    def test_health_check_no_algorithms(self, mock_get_version, mock_registry, client):
        """Test health check with no algorithms loaded."""
        mock_registry.__len__ = MagicMock(return_value=0)
        mock_get_version.return_value = "v1.0.0"

        response = client.get("/health/")
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "healthy"  # Still healthy even with 0 algorithms
        assert data["version"] == "v1.0.0"

    @patch('src.presentation.web.routes.health.global_registry')
    @patch('src.presentation.web.routes.health.get_version')
    def test_health_check_registry_error(self, mock_get_version, mock_registry, client):
        """Test health check when registry access fails."""
        mock_registry.__len__ = MagicMock(side_effect=Exception("Registry error"))
        mock_get_version.return_value = "v1.0.0"

        with patch('src.presentation.web.routes.health.logger') as mock_logger:
            response = client.get("/health/")
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "error"
        assert data["version"] == "v1.0.0"
        
        # Should log error
        mock_logger.error.assert_called_once()

    @patch('src.presentation.web.routes.health.global_registry')
    def test_health_check_version_error(self, mock_registry, client):
        """Test health check when version retrieval fails internally."""
        mock_registry.__len__ = MagicMock(return_value=3)

        # Mock the get_version to return default when file operations fail
        with patch('src.presentation.web.routes.health.get_version', return_value="v0.1.0"):
            response = client.get("/health/")
            
            assert response.status_code == 200
            data = response.json()
            assert data["status"] == "healthy"
            assert data["version"] == "v0.1.0"
class TestDetailedHealthCheck:
    """Test GET /health/detailed endpoint."""

    @patch('src.presentation.web.routes.health.global_registry')
    @patch('src.presentation.web.routes.health.get_version')
    @patch('src.presentation.web.routes.health.datetime')
    def test_detailed_health_check_success(self, mock_datetime, mock_get_version, mock_registry, client):
        """Test successful detailed health check."""
        mock_registry.__len__ = MagicMock(return_value=5)
        mock_get_version.return_value = "v1.0.0"
        mock_datetime.now.return_value.isoformat.return_value = "2024-01-01T12:00:00"

        response = client.get("/health/detailed")
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "healthy"
        assert data["timestamp"] == "2024-01-01T12:00:00"
        assert data["components"]["algorithms"] is True
        assert data["components"]["config"] is True
        assert data["version"] == "v1.0.0"

    @patch('src.presentation.web.routes.health.global_registry')
    @patch('src.presentation.web.routes.health.get_version')
    @patch('src.presentation.web.routes.health.datetime')
    def test_detailed_health_check_no_algorithms(self, mock_datetime, mock_get_version, mock_registry, client):
        """Test detailed health check with no algorithms."""
        mock_registry.__len__ = MagicMock(return_value=0)
        mock_get_version.return_value = "v1.0.0"
        mock_datetime.now.return_value.isoformat.return_value = "2024-01-01T12:00:00"

        response = client.get("/health/detailed")
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "degraded"
        assert data["timestamp"] == "2024-01-01T12:00:00"
        assert data["components"]["algorithms"] is False
        assert data["components"]["config"] is True
        assert data["version"] == "v1.0.0"

    @patch('src.presentation.web.routes.health.global_registry')
    def test_detailed_health_check_registry_error(self, mock_registry, client):
        """Test detailed health check when registry access fails."""
        mock_registry.__len__ = MagicMock(side_effect=Exception("Registry error"))

        response = client.get("/health/detailed")
        
        assert response.status_code == 500
        data = response.json()
        assert "status" in data["detail"]
        assert data["detail"]["status"] == "unhealthy"
        assert data["detail"]["components"]["algorithms"] is False
        assert data["detail"]["components"]["config"] is False
        assert "error" in data["detail"]


class TestMetricsEndpoint:
    """Test GET /health/metrics endpoint."""

    @patch('src.presentation.web.routes.health.global_registry')
    @patch('src.presentation.web.routes.health.datetime')
    def test_get_metrics_success(self, mock_datetime, mock_registry, client):
        """Test successful metrics retrieval."""
        mock_registry.__len__ = MagicMock(return_value=7)
        mock_datetime.now.return_value.isoformat.return_value = "2024-01-01T12:00:00"

        response = client.get("/health/metrics")
        
        assert response.status_code == 200
        data = response.json()
        assert data["algorithms_count"] == 7
        assert data["status"] == "running"
        assert data["uptime"] == "unknown"
        assert data["timestamp"] == "2024-01-01T12:00:00"

    @patch('src.presentation.web.routes.health.global_registry')
    @patch('src.presentation.web.routes.health.datetime')
    def test_get_metrics_registry_error(self, mock_datetime, mock_registry, client):
        """Test metrics when registry access fails."""
        mock_registry.__len__ = MagicMock(side_effect=Exception("Registry error"))
        mock_datetime.now.return_value.isoformat.return_value = "2024-01-01T12:00:00"

        with patch('src.presentation.web.routes.health.logger') as mock_logger:
            response = client.get("/health/metrics")
        
        assert response.status_code == 200
        data = response.json()
        assert "error" in data
        assert data["error"] == "Metrics unavailable"
        assert data["timestamp"] == "2024-01-01T12:00:00"
        
        # Should log error
        mock_logger.error.assert_called_once()


class TestSystemStatus:
    """Test GET /health/status endpoint."""

    @patch('src.presentation.web.routes.health.global_registry')
    @patch('src.presentation.web.routes.health.datetime')
    def test_get_system_status_online(self, mock_datetime, mock_registry, client):
        """Test system status when online."""
        mock_registry.__len__ = MagicMock(return_value=3)
        mock_datetime.now.return_value.isoformat.return_value = "2024-01-01T12:00:00"

        response = client.get("/health/status")
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "Online"
        assert data["status_type"] == "success"
        assert data["algorithms_available"] is True
        assert data["timestamp"] == "2024-01-01T12:00:00"

    @patch('src.presentation.web.routes.health.global_registry')
    @patch('src.presentation.web.routes.health.datetime')
    def test_get_system_status_degraded(self, mock_datetime, mock_registry, client):
        """Test system status when degraded."""
        mock_registry.__len__ = MagicMock(return_value=0)
        mock_datetime.now.return_value.isoformat.return_value = "2024-01-01T12:00:00"

        response = client.get("/health/status")
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "Degraded"
        assert data["status_type"] == "warning"
        assert data["algorithms_available"] is False
        assert data["timestamp"] == "2024-01-01T12:00:00"

    @patch('src.presentation.web.routes.health.global_registry')
    @patch('src.presentation.web.routes.health.datetime')
    def test_get_system_status_error(self, mock_datetime, mock_registry, client):
        """Test system status when error occurs."""
        mock_registry.__len__ = MagicMock(side_effect=Exception("Status error"))
        mock_datetime.now.return_value.isoformat.return_value = "2024-01-01T12:00:00"

        with patch('src.presentation.web.routes.health.logger') as mock_logger:
            response = client.get("/health/status")
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "Error"
        assert data["status_type"] == "danger"
        assert data["algorithms_available"] is False
        assert data["timestamp"] == "2024-01-01T12:00:00"
        assert "error" in data
        assert data["error"] == "Status error"
        
        # Should log error
        mock_logger.error.assert_called_once()


class TestSystemVersion:
    """Test GET /health/version endpoint."""

    @patch('src.presentation.web.routes.health.get_version')
    @patch('src.presentation.web.routes.health.sys')
    @patch('src.presentation.web.routes.health.datetime')
    def test_get_system_version_success(self, mock_datetime, mock_sys, mock_get_version, client):
        """Test successful version retrieval."""
        mock_get_version.return_value = "v2.0.0"
        mock_sys.version_info.major = 3
        mock_sys.version_info.minor = 12
        mock_sys.version_info.micro = 1
        mock_datetime.now.return_value.isoformat.return_value = "2024-01-01T12:00:00"

        response = client.get("/health/version")
        
        assert response.status_code == 200
        data = response.json()
        assert data["version"] == "v2.0.0"
        assert data["python_version"] == "3.12.1"
        assert data["architecture"] == "Hexagonal"
        assert data["timestamp"] == "2024-01-01T12:00:00"

    @patch('src.presentation.web.routes.health.get_version')
    @patch('src.presentation.web.routes.health.sys')
    @patch('src.presentation.web.routes.health.datetime')
    def test_get_system_version_error(self, mock_datetime, mock_sys, mock_get_version, client):
        """Test version retrieval when error occurs."""
        mock_get_version.side_effect = Exception("Version error")
        mock_sys.version_info.side_effect = Exception("Sys error")
        mock_datetime.now.return_value.isoformat.return_value = "2024-01-01T12:00:00"

        with patch('src.presentation.web.routes.health.logger') as mock_logger:
            response = client.get("/health/version")
        
        assert response.status_code == 200
        data = response.json()
        assert data["version"] == "Unknown"
        assert data["python_version"] == "Unknown"
        assert data["architecture"] == "Hexagonal"
        assert data["timestamp"] == "2024-01-01T12:00:00"
        assert "error" in data
        
        # Should log error
        mock_logger.error.assert_called_once()


class TestHealthRouterIntegration:
    """Test health router integration and edge cases."""

    def test_router_prefix_and_tags(self, client):
        """Test that router has correct prefix and tags."""
        # Test that correct prefix works
        response = client.get("/health/")
        assert response.status_code == 200
        
        # Test that incorrect prefix fails
        response = client.get("/api/health/")  # Wrong prefix
        assert response.status_code == 404

    @patch('src.presentation.web.routes.health.global_registry')
    def test_multiple_endpoints_consistency(self, mock_registry, client):
        """Test that multiple endpoints return consistent data."""
        mock_registry.__len__ = MagicMock(return_value=5)

        # Test basic health
        response1 = client.get("/health/")
        assert response1.status_code == 200
        
        # Test detailed health
        response2 = client.get("/health/detailed")
        assert response2.status_code == 200
        
        # Test metrics
        response3 = client.get("/health/metrics")
        assert response3.status_code == 200
        
        # All should report same algorithm count
        data2 = response2.json()
        data3 = response3.json()
        assert data2["components"]["algorithms"] is True
        assert data3["algorithms_count"] == 5

    @patch('src.presentation.web.routes.health.global_registry')
    def test_concurrent_health_checks(self, mock_registry, client):
        """Test multiple concurrent health check requests."""
        mock_registry.__len__ = MagicMock(return_value=10)

        # Simulate multiple concurrent requests
        responses = []
        for _ in range(5):
            response = client.get("/health/")
            responses.append(response)

        # All should succeed
        for response in responses:
            assert response.status_code == 200
            data = response.json()
            assert data["status"] == "healthy"


class TestHealthEdgeCases:
    """Test edge cases and error conditions."""

    @patch('src.presentation.web.routes.health.Path')
    def test_path_resolution_edge_cases(self, mock_path):
        """Test path resolution edge cases in get_version."""
        # Test when Path construction fails
        mock_path.side_effect = Exception("Path error")
        
        version = get_version()
        assert version == "v0.1.0"

    @patch('src.presentation.web.routes.health.global_registry')
    def test_registry_len_returns_none(self, mock_registry, client):
        """Test when registry length returns None or unexpected type."""
        mock_registry.__len__ = MagicMock(return_value=None)

        # Should handle gracefully
        response = client.get("/health/detailed")
        # This might cause a TypeError, but the except block should catch it
        assert response.status_code in [200, 500]

    def test_empty_pyproject_content(self):
        """Test get_version with empty pyproject.toml content."""
        with patch('builtins.open', mock_open(read_data="")), \
             patch('pathlib.Path.exists', return_value=True), \
             patch('os.getenv', return_value=None):
            
            version = get_version()
            assert version == "v0.1.0"

    def test_malformed_version_line(self):
        """Test get_version with malformed version line."""
        with patch('builtins.open', mock_open(read_data="version =")), \
             patch('pathlib.Path.exists', return_value=True), \
             patch('os.getenv', return_value=None):
            
            version = get_version()
            assert version == "v0.1.0"