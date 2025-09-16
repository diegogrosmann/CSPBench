"""
Unit tests for src.presentation.web.routes.algorithms module.

This module tests FastAPI algorithm endpoints including:
- GET /api/algorithms/ - list all algorithms
- GET /api/algorithms/{algorithm_name} - get specific algorithm info
- Error handling and edge cases
- Algorithm metadata extraction
- HTTPException responses
"""

import pytest
from unittest.mock import patch, MagicMock
from fastapi import HTTPException
from fastapi.testclient import TestClient
from fastapi import FastAPI

from src.presentation.web.routes.algorithms import router


@pytest.fixture
def app():
    """Create FastAPI app with algorithms router."""
    app = FastAPI()
    app.include_router(router)
    return app


@pytest.fixture
def client(app):
    """Create test client."""
    return TestClient(app)


class TestGetAlgorithms:
    """Test GET /api/algorithms/ endpoint."""

    @patch('src.presentation.web.routes.algorithms.global_registry')
    def test_get_algorithms_success(self, mock_registry, client):
        """Test successful retrieval of algorithms list."""
        # Mock algorithm classes
        mock_algorithm1 = MagicMock()
        mock_algorithm1.__doc__ = "Algorithm 1 description"
        mock_algorithm1.default_params = {"param1": "value1"}
        mock_algorithm1.is_deterministic = True
        mock_algorithm1.supports_internal_parallel = False
        mock_algorithm1.category = "Optimization"

        mock_algorithm2 = MagicMock()
        mock_algorithm2.__doc__ = "Algorithm 2 description"
        mock_algorithm2.default_params = {"param2": "value2"}
        mock_algorithm2.is_deterministic = False
        mock_algorithm2.supports_internal_parallel = True
        mock_algorithm2.category = "Heuristic"

        mock_registry.items.return_value = [
            ("algorithm1", mock_algorithm1),
            ("algorithm2", mock_algorithm2),
        ]

        response = client.get("/api/algorithms/")
        
        assert response.status_code == 200
        data = response.json()
        assert len(data) == 2
        
        # Check first algorithm
        assert data[0]["name"] == "algorithm1"
        assert data[0]["description"] == "Algorithm 1 description"
        assert data[0]["default_params"] == {"param1": "value1"}
        assert data[0]["is_deterministic"] is True
        assert data[0]["supports_internal_parallel"] is False
        assert data[0]["category"] == "Optimization"
        
        # Check second algorithm
        assert data[1]["name"] == "algorithm2"
        assert data[1]["description"] == "Algorithm 2 description"
        assert data[1]["default_params"] == {"param2": "value2"}
        assert data[1]["is_deterministic"] is False
        assert data[1]["supports_internal_parallel"] is True
        assert data[1]["category"] == "Heuristic"

    @patch('src.presentation.web.routes.algorithms.global_registry')
    def test_get_algorithms_with_missing_attributes(self, mock_registry, client):
        """Test algorithm list with missing class attributes (use defaults)."""
        mock_algorithm = MagicMock()
        mock_algorithm.__doc__ = "Test algorithm"
        # Missing default_params, is_deterministic, supports_internal_parallel, category
        del mock_algorithm.default_params
        del mock_algorithm.is_deterministic
        del mock_algorithm.supports_internal_parallel
        del mock_algorithm.category

        mock_registry.items.return_value = [("test_algo", mock_algorithm)]

        response = client.get("/api/algorithms/")
        
        assert response.status_code == 200
        data = response.json()
        assert len(data) == 1
        
        algorithm = data[0]
        assert algorithm["name"] == "test_algo"
        assert algorithm["description"] == "Test algorithm"
        assert algorithm["default_params"] == {}  # Default
        assert algorithm["is_deterministic"] is True  # Default
        assert algorithm["supports_internal_parallel"] is False  # Default
        assert algorithm["category"] == "General"  # Default

    @patch('src.presentation.web.routes.algorithms.global_registry')
    def test_get_algorithms_with_no_doc(self, mock_registry, client):
        """Test algorithm list when __doc__ is None."""
        mock_algorithm = MagicMock()
        mock_algorithm.__doc__ = None
        mock_algorithm.default_params = {}
        mock_algorithm.is_deterministic = True
        mock_algorithm.supports_internal_parallel = False
        mock_algorithm.category = "Test"

        mock_registry.items.return_value = [("no_doc_algo", mock_algorithm)]

        response = client.get("/api/algorithms/")
        
        assert response.status_code == 200
        data = response.json()
        assert len(data) == 1
        
        algorithm = data[0]
        assert algorithm["name"] == "no_doc_algo"
        assert algorithm["description"] == "no_doc_algo algorithm"  # Fallback

    @patch('src.presentation.web.routes.algorithms.global_registry')
    def test_get_algorithms_with_metadata_extraction_error(self, mock_registry, client):
        """Test algorithm list when metadata extraction fails."""
        mock_algorithm = MagicMock()
        # Make accessing attributes raise an exception
        mock_algorithm.__doc__ = "Test algorithm"
        mock_algorithm.default_params.side_effect = Exception("Attribute error")
        
        mock_registry.items.return_value = [("error_algo", mock_algorithm)]

        with patch('src.presentation.web.routes.algorithms.logger') as mock_logger:
            response = client.get("/api/algorithms/")
        
        assert response.status_code == 200
        data = response.json()
        assert len(data) == 1
        
        # Should use fallback values
        algorithm = data[0]
        assert algorithm["name"] == "error_algo"
        assert algorithm["description"] == "Test algorithm"
        assert algorithm["default_params"] == {}
        assert algorithm["is_deterministic"] is True
        assert algorithm["supports_internal_parallel"] is False
        assert algorithm["category"] == "General"
        
        # Should log warning
        mock_logger.warning.assert_called_once()

    @patch('src.presentation.web.routes.algorithms.global_registry')
    def test_get_algorithms_with_doc_access_error(self, mock_registry, client):
        """Test algorithm list when even __doc__ access fails."""
        mock_algorithm = MagicMock()
        # Make accessing __doc__ raise an exception
        type(mock_algorithm).__doc__ = PropertyMock(side_effect=Exception("Doc error"))
        mock_algorithm.default_params.side_effect = Exception("Attribute error")
        
        mock_registry.items.return_value = [("doc_error_algo", mock_algorithm)]

        with patch('src.presentation.web.routes.algorithms.logger') as mock_logger:
            response = client.get("/api/algorithms/")
        
        assert response.status_code == 200
        data = response.json()
        assert len(data) == 1
        
        # Should use algorithm name as fallback description
        algorithm = data[0]
        assert algorithm["name"] == "doc_error_algo"
        assert algorithm["description"] == "doc_error_algo algorithm"
        
        # Should log warning
        mock_logger.warning.assert_called_once()

    @patch('src.presentation.web.routes.algorithms.global_registry')
    def test_get_algorithms_empty_registry(self, mock_registry, client):
        """Test algorithm list with empty registry."""
        mock_registry.items.return_value = []

        response = client.get("/api/algorithms/")
        
        assert response.status_code == 200
        data = response.json()
        assert data == []

    @patch('src.presentation.web.routes.algorithms.global_registry')
    def test_get_algorithms_registry_error(self, mock_registry, client):
        """Test algorithm list when registry access fails."""
        mock_registry.items.side_effect = Exception("Registry error")

        with patch('src.presentation.web.routes.algorithms.logger') as mock_logger:
            response = client.get("/api/algorithms/")
        
        assert response.status_code == 500
        assert "Failed to retrieve algorithms" in response.json()["detail"]
        
        # Should log error
        mock_logger.error.assert_called_once()


class TestGetAlgorithmInfo:
    """Test GET /api/algorithms/{algorithm_name} endpoint."""

    @patch('src.presentation.web.routes.algorithms.global_registry')
    def test_get_algorithm_info_success(self, mock_registry, client):
        """Test successful retrieval of specific algorithm info."""
        mock_algorithm = MagicMock()
        mock_algorithm.__doc__ = "Detailed algorithm description"
        mock_algorithm.default_params = {"iterations": 100, "tolerance": 0.01}
        mock_algorithm.is_deterministic = True
        mock_algorithm.supports_internal_parallel = True
        mock_algorithm.category = "Optimization"
        mock_algorithm.__version__ = "2.1.0"

        mock_registry.get.return_value = mock_algorithm

        response = client.get("/api/algorithms/test_algorithm")
        
        assert response.status_code == 200
        data = response.json()
        
        assert data["name"] == "test_algorithm"
        assert data["description"] == "Detailed algorithm description"
        assert data["default_params"] == {"iterations": 100, "tolerance": 0.01}
        assert data["is_deterministic"] is True
        assert data["supports_internal_parallel"] is True
        assert data["category"] == "Optimization"
        assert data["version"] == "2.1.0"

    @patch('src.presentation.web.routes.algorithms.global_registry')
    def test_get_algorithm_info_with_defaults(self, mock_registry, client):
        """Test algorithm info with default attribute values."""
        mock_algorithm = MagicMock()
        mock_algorithm.__doc__ = "Test algorithm"
        # Missing optional attributes - should use defaults
        del mock_algorithm.default_params
        del mock_algorithm.is_deterministic
        del mock_algorithm.supports_internal_parallel
        del mock_algorithm.category
        del mock_algorithm.__version__

        mock_registry.get.return_value = mock_algorithm

        response = client.get("/api/algorithms/default_algo")
        
        assert response.status_code == 200
        data = response.json()
        
        assert data["name"] == "default_algo"
        assert data["description"] == "Test algorithm"
        assert data["default_params"] == {}
        assert data["is_deterministic"] is True
        assert data["supports_internal_parallel"] is False
        assert data["category"] == "General"
        assert data["version"] == "1.0.0"

    @patch('src.presentation.web.routes.algorithms.global_registry')
    def test_get_algorithm_info_no_doc(self, mock_registry, client):
        """Test algorithm info when __doc__ is None."""
        mock_algorithm = MagicMock()
        mock_algorithm.__doc__ = None
        mock_algorithm.default_params = {}
        mock_algorithm.is_deterministic = False
        mock_algorithm.supports_internal_parallel = False
        mock_algorithm.category = "Test"
        mock_algorithm.__version__ = "1.5.0"

        mock_registry.get.return_value = mock_algorithm

        response = client.get("/api/algorithms/no_doc_algo")
        
        assert response.status_code == 200
        data = response.json()
        
        assert data["name"] == "no_doc_algo"
        assert data["description"] == "no_doc_algo algorithm"  # Fallback

    @patch('src.presentation.web.routes.algorithms.global_registry')
    def test_get_algorithm_info_not_found(self, mock_registry, client):
        """Test algorithm info for non-existent algorithm."""
        mock_registry.get.return_value = None

        response = client.get("/api/algorithms/nonexistent")
        
        assert response.status_code == 404
        assert "Algorithm 'nonexistent' not found" in response.json()["detail"]

    @patch('src.presentation.web.routes.algorithms.global_registry')
    def test_get_algorithm_info_registry_error(self, mock_registry, client):
        """Test algorithm info when registry access fails."""
        mock_registry.get.side_effect = Exception("Registry error")

        with patch('src.presentation.web.routes.algorithms.logger') as mock_logger:
            response = client.get("/api/algorithms/error_algo")
        
        assert response.status_code == 500
        assert "Failed to retrieve algorithm information" in response.json()["detail"]
        
        # Should log error
        mock_logger.error.assert_called_once()

    @patch('src.presentation.web.routes.algorithms.global_registry')
    def test_get_algorithm_info_general_exception(self, mock_registry, client):
        """Test algorithm info when a general exception occurs."""
        # Make the registry access itself fail, not just attribute access
        mock_registry.get.side_effect = RuntimeError("General error")

        with patch('src.presentation.web.routes.algorithms.logger') as mock_logger:
            response = client.get("/api/algorithms/error_algo")
        
        assert response.status_code == 500
        assert "Failed to retrieve algorithm information" in response.json()["detail"]
        
        # Should log error
        mock_logger.error.assert_called_once()


class TestAlgorithmsRouterIntegration:
    """Test algorithms router integration and edge cases."""

    @patch('src.presentation.web.routes.algorithms.global_registry')
    def test_router_prefix_and_tags(self, mock_registry, client):
        """Test that router has correct prefix and tags."""
        mock_registry.items.return_value = []
        
        response = client.get("/api/algorithms/")
        assert response.status_code == 200
        
        # Test that incorrect prefix fails
        response = client.get("/algorithms/")  # Missing /api
        assert response.status_code == 404

    @patch('src.presentation.web.routes.algorithms.global_registry')
    def test_algorithm_name_with_special_characters(self, mock_registry, client):
        """Test algorithm retrieval with special characters in name."""
        mock_algorithm = MagicMock()
        mock_algorithm.__doc__ = "Special algorithm"
        mock_algorithm.default_params = {}
        mock_algorithm.is_deterministic = True
        mock_algorithm.supports_internal_parallel = False
        mock_algorithm.category = "Special"

        mock_registry.get.return_value = mock_algorithm

        # Test with underscores and numbers
        response = client.get("/api/algorithms/test_algorithm_v2")
        assert response.status_code == 200
        
        data = response.json()
        assert data["name"] == "test_algorithm_v2"

    @patch('src.presentation.web.routes.algorithms.global_registry')
    def test_large_algorithm_list_performance(self, mock_registry, client):
        """Test performance with large number of algorithms."""
        # Create many mock algorithms
        mock_algorithms = []
        for i in range(100):
            mock_algo = MagicMock()
            mock_algo.__doc__ = f"Algorithm {i} description"
            mock_algo.default_params = {"param": i}
            mock_algo.is_deterministic = i % 2 == 0
            mock_algo.supports_internal_parallel = i % 3 == 0
            mock_algo.category = f"Category_{i % 5}"
            mock_algorithms.append((f"algorithm_{i}", mock_algo))

        mock_registry.items.return_value = mock_algorithms

        response = client.get("/api/algorithms/")
        
        assert response.status_code == 200
        data = response.json()
        assert len(data) == 100
        
        # Verify first and last algorithms
        assert data[0]["name"] == "algorithm_0"
        assert data[99]["name"] == "algorithm_99"


# Add PropertyMock import at the top of the file for the mock that uses it
from unittest.mock import PropertyMock