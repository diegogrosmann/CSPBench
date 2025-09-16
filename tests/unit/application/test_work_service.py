"""
Tests for src.application.services.work_service module.

Tests the work service global singleton, lifecycle management,
and integration with WorkManager.
"""

from unittest.mock import MagicMock, Mock, patch

import pytest
from fastapi import FastAPI

from src.application.services import work_service


class TestWorkServiceSingleton:
    """Test work service singleton behavior."""

    def setup_method(self):
        """Reset global state before each test."""
        work_service._global_work_service = None

    def teardown_method(self):
        """Clean up global state after each test."""
        work_service._global_work_service = None

    @patch("src.application.services.work_service.WorkManager")
    def test_get_work_service_auto_initializes(self, mock_work_manager):
        """Test that get_work_service auto-initializes if not already done."""
        mock_manager_instance = Mock()
        mock_manager_instance.initialize.return_value = "initialized"
        mock_work_manager.return_value = mock_manager_instance

        service = work_service.get_work_service()

        assert service is not None
        mock_work_manager.assert_called_once()
        mock_manager_instance.initialize.assert_called_once()

    @patch("src.application.services.work_service.WorkManager")
    def test_get_work_service_returns_existing_instance(self, mock_work_manager):
        """Test that get_work_service returns existing instance if already initialized."""
        mock_manager_instance = Mock()
        mock_manager_instance.initialize.return_value = "initialized"
        mock_work_manager.return_value = mock_manager_instance

        # First call initializes
        service1 = work_service.get_work_service()
        # Second call returns same instance
        service2 = work_service.get_work_service()

        assert service1 is service2
        mock_work_manager.assert_called_once()  # Only called once

    @patch("src.application.services.work_service.WorkManager")
    def test_initialize_work_service_success(self, mock_work_manager):
        """Test successful work service initialization."""
        mock_manager_instance = Mock()
        mock_manager_instance.initialize.return_value = {"running": 0, "paused": 2}
        mock_work_manager.return_value = mock_manager_instance

        service = work_service.initialize_work_service()

        assert service is not None
        assert work_service._global_work_service is service
        mock_manager_instance.initialize.assert_called_once()

    @patch("src.application.services.work_service.WorkManager")
    def test_initialize_work_service_already_initialized(self, mock_work_manager):
        """Test initialization when service already exists."""
        mock_manager_instance = Mock()
        work_service._global_work_service = mock_manager_instance

        service = work_service.initialize_work_service()

        assert service is mock_manager_instance
        mock_work_manager.assert_not_called()

    @patch("src.application.services.work_service.WorkManager")
    def test_initialize_work_service_failure(self, mock_work_manager):
        """Test work service initialization failure."""
        mock_work_manager.side_effect = Exception("Initialization failed")

        with pytest.raises(RuntimeError, match="WorkService initialization failed"):
            work_service.initialize_work_service()

    def test_cleanup_work_service_with_existing_service(self):
        """Test cleanup when service exists."""
        mock_service = Mock()
        work_service._global_work_service = mock_service

        work_service.cleanup_work_service()

        mock_service.close.assert_called_once()
        assert work_service._global_work_service is None

    def test_cleanup_work_service_without_existing_service(self):
        """Test cleanup when no service exists."""
        work_service._global_work_service = None

        # Should not raise exception
        work_service.cleanup_work_service()

        assert work_service._global_work_service is None

    def test_cleanup_work_service_with_error(self):
        """Test cleanup when service close raises exception."""
        mock_service = Mock()
        mock_service.close.side_effect = Exception("Close failed")
        work_service._global_work_service = mock_service

        # Should not raise exception, just log error
        work_service.cleanup_work_service()

        # The service remains set because the exception prevented cleanup completion
        assert work_service._global_work_service is mock_service


class TestWorkServiceLifespan:
    """Test FastAPI lifespan integration."""

    def setup_method(self):
        """Reset global state before each test."""
        work_service._global_work_service = None

    def teardown_method(self):
        """Clean up global state after each test."""
        work_service._global_work_service = None

    @patch("src.application.services.work_service.cleanup_work_service")
    @patch("src.application.services.work_service.initialize_work_service")
    @pytest.mark.asyncio
    async def test_work_service_lifespan_success(self, mock_init, mock_cleanup):
        """Test successful lifespan management."""
        app = FastAPI()

        async with work_service.work_service_lifespan(app):
            pass  # Simulate app running

        mock_init.assert_called_once()
        mock_cleanup.assert_called_once()

    @patch("src.application.services.work_service.cleanup_work_service")
    @patch("src.application.services.work_service.initialize_work_service")
    @pytest.mark.asyncio
    async def test_work_service_lifespan_with_initialization_error(
        self, mock_init, mock_cleanup
    ):
        """Test lifespan when initialization fails."""
        mock_init.side_effect = Exception("Init failed")
        app = FastAPI()

        with pytest.raises(Exception, match="Init failed"):
            async with work_service.work_service_lifespan(app):
                pass

        mock_init.assert_called_once()
        mock_cleanup.assert_called_once()  # Cleanup should still be called

    @patch("src.application.services.work_service.cleanup_work_service")
    @patch("src.application.services.work_service.initialize_work_service")
    @pytest.mark.asyncio
    async def test_work_service_lifespan_with_cleanup_error(
        self, mock_init, mock_cleanup
    ):
        """Test lifespan when cleanup fails."""
        mock_cleanup.side_effect = Exception("Cleanup failed")
        app = FastAPI()

        # Should raise exception when cleanup fails
        with pytest.raises(Exception, match="Cleanup failed"):
            async with work_service.work_service_lifespan(app):
                pass

        mock_init.assert_called_once()
        mock_cleanup.assert_called_once()


class TestWorkServiceHelpers:
    """Test helper functions that delegate to WorkManager."""

    def setup_method(self):
        """Reset global state before each test."""
        work_service._global_work_service = None

    def teardown_method(self):
        """Clean up global state after each test."""
        work_service._global_work_service = None

    @patch("src.application.services.work_service.get_work_service")
    def test_get_work_status(self, mock_get_service):
        """Test get_work_status helper function."""
        mock_service = Mock()
        mock_service.get_status.return_value = "running"
        mock_get_service.return_value = mock_service

        status = work_service.get_work_status("work123")

        assert status == "running"
        mock_service.get_status.assert_called_once_with("work123")

    @patch("src.application.services.work_service.get_work_service")
    def test_get_work_details(self, mock_get_service):
        """Test get_work_details helper function."""
        mock_service = Mock()
        mock_work_item = Mock()
        mock_service.get.return_value = mock_work_item
        mock_get_service.return_value = mock_service

        details = work_service.get_work_details("work123")

        assert details is mock_work_item
        mock_service.get.assert_called_once_with("work123")

    @patch("src.application.services.work_service.get_work_service")
    def test_list_all_work(self, mock_get_service):
        """Test list_all_work helper function."""
        mock_service = Mock()
        mock_work_list = [Mock(), Mock()]
        mock_service.list.return_value = mock_work_list
        mock_get_service.return_value = mock_service

        work_list = work_service.list_all_work()

        assert work_list is mock_work_list
        mock_service.list.assert_called_once()

    @patch("src.application.services.work_service.get_work_service")
    def test_control_work_pause(self, mock_get_service):
        """Test control_work with pause action."""
        mock_service = Mock()
        mock_service.pause.return_value = True
        mock_get_service.return_value = mock_service

        result = work_service.control_work("work123", "pause")

        assert result is True
        mock_service.pause.assert_called_once_with("work123")

    @patch("src.application.services.work_service.get_work_service")
    def test_control_work_cancel(self, mock_get_service):
        """Test control_work with cancel action."""
        mock_service = Mock()
        mock_service.cancel.return_value = True
        mock_get_service.return_value = mock_service

        result = work_service.control_work("work123", "cancel")

        assert result is True
        mock_service.cancel.assert_called_once_with("work123")

    @patch("src.application.services.work_service.get_work_service")
    def test_control_work_restart(self, mock_get_service):
        """Test control_work with restart action."""
        mock_service = Mock()
        mock_service.restart.return_value = False
        mock_get_service.return_value = mock_service

        result = work_service.control_work("work123", "restart")

        assert result is False
        mock_service.restart.assert_called_once_with("work123")

    @patch("src.application.services.work_service.get_work_service")
    def test_control_work_invalid_action(self, mock_get_service):
        """Test control_work with invalid action."""
        mock_service = Mock()
        mock_get_service.return_value = mock_service

        result = work_service.control_work("work123", "invalid_action")

        assert result is False
        # No methods should be called on the service
        mock_service.pause.assert_not_called()
        mock_service.cancel.assert_not_called()
        mock_service.restart.assert_not_called()


class TestWorkServiceEdgeCases:
    """Test edge cases and error scenarios."""

    def setup_method(self):
        """Reset global state before each test."""
        work_service._global_work_service = None

    def teardown_method(self):
        """Clean up global state after each test."""
        work_service._global_work_service = None

    @patch("src.application.services.work_service.WorkManager")
    def test_initialization_with_manager_failure(self, mock_work_manager):
        """Test initialization when WorkManager constructor fails."""
        mock_work_manager.side_effect = RuntimeError("Database connection failed")

        with pytest.raises(RuntimeError, match="WorkService initialization failed"):
            work_service.initialize_work_service()

    @patch("src.application.services.work_service.WorkManager")
    def test_initialization_with_initialize_failure(self, mock_work_manager):
        """Test initialization when WorkManager.initialize() fails."""
        mock_manager_instance = Mock()
        mock_manager_instance.initialize.side_effect = Exception("Initialize failed")
        mock_work_manager.return_value = mock_manager_instance

        with pytest.raises(RuntimeError, match="WorkService initialization failed"):
            work_service.initialize_work_service()

    @patch("src.application.services.work_service.get_work_service")
    def test_helper_functions_with_service_error(self, mock_get_service):
        """Test helper functions when service raises errors."""
        mock_service = Mock()
        mock_service.get_status.side_effect = Exception("Service error")
        mock_get_service.return_value = mock_service

        # Should propagate the exception
        with pytest.raises(Exception, match="Service error"):
            work_service.get_work_status("work123")
