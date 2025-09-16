"""
Tests for src.application.services module.

Tests direct import functionality and service accessibility.
"""

import pytest


class TestApplicationServicesModule:
    """Test the application services module imports."""

    def test_direct_import_synthetic_dataset_generator(self):
        """Test direct import of SyntheticDatasetGenerator."""
        from src.application.services.dataset_generator import SyntheticDatasetGenerator
        
        # Should be able to import without raising an exception
        assert SyntheticDatasetGenerator is not None
        # Should be a class
        assert isinstance(SyntheticDatasetGenerator, type)

    def test_direct_import_dataset_service_function(self):
        """Test direct import of load_dataset function."""
        from src.application.services.dataset_service import load_dataset
        
        # Should be able to import without raising an exception
        assert load_dataset is not None
        # Should be a function
        assert callable(load_dataset)

    def test_getattr_invalid_service(self):
        """Test that accessing invalid service from empty __init__ raises AttributeError."""
        import src.application.services as services
        
        with pytest.raises(AttributeError, match="InvalidService"):
            _ = services.InvalidService

    def test_direct_import_work_service_functions(self):
        """Test accessing work service functions directly."""
        from src.application.services.work_service import get_work_service
        
        # Should be able to import function
        assert get_work_service is not None
        assert callable(get_work_service)

    def test_services_module_shows_submodules(self):
        """Test that services module shows imported submodules."""
        import src.application.services as services
        
        # Python automatically adds submodules to the namespace when they exist
        module_attrs = [attr for attr in dir(services) if not attr.startswith('_')]
        
        # Should show the existing submodules
        expected_modules = {'dataset_generator', 'dataset_service', 'work_service'}
        actual_modules = set(module_attrs)
        
        # All expected modules should be present
        assert expected_modules.issubset(actual_modules)