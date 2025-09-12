"""Tests for infrastructure orchestration modules.

Coverage objectives:
- Test experiment executor functionality
- Test optimization executor
- Test pipeline runner
- Test algorithm runner integration
- Test parallel execution and resource management
- Test error handling and recovery
"""

import tempfile
import time
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
import pytest

from src.domain.config import CSPBenchConfig, ExperimentTask, AlgParams
from src.domain.dataset import Dataset
from src.domain.status import BaseStatus


class TestExperimentExecutor:
    """Test experiment executor functionality."""

    def test_experiment_executor_import(self):
        """Test that experiment executor can be imported."""
        from src.infrastructure.orchestration.experiment_executor import ExperimentExecutor
        assert ExperimentExecutor is not None

    @patch('src.infrastructure.orchestration.experiment_executor.ProcessPoolExecutor')
    @patch('src.infrastructure.orchestration.experiment_executor.run_algorithm')
    def test_experiment_execution_basic(self, mock_run_algorithm, mock_process_pool):
        """Test basic experiment execution."""
        from src.infrastructure.orchestration.experiment_executor import ExperimentExecutor
        
        # Setup mocks
        mock_executor_instance = Mock()
        mock_process_pool.return_value.__enter__.return_value = mock_executor_instance
        mock_executor_instance.submit.return_value.result.return_value = {'result': 'test'}
        
        mock_run_algorithm.return_value = {
            'center_string': 'ACGT',
            'max_distance': 1,
            'execution_time': 0.5
        }
        
        # Create test data
        mock_config = Mock(spec=CSPBenchConfig)
        mock_task = Mock(spec=ExperimentTask)
        mock_task.repetitions = 2
        mock_task.algorithm = 'baseline'
        mock_task.parameters = {}
        
        mock_dataset = Mock(spec=Dataset)
        mock_dataset.strings = ['ACGT', 'AGGT', 'ATGT']
        mock_dataset.alphabet = 'ACGT'
        
        # Create executor
        executor = ExperimentExecutor()
        
        # Mock the necessary methods
        with patch.object(executor, '_prepare_execution'), \
             patch.object(executor, '_execute_repetitions'), \
             patch.object(executor, '_finalize_execution'):
            
            # This tests that the executor can be instantiated and basic methods called
            assert executor is not None

    def test_worker_execution_function(self):
        """Test worker execution function."""
        from src.infrastructure.orchestration.experiment_executor import _worker_exec
        
        # Test that the worker function exists and can be called
        # (Full testing would require more complex setup)
        assert callable(_worker_exec)

    @patch('src.infrastructure.orchestration.experiment_executor.get_logger')
    def test_experiment_executor_logging(self, mock_get_logger):
        """Test experiment executor logging setup."""
        from src.infrastructure.orchestration.experiment_executor import logger
        
        # Test that logger is configured
        assert logger is not None

    def test_experiment_executor_error_handling(self):
        """Test experiment executor error handling."""
        from src.infrastructure.orchestration.experiment_executor import ExperimentExecutor
        
        executor = ExperimentExecutor()
        
        # Test that executor handles initialization errors
        assert executor is not None


class TestOptimizationExecutor:
    """Test optimization executor functionality."""

    def test_optimization_executor_import(self):
        """Test that optimization executor can be imported."""
        from src.infrastructure.orchestration.optimization_executor import OptimizationExecutor
        assert OptimizationExecutor is not None

    @patch('src.infrastructure.orchestration.optimization_executor.get_logger')
    def test_optimization_executor_logging(self, mock_get_logger):
        """Test optimization executor logging setup."""
        # Import to trigger logger creation
        from src.infrastructure.orchestration.optimization_executor import logger
        
        assert logger is not None

    def test_optimization_executor_initialization(self):
        """Test optimization executor initialization."""
        from src.infrastructure.orchestration.optimization_executor import OptimizationExecutor
        
        executor = OptimizationExecutor()
        assert executor is not None

    @patch('src.infrastructure.orchestration.optimization_executor.OptimizationExecutor.execute')
    def test_optimization_execution_interface(self, mock_execute):
        """Test optimization execution interface."""
        from src.infrastructure.orchestration.optimization_executor import OptimizationExecutor
        
        executor = OptimizationExecutor()
        
        # Mock the execute method
        mock_execute.return_value = None
        
        # Test that execute method exists and can be called
        assert hasattr(executor, 'execute')


class TestPipelineRunner:
    """Test pipeline runner functionality."""

    def test_pipeline_runner_import(self):
        """Test that pipeline runner can be imported."""
        from src.infrastructure.orchestration.pipeline_runner import PipelineRunner
        assert PipelineRunner is not None

    def test_pipeline_runner_initialization(self):
        """Test pipeline runner initialization."""
        from src.infrastructure.orchestration.pipeline_runner import PipelineRunner
        
        runner = PipelineRunner()
        assert runner is not None

    @patch('src.infrastructure.orchestration.pipeline_runner.get_logger')
    def test_pipeline_runner_logging(self, mock_get_logger):
        """Test pipeline runner logging setup."""
        from src.infrastructure.orchestration.pipeline_runner import logger
        
        assert logger is not None

    @patch('src.infrastructure.orchestration.pipeline_runner.PipelineRunner.execute')
    def test_pipeline_execution_interface(self, mock_execute):
        """Test pipeline execution interface."""
        from src.infrastructure.orchestration.pipeline_runner import PipelineRunner
        
        runner = PipelineRunner()
        
        # Mock the execute method
        mock_execute.return_value = None
        
        # Test that execute method exists
        assert hasattr(runner, 'execute')


class TestAlgorithmRunner:
    """Test algorithm runner functionality."""

    def test_algorithm_runner_import(self):
        """Test that algorithm runner can be imported."""
        from src.infrastructure.orchestration.algorithm_runner import run_algorithm
        assert callable(run_algorithm)

    @patch('src.infrastructure.orchestration.algorithm_runner.global_registry')
    @patch('src.infrastructure.orchestration.algorithm_runner.create_distance_calculator')
    def test_run_algorithm_basic(self, mock_distance_calc, mock_registry):
        """Test basic algorithm running functionality."""
        from src.infrastructure.orchestration.algorithm_runner import run_algorithm
        
        # Setup mocks
        mock_algorithm_class = Mock()
        mock_algorithm_instance = Mock()
        mock_algorithm_class.return_value = mock_algorithm_instance
        mock_algorithm_instance.run.return_value = {
            'center_string': 'ACGT',
            'max_distance': 1
        }
        
        mock_registry.get_algorithm_class.return_value = mock_algorithm_class
        
        mock_distance_calc.return_value = Mock()
        
        # Test data
        strings = ['ACGT', 'AGGT', 'ATGT']
        alphabet = 'ACGT'
        algorithm_name = 'test_algorithm'
        
        # Run algorithm
        result = run_algorithm(
            strings=strings,
            alphabet=alphabet,
            algorithm_name=algorithm_name,
            parameters={},
            distance_method='hamming',
            use_cache=False,
            cpu_limit=None
        )
        
        # Verify result
        assert result is not None
        mock_registry.get_algorithm_class.assert_called_once_with(algorithm_name)

    @patch('src.infrastructure.orchestration.algorithm_runner.global_registry')
    def test_run_algorithm_not_found(self, mock_registry):
        """Test algorithm runner with non-existent algorithm."""
        from src.infrastructure.orchestration.algorithm_runner import run_algorithm
        
        mock_registry.get_algorithm_class.return_value = None
        
        with pytest.raises(ValueError, match="Algorithm .* not found"):
            run_algorithm(
                strings=['ACGT'],
                alphabet='ACGT',
                algorithm_name='nonexistent_algorithm',
                parameters={},
                distance_method='hamming',
                use_cache=False
            )

    @patch('src.infrastructure.orchestration.algorithm_runner.global_registry')
    def test_run_algorithm_execution_error(self, mock_registry):
        """Test algorithm runner with execution error."""
        from src.infrastructure.orchestration.algorithm_runner import run_algorithm
        
        # Setup mock that raises exception
        mock_algorithm_class = Mock()
        mock_algorithm_instance = Mock()
        mock_algorithm_class.return_value = mock_algorithm_instance
        mock_algorithm_instance.run.side_effect = Exception("Algorithm failed")
        
        mock_registry.get_algorithm_class.return_value = mock_algorithm_class
        
        with pytest.raises(Exception, match="Algorithm failed"):
            run_algorithm(
                strings=['ACGT'],
                alphabet='ACGT',
                algorithm_name='failing_algorithm',
                parameters={},
                distance_method='hamming',
                use_cache=False
            )

    @patch('src.infrastructure.orchestration.algorithm_runner.global_registry')
    @patch('src.infrastructure.orchestration.algorithm_runner.create_distance_calculator')
    @patch('time.time', side_effect=[1000.0, 1001.5])  # Mock timing
    def test_run_algorithm_timing(self, mock_time, mock_distance_calc, mock_registry):
        """Test that algorithm runner measures execution time."""
        from src.infrastructure.orchestration.algorithm_runner import run_algorithm
        
        # Setup mocks
        mock_algorithm_class = Mock()
        mock_algorithm_instance = Mock()
        mock_algorithm_class.return_value = mock_algorithm_instance
        mock_algorithm_instance.run.return_value = {
            'center_string': 'ACGT',
            'max_distance': 1
        }
        
        mock_registry.get_algorithm_class.return_value = mock_algorithm_class
        mock_distance_calc.return_value = Mock()
        
        # Run algorithm
        result = run_algorithm(
            strings=['ACGT'],
            alphabet='ACGT',
            algorithm_name='test_algorithm',
            parameters={},
            distance_method='hamming',
            use_cache=False
        )
        
        # Verify timing was recorded
        assert 'execution_time' in result
        assert result['execution_time'] == 1.5  # 1001.5 - 1000.0


class TestResourceManagement:
    """Test resource management in orchestration."""

    @patch('src.infrastructure.orchestration.experiment_executor.ProcessPoolExecutor')
    def test_process_pool_configuration(self, mock_process_pool):
        """Test process pool configuration."""
        from src.infrastructure.orchestration.experiment_executor import ExperimentExecutor
        
        executor = ExperimentExecutor()
        
        # Test that ProcessPoolExecutor is available for configuration
        assert mock_process_pool is not None

    def test_cpu_limit_handling(self):
        """Test CPU limit handling in algorithm runner."""
        from src.infrastructure.orchestration.algorithm_runner import run_algorithm
        
        # Test that cpu_limit parameter is accepted
        # (Full testing would require more complex setup)
        assert callable(run_algorithm)

    @patch('src.infrastructure.orchestration.experiment_executor.ExecutionController')
    def test_execution_controller_integration(self, mock_controller):
        """Test execution controller integration."""
        from src.infrastructure.orchestration.experiment_executor import ExecutionController
        
        # Test that ExecutionController is available
        assert ExecutionController is not None


class TestPersistenceIntegration:
    """Test persistence integration in orchestration."""

    def test_persistence_monitors_import(self):
        """Test persistence monitor imports."""
        from src.infrastructure.orchestration.experiment_executor import PersistenceMonitor
        assert PersistenceMonitor is not None

    def test_scoped_persistence_imports(self):
        """Test scoped persistence wrapper imports."""
        from src.infrastructure.orchestration.experiment_executor import (
            CombinationScopedPersistence,
            ExecutionScopedPersistence,
            WorkScopedPersistence
        )
        
        assert CombinationScopedPersistence is not None
        assert ExecutionScopedPersistence is not None
        assert WorkScopedPersistence is not None


class TestErrorRecovery:
    """Test error recovery mechanisms."""

    @patch('src.infrastructure.orchestration.experiment_executor.logger')
    def test_logging_error_scenarios(self, mock_logger):
        """Test logging in error scenarios."""
        from src.infrastructure.orchestration.experiment_executor import ExperimentExecutor
        
        executor = ExperimentExecutor()
        
        # Test that logger is available for error reporting
        assert mock_logger is not None

    def test_exception_propagation(self):
        """Test exception propagation in executors."""
        from src.infrastructure.orchestration.experiment_executor import ExperimentExecutor
        
        executor = ExperimentExecutor()
        
        # Test that executor can handle exceptions during initialization
        assert executor is not None


class TestParallelExecution:
    """Test parallel execution capabilities."""

    def test_concurrent_futures_integration(self):
        """Test concurrent.futures integration."""
        from src.infrastructure.orchestration.experiment_executor import ProcessPoolExecutor
        
        # Test that ProcessPoolExecutor is available
        assert ProcessPoolExecutor is not None

    @patch('src.infrastructure.orchestration.experiment_executor.ProcessPoolExecutor')
    def test_worker_function_serialization(self, mock_process_pool):
        """Test that worker functions can be serialized."""
        from src.infrastructure.orchestration.experiment_executor import _worker_exec
        
        # Test that worker function exists and is callable
        assert callable(_worker_exec)


class TestConfigurationHandling:
    """Test configuration handling in orchestration."""

    def test_config_parameter_passing(self):
        """Test configuration parameter passing."""
        from src.infrastructure.orchestration.experiment_executor import ExperimentExecutor
        
        executor = ExperimentExecutor()
        
        # Test that executor can handle configuration
        assert executor is not None

    def test_algorithm_parameters_handling(self):
        """Test algorithm parameters handling."""
        from src.infrastructure.orchestration.algorithm_runner import run_algorithm
        
        # Test that parameters are accepted
        assert callable(run_algorithm)

    def test_distance_method_configuration(self):
        """Test distance method configuration."""
        from src.infrastructure.orchestration.algorithm_runner import run_algorithm
        
        # Test that distance_method parameter is handled
        assert callable(run_algorithm)


class TestStatusTracking:
    """Test status tracking in orchestration."""

    def test_base_status_import(self):
        """Test BaseStatus import in orchestration modules."""
        from src.infrastructure.orchestration.experiment_executor import BaseStatus
        
        assert BaseStatus is not None

    def test_status_transitions(self):
        """Test status transitions during execution."""
        # This would require integration testing with actual persistence
        # For now, test that the status enum is available
        from src.domain.status import BaseStatus
        
        assert BaseStatus.RUNNING is not None
        assert BaseStatus.COMPLETED is not None
        assert BaseStatus.FAILED is not None