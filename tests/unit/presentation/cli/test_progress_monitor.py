"""
Test suite for ProgressMonitor curses-based monitoring interface.
"""

import curses
import signal
import time
from unittest.mock import Mock, patch

import pytest

from src.infrastructure.persistence.work_state.core import (
    ExecutionDetail,
    ProgressSummary,
)
from src.presentation.cli.progress_monitor import ProgressMonitor, main


class TestProgressMonitor:
    """Test ProgressMonitor class."""

    @pytest.fixture
    def mock_persistence(self):
        """Create mock persistence."""
        mock = Mock()
        mock.work_get.return_value = {"status": "running", "id": "test-work-id"}

        progress_summary = Mock(spec=ProgressSummary)
        progress_summary.global_progress = 0.5
        progress_summary.global_execution = {"Finished": 50, "Total": 100}
        progress_summary.tasks = {"Running": "Task-1", "Finished": 0, "Total": 1}
        progress_summary.datasets = {"Running": "Dataset-1", "Finished": 0, "Total": 1}
        progress_summary.configs = {"Running": "Config-1", "Finished": 0, "Total": 1}
        progress_summary.algorithms = {
            "Running": "Algorithm-1",
            "Finished": 0,
            "Total": 1,
        }
        progress_summary.execution = {"Finished": 25, "Running": 1, "Total": 50}
        progress_summary.current_combination_details = {
            "task_id": "Task-1",
            "dataset_id": "Dataset-1",
            "preset_id": "Config-1",
        }
        mock.get_work_progress_summary.return_value = progress_summary

        exec_detail = Mock(spec=ExecutionDetail)
        exec_detail.sequencia = 1
        exec_detail.algorithm_id = "TestAlgorithm"
        exec_detail.status = "running"
        exec_detail.progress = 0.75
        exec_detail.progress_message = "Processing..."
        exec_detail.started_at = time.time() - 10
        exec_detail.finished_at = None
        mock.get_combination_executions_detail.return_value = [exec_detail]

        mock.get_events.return_value = []
        mock.close.return_value = None
        return mock

    def test_init_successful(self, mock_persistence):
        """Test successful initialization."""
        with patch("src.presentation.cli.progress_monitor.WorkPersistence") as mock_wp:
            mock_wp.return_value = mock_persistence
            with patch.object(
                ProgressMonitor, "_wait_for_work_ready", return_value=True
            ):
                monitor = ProgressMonitor("test-work-id")

                assert monitor.work_id == "test-work-id"
                assert monitor.running is True
                assert monitor.refresh_interval == 0.1
                assert monitor.scroll_pos == 0

    def test_init_work_not_found(self, mock_persistence):
        """Test initialization failure."""
        mock_persistence.work_get.return_value = None

        with patch("src.presentation.cli.progress_monitor.WorkPersistence") as mock_wp:
            mock_wp.return_value = mock_persistence
            with patch.object(
                ProgressMonitor, "_wait_for_work_ready", return_value=False
            ):
                with pytest.raises(ValueError, match="Timeout.*not found"):
                    ProgressMonitor("nonexistent-work")

    def test_signal_handler_pause(self, mock_persistence):
        """Test signal handler."""
        with patch("src.presentation.cli.progress_monitor.WorkPersistence") as mock_wp:
            mock_wp.return_value = mock_persistence
            with patch.object(
                ProgressMonitor, "_wait_for_work_ready", return_value=True
            ):
                mock_work_service = Mock()
                mock_work_service.pause.return_value = True

                with patch(
                    "src.application.services.work_service.get_work_service"
                ) as mock_gws:
                    mock_gws.return_value = mock_work_service

                    monitor = ProgressMonitor("test-work-id")
                    monitor.signal_handler(signal.SIGINT, None)

                    mock_work_service.pause.assert_called_once_with("test-work-id")
                    assert monitor.running is False

    def test_format_progress_bar(self, mock_persistence):
        """Test progress bar formatting."""
        with patch("src.presentation.cli.progress_monitor.WorkPersistence") as mock_wp:
            mock_wp.return_value = mock_persistence
            with patch.object(
                ProgressMonitor, "_wait_for_work_ready", return_value=True
            ):
                monitor = ProgressMonitor("test-work-id")

                bar = monitor.format_progress_bar(0.5, 10, show_percent=True)
                assert "█████─────" in bar
                assert "50.0%" in bar

    def test_format_execution_time(self, mock_persistence):
        """Test execution time formatting."""
        with patch("src.presentation.cli.progress_monitor.WorkPersistence") as mock_wp:
            mock_wp.return_value = mock_persistence
            with patch.object(
                ProgressMonitor, "_wait_for_work_ready", return_value=True
            ):
                monitor = ProgressMonitor("test-work-id")

                time_str = monitor.format_execution_time(125.456)
                assert time_str == "02:05.456"

                time_str = monitor.format_execution_time(45.123)
                assert time_str == "45.123s"

    def test_draw_methods_exist(self, mock_persistence):
        """Test that drawing methods exist and are callable."""
        with patch("src.presentation.cli.progress_monitor.WorkPersistence") as mock_wp:
            mock_wp.return_value = mock_persistence
            with patch.object(
                ProgressMonitor, "_wait_for_work_ready", return_value=True
            ):
                monitor = ProgressMonitor("test-work-id")

                # Test that drawing methods exist
                assert hasattr(monitor, "draw_header")
                assert hasattr(monitor, "draw_executions_list")
                assert hasattr(monitor, "draw_footer")
                assert hasattr(monitor, "draw_log_panel")
                assert callable(monitor.draw_header)
                assert callable(monitor.draw_executions_list)

    @patch("curses.curs_set")
    @patch("curses.start_color")
    @patch("curses.init_pair")
    @patch("curses.color_pair", return_value=1)
    def test_run_basic_loop(
        self,
        mock_color_pair,
        mock_init_pair,
        mock_start_color,
        mock_curs_set,
        mock_persistence,
    ):
        """Test basic run loop."""
        with patch("src.presentation.cli.progress_monitor.WorkPersistence") as mock_wp:
            mock_wp.return_value = mock_persistence
            with patch.object(
                ProgressMonitor, "_wait_for_work_ready", return_value=True
            ):
                monitor = ProgressMonitor("test-work-id")

                mock_stdscr = Mock()
                mock_stdscr.getmaxyx.return_value = (24, 80)
                mock_stdscr.getch.return_value = ord("q")
                mock_stdscr.addstr = Mock()
                mock_stdscr.refresh = Mock()
                mock_stdscr.clear = Mock()
                mock_stdscr.nodelay = Mock()
                mock_stdscr.timeout = Mock()

                monitor.run(mock_stdscr)

                mock_curs_set.assert_called_once_with(0)
                mock_stdscr.nodelay.assert_called_once_with(1)
                mock_start_color.assert_called_once()

    @patch("curses.wrapper")
    def test_start_success(self, mock_wrapper, mock_persistence):
        """Test successful monitor start."""
        with patch("src.presentation.cli.progress_monitor.WorkPersistence") as mock_wp:
            mock_wp.return_value = mock_persistence
            with patch.object(
                ProgressMonitor, "_wait_for_work_ready", return_value=True
            ):
                with patch("signal.signal") as mock_signal:
                    monitor = ProgressMonitor("test-work-id")

                    result = monitor.start()

                    assert result == 0
                    mock_signal.assert_called_once()
                    mock_wrapper.assert_called_once()


class TestMainFunction:
    """Test main function."""

    @patch("argparse.ArgumentParser.parse_args")
    def test_main_successful_execution(self, mock_parse_args):
        """Test successful execution."""
        mock_args = Mock()
        mock_args.work_id = "test-work-id"
        mock_parse_args.return_value = mock_args

        with patch("src.presentation.cli.progress_monitor.ProgressMonitor") as mock_pm:
            mock_monitor = Mock()
            mock_monitor.start.return_value = 0
            mock_pm.return_value = mock_monitor

            result = main()

            assert result == 0
            mock_pm.assert_called_once_with("test-work-id")

    @patch("argparse.ArgumentParser.parse_args")
    def test_main_value_error(self, mock_parse_args):
        """Test main with value error."""
        mock_args = Mock()
        mock_args.work_id = "test-work-id"
        mock_parse_args.return_value = mock_args

        with patch("src.presentation.cli.progress_monitor.ProgressMonitor") as mock_pm:
            mock_pm.side_effect = ValueError("Work not found")

            result = main()

            assert result == 1
