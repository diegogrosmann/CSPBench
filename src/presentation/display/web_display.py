"""Web display adapter for progress events - Refactored for structured batch tracking."""

import logging
from datetime import datetime
from typing import Any, Dict, Optional


class WebDisplay:
    """
    Web display adapter for generating structured batch tracking data.

    Generates three key JSON structures:
    1. batch_structure.json - Static batch structure (generated once)
    2. execution_status.json - Dynamic execution status (updated continuously)
    3. callbacks_history.json - Dynamic callbacks history (updated per repetition)
    """

    def __init__(self, session_manager, session_id: str):
        """
        Initialize web display.

        Args:
            session_manager: Web session manager
            session_id: Session ID to update
        """
        self._session_manager = session_manager
        self._session_id = session_id
        self._logger = logging.getLogger(__name__)

        # Batch structure data (static, generated once)
        self._batch_structure: Optional[Dict[str, Any]] = None
        self._batch_metadata: Optional[Dict[str, Any]] = None

        # Execution status data (dynamic, updated continuously)
        self._execution_status: Dict[str, Any] = {
            "status": "not_started",
            "start_time": None,
            "end_time": None,
            "total_repetitions": 0,
            "completed_repetitions_total": 0,
            "batch_progress": 0.0,
            "current_execution": {},
        }

        # Callbacks history data (dynamic, updated per repetition)
        self._callbacks_history: Dict[str, Any] = {
            "batch_name": "Unknown Batch",
            "history": [],
        }

        # Tracking current state
        self._current_execution_name = ""
        self._current_dataset_id = ""
        self._current_config_id = ""
        self._current_algorithm = ""
        self._current_repetition = 0
        self._completed_combinations_count = 0
        self._current_repetition_entry = None  # Current repetition entry for callbacks

        # Track completed datasets and configs
        self._completed_datasets = set()
        self._completed_configs = set()

        # Track algorithm completion per dataset/config combination
        self._completed_algorithms_per_dataset_config = {}

        # Track completed repetitions per execution context
        self._completed_repetitions_count = 0

        self._logger.info(f"[WEB_DISPLAY] Initialized for session {session_id}")

    def handle_event(self, event: Any) -> None:
        """Handle a progress event and update appropriate data structures."""
        try:
            event_type = type(event).__name__
            self._logger.info(f"[WEB_DISPLAY] 🔥 HANDLING EVENT: {event_type}")

            # Log detailed event attributes for debugging
            self._logger.info(
                f"[WEB_DISPLAY] 📊 Event attributes: {[attr for attr in dir(event) if not attr.startswith('_')]}"
            )

            # Log specific key attributes if they exist
            if hasattr(event, "execution_name"):
                self._logger.info(
                    f"[WEB_DISPLAY] 🎯 execution_name: {event.execution_name}"
                )
            if hasattr(event, "context"):
                self._logger.info(f"[WEB_DISPLAY] 🔍 context: {event.context}")
            if hasattr(event, "metadata"):
                self._logger.info(f"[WEB_DISPLAY] 🏷️ metadata: {event.metadata}")
            if hasattr(event, "algorithm_name"):
                self._logger.info(
                    f"[WEB_DISPLAY] ⚙️ algorithm_name: {event.algorithm_name}"
                )
            if hasattr(event, "message"):
                self._logger.info(f"[WEB_DISPLAY] 💬 message: {event.message}")

            # Route to specific handlers (including unified DisplayEvent)
            if event_type == "DisplayEvent":
                self._handle_display(event)
            elif hasattr(self, f'_handle_{event_type.lower().replace("event", "")}'):
                handler = getattr(
                    self, f'_handle_{event_type.lower().replace("event", "")}'
                )
                handler(event)
            else:
                self._handle_generic(event)

            # Always update session after handling event
            self._update_web_session()

        except Exception as e:
            self._logger.error(
                f"[WEB_DISPLAY] ❌ Error handling event {type(event).__name__}: {e}",
                exc_info=True,
            )

    def _handle_display(self, event) -> None:
        """Handle unified DisplayEvent by appending to callbacks history and updating current context."""
        try:
            phase = getattr(event, "phase", None)
            dataset_id = getattr(event, "dataset_id", None)
            algorithm_name = getattr(event, "algorithm_name", None)
            task_id = getattr(event, "task_id", None)
            trial_no = getattr(event, "trial_no", None)
            rep_idx = getattr(event, "rep_idx", None)
            progress = getattr(event, "progress", 0.0)
            message = getattr(event, "message", "")
            payload = getattr(event, "payload", {}) or {}

            # Update current context
            if dataset_id:
                self._current_dataset_id = dataset_id
            if algorithm_name:
                self._current_algorithm = algorithm_name

            # Append history entry
            entry = {
                "ts": datetime.now().isoformat(),
                "phase": getattr(phase, "value", str(phase)),
                "dataset": dataset_id,
                "algorithm": algorithm_name,
                "task_id": task_id,
                "trial": trial_no,
                "rep": rep_idx,
                "progress": progress,
                "message": message,
                "payload": payload,
            }
            self._callbacks_history.setdefault("history", []).append(entry)

            # Update execution status
            self._execution_status["current_execution"] = {
                "dataset": self._current_dataset_id,
                "algorithm": self._current_algorithm,
                "phase": getattr(phase, "value", str(phase)),
                "progress": progress,
                "message": message,
            }
        except Exception as e:
            self._logger.warning(f"[WEB_DISPLAY] Failed to handle DisplayEvent: {e}")

    def _handle_taskstarted(self, event) -> None:
        """Handle TaskStartedEvent - Initialize batch structure."""
        self._logger.info(f"[WEB_DISPLAY] 🚀 TASK STARTED: {event.task_name}")
        self._logger.info(
            f"[WEB_DISPLAY] 📋 Task type: {getattr(event, 'task_type', 'unknown')}"
        )
        self._logger.info(
            f"[WEB_DISPLAY] 🕒 Timestamp: {getattr(event, 'timestamp', 'unknown')}"
        )

        # Deep metadata logging
        if hasattr(event, "metadata"):
            self._logger.info(
                f"[WEB_DISPLAY] 🗂️ RAW METADATA: {len(str(event.metadata))} characters"
            )

            # Check for batch_config in metadata
            if "batch_config" in event.metadata:
                batch_config = event.metadata["batch_config"]
                self._logger.info(
                    f"[WEB_DISPLAY] 🔧 BATCH CONFIG FOUND: {type(batch_config)}"
                )
                self._initialize_batch_structure({"batch_config": batch_config})
            else:
                # Check if metadata contains direct configuration
                if all(
                    key in event.metadata
                    for key in ["metadata", "execution", "datasets", "algorithms"]
                ):
                    self._logger.info(
                        "[WEB_DISPLAY] 🔧 FULL CONFIG IN METADATA - Using as batch_config"
                    )
                    self._initialize_batch_structure({"batch_config": event.metadata})
                else:
                    self._logger.warning("[WEB_DISPLAY] ⚠️ NO BATCH_CONFIG IN METADATA")
        else:
            self._logger.warning("[WEB_DISPLAY] ⚠️ NO METADATA IN TaskStartedEvent")

        # Update execution status
        if self._execution_status["status"] == "not_started":
            self._execution_status["status"] = "running"
            self._execution_status["start_time"] = datetime.now().isoformat()
            self._logger.info(
                f"[WEB_DISPLAY] ✅ Execution status updated to running at {self._execution_status['start_time']}"
            )

    def _handle_taskfinished(self, event) -> None:
        """Handle TaskFinishedEvent - Finalize execution status."""
        self._logger.info(f"[WEB_DISPLAY] Task finished: {event.task_name}")

        # Update execution status
        if event.success:
            self._execution_status["status"] = "completed"
        else:
            self._execution_status["status"] = "failed"

        self._execution_status["end_time"] = datetime.now().isoformat()
        self._execution_status["batch_progress"] = 1.0

    def _handle_executionstarted(self, event) -> None:
        """Handle ExecutionStartedEvent - Set current execution context and enhance batch structure."""
        self._current_execution_name = event.execution_name
        self._logger.info(
            f"[WEB_DISPLAY] Execution started: {self._current_execution_name}"
        )
        self._logger.info(f"[WEB_DISPLAY] Total items: {event.total_items}")

        # If batch structure wasn't properly initialized, create a basic one
        if not self._batch_structure:
            self._logger.info("[WEB_DISPLAY] Creating fallback batch structure")
            self._create_fallback_batch_structure()

        # Update total repetitions estimate if available, but don't overwrite calculated total
        if hasattr(event, "total_items") and event.total_items > 0:
            # Only update if we don't have a calculated total yet
            if self._execution_status["total_repetitions"] == 0:
                self._execution_status["total_repetitions"] = event.total_items
                if self._batch_structure:
                    self._batch_structure["total_repetitions"] = event.total_items
                self._logger.info(
                    f"[WEB_DISPLAY] ✅ Using event total_items as total_repetitions: {event.total_items}"
                )
            else:
                self._logger.info(
                    f"[WEB_DISPLAY] 🔒 Keeping calculated total_repetitions: {self._execution_status['total_repetitions']} (ignoring event.total_items: {event.total_items})"
                )

        # Update execution status current execution
        self._execution_status["current_execution"] = {
            "name": self._current_execution_name,
            "dataset": {"id": "unknown", "total": 0, "completed": 0},
            "config_id": {"id": "unknown", "total": 0, "completed": 0},
            "algorithm": {
                "name": "unknown",
                "total": 0,
                "completed": 0,
                "repetitions": {"total": 0, "completed": 0},
            },
        }

    def _handle_executionprogress(self, event) -> None:
        """Handle ExecutionProgressEvent - Track individual repetitions."""
        context = getattr(event, "context", {})
        dataset_name = context.get("dataset_name", "Unknown")
        algorithm_config = context.get("algorithm_config_name", "Unknown")
        algorithm_name = context.get("algorithm_name", "Unknown")

        self._logger.info("[WEB_DISPLAY] 📈 EXECUTION PROGRESS EVENT")
        self._logger.info(f"[WEB_DISPLAY] 🎯 execution_name: {event.execution_name}")
        self._logger.info(
            f"[WEB_DISPLAY] 🗂️ dataset_name (from context): {dataset_name}"
        )
        self._logger.info(
            f"[WEB_DISPLAY] ⚙️ algorithm_config_name (from context): {algorithm_config}"
        )
        self._logger.info(
            f"[WEB_DISPLAY] 🔧 algorithm_name (from context): {algorithm_name}"
        )
        self._logger.info(
            f"[WEB_DISPLAY] 📊 current_item: {getattr(event, 'current_item', 'not_found')}"
        )
        self._logger.info(
            f"[WEB_DISPLAY] 📊 total_items: {getattr(event, 'total_items', 'not_found')}"
        )
        self._logger.info(
            f"[WEB_DISPLAY] 💬 message: {getattr(event, 'message', 'not_found')}"
        )
        self._logger.info(f"[WEB_DISPLAY] 🔍 FULL CONTEXT: {context}")

        # Update current context
        old_execution = self._current_execution_name
        old_dataset = self._current_dataset_id
        old_config = self._current_config_id
        old_algorithm = self._current_algorithm
        old_repetition = self._current_repetition

        self._current_execution_name = event.execution_name
        self._current_dataset_id = self._get_dataset_id_from_name(dataset_name)
        self._current_config_id = self._get_config_id_from_name(algorithm_config)
        self._current_algorithm = algorithm_name
        self._current_repetition += 1

        self._logger.info("[WEB_DISPLAY] 🔄 CONTEXT CHANGES:")
        self._logger.info(
            f"[WEB_DISPLAY]    execution: {old_execution} → {self._current_execution_name}"
        )
        self._logger.info(
            f"[WEB_DISPLAY]    dataset: {old_dataset} → {self._current_dataset_id}"
        )
        self._logger.info(
            f"[WEB_DISPLAY]    config: {old_config} → {self._current_config_id}"
        )
        self._logger.info(
            f"[WEB_DISPLAY]    algorithm: {old_algorithm} → {self._current_algorithm}"
        )
        self._logger.info(
            f"[WEB_DISPLAY]    repetition: {old_repetition} → {self._current_repetition}"
        )

        # Update execution status
        self._update_execution_status_progress()

        # Initialize new repetition in callbacks history
        self._start_new_repetition()

    def _handle_executionfinished(self, event) -> None:
        """Handle ExecutionFinishedEvent."""
        self._logger.info(f"[WEB_DISPLAY] Execution finished: {event.execution_name}")

        # Finalize current repetition if exists
        if self._current_repetition_entry:
            self._current_repetition_entry["status"] = (
                "completed" if event.success else "failed"
            )

    def _handle_algorithmprogress(self, event) -> None:
        """Handle AlgorithmProgressEvent - Add to current repetition callbacks."""
        algorithm_name = getattr(event, "algorithm_name", "Unknown")
        message = getattr(event, "message", "")
        progress_percent = getattr(event, "progress_percent", 0.0)

        self._logger.info("[WEB_DISPLAY] 🔧 ALGORITHM PROGRESS:")
        self._logger.info(f"[WEB_DISPLAY]    algorithm_name: {algorithm_name}")
        self._logger.info(f"[WEB_DISPLAY]    message: {message}")
        self._logger.info(f"[WEB_DISPLAY]    progress_percent: {progress_percent}")

        if self._current_repetition_entry and message.strip():
            # Filter out generic "Completed run X/Y" messages as they are not from algorithm internals
            raw_message = message.strip()

            # Skip generic execution status messages - these are from the orchestrator, not the algorithm
            if "Completed run" in raw_message and "/" in raw_message:
                self._logger.debug(
                    f"[WEB_DISPLAY] 🚫 Skipping orchestrator message: {raw_message}"
                )
                return

            # Only add meaningful algorithm progress messages
            if raw_message and raw_message not in ["algorithm_progress", "start"]:
                self._current_repetition_entry["callbacks"].append(raw_message)
                self._logger.info(
                    f"[WEB_DISPLAY] ✅ Added RAW algorithm callback: {raw_message}"
                )

            # Update algorithm name in current repetition if it was unknown
            if (
                self._current_repetition_entry["algorithm"] == "Unknown"
                and algorithm_name != "Unknown"
            ):
                self._current_repetition_entry["algorithm"] = algorithm_name
                self._logger.info(
                    f"[WEB_DISPLAY] 🔄 Updated repetition algorithm: Unknown → {algorithm_name}"
                )
        else:
            if not self._current_repetition_entry:
                self._logger.warning(
                    "[WEB_DISPLAY] ⚠️ No current repetition entry to add callback to"
                )
            if not message.strip():
                self._logger.debug("[WEB_DISPLAY] 🚫 Skipping empty message")

    def _handle_algorithmfinished(self, event) -> None:
        """Handle AlgorithmFinishedEvent - Finalize current repetition."""
        success = getattr(event, "success", True)
        algorithm_name = getattr(event, "algorithm_name", "Algorithm")

        self._logger.info("[WEB_DISPLAY] 🏁 ALGORITHM FINISHED:")
        self._logger.info(f"[WEB_DISPLAY]    algorithm_name: {algorithm_name}")
        self._logger.info(f"[WEB_DISPLAY]    success: {success}")

        # Update current algorithm name from the event (more reliable than context)
        if algorithm_name != "Algorithm" and algorithm_name != "Unknown":
            self._current_algorithm = algorithm_name
            self._logger.info(
                f"[WEB_DISPLAY] 🔄 Updated current algorithm: {self._current_algorithm}"
            )

        # Track algorithm completion for current dataset/config combination
        if success and self._current_dataset_id and self._current_config_id:
            key = f"{self._current_dataset_id}_{self._current_config_id}"
            # Note: This assumes one repetition = one algorithm completion
            # For algorithms that run multiple times per repetition, this might need adjustment
            if key not in self._completed_algorithms_per_dataset_config:
                self._completed_algorithms_per_dataset_config[key] = 0
            self._completed_algorithms_per_dataset_config[key] += 1
            self._logger.info(
                f"[WEB_DISPLAY] 🎯 Incremented algorithm completion for {key}: {self._completed_algorithms_per_dataset_config[key]}"
            )

            # Check if this completes a dataset/config combination
            total_algorithms_for_config = self._get_algorithm_count_for_config(
                self._current_config_id
            )
            repetitions = self._get_current_execution_repetitions()
            expected_total_for_combination = total_algorithms_for_config * repetitions

            if (
                self._completed_algorithms_per_dataset_config[key]
                >= expected_total_for_combination
            ):
                # This dataset/config combination is complete
                self._logger.info(
                    f"[WEB_DISPLAY] ✅ Dataset/Config combination {key} is complete!"
                )

                # Check if this config is now complete across all datasets for current execution
                current_exec = self._get_current_execution_info()
                if current_exec:
                    config_complete_across_datasets = True
                    for dataset in current_exec["datasets"]:
                        dataset_config_key = (
                            f"{dataset['id']}_{self._current_config_id}"
                        )
                        if (
                            dataset_config_key
                            not in self._completed_algorithms_per_dataset_config
                            or self._completed_algorithms_per_dataset_config[
                                dataset_config_key
                            ]
                            < expected_total_for_combination
                        ):
                            config_complete_across_datasets = False
                            break

                    if config_complete_across_datasets:
                        self._completed_configs.add(self._current_config_id)
                        self._logger.info(
                            f"[WEB_DISPLAY] ✅ Config {self._current_config_id} is complete across all datasets!"
                        )

                    # Check if this dataset is now complete across all configs for current execution
                    dataset_complete_across_configs = True
                    for config in current_exec["config_algorithms"]:
                        dataset_config_key = (
                            f"{self._current_dataset_id}_{config['id']}"
                        )
                        config_algo_count = self._get_algorithm_count_for_config(
                            config["id"]
                        )
                        expected_for_this_config = config_algo_count * repetitions
                        if (
                            dataset_config_key
                            not in self._completed_algorithms_per_dataset_config
                            or self._completed_algorithms_per_dataset_config[
                                dataset_config_key
                            ]
                            < expected_for_this_config
                        ):
                            dataset_complete_across_configs = False
                            break

                    if dataset_complete_across_configs:
                        self._completed_datasets.add(self._current_dataset_id)
                        self._logger.info(
                            f"[WEB_DISPLAY] ✅ Dataset {self._current_dataset_id} is complete across all configs!"
                        )

        if self._current_repetition_entry:
            # Update algorithm name in current repetition if it was unknown
            if (
                self._current_repetition_entry["algorithm"] == "Unknown"
                and algorithm_name != "Unknown"
            ):
                self._current_repetition_entry["algorithm"] = algorithm_name
                self._logger.info(
                    f"[WEB_DISPLAY] 🔄 Updated repetition algorithm: Unknown → {algorithm_name}"
                )

            # Add algorithm completion callback with result information
            completion_message = (
                f"FINISHED: {algorithm_name} - {'SUCCESS' if success else 'FAILED'}"
            )
            self._current_repetition_entry["callbacks"].append(completion_message)

            # Mark this repetition as completed and increment counter
            self._current_repetition_entry["status"] = (
                "completed" if success else "failed"
            )
            self._completed_combinations_count += 1
            self._completed_repetitions_count += 1  # Track actual completed repetitions
            self._execution_status["completed_repetitions_total"] = (
                self._completed_combinations_count
            )

            # Update batch progress
            if self._execution_status["total_repetitions"] > 0:
                old_progress = self._execution_status["batch_progress"]
                self._execution_status["batch_progress"] = (
                    self._completed_combinations_count
                    / self._execution_status["total_repetitions"]
                )
                self._logger.info(
                    f"[WEB_DISPLAY] 🎉 REPETITION #{self._completed_combinations_count} COMPLETED!"
                )
                self._logger.info(
                    f"[WEB_DISPLAY] 📊 Total completed repetitions: {self._completed_combinations_count}/{self._execution_status['total_repetitions']}"
                )
                self._logger.info(
                    f"[WEB_DISPLAY] 📈 Batch progress: {old_progress:.3f} → {self._execution_status['batch_progress']:.3f}"
                )
            else:
                self._logger.warning(
                    "[WEB_DISPLAY] ⚠️ Cannot calculate progress: total_repetitions is 0"
                )
        else:
            self._logger.warning(
                "[WEB_DISPLAY] ⚠️ No current repetition entry to finalize"
            )

        # Update execution status with correct algorithm name
        self._update_execution_status_progress()

    def _handle_warning(self, event) -> None:
        """Handle WarningEvent - Add warning to current repetition."""
        if self._current_repetition_entry:
            # Use the RAW warning message from the algorithm
            raw_warning = getattr(event, "warning_message", "Unknown warning")
            warning_text = f"WARNING: {raw_warning}"
            self._current_repetition_entry["callbacks"].append(warning_text)
            self._logger.warning(f"[WEB_DISPLAY] Warning: {warning_text}")

    def _handle_error(self, event) -> None:
        """Handle ErrorEvent - Add error to current repetition."""
        if self._current_repetition_entry:
            # Use the RAW error message and type from the event
            error_message = getattr(event, "error_message", "Unknown error")
            error_type = getattr(event, "error_type", "generic")
            error_text = f"ERROR ({error_type}): {error_message}"
            self._current_repetition_entry["callbacks"].append(error_text)
            self._current_repetition_entry["status"] = "failed"
            self._logger.error(f"[WEB_DISPLAY] Error: {error_text}")

    def _handle_generic(self, event) -> None:
        """Handle generic events."""
        event_type = type(event).__name__
        self._logger.info(f"[WEB_DISPLAY] 🔶 GENERIC EVENT: {event_type}")

        # Log all available attributes
        attrs = [attr for attr in dir(event) if not attr.startswith("_")]
        self._logger.info(f"[WEB_DISPLAY] 📋 Available attributes: {attrs}")

        # Log key attributes if they exist
        for attr in [
            "message",
            "algorithm_name",
            "execution_name",
            "success",
            "error_message",
        ]:
            if hasattr(event, attr):
                value = getattr(event, attr)
                self._logger.info(f"[WEB_DISPLAY] 🔍 {attr}: {value}")

        # Special handling for common progress patterns
        if hasattr(event, "message") and hasattr(event, "algorithm_name"):
            message = getattr(event, "message", "")
            algorithm_name = getattr(event, "algorithm_name", "Unknown")

            # Filter out orchestrator messages that are not real algorithm callbacks
            if message.strip() and self._current_repetition_entry:
                raw_message = message.strip()

                # Skip generic execution status messages - these are from the orchestrator
                if "Completed run" in raw_message and "/" in raw_message:
                    self._logger.debug(
                        f"[WEB_DISPLAY] 🚫 Skipping orchestrator message in generic handler: {raw_message}"
                    )
                    return

                # Only add if it's a meaningful message not already added
                if raw_message not in self._current_repetition_entry["callbacks"]:
                    self._current_repetition_entry["callbacks"].append(raw_message)
                    self._logger.info(
                        f"[WEB_DISPLAY] ✅ Added RAW callback from generic handler: {raw_message}"
                    )

                # Update algorithm name if needed
                if (
                    self._current_repetition_entry["algorithm"] == "Unknown"
                    and algorithm_name != "Unknown"
                ):
                    self._current_repetition_entry["algorithm"] = algorithm_name
                    self._logger.info(
                        f"[WEB_DISPLAY] 🔄 Updated algorithm from generic event: {algorithm_name}"
                    )

        self._logger.debug(f"[WEB_DISPLAY] 🔶 Generic event processed: {event_type}")

    def _initialize_batch_structure(self, metadata: Dict[str, Any]) -> None:
        """Initialize batch structure from task metadata."""
        try:
            self._logger.info("[WEB_DISPLAY] 🚀 INITIALIZING BATCH STRUCTURE")
            self._logger.info(f"[WEB_DISPLAY] 📥 Raw metadata: {metadata}")

            # Extract batch configuration from metadata
            batch_config = metadata.get("batch_config", {})
            self._logger.info(f"[WEB_DISPLAY] 🔧 Batch config: {batch_config}")

            # Initialize default structure
            self._batch_structure = {
                "batch": {
                    "name": "Unknown Batch",
                    "description": "No description",
                    "author": "Unknown",
                    "version": "1.0",
                    "creation_date": datetime.now().strftime("%Y-%m-%d"),
                    "tags": [],
                },
                "executions": [],
                "total_repetitions": 0,
            }

            # Extract metadata section
            if "metadata" in batch_config:
                meta = batch_config["metadata"]
                self._logger.info(f"[WEB_DISPLAY] 📝 Batch metadata found: {meta}")

                self._batch_structure["batch"].update(
                    {
                        "name": meta.get("name", "Unknown Batch"),
                        "description": meta.get("description", "No description"),
                        "author": meta.get("author", "Unknown"),
                        "version": meta.get("version", "1.0"),
                        "creation_date": meta.get(
                            "creation_date", datetime.now().strftime("%Y-%m-%d")
                        ),
                        "tags": meta.get("tags", []),
                    }
                )
                self._logger.info(
                    f"[WEB_DISPLAY] ✅ Batch metadata updated: {self._batch_structure['batch']}"
                )
            else:
                self._logger.warning(
                    "[WEB_DISPLAY] ⚠️ No metadata section in batch config"
                )

            # Process executions from batch config
            total_repetitions = 0

            if (
                "execution" in batch_config
                and "executions" in batch_config["execution"]
            ):
                executions_config = batch_config["execution"]["executions"]
                self._logger.info(
                    f"[WEB_DISPLAY] 🎯 Processing {len(executions_config)} executions"
                )

                for i, exec_config in enumerate(executions_config):
                    exec_name = exec_config.get("name", f"Execution_{i}")
                    datasets = exec_config.get("datasets", [])
                    algorithm_configs = exec_config.get("algorithms", [])
                    repetitions = exec_config.get("repetitions", 1)

                    self._logger.info(
                        f"[WEB_DISPLAY] 🎯 Processing execution '{exec_name}':"
                    )
                    self._logger.info(f"[WEB_DISPLAY]    datasets: {datasets}")
                    self._logger.info(
                        f"[WEB_DISPLAY]    algorithm_configs: {algorithm_configs}"
                    )
                    self._logger.info(f"[WEB_DISPLAY]    repetitions: {repetitions}")

                    # Calculate total repetitions for this execution
                    # Need to count actual algorithms in each config, not just config count
                    exec_total_repetitions = 0

                    # Get algorithms section to count actual algorithms per config
                    algorithms_section = batch_config.get("algorithms", [])

                    for algo_config_id in algorithm_configs:
                        # Find the algorithm config
                        algo_config = None
                        for algo_def in algorithms_section:
                            if algo_def.get("id") == algo_config_id:
                                algo_config = algo_def
                                break

                        if algo_config:
                            actual_algorithms = algo_config.get("algorithms", [])
                            num_algorithms = len(actual_algorithms)
                            self._logger.info(
                                f"[WEB_DISPLAY]    config '{algo_config_id}' has {num_algorithms} algorithms: {actual_algorithms}"
                            )
                        else:
                            # Fallback if config not found
                            num_algorithms = 1
                            self._logger.warning(
                                f"[WEB_DISPLAY]    config '{algo_config_id}' not found, assuming 1 algorithm"
                            )

                        # Add combinations for this config
                        config_total = len(datasets) * num_algorithms * repetitions
                        exec_total_repetitions += config_total
                        self._logger.info(
                            f"[WEB_DISPLAY]    config '{algo_config_id}': {len(datasets)} datasets × {num_algorithms} algorithms × {repetitions} repetitions = {config_total}"
                        )

                    total_repetitions += exec_total_repetitions

                    self._logger.info(
                        f"[WEB_DISPLAY]    execution total: {exec_total_repetitions}"
                    )

                    # Build datasets info
                    datasets_info = []
                    for dataset_id in datasets:
                        dataset_info = {
                            "id": dataset_id,
                            "name": dataset_id.replace("_", " ").title(),
                        }
                        datasets_info.append(dataset_info)
                        self._logger.info(f"[WEB_DISPLAY]    dataset: {dataset_info}")

                    # Build config_algorithms info
                    config_algorithms_info = []
                    for config_id in algorithm_configs:
                        # Find the algorithm config to get actual algorithm names
                        algo_config = None
                        for algo_def in algorithms_section:
                            if algo_def.get("id") == config_id:
                                algo_config = algo_def
                                break

                        actual_algorithms = ["Unknown"]
                        if algo_config:
                            actual_algorithms = algo_config.get(
                                "algorithms", ["Unknown"]
                            )

                        config_info = {
                            "id": config_id,
                            "name": config_id.replace("_", " ").title(),
                            "description": f"Configuration {config_id}",
                            "algorithms": actual_algorithms,
                        }
                        config_algorithms_info.append(config_info)
                        self._logger.info(
                            f"[WEB_DISPLAY]    algorithm config: {config_info}"
                        )

                    execution_info = {
                        "name": exec_name,
                        "datasets": datasets_info,
                        "config_algorithms": config_algorithms_info,
                        "repetitions": repetitions,
                        "total_repetitions": exec_total_repetitions,
                    }

                    self._batch_structure["executions"].append(execution_info)
                    self._logger.info(
                        f"[WEB_DISPLAY] ✅ Added execution: {execution_info}"
                    )

                self._logger.info(
                    f"[WEB_DISPLAY] 🔢 TOTAL REPETITIONS CALCULATED: {total_repetitions}"
                )
            else:
                self._logger.warning(
                    "[WEB_DISPLAY] ⚠️ No execution/executions in batch config"
                )

            # Update totals
            self._batch_structure["total_repetitions"] = total_repetitions
            self._execution_status["total_repetitions"] = total_repetitions

            # Update callbacks history batch name
            self._callbacks_history["batch_name"] = self._batch_structure["batch"][
                "name"
            ]

            self._logger.info(
                "[WEB_DISPLAY] 🎉 BATCH STRUCTURE INITIALIZED SUCCESSFULLY:"
            )
            self._logger.info(
                f"[WEB_DISPLAY]    Batch name: {self._batch_structure['batch']['name']}"
            )
            self._logger.info(
                f"[WEB_DISPLAY]    Total executions: {len(self._batch_structure['executions'])}"
            )
            self._logger.info(
                f"[WEB_DISPLAY]    Total repetitions: {total_repetitions}"
            )

        except Exception as e:
            self._logger.error(
                f"[WEB_DISPLAY] ❌ Error initializing batch structure: {e}",
                exc_info=True,
            )

            # Fallback structure
            self._batch_structure = {
                "batch": {
                    "name": "Unknown Batch",
                    "description": "Batch structure could not be parsed",
                    "author": "Unknown",
                    "version": "1.0",
                    "creation_date": datetime.now().strftime("%Y-%m-%d"),
                    "tags": [],
                },
                "executions": [],
                "total_repetitions": 0,
            }
            self._logger.info("[WEB_DISPLAY] 🔄 Created fallback batch structure")

    def _create_fallback_batch_structure(self) -> None:
        """Create a fallback batch structure when metadata is not available."""
        self._batch_structure = {
            "batch": {
                "name": "Batch Execution",
                "description": "Batch execution in progress",
                "author": "CSPBench",
                "version": "1.0",
                "creation_date": datetime.now().strftime("%Y-%m-%d"),
                "tags": ["execution"],
            },
            "executions": [
                {
                    "name": self._current_execution_name or "Unknown Execution",
                    "datasets": [],
                    "config_algorithms": [],
                    "repetitions": 1,
                    "total_repetitions": 0,
                }
            ],
            "total_repetitions": 0,
        }

        # Update callbacks history batch name
        self._callbacks_history["batch_name"] = self._batch_structure["batch"]["name"]

        self._logger.info("[WEB_DISPLAY] Created fallback batch structure")

    def _update_execution_status_progress(self) -> None:
        """Update execution status with current progress."""
        if not self._current_execution_name:
            return

        # Find the current execution in batch structure
        current_exec_info = None
        if self._batch_structure:
            for exec_info in self._batch_structure["executions"]:
                if exec_info["name"] == self._current_execution_name:
                    current_exec_info = exec_info
                    break

        if current_exec_info:
            # Count datasets and configs
            total_datasets = len(current_exec_info["datasets"])
            total_configs = len(current_exec_info["config_algorithms"])

            # Calculate completed counts more intelligently
            dataset_completed = self._calculate_completed_datasets(current_exec_info)
            config_completed = self._calculate_completed_configs(current_exec_info)

            # Calculate algorithm counts for current config
            total_algorithms_in_config = self._get_algorithm_count_for_config(
                self._current_config_id
            )
            completed_algorithms_in_config = (
                self._calculate_completed_algorithms_in_config()
            )

            # Calculate completed repetitions for current context
            completed_repetitions = self._calculate_completed_repetitions_for_context()
            total_repetitions_for_context = current_exec_info.get("repetitions", 1)

            # Update execution status
            self._execution_status["current_execution"] = {
                "name": self._current_execution_name,
                "dataset": {
                    "id": self._current_dataset_id or "unknown",
                    "total": total_datasets,
                    "completed": dataset_completed,
                },
                "config_id": {
                    "id": self._current_config_id or "unknown",
                    "total": total_configs,
                    "completed": config_completed,
                },
                "algorithm": {
                    "name": self._current_algorithm or "unknown",
                    "total": total_algorithms_in_config,
                    "completed": completed_algorithms_in_config,
                    "repetitions": {
                        "total": total_repetitions_for_context,
                        "completed": completed_repetitions,
                    },
                },
            }

    def _get_algorithm_count_for_config(self, config_id: Optional[str]) -> int:
        """Get the number of algorithms in a specific configuration."""
        if not config_id or not self._batch_structure:
            return 1

        # Find the config in executions -> config_algorithms section
        if "executions" in self._batch_structure:
            for execution in self._batch_structure["executions"]:
                if execution.get("name") == self._current_execution_name:
                    for config_algo in execution.get("config_algorithms", []):
                        if config_algo.get("id") == config_id:
                            return len(config_algo.get("algorithms", []))

        return 1  # fallback

    def _get_algorithm_completed_count(self) -> int:
        """Get the number of algorithms completed in current dataset/config combination."""
        if not self._current_dataset_id or not self._current_config_id:
            return 0

        # Create key for current dataset/config combination
        key = f"{self._current_dataset_id}_{self._current_config_id}"
        return self._completed_algorithms_per_dataset_config.get(key, 0)

    def _get_current_execution_info(self) -> Optional[dict]:
        """Get the current execution information from batch structure."""
        if not self._current_execution_name or not self._batch_structure:
            return None

        for execution in self._batch_structure.get("executions", []):
            if execution["name"] == self._current_execution_name:
                return execution
        return None

    def _get_current_execution_repetitions(self) -> int:
        """Get the number of repetitions for current execution."""
        current_exec = self._get_current_execution_info()
        if current_exec:
            return current_exec.get("repetitions", 1)
        return 1

    def _calculate_completed_datasets(self, current_exec_info: dict) -> int:
        """Calculate how many datasets have been completed in current execution."""
        if not current_exec_info or not self._current_dataset_id:
            return 0

        completed_count = 0
        repetitions = current_exec_info.get("repetitions", 1)

        self._logger.info(
            f"[WEB_DISPLAY] 🔍 Calculating completed datasets for execution: {current_exec_info.get('name', 'unknown')}"
        )
        self._logger.info(
            f"[WEB_DISPLAY] 🔍 Available completion data: {self._completed_algorithms_per_dataset_config}"
        )

        for dataset in current_exec_info.get("datasets", []):
            dataset_id = dataset.get("id", "")
            dataset_is_complete = True

            self._logger.info(f"[WEB_DISPLAY] 🔍 Checking dataset: {dataset_id}")

            # Check if this dataset is complete across all configs for THIS execution
            for config in current_exec_info.get("config_algorithms", []):
                config_id = config.get("id", "")
                key = f"{dataset_id}_{config_id}"

                total_algorithms_for_config = len(config.get("algorithms", []))
                expected_total = total_algorithms_for_config * repetitions
                actual_completed = self._completed_algorithms_per_dataset_config.get(
                    key, 0
                )

                self._logger.info(
                    f"[WEB_DISPLAY] 🔍   Config {config_id}: {actual_completed}/{expected_total} (key: {key})"
                )

                if actual_completed < expected_total:
                    dataset_is_complete = False
                    break

            if dataset_is_complete:
                completed_count += 1
                self._logger.info(f"[WEB_DISPLAY] ✅ Dataset {dataset_id} is complete!")
            else:
                self._logger.info(
                    f"[WEB_DISPLAY] ❌ Dataset {dataset_id} is incomplete"
                )

        self._logger.info(
            f"[WEB_DISPLAY] 📊 Total completed datasets: {completed_count}/{len(current_exec_info.get('datasets', []))}"
        )
        return completed_count

    def _calculate_completed_configs(self, current_exec_info: dict) -> int:
        """Calculate how many configs have been completed in current execution."""
        if not current_exec_info or not self._current_config_id:
            return 0

        completed_count = 0
        repetitions = current_exec_info.get("repetitions", 1)

        for config in current_exec_info.get("config_algorithms", []):
            config_id = config.get("id", "")
            config_is_complete = True

            # Check if this config is complete across all datasets for THIS execution
            for dataset in current_exec_info.get("datasets", []):
                dataset_id = dataset.get("id", "")
                key = f"{dataset_id}_{config_id}"

                total_algorithms_for_config = len(config.get("algorithms", []))
                expected_total = total_algorithms_for_config * repetitions
                actual_completed = self._completed_algorithms_per_dataset_config.get(
                    key, 0
                )

                if actual_completed < expected_total:
                    config_is_complete = False
                    break

            if config_is_complete:
                completed_count += 1

        return completed_count

    def _calculate_completed_algorithms_in_config(self) -> int:
        """Calculate how many unique algorithms have been completed in current config."""
        if not self._current_config_id or not self._current_dataset_id:
            return 0

        # Get algorithms for current config
        total_algorithms_in_config = self._get_algorithm_count_for_config(
            self._current_config_id
        )
        if total_algorithms_in_config == 0:
            return 0

        # Get total executions for current dataset/config combination
        key = f"{self._current_dataset_id}_{self._current_config_id}"
        total_executions = self._completed_algorithms_per_dataset_config.get(key, 0)

        # Since each algorithm runs once per repetition, divide by repetitions to get unique algorithms
        current_exec = self._get_current_execution_info()
        repetitions = current_exec.get("repetitions", 1) if current_exec else 1

        if repetitions > 0:
            completed_unique_algorithms = min(
                total_executions // repetitions, total_algorithms_in_config
            )
            return completed_unique_algorithms

        return 0

    def _calculate_completed_repetitions_for_context(self) -> int:
        """Calculate completed repetitions for current algorithm in current dataset/config."""
        if not self._current_config_id or not self._current_dataset_id:
            return 0

        # Get total algorithm executions for current dataset/config combination
        key = f"{self._current_dataset_id}_{self._current_config_id}"
        total_executions = self._completed_algorithms_per_dataset_config.get(key, 0)

        # Get total algorithms in current config
        total_algorithms_in_config = self._get_algorithm_count_for_config(
            self._current_config_id
        )

        if total_algorithms_in_config > 0:
            # Each completed repetition means all algorithms in config were executed once
            completed_repetitions = total_executions // total_algorithms_in_config
            return completed_repetitions

        return 0

    def _get_dataset_id_from_name(self, dataset_name: str) -> str:
        """Convert dataset name to dataset ID."""
        if not self._batch_structure:
            return dataset_name.lower()

        # First check in all executions' datasets
        for execution in self._batch_structure.get("executions", []):
            for dataset in execution.get("datasets", []):
                if dataset.get("name", "").lower() == dataset_name.lower():
                    return dataset.get("id", dataset_name.lower())

        # Fallback: convert name to ID format
        return dataset_name.lower().replace(" ", "_")

    def _get_config_id_from_name(self, config_name: str) -> str:
        """Convert config name to config ID."""
        if not self._batch_structure:
            return config_name.lower()

        # First check in all executions' config_algorithms
        for execution in self._batch_structure.get("executions", []):
            for config in execution.get("config_algorithms", []):
                if config.get("name", "").lower() == config_name.lower():
                    return config.get("id", config_name.lower())

        # Fallback: convert name to ID format
        return config_name.lower().replace(" ", "_")

    def _start_new_repetition(self) -> None:
        """Start a new repetition entry in callbacks history."""
        self._logger.info(
            f"[WEB_DISPLAY] 🚀 STARTING NEW REPETITION {self._current_repetition}"
        )
        self._logger.info("[WEB_DISPLAY] 🎯 Current context:")
        self._logger.info(
            f"[WEB_DISPLAY]    execution_name: {self._current_execution_name}"
        )
        self._logger.info(f"[WEB_DISPLAY]    dataset_id: {self._current_dataset_id}")
        self._logger.info(f"[WEB_DISPLAY]    config_id: {self._current_config_id}")
        self._logger.info(f"[WEB_DISPLAY]    algorithm: {self._current_algorithm}")
        self._logger.info(f"[WEB_DISPLAY]    repetition: {self._current_repetition}")

        # Finalize previous repetition if it exists and is still running
        if (
            self._current_repetition_entry
            and self._current_repetition_entry["status"] == "running"
        ):
            self._current_repetition_entry["status"] = "completed"
            self._logger.info(
                f"[WEB_DISPLAY] 🏁 Finalized previous repetition with {len(self._current_repetition_entry['callbacks'])} callbacks"
            )

        # Create new repetition entry
        self._current_repetition_entry = {
            "execution_name": self._current_execution_name or "unknown",
            "dataset_id": self._current_dataset_id or "unknown",
            "config_id": self._current_config_id or "unknown",
            "algorithm": self._current_algorithm or "unknown",
            "repetition": self._current_repetition,
            "status": "running",
            "callbacks": ["start"],
        }

        self._callbacks_history["history"].append(self._current_repetition_entry)
        self._logger.info(
            f"[WEB_DISPLAY] ✅ Created repetition entry: {self._current_repetition_entry}"
        )
        self._logger.info(
            f"[WEB_DISPLAY] 📊 Total repetitions in history: {len(self._callbacks_history['history'])}"
        )

    def _update_web_session(self) -> None:
        """Update web session with current data structures."""
        try:
            # Prepare structured data for web interface
            session_data = {
                "batch_structure": self._batch_structure,
                "execution_status": self._execution_status,
                "callbacks_history": self._callbacks_history,
                "last_updated": datetime.now().isoformat(),
            }

            self._logger.debug("[WEB_DISPLAY] 🔄 UPDATING WEB SESSION:")
            self._logger.debug(
                f"[WEB_DISPLAY]    batch_structure: {self._batch_structure is not None}"
            )
            if self._batch_structure:
                self._logger.debug(
                    f"[WEB_DISPLAY]    batch_name: {self._batch_structure.get('batch', {}).get('name', 'None')}"
                )
                self._logger.debug(
                    f"[WEB_DISPLAY]    total_repetitions: {self._batch_structure.get('total_repetitions', 0)}"
                )
            self._logger.debug(
                f"[WEB_DISPLAY]    execution_status.total_repetitions: {self._execution_status.get('total_repetitions', 0)}"
            )
            self._logger.debug(
                f"[WEB_DISPLAY]    execution_status.completed_repetitions_total: {self._execution_status.get('completed_repetitions_total', 0)}"
            )
            self._logger.debug(
                f"[WEB_DISPLAY]    execution_status.batch_progress: {self._execution_status.get('batch_progress', 0.0)}"
            )
            self._logger.debug(
                f"[WEB_DISPLAY]    callbacks_history entries: {len(self._callbacks_history.get('history', []))}"
            )

            # Update session with structured data
            self._session_manager.update_session(
                self._session_id,
                {
                    "progress_state": session_data,
                    "last_updated": session_data["last_updated"],
                },
            )

            self._logger.debug("[WEB_DISPLAY] ✅ Web session updated successfully")

        except Exception as e:
            self._logger.error(
                f"[WEB_DISPLAY] ❌ Error updating web session: {e}", exc_info=True
            )

    def get_batch_structure(self) -> Optional[Dict[str, Any]]:
        """Get the current batch structure."""
        return self._batch_structure

    def get_execution_status(self) -> Dict[str, Any]:
        """Get the current execution status."""
        return self._execution_status

    def get_callbacks_history(self) -> Dict[str, Any]:
        """Get the current callbacks history."""
        return self._callbacks_history
