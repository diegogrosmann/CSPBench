"""Advanced database queries mixin for complex operations.

This module provides sophisticated query operations for monitoring, reporting,
and analytics on work execution data. It includes dataclasses for structured
results and methods for complex multi-table queries.
"""

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple, Union

from sqlalchemy import case, desc, func


@dataclass
class ExecutionDetail:
    """Detailed execution information for monitoring and reporting.

    This dataclass contains comprehensive information about a single execution,
    including progress tracking, timing information, and hierarchical context
    (task, dataset, preset, algorithm).

    Attributes:
        unit_id: Unique identifier for the execution unit
        combination_id: ID of the combination this execution belongs to
        sequencia: Sequence number within the combination
        status: Current execution status
        progress: Progress percentage (0.0 to 1.0)
        progress_message: Human-readable progress description
        started_at: Timestamp when execution started
        finished_at: Timestamp when execution finished
        objective: Objective function value (if completed)
        task_id: Task identifier
        dataset_id: Dataset identifier
        preset_id: Configuration preset identifier
        algorithm_id: Algorithm identifier
        mode: Execution mode
        total_sequences: Total number of sequences in the combination
    """

    unit_id: str
    combination_id: int
    sequencia: int
    status: str
    progress: float
    progress_message: Optional[str]
    started_at: Optional[float]
    finished_at: Optional[float]
    objective: Optional[float]
    task_id: str
    dataset_id: str
    preset_id: str
    algorithm_id: str
    mode: str
    total_sequences: int


@dataclass
class ProgressSummary:
    """Comprehensive progress summary for a work item.

    Provides hierarchical progress information across tasks, datasets,
    configurations, algorithms, and executions for monitoring purposes.

    Attributes:
        work_id: Work identifier this summary belongs to
        tasks: Task-level progress information
        datasets: Dataset-level progress information
        configs: Configuration-level progress information
        algorithms: Algorithm-level progress information
        execution: Current combination execution progress
        global_execution: Overall execution statistics
        global_progress: Global progress percentage (0.0 to 1.0)
        current_combination_details: Details of currently active combination
    """

    work_id: str
    tasks: Dict[str, Union[int, str]]
    datasets: Dict[str, Union[int, str]]
    configs: Dict[str, Union[int, str]]
    algorithms: Dict[str, Union[int, str]]
    execution: Dict[str, int]
    global_execution: Dict[str, int]
    global_progress: float
    current_combination_details: Optional[Dict[str, Any]] = None


@dataclass
class ErrorSummary:
    """Summary of execution errors for troubleshooting.

    Attributes:
        unit_id: Execution unit where the error occurred
        error_type: Type/category of the error
        error_message: Detailed error message
        timestamp: When the error occurred
    """

    unit_id: str
    error_type: str
    error_message: str
    timestamp: float


class QueriesMixin:
    """Advanced database queries for complex operations and analytics.

    This mixin provides sophisticated query methods for monitoring work progress,
    analyzing execution statistics, retrieving detailed execution information,
    and supporting data export operations. It complements the basic CRUD mixins
    with complex multi-table queries and aggregations.
    """

    def get_running_combination(self, work_id: str) -> Optional[Dict[str, Any]]:
        """Get the currently running combination for a work.

        Args:
            work_id: Work identifier to search for running combinations

        Returns:
            Dictionary with combination details if a running combination exists,
            None otherwise. Contains combination_id, task_id, dataset_id,
            preset_id, algorithm_id, and total_sequences.
        """
        from ..models import Combination

        with self.session_scope() as session:
            combination = (
                session.query(Combination)
                .filter(Combination.work_id == work_id, Combination.status == "running")
                .order_by(
                    Combination.task_id,
                    Combination.dataset_id,
                    Combination.preset_id,
                    Combination.algorithm_id,
                )
                .first()
            )

            if combination:
                return {
                    "id": combination.id,
                    "combination_id": combination.id,
                    "task_id": combination.task_id,
                    "dataset_id": combination.dataset_id,
                    "preset_id": combination.preset_id,
                    "algorithm_id": combination.algorithm_id,
                    "total_sequences": combination.total_sequences,
                }
            return None

    def get_last_paused_or_incomplete_combination(
        self, work_id: str
    ) -> Optional[Dict[str, Any]]:
        """Get the most recent paused or incomplete combination for context.

        When no running combination exists, this method helps provide context
        by finding the last combination that was paused, queued, failed, or canceled.
        Prioritizes paused combinations, then queued, then failed/error, then canceled.

        Args:
            work_id: Work identifier to search combinations for

        Returns:
            Dictionary with combination details if found, None otherwise.
            Same structure as get_running_combination.
        """
        from ..models import Combination

        with self.session_scope() as session:
            # Try paused first
            combination = (
                session.query(Combination)
                .filter(Combination.work_id == work_id, Combination.status == "paused")
                .order_by(desc(Combination.id))
                .first()
            )

            if combination:
                return self._combination_to_dict(combination)

            # Try queued
            combination = (
                session.query(Combination)
                .filter(Combination.work_id == work_id, Combination.status == "queued")
                .order_by(Combination.id)
                .first()
            )

            if combination:
                return self._combination_to_dict(combination)

            # Try failed/error
            combination = (
                session.query(Combination)
                .filter(
                    Combination.work_id == work_id,
                    Combination.status.in_(["failed", "error"]),
                )
                .order_by(desc(Combination.id))
                .first()
            )

            if combination:
                return self._combination_to_dict(combination)

            # Try cancelled
            combination = (
                session.query(Combination)
                .filter(
                    Combination.work_id == work_id, Combination.status == "canceled"
                )
                .order_by(desc(Combination.id))
                .first()
            )

            if combination:
                return self._combination_to_dict(combination)

            return None

    def _combination_to_dict(self, combination) -> Dict[str, Any]:
        """Convert combination model to dictionary format.

        Internal helper method to standardize combination data format
        across different query methods.

        Args:
            combination: SQLAlchemy combination model instance

        Returns:
            Dictionary with standardized combination fields
        """
        return {
            "id": combination.id,
            "combination_id": combination.id,
            "task_id": combination.task_id,
            "dataset_id": combination.dataset_id,
            "preset_id": combination.preset_id,
            "algorithm_id": combination.algorithm_id,
            "total_sequences": combination.total_sequences,
        }

    def get_task_status_lists(
        self, work_id: str
    ) -> Tuple[List[str], List[str], List[str]]:
        """Get task status categorization for a work.

        Analyzes all combinations within a work to categorize tasks by their
        completion status based on their constituent combinations.

        Args:
            work_id: Work identifier to analyze

        Returns:
            Tuple containing three lists:
                - finished_tasks: Task IDs where all combinations are completed
                - running_tasks: Task IDs with mixed completion status
                - queued_tasks: Task IDs where all combinations are queued
        """
        from ..models import Combination

        with self.session_scope() as session:
            # Get task status aggregation
            result = (
                session.query(
                    Combination.task_id,
                    func.sum(
                        case(
                            (
                                Combination.status.in_(
                                    ["completed", "failed", "error"]
                                ),
                                1,
                            ),
                            else_=0,
                        )
                    ).label("completed_cnt"),
                    func.sum(case((Combination.status == "queued", 1), else_=0)).label(
                        "queued_cnt"
                    ),
                    func.count().label("total_cnt"),
                )
                .filter(Combination.work_id == work_id)
                .group_by(Combination.task_id)
                .all()
            )

            finished, queued, running = [], [], []
            for row in result:
                if row.completed_cnt == row.total_cnt:
                    finished.append(row.task_id)
                elif row.queued_cnt == row.total_cnt:
                    queued.append(row.task_id)
                else:
                    running.append(row.task_id)

            return finished, running, queued

    def get_dataset_status_lists(
        self, work_id: str, task_id: str
    ) -> Tuple[List[str], List[str], List[str]]:
        """Get dataset status categorization within a specific task.

        Args:
            work_id: Work identifier
            task_id: Task identifier to analyze datasets for

        Returns:
            Tuple containing three lists:
                - finished_datasets: Dataset IDs with all combinations completed
                - running_datasets: Dataset IDs with mixed completion status
                - queued_datasets: Dataset IDs with all combinations queued
        """
        from ..models import Combination

        with self.session_scope() as session:
            result = (
                session.query(
                    Combination.dataset_id,
                    func.sum(
                        case(
                            (
                                Combination.status.in_(
                                    ["completed", "failed", "error"]
                                ),
                                1,
                            ),
                            else_=0,
                        )
                    ).label("completed_cnt"),
                    func.sum(case((Combination.status == "queued", 1), else_=0)).label(
                        "queued_cnt"
                    ),
                    func.count().label("total_cnt"),
                )
                .filter(Combination.work_id == work_id, Combination.task_id == task_id)
                .group_by(Combination.dataset_id)
                .all()
            )

            finished, queued, running = [], [], []
            for row in result:
                if row.completed_cnt == row.total_cnt:
                    finished.append(row.dataset_id)
                elif row.queued_cnt == row.total_cnt:
                    queued.append(row.dataset_id)
                else:
                    running.append(row.dataset_id)

            return finished, running, queued

    def get_config_status_lists(
        self, work_id: str, task_id: str, dataset_id: str
    ) -> Tuple[List[str], List[str], List[str]]:
        """Get configuration preset status categorization within a dataset.

        Args:
            work_id: Work identifier
            task_id: Task identifier
            dataset_id: Dataset identifier to analyze configurations for

        Returns:
            Tuple containing three lists:
                - finished_configs: Preset IDs with all combinations completed
                - running_configs: Preset IDs with mixed completion status
                - queued_configs: Preset IDs with all combinations queued
        """
        from ..models import Combination

        with self.session_scope() as session:
            result = (
                session.query(
                    Combination.preset_id,
                    func.sum(
                        case(
                            (
                                Combination.status.in_(
                                    ["completed", "failed", "error"]
                                ),
                                1,
                            ),
                            else_=0,
                        )
                    ).label("completed_cnt"),
                    func.sum(case((Combination.status == "queued", 1), else_=0)).label(
                        "queued_cnt"
                    ),
                    func.count().label("total_cnt"),
                )
                .filter(
                    Combination.work_id == work_id,
                    Combination.task_id == task_id,
                    Combination.dataset_id == dataset_id,
                )
                .group_by(Combination.preset_id)
                .all()
            )

            finished, queued, running = [], [], []
            for row in result:
                if row.completed_cnt == row.total_cnt:
                    finished.append(row.preset_id)
                elif row.queued_cnt == row.total_cnt:
                    queued.append(row.preset_id)
                else:
                    running.append(row.preset_id)

            return finished, running, queued

    def get_algorithm_status_lists(
        self, work_id: str, task_id: str, dataset_id: str, preset_id: str
    ) -> Tuple[List[str], List[str], List[str]]:
        """Get algorithm status categorization within a configuration.

        Args:
            work_id: Work identifier
            task_id: Task identifier
            dataset_id: Dataset identifier
            preset_id: Configuration preset identifier

        Returns:
            Tuple containing three lists:
                - finished_algorithms: Algorithm IDs with all combinations completed
                - running_algorithms: Algorithm IDs with mixed completion status
                - queued_algorithms: Algorithm IDs with all combinations queued
        """
        from ..models import Combination

        with self.session_scope() as session:
            result = (
                session.query(
                    Combination.algorithm_id,
                    func.sum(
                        case(
                            (
                                Combination.status.in_(
                                    ["completed", "failed", "error"]
                                ),
                                1,
                            ),
                            else_=0,
                        )
                    ).label("completed_cnt"),
                    func.sum(case((Combination.status == "queued", 1), else_=0)).label(
                        "queued_cnt"
                    ),
                    func.count().label("total_cnt"),
                )
                .filter(
                    Combination.work_id == work_id,
                    Combination.task_id == task_id,
                    Combination.dataset_id == dataset_id,
                    Combination.preset_id == preset_id,
                )
                .group_by(Combination.algorithm_id)
                .all()
            )

            finished, queued, running = [], [], []
            for row in result:
                if row.completed_cnt == row.total_cnt:
                    finished.append(row.algorithm_id)
                elif row.queued_cnt == row.total_cnt:
                    queued.append(row.algorithm_id)
                else:
                    running.append(row.algorithm_id)

            return finished, running, queued

    def get_execution_status_lists(
        self, combination_id: int
    ) -> Tuple[List[int], List[int], List[int]]:
        """Get execution status categorization for a specific combination.

        Args:
            combination_id: Combination identifier to analyze executions for

        Returns:
            Tuple containing three lists:
                - finished_executions: Execution IDs that are completed/failed/error
                - running_executions: Execution IDs that are currently running
                - queued_executions: Execution IDs that are queued
        """
        from ..models import Execution

        with self.session_scope() as session:
            executions = (
                session.query(Execution)
                .filter(Execution.combination_id == combination_id)
                .all()
            )

            finished, running, queued = [], [], []
            for execution in executions:
                if execution.status in ("completed", "failed", "error"):
                    finished.append(execution.id)
                elif execution.status == "running":
                    running.append(execution.id)
                elif execution.status == "queued":
                    queued.append(execution.id)

            return finished, running, queued

    def get_global_execution_stats(self, work_id: str) -> Dict[str, Union[int, float]]:
        """Get global execution statistics across all combinations in a work.

        Args:
            work_id: Work identifier to calculate statistics for

        Returns:
            Dictionary containing:
                - Finished: Number of completed executions
                - Total: Total number of sequences across all combinations
                - Progress: Overall progress ratio (0.0 to 1.0)
        """
        from ..models import Combination, Execution

        with self.session_scope() as session:
            # Total sequences from combinations
            total_sequences = (
                session.query(func.coalesce(func.sum(Combination.total_sequences), 0))
                .filter(Combination.work_id == work_id)
                .scalar()
                or 0
            )

            # Finished executions
            finished = (
                session.query(func.count(Execution.id))
                .join(Combination, Execution.combination_id == Combination.id)
                .filter(
                    Combination.work_id == work_id,
                    Execution.status.in_(["completed", "failed", "error"]),
                )
                .scalar()
                or 0
            )

            progress = finished / total_sequences if total_sequences else 0.0

            return {
                "Finished": finished,
                "Total": total_sequences,
                "Progress": progress,
            }

    def get_work_progress_summary(self, work_id: str) -> Optional[ProgressSummary]:
        """Get comprehensive progress summary for a work item.

        Builds a complete hierarchical progress summary including status
        information at task, dataset, configuration, algorithm, and execution
        levels, plus global statistics and current combination context.

        Args:
            work_id: Work identifier to generate summary for

        Returns:
            ProgressSummary dataclass instance with complete progress information,
            or None if the work doesn't exist.
        """
        work = self.work_get(work_id)
        if not work:
            return None

        current_combo = self.get_running_combination(work_id)
        if not current_combo:
            current_combo = self.get_last_paused_or_incomplete_combination(work_id)

        # TASKS
        task_finished, task_running, task_queued = self.get_task_status_lists(work_id)
        running_task_id = current_combo["task_id"] if current_combo else ""

        # DATASETS (only for current task)
        if current_combo:
            dataset_finished, dataset_running, dataset_queued = (
                self.get_dataset_status_lists(work_id, current_combo["task_id"])
            )
            running_dataset_id = current_combo["dataset_id"]
        else:
            dataset_finished = dataset_running = dataset_queued = []
            running_dataset_id = ""

        # CONFIGS
        if current_combo and current_combo["dataset_id"]:
            config_finished, config_running, config_queued = (
                self.get_config_status_lists(
                    work_id, current_combo["task_id"], current_combo["dataset_id"]
                )
            )
            running_config_id = current_combo["preset_id"]
        else:
            config_finished = config_running = config_queued = []
            running_config_id = ""

        # ALGORITHMS
        if current_combo and current_combo["preset_id"]:
            alg_finished, alg_running, alg_queued = self.get_algorithm_status_lists(
                work_id,
                current_combo["task_id"],
                current_combo["dataset_id"],
                current_combo["preset_id"],
            )
            running_alg_id = current_combo["algorithm_id"]
        else:
            alg_finished = alg_running = alg_queued = []
            running_alg_id = ""

        # EXECUTIONS
        if current_combo:
            exec_finished_list, exec_running_list, exec_queued_list = (
                self.get_execution_status_lists(current_combo["combination_id"])
            )
            execution_total = current_combo.get("total_sequences") or 0
        else:
            exec_finished_list = exec_running_list = exec_queued_list = []
            execution_total = 0

        # GLOBAL EXECUTIONS
        global_stats = self.get_global_execution_stats(work_id)
        global_finished = global_stats["Finished"]
        global_total = global_stats["Total"]
        global_progress = global_stats["Progress"]

        return ProgressSummary(
            work_id=work_id,
            tasks={
                "Finished": len(task_finished),
                "Running": running_task_id,
                "Queued": len(task_queued),
                "Total": len(task_finished) + len(task_running) + len(task_queued),
            },
            datasets={
                "Finished": len(dataset_finished),
                "Running": running_dataset_id,
                "Queued": len(dataset_queued),
                "Total": len(dataset_finished)
                + len(dataset_running)
                + len(dataset_queued),
            },
            configs={
                "Finished": len(config_finished),
                "Running": running_config_id,
                "Queued": len(config_queued),
                "Total": len(config_finished)
                + len(config_running)
                + len(config_queued),
            },
            algorithms={
                "Finished": len(alg_finished),
                "Running": running_alg_id,
                "Queued": len(alg_queued),
                "Total": len(alg_finished) + len(alg_running) + len(alg_queued),
            },
            execution={
                "Finished": len(exec_finished_list),
                "Running": len(exec_running_list),
                "Queued": len(exec_queued_list),
                "Total": execution_total,
            },
            global_execution={
                "Finished": global_finished,
                "Total": global_total,
            },
            global_progress=global_progress,
            current_combination_details=dict(current_combo) if current_combo else None,
        )

    def get_combination_executions_detail(
        self, combination_id: int
    ) -> List[ExecutionDetail]:
        """Get detailed information about all executions in a combination.

        Retrieves comprehensive execution details including progress tracking,
        timing information, and results for monitoring and analysis purposes.

        Args:
            combination_id: Combination identifier to get execution details for

        Returns:
            List of ExecutionDetail dataclass instances, ordered by sequence number.
        """
        from ..models import Combination, Execution, ExecutionProgress

        with self.session_scope() as session:
            # Query executions with combination data
            query = (
                session.query(Execution, Combination)
                .join(Combination, Execution.combination_id == Combination.id)
                .filter(Execution.combination_id == combination_id)
                .order_by(Execution.sequencia)
            )

            results = []
            for execution, combination in query.all():
                # Get latest progress
                latest_progress = (
                    session.query(ExecutionProgress)
                    .filter(ExecutionProgress.execution_id == execution.id)
                    .order_by(desc(ExecutionProgress.timestamp))
                    .first()
                )

                progress = latest_progress.progress if latest_progress else 0.0
                progress_message = latest_progress.message if latest_progress else None

                # Format progress message based on status
                if execution.status == "completed":
                    progress_message = f"Objective: {execution.objective}"
                elif execution.status in ("failed", "error") and execution.result_json:
                    progress_message = execution.result_json.get(
                        "error_message", "Unknown error"
                    )

                results.append(
                    ExecutionDetail(
                        unit_id=execution.unit_id,
                        combination_id=execution.combination_id,
                        sequencia=execution.sequencia,
                        status=execution.status,
                        progress=progress,
                        progress_message=progress_message,
                        started_at=execution.started_at,
                        finished_at=execution.finished_at,
                        objective=execution.objective,
                        task_id=combination.task_id,
                        dataset_id=combination.dataset_id,
                        preset_id=combination.preset_id,
                        algorithm_id=combination.algorithm_id,
                        mode=combination.mode,
                        total_sequences=combination.total_sequences,
                    )
                )

        return results

    def get_running_executions_detail(
        self, work_id: str, limit: int = 20
    ) -> List[ExecutionDetail]:
        """Get detailed information about currently running executions.

        Args:
            work_id: Work identifier to find running executions for
            limit: Maximum number of executions to return

        Returns:
            List of ExecutionDetail instances for running executions,
            ordered by start time (most recent first).
        """
        from ..models import Combination, Execution, ExecutionProgress

        with self.session_scope() as session:
            query = (
                session.query(Execution, Combination)
                .join(Combination, Execution.combination_id == Combination.id)
                .filter(Combination.work_id == work_id, Execution.status == "running")
                .order_by(desc(Execution.started_at))
                .limit(limit)
            )

            results = []
            for execution, combination in query.all():
                # Get latest progress
                latest_progress = (
                    session.query(ExecutionProgress)
                    .filter(ExecutionProgress.execution_id == execution.id)
                    .order_by(desc(ExecutionProgress.timestamp))
                    .first()
                )

                progress = latest_progress.progress if latest_progress else 0.0
                progress_message = latest_progress.message if latest_progress else None

                results.append(
                    ExecutionDetail(
                        unit_id=execution.unit_id,
                        combination_id=execution.combination_id,
                        sequencia=execution.sequencia or 0,
                        status=execution.status,
                        progress=progress,
                        progress_message=progress_message,
                        started_at=execution.started_at,
                        finished_at=execution.finished_at,
                        objective=execution.objective,
                        task_id=combination.task_id,
                        dataset_id=combination.dataset_id,
                        preset_id=combination.preset_id,
                        algorithm_id=combination.algorithm_id,
                        mode=combination.mode,
                        total_sequences=combination.total_sequences or 0,
                    )
                )

        return results

    def get_error_summary(self, work_id: str, limit: int = 10) -> List[ErrorSummary]:
        """Get summary of recent errors for troubleshooting.

        Args:
            work_id: Work identifier to get error summary for
            limit: Maximum number of error entries to return

        Returns:
            List of ErrorSummary instances for failed/error executions,
            ordered by occurrence time (most recent first).
        """
        from ..models import Combination, Execution

        with self.session_scope() as session:
            query = (
                session.query(Execution)
                .join(Combination, Execution.combination_id == Combination.id)
                .filter(
                    Combination.work_id == work_id,
                    Execution.status.in_(["failed", "error"]),
                )
                .order_by(desc(Execution.finished_at))
                .limit(limit)
            )

            results = []
            for execution in query.all():
                error_message = "Unknown error"
                if execution.result_json and "error_message" in execution.result_json:
                    error_message = execution.result_json["error_message"]

                results.append(
                    ErrorSummary(
                        unit_id=execution.unit_id,
                        error_type=execution.status,
                        error_message=error_message,
                        timestamp=execution.finished_at or 0.0,
                    )
                )

        return results

    def get_execution_warnings(
        self, work_id: str, limit: int = 10
    ) -> List[Dict[str, Any]]:
        """Get recent warning events from the events table.

        Args:
            work_id: Work identifier to get warnings for
            limit: Maximum number of warning entries to return

        Returns:
            List of dictionaries containing warning event details,
            ordered by timestamp (most recent first).
        """
        from ..models import Event

        with self.session_scope() as session:
            query = (
                session.query(Event)
                .filter(Event.work_id == work_id, Event.event_type == "warning")
                .order_by(desc(Event.timestamp))
                .limit(limit)
            )

            results = []
            for event in query.all():
                entity_data = event.entity_data_json or {}
                results.append(
                    {
                        "event_type": event.event_type,
                        "event_category": event.event_category,
                        "message": entity_data.get("message", "Unknown warning"),
                        "unit_id": entity_data.get("unit_id"),
                        "combination_id": entity_data.get("combination_id"),
                        "timestamp": event.timestamp,
                    }
                )

        return results

    def get_events(
        self, work_id: str, limit: int = 100, event_types: Optional[List[str]] = None
    ) -> List[Dict[str, Any]]:
        """Get events for a work item with optional filtering by type.

        Args:
            work_id: Work identifier to get events for
            limit: Maximum number of events to return
            event_types: Optional list of event types to filter by

        Returns:
            List of dictionaries containing event details,
            ordered by timestamp (most recent first).
        """
        from ..models import Event

        results = []
        with self.session_scope() as session:
            query = session.query(Event).filter(Event.work_id == work_id)

            if event_types:
                query = query.filter(Event.event_type.in_(event_types))

            query = query.order_by(desc(Event.timestamp)).limit(limit)

            for event in query.all():
                entity_data = event.entity_data_json or {}
                results.append(
                    {
                        "id": event.id,
                        "work_id": event.work_id,
                        "event_type": event.event_type,
                        "event_category": event.event_category,
                        "entity_data": entity_data,
                        "timestamp": event.timestamp,
                        "message": entity_data.get("message", ""),
                    }
                )

        return results

    def get_combination_status_counts(self, work_id: str) -> Dict[str, int]:
        """Get count statistics of combinations grouped by status.

        Args:
            work_id: Work identifier to get statistics for

        Returns:
            Dictionary mapping status names to their counts.
        """
        from ..models import Combination

        with self.session_scope() as session:
            result = (
                session.query(Combination.status, func.count().label("count"))
                .filter(Combination.work_id == work_id)
                .group_by(Combination.status)
                .all()
            )

            return {row.status: row.count for row in result}

    def get_execution_status_counts(self, work_id: str) -> Dict[str, int]:
        """Get count statistics of executions grouped by status.

        Args:
            work_id: Work identifier to get statistics for

        Returns:
            Dictionary mapping status names to their counts.
        """
        from ..models import Combination, Execution

        with self.session_scope() as session:
            result = (
                session.query(Execution.status, func.count().label("count"))
                .join(Combination, Execution.combination_id == Combination.id)
                .filter(Combination.work_id == work_id)
                .group_by(Execution.status)
                .all()
            )

            return {row.status: row.count for row in result}

    def get_execution_stats(self, work_id: str) -> Dict[str, Any]:
        """Get aggregated execution statistics across all combinations.

        Args:
            work_id: Work identifier to get statistics for

        Returns:
            Dictionary containing execution counts by status and totals:
                - completed: Number of completed executions
                - running: Number of running executions
                - queued: Number of queued executions
                - failed: Number of failed/error executions
                - total_sequences: Total sequences across all combinations
                - raw_status_counts: Raw status counts dictionary
        """
        from ..models import Combination, Execution

        with self.session_scope() as session:
            # Status counts at execution level
            status_counts = (
                session.query(Execution.status, func.count().label("cnt"))
                .join(Combination, Execution.combination_id == Combination.id)
                .filter(Combination.work_id == work_id)
                .group_by(Execution.status)
                .all()
            )

            status_dict = {row.status: row.cnt for row in status_counts}

            # Total sequences across combinations
            total_sequences = (
                session.query(func.coalesce(func.sum(Combination.total_sequences), 0))
                .filter(Combination.work_id == work_id)
                .scalar()
                or 0
            )

            completed = status_dict.get("completed", 0)
            running = status_dict.get("running", 0)
            queued = status_dict.get("queued", 0)
            failed = status_dict.get("failed", 0) + status_dict.get("error", 0)

            return {
                "completed": completed,
                "running": running,
                "queued": queued,
                "failed": failed,
                "total_sequences": total_sequences,
                "raw_status_counts": status_dict,
            }

    # -------------- EXPORT SPECIFIC QUERIES --------------

    def get_work_export_data(self, work_id: str) -> Dict[str, Any]:
        """Get work data for export operations.

        Args:
            work_id: Work identifier to export

        Returns:
            Dictionary containing complete work data for export.
        """
        from ..models import Work

        with self.session_scope() as session:
            work = session.query(Work).filter(Work.id == work_id).first()
            return work.to_dict() if work else {}

    def list_combinations(self, work_id: str) -> List[Dict[str, Any]]:
        """List all combinations for a work with combination_id field.

        Args:
            work_id: Work identifier to list combinations for

        Returns:
            List of combination dictionaries, each including a combination_id
            field for compatibility with legacy code.
        """
        from ..models import Combination

        with self.session_scope() as session:
            combinations = (
                session.query(Combination)
                .filter(Combination.work_id == work_id)
                .order_by(Combination.id)
                .all()
            )

            result = []
            for combo in combinations:
                combo_dict = combo.to_dict()
                # Ensure combination_id field is present (some code expects this field name)
                combo_dict["combination_id"] = combo_dict.get("id")
                result.append(combo_dict)

            return result

    def get_combinations_for_export(self, work_id: str) -> List[Dict[str, Any]]:
        """Get all combinations for a work for export operations.

        Args:
            work_id: Work identifier to export combinations for

        Returns:
            List of combination dictionaries with all fields.
        """
        from ..models import Combination

        with self.session_scope() as session:
            combinations = (
                session.query(Combination).filter(Combination.work_id == work_id).all()
            )
            return [combo.to_dict() for combo in combinations]

    def get_executions_for_export(self, work_id: str) -> List[Dict[str, Any]]:
        """Get all executions for a work for export operations.

        Args:
            work_id: Work identifier to export executions for

        Returns:
            List of execution dictionaries with all fields.
        """
        from ..models import Combination, Execution

        with self.session_scope() as session:
            executions = (
                session.query(Execution)
                .join(Combination, Execution.combination_id == Combination.id)
                .filter(Combination.work_id == work_id)
                .all()
            )
            return [exec.to_dict() for exec in executions]

    def get_execution_progress_for_export(self, work_id: str) -> List[Dict[str, Any]]:
        """Get all execution progress records for export operations.

        Args:
            work_id: Work identifier to export progress data for

        Returns:
            List of execution progress dictionaries with all fields.
        """
        from ..models import Combination, Execution, ExecutionProgress

        with self.session_scope() as session:
            progress_records = (
                session.query(ExecutionProgress)
                .join(Execution, ExecutionProgress.execution_id == Execution.id)
                .join(Combination, Execution.combination_id == Combination.id)
                .filter(Combination.work_id == work_id)
                .all()
            )
            return [progress.to_dict() for progress in progress_records]

    def get_events_for_export(self, work_id: str) -> List[Dict[str, Any]]:
        """Get all events for a work for export operations.

        Args:
            work_id: Work identifier to export events for

        Returns:
            List of event dictionaries. Falls back to all events if
            work_id column doesn't exist in the events table.
        """
        from ..models import Event

        with self.session_scope() as session:
            # Check if work_id column exists
            try:
                events = session.query(Event).filter(Event.work_id == work_id).all()
                return [event.to_dict() for event in events]
            except Exception:
                # Fallback: return all events if work_id column doesn't exist
                events = session.query(Event).all()
                return [event.to_dict() for event in events]

    def get_datasets_for_export(self, work_id: str) -> List[Dict[str, Any]]:
        """Get datasets used by a work for export operations.

        Args:
            work_id: Work identifier to export datasets for

        Returns:
            List of dataset dictionaries with all fields.
        """
        from ..models import Dataset

        with self.session_scope() as session:
            datasets = session.query(Dataset).filter(Dataset.work_id == work_id).all()
            return [dataset.to_dict() for dataset in datasets]

    def get_dataset_sequences_for_export(self, work_id: str) -> List[Dict[str, Any]]:
        """Get dataset sequences used by a work for export operations.

        Args:
            work_id: Work identifier to export dataset sequences for

        Returns:
            List of dataset sequence dictionaries with all fields.
        """
        from ..models import Dataset, DatasetSequence

        with self.session_scope() as session:
            sequences = (
                session.query(DatasetSequence)
                .join(Dataset, DatasetSequence.dataset_id == Dataset.id)
                .filter(Dataset.work_id == work_id)
                .all()
            )
            return [seq.to_dict() for seq in sequences]

    def get_optimization_executions_for_export(
        self, work_id: str
    ) -> List[Dict[str, Any]]:
        """Get optimization-specific executions for export operations.

        Filters executions that have unit_id starting with 'optimization:'.

        Args:
            work_id: Work identifier to export optimization executions for

        Returns:
            List of optimization execution dictionaries ordered by unit_id and sequence.
        """
        from ..models import Combination, Execution

        with self.session_scope() as session:
            executions = (
                session.query(Execution)
                .join(Combination, Execution.combination_id == Combination.id)
                .filter(
                    Combination.work_id == work_id,
                    Execution.unit_id.like("optimization:%"),
                )
                .order_by(Execution.unit_id, Execution.sequencia)
                .all()
            )
            return [exec.to_dict() for exec in executions]

    def get_sensitivity_events_for_export(self, work_id: str) -> List[Dict[str, Any]]:
        """Get sensitivity analysis events for export operations.

        Filters progress events related to sensitivity analysis by looking
        for unit_id containing 'sensitivity_analysis'.

        Args:
            work_id: Work identifier to export sensitivity events for

        Returns:
            List of sensitivity analysis event dictionaries ordered by timestamp.
            Falls back to all sensitivity events if work_id column doesn't exist.
        """
        from ..models import Event

        with self.session_scope() as session:
            try:
                # Try with work_id filter first - look for progress events with sensitivity_analysis unit_id
                events = (
                    session.query(Event)
                    .filter(Event.event_type == "progress", Event.work_id == work_id)
                    .filter(
                        Event.entity_data_json.like(
                            '%"unit_id": "sensitivity_analysis"%'
                        )
                    )
                    .order_by(Event.timestamp)
                    .all()
                )
                return [event.to_dict() for event in events]
            except Exception:
                # Fallback: return all sensitivity events if work_id column doesn't exist
                events = (
                    session.query(Event)
                    .filter(Event.event_type == "progress")
                    .filter(
                        Event.entity_data_json.like(
                            '%"unit_id": "sensitivity_analysis"%'
                        )
                    )
                    .order_by(Event.timestamp)
                    .all()
                )
                return [event.to_dict() for event in events]
