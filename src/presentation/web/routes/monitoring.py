"""
API routes for monitoring functionality.
"""

import json
from pathlib import Path
from typing import Any, Dict

from fastapi import APIRouter, HTTPException

from src.application.services.work_service import get_work_service
from src.infrastructure.persistence.work_state.core import WorkPersistence

router = APIRouter(prefix="/api/monitor", tags=["monitoring"])


def _get_persistence() -> WorkPersistence:
    """Get WorkPersistence instance following console pattern."""
    return WorkPersistence()


@router.get("/works")
async def get_works(
    page: int = 1,
    per_page: int = 20,
    status: str = None,
    search: str = None,
    author: str = None,
    tags: str = None,
) -> Dict[str, Any]:
    """Get list of all works with status information, filtering and pagination."""
    try:
        persistence = _get_persistence()
        
        # Get all works using WorkPersistence
        work_items, _ = persistence.work_list()

        # Enrich with status information and batch metadata
        enriched_works = []
        for work_item in work_items:
            # Extract config from config_json field
            config = work_item.get("config_json", {})
            # If config_json is already a dict (SQLAlchemy JSONType), use it directly
            # If it's a string, parse it as JSON
            if isinstance(config, str):
                try:
                    config = json.loads(config)
                except json.JSONDecodeError:
                    config = {}

            metadata = config.get("metadata", {})
            work_data = {
                "id": work_item["id"],
                "status": work_item.get("status", "unknown"),
                "created_at": work_item.get("created_at", 0),
                "updated_at": work_item.get("updated_at", 0),
                "output_path": work_item.get("output_path", ""),
                "error": work_item.get("error", None),
                "config_name": metadata.get("name", "Unknown"),
                "config_description": metadata.get("description", ""),
                "config_author": metadata.get("author", ""),
                "config_version": metadata.get("version", ""),
                "config_creation_date": metadata.get("creation_date", ""),
                "config_tags": metadata.get("tags", []),
            }
            enriched_works.append(work_data)

        # Apply filters
        filtered_works = enriched_works

        # Filter by status
        if status:
            filtered_works = [w for w in filtered_works if w["status"] == status]

        # Filter by author
        if author:
            filtered_works = [
                w
                for w in filtered_works
                if author.lower() in w["config_author"].lower()
            ]

        # Filter by tags
        if tags:
            tag_list = [t.strip().lower() for t in tags.split(",")]
            filtered_works = [
                w
                for w in filtered_works
                if any(tag in [t.lower() for t in w["config_tags"]] for tag in tag_list)
            ]

        # Filter by search term (searches in multiple fields)
        if search:
            search_lower = search.lower()
            filtered_works = [
                w
                for w in filtered_works
                if (
                    search_lower in w["id"].lower()
                    or search_lower in w["config_name"].lower()
                    or search_lower in w["config_description"].lower()
                    or search_lower in w["config_author"].lower()
                    or any(search_lower in tag.lower() for tag in w["config_tags"])
                )
            ]

        # Sort by creation date (newest first)
        filtered_works.sort(key=lambda x: x["created_at"], reverse=True)

        # Calculate pagination
        total_items = len(filtered_works)
        total_pages = (total_items + per_page - 1) // per_page
        start_idx = (page - 1) * per_page
        end_idx = start_idx + per_page
        paginated_works = filtered_works[start_idx:end_idx]

        return {
            "works": paginated_works,
            "total": total_items,
            "page": page,
            "per_page": per_page,
            "total_pages": total_pages,
            "has_next": page < total_pages,
            "has_prev": page > 1,
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to get works: {str(e)}")


@router.get("/work/{work_id}/status")
async def get_work_status(work_id: str) -> Dict[str, Any]:
    """Get detailed status of a specific work."""
    try:
        persistence = _get_persistence()
        work = persistence.work_get(work_id)

        if not work:
            raise HTTPException(status_code=404, detail=f"Work {work_id} not found")

        # Extract config from work data
        config = work.get("config_json", {})
        # If config_json is already a dict (SQLAlchemy JSONType), use it directly
        # If it's a string, parse it as JSON
        if isinstance(config, str):
            try:
                config = json.loads(config)
            except json.JSONDecodeError:
                config = {}

        return {
            "work_id": work_id,
            "status": work.get("status", "unknown"),
            "created_at": work.get("created_at"),
            "updated_at": work.get("updated_at"),
            "output_path": work.get("output_path"),
            "error": work.get("error"),
            "config": config,
        }

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(
            status_code=500, detail=f"Failed to get work status: {str(e)}"
        )


@router.get("/work/{work_id}/database-status")
async def get_work_database_status(work_id: str) -> Dict[str, Any]:
    """Check if work database is ready for monitoring."""
    try:
        persistence = _get_persistence()
        work = persistence.work_get(work_id)

        if not work:
            raise HTTPException(status_code=404, detail=f"Work {work_id} not found")

        # Check if there's meaningful progress data
        try:
            progress = persistence.get_work_progress_summary(work_id)
            if not progress:
                return {
                    "ready": False,
                    "reason": "no_progress_data",
                    "message": "No progress data available yet",
                }

            return {
                "ready": True,
                "reason": "ready",
                "message": "Database ready for monitoring",
            }

        except Exception as db_error:
            return {
                "ready": False,
                "reason": "database_error",
                "message": f"Database error: {str(db_error)}",
            }

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(
            status_code=500, detail=f"Failed to check database status: {str(e)}"
        )


@router.post("/work/{work_id}/action/{action}")
async def work_action(work_id: str, action: str) -> Dict[str, Any]:
    """Perform action on work (pause, resume, cancel)."""
    try:
        work_service = get_work_service()

        if action == "pause":
            success = work_service.pause(work_id)
        elif action == "resume":
            # Use WorkManager para retomar execução real
            success = work_service.resume_execution(work_id)
        elif action == "cancel":
            success = work_service.cancel(work_id)
        else:
            raise HTTPException(status_code=400, detail=f"Invalid action: {action}")

        if not success:
            raise HTTPException(
                status_code=400, detail=f"Failed to {action} work {work_id}"
            )

        return {
            "success": True,
            "action": action,
            "work_id": work_id,
            "message": f"Work {action} executed successfully",
        }

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(
            status_code=500, detail=f"Failed to {action} work: {str(e)}"
        )


@router.get("/active-works")
async def get_active_works() -> Dict[str, Any]:
    """Get only active (running/queued) works."""
    try:
        persistence = _get_persistence()
        work_items, _ = persistence.work_list()  # Unpack the tuple

        active_statuses = ["queued", "running"]
        active_works = [
            work_item
            for work_item in work_items
            if work_item.get("status") in active_statuses
        ]

        return {"works": active_works, "count": len(active_works)}

    except Exception as e:
        raise HTTPException(
            status_code=500, detail=f"Failed to get active works: {str(e)}"
        )


@router.get("/work/{work_id}/combinations")
async def list_work_combinations(work_id: str) -> Dict[str, Any]:
    """Lista combinações de um work para popular filtros de execução."""
    try:
        persistence = _get_persistence()
        work = persistence.work_get(work_id)
        if not work:
            raise HTTPException(status_code=404, detail=f"Work {work_id} not found")
        
        combos = persistence.list_combinations(work_id)
        return {"combinations": combos, "total": len(combos)}
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(
            status_code=500, detail=f"Failed to list combinations: {str(e)}"
        )


@router.get("/work/{work_id}/combinations/{combination_id}/executions")
async def get_combination_executions(
    work_id: str, combination_id: int
) -> Dict[str, Any]:
    """Detalhes e estatísticas das execuções de uma combinação específica."""
    try:
        persistence = _get_persistence()
        work = persistence.work_get(work_id)
        if not work:
            raise HTTPException(status_code=404, detail=f"Work {work_id} not found")
        
        # Verificar se combinação existe
        combo_list = persistence.list_combinations(work_id)
        combo_map = {c["combination_id"]: c for c in combo_list}
        combo = combo_map.get(combination_id)
        if not combo:
            raise HTTPException(
                status_code=404,
                detail=f"Combination {combination_id} not found for work {work_id}",
            )
        
        executions = persistence.get_combination_executions_detail(combination_id)
        
        # Agregar estatísticas
        stats = {
            "Running": 0,
            "Completed": 0,
            "Queued": 0,
            "Failed": 0,
            "Total": combo.get("total_sequences") or 0,
        }
        exec_rows = []
        for ex in executions:
            status = ex.status
            if status == "running":
                stats["Running"] += 1
            elif status == "completed":
                stats["Completed"] += 1
            elif status == "queued":
                stats["Queued"] += 1
            elif status in ("failed", "error"):
                stats["Failed"] += 1
            exec_rows.append(
                {
                    "unit_id": ex.unit_id,
                    "combination_id": ex.combination_id,
                    "sequencia": ex.sequencia,
                    "status": ex.status,
                    "progress": ex.progress,
                    "progress_message": ex.progress_message,
                    "started_at": ex.started_at,
                    "finished_at": ex.finished_at,
                    "objective": ex.objective,
                    "task_id": ex.task_id,
                    "dataset_id": ex.dataset_id,
                    "preset_id": ex.preset_id,
                    "algorithm_id": ex.algorithm_id,
                    "mode": ex.mode,
                    "total_sequences": ex.total_sequences,
                }
            )
        
        return {
            "work_id": work_id,
            "combination_id": combination_id,
            "combination": combo,
            "stats": stats,
            "executions": exec_rows,
            "count": len(exec_rows),
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(
            status_code=500, detail=f"Failed to get executions: {str(e)}"
        )
