"""
API routes for monitoring functionality.
"""

from typing import List, Dict, Any
from fastapi import APIRouter, HTTPException
from src.application.services.work_service import get_work_service
from pathlib import Path
from src.infrastructure.persistence.work_state.queries import WorkStateQueries

router = APIRouter(prefix="/api/monitor", tags=["monitoring"])


@router.get("/works")
async def get_works() -> Dict[str, Any]:
    """Get list of all works with status information."""
    try:
        work_service = get_work_service()
        works = work_service.list()
        
        # Enrich with status information
        enriched_works = []
        for work in works:
            work_data = {
                "id": work["id"],
                "status": work.get("status", "unknown"),
                "created_at": work.get("created_at", 0),
                "updated_at": work.get("updated_at", 0),
                "output_path": work.get("output_path", ""),
                "error": work.get("error", None),
                "config_name": work.get("config", {}).get("metadata", {}).get("name", "Unknown"),
                "config_description": work.get("config", {}).get("metadata", {}).get("description", ""),
                "config_author": work.get("config", {}).get("metadata", {}).get("author", ""),
            }
            enriched_works.append(work_data)
        
        # Sort by creation date (newest first)
        enriched_works.sort(key=lambda x: x["created_at"], reverse=True)
        
        return {
            "works": enriched_works,
            "total": len(enriched_works)
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to get works: {str(e)}")


@router.get("/work/{work_id}/status")
async def get_work_status(work_id: str) -> Dict[str, Any]:
    """Get detailed status of a specific work."""
    try:
        work_service = get_work_service()
        work = work_service.get(work_id)
        
        if not work:
            raise HTTPException(status_code=404, detail=f"Work {work_id} not found")
        
        return {
            "work_id": work_id,
            "status": work.get("status", "unknown"),
            "created_at": work.get("created_at", 0),
            "updated_at": work.get("updated_at", 0),
            "output_path": work.get("output_path", ""),
            "error": work.get("error", None),
            "config": work.get("config", {}),
        }
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to get work status: {str(e)}")


@router.get("/work/{work_id}/database-status")
async def get_work_database_status(work_id: str) -> Dict[str, Any]:
    """Check if work database is ready for monitoring."""
    try:
    # imports movidos para topo
        
        work_service = get_work_service()
        work = work_service.get(work_id)
        
        if not work:
            raise HTTPException(status_code=404, detail=f"Work {work_id} not found")
        
        db_path = Path(work["output_path"]) / "state.db"
        
        if not db_path.exists():
            return {
                "ready": False,
                "reason": "database_not_created",
                "message": "Database file not yet created"
            }
        
        try:
            with WorkStateQueries(db_path) as queries:
                if not queries.work_exists(work_id):
                    return {
                        "ready": False,
                        "reason": "work_not_in_database",
                        "message": "Work not yet written to database"
                    }
                
                # Check if there's meaningful progress data
                progress = queries.get_work_progress_summary(work_id)
                if not progress:
                    return {
                        "ready": False,
                        "reason": "no_progress_data",
                        "message": "No progress data available yet"
                    }
                
                return {
                    "ready": True,
                    "reason": "ready",
                    "message": "Database ready for monitoring"
                }
                
        except Exception as db_error:
            return {
                "ready": False,
                "reason": "database_error",
                "message": f"Database error: {str(db_error)}"
            }
            
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to check database status: {str(e)}")


@router.post("/work/{work_id}/action/{action}")
async def work_action(work_id: str, action: str) -> Dict[str, Any]:
    """Perform action on work (pause, resume, cancel)."""
    try:
        work_service = get_work_service()
        
        if action == "pause":
            success = work_service.pause(work_id)
        elif action == "resume":
            # Use ExecutionManager para retomar execução real
            from src.application.services.execution_manager import ExecutionManager
            exec_manager = ExecutionManager(work_service=work_service)
            success = exec_manager.resume(work_id)
        elif action == "cancel":
            success = work_service.cancel(work_id)
        else:
            raise HTTPException(status_code=400, detail=f"Invalid action: {action}")
        
        if not success:
            raise HTTPException(status_code=400, detail=f"Failed to {action} work {work_id}")
        
        return {
            "success": True,
            "action": action,
            "work_id": work_id,
            "message": f"Work {action} executed successfully"
        }
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to {action} work: {str(e)}")


@router.get("/active-works")
async def get_active_works() -> Dict[str, Any]:
    """Get only active (running/queued) works."""
    try:
        work_service = get_work_service()
        all_works = work_service.list()
        
        active_statuses = ["queued", "running"]
        active_works = [
            work for work in all_works 
            if work.get("status", "unknown") in active_statuses
        ]
        
        return {
            "works": active_works,
            "count": len(active_works)
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to get active works: {str(e)}")


@router.get("/work/{work_id}/combinations")
async def list_work_combinations(work_id: str) -> Dict[str, Any]:
    """Lista combinações de um work para popular filtros de execução."""
    try:
        work_service = get_work_service()
        work = work_service.get(work_id)
        if not work:
            raise HTTPException(status_code=404, detail=f"Work {work_id} not found")
        db_path = Path(work["output_path"]) / "state.db"
        if not db_path.exists():
            return {"combinations": [], "total": 0}
        with WorkStateQueries(db_path) as queries:
            if not queries.work_exists(work_id):
                return {"combinations": [], "total": 0}
            combos = queries.list_combinations(work_id)
        return {"combinations": combos, "total": len(combos)}
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to list combinations: {str(e)}")


@router.get("/work/{work_id}/combinations/{combination_id}/executions")
async def get_combination_executions(work_id: str, combination_id: int) -> Dict[str, Any]:
    """Detalhes e estatísticas das execuções de uma combinação específica."""
    try:
        work_service = get_work_service()
        work = work_service.get(work_id)
        if not work:
            raise HTTPException(status_code=404, detail=f"Work {work_id} not found")
        db_path = Path(work["output_path"]) / "state.db"
        if not db_path.exists():
            raise HTTPException(status_code=404, detail="State database not found")
        with WorkStateQueries(db_path) as queries:
            if not queries.work_exists(work_id):
                raise HTTPException(status_code=404, detail="Work not registered in state database")
            # Verificar se combinação existe
            combo_list = queries.list_combinations(work_id)
            combo_map = {c["combination_id"]: c for c in combo_list}
            combo = combo_map.get(combination_id)
            if not combo:
                raise HTTPException(status_code=404, detail=f"Combination {combination_id} not found for work {work_id}")
            executions = queries.get_combination_executions_detail(combination_id)
            # Agregar estatísticas
            stats = {"Running": 0, "Completed": 0, "Queued": 0, "Failed": 0, "Total": combo.get("total_sequences") or 0}
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
                exec_rows.append({
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
                })
        return {
            "work_id": work_id,
            "combination_id": combination_id,
            "combination": combo,
            "stats": stats,
            "executions": exec_rows,
            "count": len(exec_rows)
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to get executions: {str(e)}")
