"""
Hierarchical Monitoring API Routes

Provides endpoints for accessing hierarchical monitoring data.
"""

import logging
from typing import Dict, Any, List
from fastapi import APIRouter, HTTPException
from fastapi.responses import JSONResponse

# Import hierarchical monitoring components
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../../'))

from hierarchical_web_monitor import HierarchicalWebMonitor

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/monitoring", tags=["hierarchical-monitoring"])

# Global monitoring instance
_monitor_instance = None

def get_monitor() -> HierarchicalWebMonitor:
    """Get or create the hierarchical monitoring instance."""
    global _monitor_instance
    if _monitor_instance is None:
        _monitor_instance = HierarchicalWebMonitor()
    return _monitor_instance


@router.get("/health")
async def monitoring_health():
    """Health check for monitoring system."""
    try:
        monitor = get_monitor()
        return {
            "status": "healthy",
            "monitor_type": "hierarchical",
            "active_runs": len(monitor.active_runs),
            "hierarchy_initialized": monitor.hierarchy is not None
        }
    except Exception as e:
        logger.error(f"Monitoring health check failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/status")
async def get_full_status():
    """Get complete hierarchical status."""
    try:
        monitor = get_monitor()
        return monitor.get_full_status()
    except Exception as e:
        logger.error(f"Failed to get full status: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/hierarchy")
async def get_hierarchy():
    """Get execution hierarchy structure."""
    try:
        monitor = get_monitor()
        return monitor.get_hierarchy_status()
    except Exception as e:
        logger.error(f"Failed to get hierarchy: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/active-runs")
async def get_active_runs():
    """Get all active runs."""
    try:
        monitor = get_monitor()
        active_runs = monitor.get_active_runs()
        return [
            {
                "run_id": run.run_id,
                "algorithm_name": run.algorithm_name,
                "run_info": run.run_info,
                "status": run.status.value if hasattr(run.status, 'value') else str(run.status),
                "progress": run.progress,
                "start_time": run.start_time.isoformat() if run.start_time else None,
                "current_generation": run.current_generation,
                "total_generations": run.total_generations,
                "current_message": run.current_message
            }
            for run in active_runs
        ]
    except Exception as e:
        logger.error(f"Failed to get active runs: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/statistics")
async def get_statistics():
    """Get algorithm statistics."""
    try:
        monitor = get_monitor()
        return monitor.get_algorithm_statistics()
    except Exception as e:
        logger.error(f"Failed to get statistics: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/real-time")
async def get_real_time_data():
    """Get real-time monitoring data."""
    try:
        monitor = get_monitor()
        return monitor.get_real_time_data()
    except Exception as e:
        logger.error(f"Failed to get real-time data: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/callbacks/recent")
async def get_recent_callbacks():
    """Get recent callback entries."""
    try:
        monitor = get_monitor()
        callbacks = monitor.get_recent_callbacks()
        return [
            {
                "timestamp": callback.timestamp.isoformat() if callback.timestamp else None,
                "run_id": callback.run_id,
                "algorithm_name": callback.algorithm_name,
                "progress": callback.progress,
                "message": callback.message,
                "generation": callback.generation,
                "context": callback.context.get_path() if callback.context else None
            }
            for callback in callbacks
        ]
    except Exception as e:
        logger.error(f"Failed to get recent callbacks: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/performance")
async def get_performance_metrics():
    """Get performance metrics."""
    try:
        monitor = get_monitor()
        return monitor.get_performance_metrics()
    except Exception as e:
        logger.error(f"Failed to get performance metrics: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/history")
async def get_completed_runs_history():
    """Get completed runs history."""
    try:
        monitor = get_monitor()
        return monitor.get_completed_runs_history()
    except Exception as e:
        logger.error(f"Failed to get completed runs history: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/algorithm/{algorithm_name}")
async def get_algorithm_details(algorithm_name: str):
    """Get detailed data for specific algorithm."""
    try:
        monitor = get_monitor()
        return monitor.get_detailed_algorithm_data(algorithm_name)
    except Exception as e:
        logger.error(f"Failed to get algorithm details: {e}")
        raise HTTPException(status_code=500, detail=str(e))
