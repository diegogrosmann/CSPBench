"""
WebSocket routes for real-time monitoring.
"""

import logging
from fastapi import APIRouter, WebSocket
from .server import websocket_server

logger = logging.getLogger(__name__)

router = APIRouter()


@router.websocket("/ws/{client_id}")
async def websocket_endpoint(websocket: WebSocket, client_id: str):
    """General WebSocket endpoint for testing and legacy support."""
    await websocket_server.handle_general_websocket(websocket, client_id)


@router.websocket("/ws/work/{work_id}")
async def websocket_work_monitor(websocket: WebSocket, work_id: str):
    """WebSocket endpoint for monitoring specific work progress."""
    await websocket_server.handle_work_monitor(websocket, work_id)
