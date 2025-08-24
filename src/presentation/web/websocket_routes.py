"""WebSocket routes for CSPBench web interface.

Provides a minimal `/ws/{client_id}` endpoint that registers connections in the
`connection_manager` so monitoring and other parts of the app can report stats.

This implementation is intentionally small and robust for environments where a
full-featured WebSocket frontend is optional.
"""

from fastapi import APIRouter, WebSocket, WebSocketDisconnect
from typing import Optional
from datetime import datetime

from .websocket_manager import connection_manager

router = APIRouter()


@router.websocket("/ws/{client_id}")
async def websocket_endpoint(
    websocket: WebSocket, client_id: str, token: Optional[str] = None
):
    """Accept a WebSocket connection and register it in connection_manager.

    The optional `token` query parameter is accepted for future auth checks but
    currently ignored by this minimal implementation.
    """
    await websocket.accept()
    # Register connection with metadata
    await connection_manager.register(
        client_id, {"connected_at": datetime.utcnow().isoformat(), "token": token}
    )

    try:
        while True:
            # Keep the connection open; handle incoming messages if any.
            data = await websocket.receive_text()
            # Echo received messages back to the client
            await websocket.send_text(f"echo: {data}")
    except WebSocketDisconnect:
        await connection_manager.unregister(client_id)
    except Exception:
        # Ensure cleanup on unexpected errors
        await connection_manager.unregister(client_id)
        try:
            await websocket.close()
        except Exception:
            pass
