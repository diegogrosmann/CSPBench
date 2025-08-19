"""Minimal WebSocket connection manager used by the web presentation layer.

This file provides a lightweight ConnectionManager implementation exposing a
`connection_manager` instance with an async `get_connection_stats()` method.

The implementation is intentionally minimal: it does not start a server or
manage real socket objects, but keeps an in-memory registry that other parts
of the application can use to register/unregister connections if needed.
"""
from typing import Dict, List, Any
import asyncio


class ConnectionManager:
    """A minimal connection manager for WebSocket endpoints.

    Methods are async-friendly to match FastAPI WebSocket handlers.
    """

    def __init__(self):
        # Protect internal state with a lock for concurrency
        self._lock = asyncio.Lock()
        # connections: mapping from connection id -> metadata
        self._connections: Dict[str, Dict[str, Any]] = {}

    async def register(self, conn_id: str, metadata: Dict[str, Any] = None) -> None:
        """Register a connection with optional metadata."""
        async with self._lock:
            self._connections[conn_id] = metadata or {}

    async def unregister(self, conn_id: str) -> None:
        """Unregister/close a connection."""
        async with self._lock:
            self._connections.pop(conn_id, None)

    async def get_connection(self, conn_id: str) -> Dict[str, Any]:
        """Return metadata for a connection or empty dict if missing."""
        async with self._lock:
            return dict(self._connections.get(conn_id, {}))

    async def get_connection_stats(self) -> Dict[str, Any]:
        """Return summary stats about current connections.

        This mirrors the small API expected by monitoring routes. Example response:
        {
            "total": 0,
            "connections": []
        }
        """
        async with self._lock:
            connections = [
                {"id": cid, **(meta if meta is not None else {})}
                for cid, meta in self._connections.items()
            ]
            return {"total": len(connections), "connections": connections}


# Single global instance imported by other modules
connection_manager = ConnectionManager()
