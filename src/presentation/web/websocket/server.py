"""
WebSocket server for real-time progress monitoring.
"""

import asyncio
import logging
from pathlib import Path
from typing import Any, Dict, Optional

from fastapi import WebSocket, WebSocketDisconnect

from src.application.services.work_service import get_work_service

from .manager import work_monitor_manager

logger = logging.getLogger(__name__)


class WebSocketServer:
    """WebSocket server for progress monitoring."""

    def __init__(self):
        self.work_service = get_work_service()

    async def handle_work_monitor(self, websocket: WebSocket, work_id: str) -> None:
        """Handle WebSocket connection for work monitoring."""
        await websocket.accept()
        logger.info(f"[WEBSOCKET] 🟢 WebSocket conectado para monitoramento do work: {work_id}")

        # Create message queue for this connection
        message_queue = asyncio.Queue()

        try:
            # Get or create monitoring session
            session = await work_monitor_manager.get_or_create_session(work_id)

            # Add this connection as subscriber
            await session.add_subscriber(message_queue)

            # Start monitoring if not already running
            if not session.running:
                await session.start()

            # Garantir que primeiro subscriber receba snapshot imediatamente
            if not session._initial_snapshot_sent:
                try:
                    # Construir e enviar snapshot apenas para este subscriber
                    snapshot_message = await session._build_snapshot_message(
                        update_tracking=True
                    )
                    if snapshot_message:
                        ws_message = snapshot_message.to_websocket_message().to_json()
                        await message_queue.put(ws_message)
                        session._initial_snapshot_sent = True
                        logger.debug(
                            f"Immediate snapshot sent on connection for work {work_id}"
                        )
                except Exception as e:
                    logger.warning(
                        f"Failed to send immediate snapshot for work {work_id}: {e}"
                    )

            # Send/receive loop
            await self._handle_websocket_communication(
                websocket, message_queue, session, work_id
            )

        except WebSocketDisconnect as e:
            logger.info(f"[WEBSOCKET] 🔴 WebSocket desconectado para work {work_id}: code={e.code}, reason='{e.reason}'")
        except Exception as e:
            logger.error(f"[WEBSOCKET] ❌ Erro no WebSocket para work {work_id}: {e}")
            try:
                await websocket.send_json(
                    {
                        "type": "error",
                        "work_id": work_id,
                        "timestamp": asyncio.get_event_loop().time(),
                        "payload": {"code": "INTERNAL_ERROR", "message": str(e)},
                    }
                )
            except:
                pass
        finally:
            # Clean up
            try:
                if "session" in locals():
                    logger.info(f"[WEBSOCKET] 🧹 Limpando subscriber para work {work_id}")
                    await session.remove_subscriber(message_queue)
                else:
                    logger.warning(f"[WEBSOCKET] ⚠️ Session não encontrada durante cleanup para work {work_id}")
            except Exception as e:
                logger.error(f"[WEBSOCKET] ❌ Erro durante cleanup do WebSocket para work {work_id}: {e}")

    async def _handle_websocket_communication(
        self, websocket: WebSocket, message_queue: asyncio.Queue, session, work_id: str
    ) -> None:
        """Handle bidirectional WebSocket communication."""

        async def send_messages():
            """Send messages from queue to WebSocket."""
            import json
            message_count = 0
            try:
                logger.debug(f"[WEBSOCKET] 📤 Iniciando loop de envio de mensagens para work {work_id}")
                while True:
                    message = await message_queue.get()
                    message_count += 1
                    
                    # Log detalhado para mensagens importantes
                    try:
                        parsed_msg = json.loads(message)
                        msg_type = parsed_msg.get('type', 'unknown')
                        if msg_type in ['event', 'error'] or message_count % 10 == 0:  # Log eventos importantes ou a cada 10 mensagens
                            logger.debug(f"[WEBSOCKET] 📨 Enviando mensagem #{message_count} tipo '{msg_type}' para work {work_id}")
                    except:
                        pass
                    
                    await websocket.send_text(message)
            except WebSocketDisconnect as e:
                logger.info(f"[WEBSOCKET] 🔌 Cliente desconectou durante envio para work {work_id} (msg #{message_count}): code={e.code}, reason='{e.reason}'")
                raise  # Re-raise to be handled by the main task handler
            except Exception as e:
                logger.info(f"[WEBSOCKET] 📤 Loop de envio finalizado para work {work_id} após {message_count} mensagens: {e}")

        async def receive_messages():
            """Receive messages from client (heartbeat / control)."""
            import json

            from .schemas import MessageType, ProgressMessage, ProgressUpdate

            try:
                while True:
                    raw = await websocket.receive_text()
                    logger.debug(
                        f"Received WebSocket message for work {work_id}: {raw}"
                    )
                    try:
                        msg = json.loads(raw)
                    except Exception:
                        continue
                    action = msg.get("action")
                    if action in {"refresh_executions", "set_combination"}:
                        # Permite refresh para combinação arbitrária (se fornecida) ou corrente
                        try:
                            combination_id = msg.get("combination_id")
                            executions = None
                            if combination_id is not None:
                                try:
                                    combination_id = int(combination_id)
                                except ValueError:
                                    combination_id = None
                            if combination_id is not None:
                                # Busca lista completa para a combinação solicitada
                                executions = (
                                    await session.get_full_executions_for_combination(
                                        combination_id
                                    )
                                )
                                update_payload = {
                                    "executions": executions or [],
                                    "executions_combination_id": combination_id,
                                }
                            else:
                                # Usa snapshot parcial da combinação corrente
                                snapshot_message = (
                                    await session._build_snapshot_message(
                                        update_tracking=False
                                    )
                                )
                                if snapshot_message and snapshot_message.snapshot:
                                    # Current combination (se existir) já vem em snapshot.executions
                                    update_payload = {
                                        "executions": snapshot_message.snapshot.executions,
                                    }
                            if update_payload:
                                message = ProgressMessage(
                                    type=MessageType.UPDATE,
                                    work_id=work_id,
                                    timestamp=asyncio.get_event_loop().time(),
                                    update=ProgressUpdate(**update_payload),
                                )
                                await message_queue.put(
                                    message.to_websocket_message().to_json()
                                )
                        except Exception as e:
                            logger.warning(f"Failed {action} action: {e}")
                    # futuro: ping, pause stream, etc.
            except WebSocketDisconnect as e:
                logger.debug(f"Client disconnected from work {work_id}: code={e.code}, reason='{e.reason}'")
                raise  # Re-raise to be handled by the main task handler
            except Exception as e:
                logger.debug(f"Receive messages loop ended for work {work_id}: {e}")

        # Run both send and receive tasks concurrently
        send_task = asyncio.create_task(send_messages())
        receive_task = asyncio.create_task(receive_messages())

        try:
            # Wait for either task to complete/fail
            done, pending = await asyncio.wait(
                [send_task, receive_task], return_when=asyncio.FIRST_COMPLETED
            )

            # Cancel pending tasks
            for task in pending:
                task.cancel()
                try:
                    await task
                except asyncio.CancelledError:
                    pass

            # Check for exceptions in completed tasks
            for task in done:
                try:
                    task.result()  # This will re-raise any exception that occurred
                except WebSocketDisconnect as e:
                    logger.debug(f"WebSocket disconnected for work {work_id}: code={e.code}, reason='{e.reason}'")
                    # This is expected behavior, don't re-raise
                    break
                except Exception as e:
                    logger.error(f"Task exception for work {work_id}: {e}")
                    raise

        except WebSocketDisconnect as e:
            logger.debug(f"WebSocket disconnected for work {work_id}: code={e.code}, reason='{e.reason}'")
            # Cancel both tasks gracefully
            send_task.cancel()
            receive_task.cancel()
            try:
                await asyncio.gather(send_task, receive_task, return_exceptions=True)
            except Exception:
                pass  # Ignore cancellation exceptions
        except Exception as e:
            logger.error(f"Error in WebSocket communication for work {work_id}: {e}")
            # Cancel both tasks
            send_task.cancel()
            receive_task.cancel()
            try:
                await asyncio.gather(send_task, receive_task, return_exceptions=True)
            except Exception:
                pass  # Ignore cancellation exceptions
            raise

    def _get_work_details(self, work_id: str) -> Optional[Any]:
        """Get work details from work service."""
        try:
            work_details = self.work_service.get(work_id)
            return work_details
        except Exception as e:
            logger.error(f"Failed to get work details for {work_id}: {e}")
            return None

    async def handle_general_websocket(
        self, websocket: WebSocket, client_id: str
    ) -> None:
        """Handle general WebSocket connection (legacy support)."""
        await websocket.accept()
        logger.info(f"General WebSocket connected: {client_id}")

        try:
            # Send welcome message
            await websocket.send_json(
                {
                    "type": "welcome",
                    "client_id": client_id,
                    "timestamp": asyncio.get_event_loop().time(),
                    "payload": {
                        "message": "Connected to CSPBench WebSocket",
                        "available_endpoints": [
                            "/ws/work/{work_id} - Monitor specific work progress"
                        ],
                    },
                }
            )

            # Simple echo loop for testing
            while True:
                data = await websocket.receive_text()
                await websocket.send_json(
                    {
                        "type": "echo",
                        "client_id": client_id,
                        "timestamp": asyncio.get_event_loop().time(),
                        "payload": {"received": data},
                    }
                )

        except WebSocketDisconnect:
            logger.info(f"General WebSocket disconnected: {client_id}")
        except Exception as e:
            logger.error(f"General WebSocket error for {client_id}: {e}")


# Global server instance
websocket_server = WebSocketServer()
