"""
WebSocket server for real-time progress monitoring and communication.

This module implements a comprehensive WebSocket server that provides real-time
monitoring capabilities for CSPBench work executions. It handles bidirectional
communication, connection management, and integrates with the work monitoring
system to deliver live progress updates.
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
    """WebSocket server for comprehensive progress monitoring and communication.
    
    This class handles WebSocket connections for real-time monitoring of work
    executions, providing live progress updates, status changes, and bidirectional
    communication capabilities between the server and connected clients.
    
    Features:
        - Real-time work progress monitoring with sub-second updates
        - Bidirectional communication for client commands and server responses
        - Automatic connection management and cleanup
        - Message queuing and delivery optimization
        - Error handling and recovery mechanisms
        - Support for multiple concurrent connections per work
    
    Attributes:
        work_service: Service instance for accessing work data and operations
    """

    def __init__(self):
        """Initialize WebSocket server with work service integration.
        
        Sets up the server with necessary service dependencies and
        prepares for handling WebSocket connections.
        """
        self.work_service = get_work_service()

    async def handle_work_monitor(self, websocket: WebSocket, work_id: str) -> None:
        """Handle WebSocket connection for comprehensive work monitoring.
        
        This method manages the complete lifecycle of a WebSocket connection
        dedicated to monitoring a specific work execution. It handles connection
        establishment, message routing, real-time updates, and proper cleanup.
        
        Args:
            websocket (WebSocket): FastAPI WebSocket connection instance
            work_id (str): Unique identifier of the work to monitor
            
        Features:
            - Immediate connection acceptance and logging
            - Real-time progress updates via message queues
            - Automatic session management and subscriber registration
            - Initial snapshot delivery for immediate state synchronization
            - Bidirectional communication support
            - Graceful error handling and connection cleanup
            
        Raises:
            WebSocketDisconnect: When client disconnects (handled gracefully)
            Exception: For other connection or processing errors
            
        Note:
            Connection lifecycle is fully managed, including automatic cleanup
            of resources when the connection terminates.
        """
        await websocket.accept()
        logger.info(f"[WEBSOCKET] 🟢 WebSocket connected for work monitoring: {work_id}")

        # Create dedicated message queue for this connection
        message_queue = asyncio.Queue()

        try:
            # Get or create monitoring session for this work
            session = await work_monitor_manager.get_or_create_session(work_id)

            # Register this connection as a subscriber
            await session.add_subscriber(message_queue)

            # Start monitoring session if not already running
            if not session.running:
                await session.start()

            # Send immediate snapshot for new connections (first subscriber priority)
            if not session._initial_snapshot_sent:
                try:
                    # Build and send initial snapshot only for this subscriber
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

            # Handle bidirectional communication
            await self._handle_websocket_communication(
                websocket, message_queue, session, work_id
            )

        except WebSocketDisconnect as e:
            logger.info(f"[WEBSOCKET] 🔴 WebSocket disconnected for work {work_id}: code={e.code}, reason='{e.reason}'")
        except Exception as e:
            logger.error(f"[WEBSOCKET] ❌ WebSocket error for work {work_id}: {e}")
            try:
                # Send error notification to client before closing
                await websocket.send_json(
                    {
                        "type": "error",
                        "work_id": work_id,
                        "timestamp": asyncio.get_event_loop().time(),
                        "payload": {"code": "INTERNAL_ERROR", "message": str(e)},
                    }
                )
            except:
                pass  # Connection may already be closed
        finally:
            # Comprehensive cleanup of resources
            try:
                if "session" in locals():
                    logger.info(f"[WEBSOCKET] 🧹 Cleaning up subscriber for work {work_id}")
                    await session.remove_subscriber(message_queue)
                else:
                    logger.warning(f"[WEBSOCKET] ⚠️ Session not found during cleanup for work {work_id}")
            except Exception as e:
                logger.error(f"[WEBSOCKET] ❌ Error during WebSocket cleanup for work {work_id}: {e}")

    async def _handle_websocket_communication(
        self, websocket: WebSocket, message_queue: asyncio.Queue, session, work_id: str
    ) -> None:
        """Handle bidirectional WebSocket communication with optimized message flow.
        
        This method manages the concurrent send/receive operations for WebSocket
        communication, ensuring optimal performance and proper error handling
        for both directions of data flow.
        
        Args:
            websocket (WebSocket): Active WebSocket connection
            message_queue (asyncio.Queue): Queue for outbound messages
            session: Monitoring session instance
            work_id (str): Work identifier for logging and context
            
        Features:
            - Concurrent send and receive message handling
            - Message counting and performance monitoring
            - Client command processing (refresh, combination selection)
            - Heartbeat and connection health monitoring
            - Graceful task cancellation on disconnection
            - Detailed logging for debugging and monitoring
            
        Note:
            Uses asyncio tasks for concurrent operation handling and proper
            cancellation when either direction encounters an error.
        """

        async def send_messages():
            """Send messages from queue to WebSocket with performance monitoring.
            
            Continuously processes the message queue and sends messages to the
            connected WebSocket client. Includes performance monitoring and
            detailed logging for debugging purposes.
            
            Raises:
                WebSocketDisconnect: When client disconnects
                Exception: For other send-related errors
            """
            import json
            message_count = 0
            try:
                logger.debug(f"[WEBSOCKET] 📤 Starting message send loop for work {work_id}")
                while True:
                    message = await message_queue.get()
                    message_count += 1
                    
                    # Detailed logging for important messages and periodic updates
                    try:
                        parsed_msg = json.loads(message)
                        msg_type = parsed_msg.get('type', 'unknown')
                        # Log important events or every 10th message for monitoring
                        if msg_type in ['event', 'error'] or message_count % 10 == 0:
                            logger.debug(f"[WEBSOCKET] 📨 Sending message #{message_count} type '{msg_type}' for work {work_id}")
                    except:
                        pass  # JSON parsing is optional for logging
                    
                    await websocket.send_text(message)
            except WebSocketDisconnect as e:
                logger.info(f"[WEBSOCKET] 🔌 Client disconnected during send for work {work_id} (msg #{message_count}): code={e.code}, reason='{e.reason}'")
                raise  # Re-raise for main task handler
            except Exception as e:
                logger.info(f"[WEBSOCKET] 📤 Send loop ended for work {work_id} after {message_count} messages: {e}")

        async def receive_messages():
            """Receive and process messages from client with command handling.
            
            Handles incoming messages from WebSocket clients, including heartbeat
            responses, refresh commands, and combination selection requests.
            Processes commands and triggers appropriate server responses.
            
            Commands Supported:
                - refresh_executions: Request updated execution status
                - set_combination: Switch to specific parameter combination
                - heartbeat: Connection health maintenance (future)
                
            Raises:
                WebSocketDisconnect: When client disconnects
                Exception: For other receive-related errors
            """
            import json

            from .schemas import MessageType, ProgressMessage, ProgressUpdate

            try:
                while True:
                    raw = await websocket.receive_text()
                    logger.debug(
                        f"Received WebSocket message for work {work_id}: {raw}"
                    )
                    
                    # Parse and process client commands
                    try:
                        msg = json.loads(raw)
                    except Exception:
                        continue  # Skip invalid JSON messages
                        
                    action = msg.get("action")
                    if action in {"refresh_executions", "set_combination"}:
                        # Handle execution refresh for specific or current combination
                        try:
                            combination_id = msg.get("combination_id")
                            executions = None
                            
                            # Validate and process combination ID
                            if combination_id is not None:
                                try:
                                    combination_id = int(combination_id)
                                except ValueError:
                                    combination_id = None
                                    
                            if combination_id is not None:
                                # Fetch complete execution list for requested combination
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
                                # Use partial snapshot from current combination
                                snapshot_message = (
                                    await session._build_snapshot_message(
                                        update_tracking=False
                                    )
                                )
                                if snapshot_message and snapshot_message.snapshot:
                                    # Current combination executions from existing snapshot
                                    update_payload = {
                                        "executions": snapshot_message.snapshot.executions,
                                    }
                                    
                            # Send update response to client
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
                    
                    # Future commands: ping/pong, pause stream, connection control
                    
            except WebSocketDisconnect as e:
                logger.debug(f"Client disconnected from work {work_id}: code={e.code}, reason='{e.reason}'")
                raise  # Re-raise for main task handler
            except Exception as e:
                logger.debug(f"Receive messages loop ended for work {work_id}: {e}")

        # Execute concurrent send and receive tasks
        send_task = asyncio.create_task(send_messages())
        receive_task = asyncio.create_task(receive_messages())

        try:
            # Wait for either task to complete or fail
            done, pending = await asyncio.wait(
                [send_task, receive_task], return_when=asyncio.FIRST_COMPLETED
            )

            # Cancel remaining pending tasks
            for task in pending:
                task.cancel()
                try:
                    await task
                except asyncio.CancelledError:
                    pass  # Expected for cancelled tasks

            # Process completed task results and handle exceptions
            for task in done:
                try:
                    task.result()  # Re-raise any exceptions that occurred
                except WebSocketDisconnect as e:
                    logger.debug(f"WebSocket disconnected for work {work_id}: code={e.code}, reason='{e.reason}'")
                    # Expected behavior, don't re-raise
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
            # Cancel both tasks on error
            send_task.cancel()
            receive_task.cancel()
            try:
                await asyncio.gather(send_task, receive_task, return_exceptions=True)
            except Exception:
                pass  # Ignore cancellation exceptions
            raise

    def _get_work_details(self, work_id: str) -> Optional[Any]:
        """Get comprehensive work details from work service.
        
        Retrieves detailed work information including status, configuration,
        and metadata through the work service interface.
        
        Args:
            work_id (str): Unique identifier of the work
            
        Returns:
            Optional[Any]: Work details object or None if not found/error
            
        Note:
            Includes error handling to prevent service failures from
            affecting WebSocket operations.
        """
        try:
            work_details = self.work_service.get(work_id)
            return work_details
        except Exception as e:
            logger.error(f"Failed to get work details for {work_id}: {e}")
            return None

    async def handle_general_websocket(
        self, websocket: WebSocket, client_id: str
    ) -> None:
        """Handle general WebSocket connection for testing and legacy support.
        
        Provides a general-purpose WebSocket endpoint for testing, debugging,
        and legacy client support. Includes welcome messages and echo functionality
        for connection validation and basic communication testing.
        
        Args:
            websocket (WebSocket): WebSocket connection instance
            client_id (str): Unique client identifier for connection tracking
            
        Features:
            - Welcome message with endpoint information
            - Echo functionality for connection testing
            - Basic heartbeat and communication validation
            - Legacy client support for older implementations
            
        Note:
            This endpoint is primarily for development and testing purposes.
            Production clients should use specific work monitoring endpoints.
        """
        await websocket.accept()
        logger.info(f"General WebSocket connected: {client_id}")

        try:
            # Send comprehensive welcome message with endpoint information
            await websocket.send_json(
                {
                    "type": "welcome",
                    "client_id": client_id,
                    "timestamp": asyncio.get_event_loop().time(),
                    "payload": {
                        "message": "Connected to CSPBench WebSocket Server",
                        "server_version": "1.0.0",
                        "available_endpoints": [
                            "/ws/work/{work_id} - Monitor specific work progress with real-time updates"
                        ],
                        "supported_features": [
                            "Real-time progress monitoring",
                            "Bidirectional communication", 
                            "Automatic reconnection handling",
                            "Multi-combination execution tracking"
                        ]
                    },
                }
            )

            # Simple echo loop for connection testing and validation
            while True:
                data = await websocket.receive_text()
                await websocket.send_json(
                    {
                        "type": "echo",
                        "client_id": client_id,
                        "timestamp": asyncio.get_event_loop().time(),
                        "payload": {"received": data, "echo": f"Server received: {data}"},
                    }
                )

        except WebSocketDisconnect:
            logger.info(f"General WebSocket disconnected: {client_id}")
        except Exception as e:
            logger.error(f"General WebSocket error for {client_id}: {e}")


# Global server instance for application-wide access
websocket_server = WebSocketServer()
