#!/usr/bin/env python3
"""
WebSocket Monitor Debug Script

Este script mostra em detalhes o que o WebSocket está enviando.
"""

import asyncio
import json
import logging
import sys
import time
import websockets
from pathlib import Path
import pprint

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

async def debug_websocket_messages():
    """Debug WebSocket messages in detail."""
    
    print("🔍 WebSocket Message Debug")
    print("=" * 50)
    
    # Test general WebSocket first
    print("1. Testing General WebSocket...")
    try:
        websocket_url = "ws://localhost:8000/ws/debug_client"
        print(f"🔗 Connecting to: {websocket_url}")
        
        async with websockets.connect(websocket_url) as websocket:
            print("✅ Connected!")
            
            # Receive welcome message
            message = await asyncio.wait_for(websocket.recv(), timeout=5.0)
            data = json.loads(message)
            
            print("📨 Welcome Message:")
            pprint.pprint(data, indent=2)
            print()
            
    except Exception as e:
        print(f"❌ General WebSocket test failed: {e}")
        return False

    # Test work monitoring with an existing work if available
    print("2. Testing Work Monitor WebSocket...")
    
    # First, let's see what works are available
    try:
        from src.application.services.work_service import get_work_service
        work_service = get_work_service()
        works = work_service.list()
        
        print(f"📋 Available works: {len(works)}")
        
        if not works:
            print("ℹ️  No works available, creating a test work...")
            return await create_and_monitor_test_work()
        
        # Use the most recent work
        latest_work = max(works, key=lambda w: w.get('created_at', 0))
        work_id = latest_work['id']
        print(f"🎯 Using work: {work_id}")
        print(f"   Status: {latest_work.get('status', 'unknown')}")
        print(f"   Created: {latest_work.get('created_at', 'unknown')}")
        
        # Try to connect to this work
        await monitor_work_websocket(work_id)
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    return True


async def create_and_monitor_test_work():
    """Create a minimal test work and monitor it."""
    print("🛠️  Creating minimal test work...")
    
    try:
        from src.application.services.execution_manager import ExecutionManager
        from src.domain.config import load_cspbench_config
        import tempfile
        import yaml
        
        # Minimal test configuration
        test_config = {
            "metadata": {
                "name": "WebSocket Debug Test",
                "description": "Minimal test for WebSocket debugging",
                "author": "Debug",
                "version": "1.0",
                "creation_date": "2025-08-26",
                "tags": ["debug", "websocket"]
            },
            "datasets": [
                {
                    "id": "debug_dataset",
                    "name": "Debug Dataset",
                    "type": "synthetic",
                    "sequences": 2,
                    "length": 8,
                    "alphabet": "ATCG",
                    "seed": 42
                }
            ],
            "algorithms": [
                {
                    "id": "debug_config",
                    "name": "Debug Config",
                    "algorithms": ["Baseline"],
                    "algorithm_params": {
                        "Baseline": {
                            "tie_break": "first"
                        }
                    }
                }
            ],
            "task": {"type": "experiment"},
            "experiment": {
                "tasks": [
                    {
                        "id": "debug_task",
                        "name": "Debug Task",
                        "datasets": ["debug_dataset"],
                        "algorithms": ["debug_config"],
                        "repetitions": 1,
                        "timeout": 30
                    }
                ]
            },
            "output": {
                "logging": True,
                "results": {"partial_results": True}
            }
        }
        
        # Save to temp file and load
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(test_config, f, default_flow_style=False)
            temp_file = f.name
        
        try:
            # Load and execute
            config = load_cspbench_config(Path(temp_file))
            execution_manager = ExecutionManager()
            work_id = execution_manager.execute(config=config, extra={"origin": "debug_test"})
            
            print(f"✅ Created work: {work_id}")
            
            # Wait a moment
            await asyncio.sleep(2)
            
            # Monitor it
            await monitor_work_websocket(work_id)
            
        finally:
            Path(temp_file).unlink(missing_ok=True)
            
    except Exception as e:
        print(f"❌ Failed to create test work: {e}")
        import traceback
        traceback.print_exc()
        return False


async def monitor_work_websocket(work_id: str):
    """Monitor a specific work via WebSocket."""
    print(f"\n3. Monitoring Work {work_id} via WebSocket...")
    
    websocket_url = f"ws://localhost:8000/ws/work/{work_id}"
    print(f"🔗 Connecting to: {websocket_url}")
    
    try:
        async with websockets.connect(websocket_url) as websocket:
            print("✅ WebSocket connected!")
            
            message_count = 0
            start_time = time.time()
            
            # Receive messages for up to 30 seconds
            while time.time() - start_time < 30:
                try:
                    message = await asyncio.wait_for(websocket.recv(), timeout=5.0)
                    data = json.loads(message)
                    message_count += 1
                    
                    print(f"\n📨 Message #{message_count}:")
                    print(f"   Type: {data.get('type', 'unknown')}")
                    print(f"   Work ID: {data.get('work_id', 'unknown')}")
                    print(f"   Timestamp: {data.get('timestamp', 'unknown')}")
                    
                    payload = data.get('payload', {})
                    
                    if data.get('type') == 'snapshot':
                        print("   📊 SNAPSHOT Content:")
                        progress = payload.get('progress', {})
                        executions = payload.get('executions', [])
                        logs = payload.get('logs', {})
                        
                        print(f"      Progress:")
                        print(f"         Global Progress: {progress.get('global_progress', 'N/A')}")
                        print(f"         Global Execution: {progress.get('global_execution', {})}")
                        print(f"         Tasks: {progress.get('tasks', {})}")
                        print(f"         Current Combination: {progress.get('current_combination_details', {})}")
                        
                        print(f"      Executions: {len(executions)} items")
                        for i, exec_data in enumerate(executions[:3]):  # Show first 3
                            print(f"         [{i}] {exec_data.get('unit_id', 'N/A')}: {exec_data.get('status', 'N/A')} ({exec_data.get('progress', 0):.1%})")
                        
                        print(f"      Logs:")
                        print(f"         Errors: {len(logs.get('errors', []))}")
                        print(f"         Warnings: {len(logs.get('warnings', []))}")
                        
                    elif data.get('type') == 'update':
                        print("   📈 UPDATE Content:")
                        if 'progress' in payload:
                            print(f"      Progress Update: {payload['progress']}")
                        if 'executions_changed' in payload:
                            print(f"      Executions Changed: {payload['executions_changed']}")
                        if 'logs_appended' in payload:
                            print(f"      New Logs: {payload['logs_appended']}")
                            
                    elif data.get('type') == 'event':
                        print("   🎯 EVENT Content:")
                        event_type = payload.get('event_type', 'unknown')
                        print(f"      Event Type: {event_type}")
                        print(f"      Event Data: {payload}")
                        
                    elif data.get('type') == 'error':
                        print("   ❌ ERROR Content:")
                        print(f"      Code: {payload.get('code', 'unknown')}")
                        print(f"      Message: {payload.get('message', 'unknown')}")
                        
                    elif data.get('type') == 'heartbeat':
                        print("   💓 HEARTBEAT")
                    
                    else:
                        print("   ❓ UNKNOWN MESSAGE TYPE")
                        print("   Raw payload:")
                        pprint.pprint(payload, indent=6)
                    
                    # If we got an error, break
                    if data.get('type') == 'error':
                        break
                        
                except asyncio.TimeoutError:
                    print("   ⏰ No message in 5 seconds...")
                    # Check if we should continue
                    if message_count == 0:
                        print("   No messages received, continuing to wait...")
                        continue
                    else:
                        print("   Stopping after timeout with existing messages")
                        break
                        
            print(f"\n📊 Debug Summary:")
            print(f"   Total messages: {message_count}")
            print(f"   Duration: {time.time() - start_time:.1f}s")
            
    except Exception as e:
        print(f"❌ WebSocket monitor failed: {e}")
        import traceback
        traceback.print_exc()


async def main():
    """Main debug function."""
    print("🔍 CSPBench WebSocket Debug Tool")
    print("=" * 60)
    
    # Check if web interface is running
    try:
        import requests
        response = requests.get("http://localhost:8000/health", timeout=5)
        if response.status_code != 200:
            print("❌ Web interface not responding")
            return False
    except Exception as e:
        print("❌ Cannot connect to web interface")
        print("   Please start the web interface first:")
        print("   python main.py web")
        return False
        
    print("✅ Web interface is running")
    
    # Run debug
    success = await debug_websocket_messages()
    
    if success:
        print("\n🎉 Debug completed successfully!")
    else:
        print("\n❌ Debug failed - check logs above")
    
    return success


if __name__ == "__main__":
    # Setup detailed logging
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Run debug
    success = asyncio.run(main())
    sys.exit(0 if success else 1)
