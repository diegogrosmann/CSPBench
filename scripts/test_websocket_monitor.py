#!/usr/bin/env python3
"""
Manual WebSocket Monitor Test Script

This script tests the WebSocket monitor with real batch executions.
It starts a simple batch and monitors it via WebSocket to validate the implementation.
"""

import asyncio
import json
import logging
import sys
import time
import websockets
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from src.application.services.execution_manager import ExecutionManager
from src.domain.config import CSPBenchConfig


async def test_websocket_monitor():
    """Test WebSocket monitoring with real execution."""
    
    print("🚀 Starting WebSocket Monitor Test")
    print("=" * 50)
    
    # Step 1: Start a test batch
    print("1. Creating test batch...")
    
    execution_manager = ExecutionManager()
    
    # Simple test batch configuration dict
    test_batch_dict = {
        "metadata": {
            "name": "WebSocket Test Batch",
            "description": "Simple batch for WebSocket testing",
            "author": "Test",
            "version": "1.0",
            "creation_date": "2025-08-26",
            "tags": ["test", "websocket"]
        },
        "datasets": [
            {
                "id": "test_dataset",
                "name": "Test Dataset",
                "type": "synthetic",
                "sequences": 3,  # Very small for quick testing
                "length": 10,
                "alphabet": "ATCG",
                "seed": 42
            }
        ],
        "algorithms": [
            {
                "id": "quick_test",
                "name": "Quick Test Config",
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
                    "id": "websocket_test",
                    "name": "WebSocket Test Task",
                    "datasets": ["test_dataset"], 
                    "algorithms": ["quick_test"],
                    "repetitions": 2,
                    "timeout": 30
                }
            ]
        },
        "output": {
            "logging": True,
            "results": {"partial_results": True}
        },
        "resources": {
            "cpu": {"max_cores": 1},
            "timeouts": {"timeout_per_item": 30}
        }
    }
    
    try:
        # Save test batch to temporary file and load with CSPBench
        import tempfile
        import yaml
        from src.domain.config import load_cspbench_config
        
        # Create temporary file with test batch
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(test_batch_dict, f, default_flow_style=False)
            temp_batch_file = f.name
        
        try:
            # Load configuration
            config = load_cspbench_config(Path(temp_batch_file))
            
            # Execute the batch
            work_id = execution_manager.execute(config=config, extra={"origin": "websocket_test"})
            print(f"✅ Created work: {work_id}")
        finally:
            # Clean up temp file
            Path(temp_batch_file).unlink(missing_ok=True)
        
        # Wait a moment for work to initialize
        await asyncio.sleep(2)
        
        # Get work details to find database path
        from src.application.services.work_service import get_work_service
        work_service = get_work_service()
        work_details = work_service.get(work_id)
        if not work_details:
            print("❌ Failed to get work details")
            return False
            
        db_path = Path(work_details["output_path"]) / "state.db"
        print(f"📁 Database path: {db_path}")
        
        # Step 2: Work execution starts automatically
        print("\n2. Work execution started automatically...")
        print("   (ExecutionManager runs work in background thread)")
        
        # Wait for database to be created
        max_wait = 30
        wait_time = 0
        while not db_path.exists() and wait_time < max_wait:
            await asyncio.sleep(1)
            wait_time += 1
            
        if not db_path.exists():
            print("❌ Database not created within timeout")
            return False
            
        print("✅ Database created, starting WebSocket test...")
        
        # Step 3: Connect to WebSocket and monitor
        print("\n3. Testing WebSocket connection...")
        
        websocket_url = f"ws://localhost:8000/ws/work/{work_id}"
        print(f"🔗 Connecting to: {websocket_url}")
        
        async with websockets.connect(websocket_url) as websocket:
            print("✅ WebSocket connected!")
            
            message_count = 0
            start_time = time.time()
            last_progress = 0
            
            try:
                while time.time() - start_time < 60:  # Max 1 minute test
                    try:
                        # Wait for message with timeout
                        message = await asyncio.wait_for(
                            websocket.recv(), 
                            timeout=10.0
                        )
                        
                        data = json.loads(message)
                        message_count += 1
                        
                        print(f"\n📨 Message #{message_count} ({data['type']}):")
                        
                        if data['type'] == 'snapshot':
                            payload = data['payload']
                            progress = payload['progress']
                            print(f"   📊 Global Progress: {progress['global_progress']:.1%}")
                            print(f"   📈 Execution: {progress['global_execution']['Finished']}/{progress['global_execution']['Total']}")
                            
                            if payload['executions']:
                                print(f"   🔄 Running Executions: {len(payload['executions'])}")
                                for exec_data in payload['executions'][:3]:  # Show first 3
                                    print(f"      - {exec_data['unit_id']}: {exec_data['status']} ({exec_data['progress']:.1%})")
                                    
                            last_progress = progress['global_progress']
                            
                        elif data['type'] == 'update':
                            payload = data['payload']
                            if 'progress' in payload:
                                progress = payload['progress']['global_progress']
                                print(f"   📊 Progress Update: {progress:.1%}")
                                last_progress = progress
                                
                            if 'executions_changed' in payload:
                                changes = payload['executions_changed']
                                if changes.get('completed'):
                                    print(f"   ✅ Completed: {changes['completed']}")
                                if changes.get('new'):
                                    print(f"   🆕 New: {changes['new']}")
                                    
                        elif data['type'] == 'event':
                            event = data['payload']
                            print(f"   🎯 Event: {event['event_type']}")
                            if event['event_type'] == 'work_status_changed':
                                print(f"      Status: {event['old_status']} → {event['new_status']}")
                                
                        elif data['type'] == 'error':
                            error = data['payload']
                            print(f"   ❌ Error: {error['code']} - {error['message']}")
                            
                        # Check if work completed
                        if last_progress >= 1.0:
                            print("\n🎉 Work completed, monitoring for final messages...")
                            await asyncio.sleep(5)  # Wait for final messages
                            break
                            
                    except asyncio.TimeoutError:
                        print("⏰ No message received in 10 seconds, checking work status...")
                        
                        # Check if work is still running
                        try:
                            work_status = work_service.get_status(work_id)
                            print(f"   Work status: {work_status}")
                            
                            if work_status in ["completed", "failed", "cancelled"]:
                                print("   Work finished, ending test")
                                break
                        except Exception as e:
                            print(f"   Error checking work status: {e}")
                            
            except Exception as e:
                print(f"❌ WebSocket error: {e}")
                return False
                
        print(f"\n📊 Test Summary:")
        print(f"   Messages received: {message_count}")
        print(f"   Final progress: {last_progress:.1%}")
        print(f"   Duration: {time.time() - start_time:.1f}s")
        
        # Step 4: Cleanup
        print("\n4. Cleaning up...")
        try:
            work_service.cancel(work_id)
            print("✅ Work cancelled")
        except Exception as e:
            print(f"⚠️  Cleanup warning: {e}")
            
        return True
        
    except Exception as e:
        print(f"❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


async def test_general_websocket():
    """Test general WebSocket endpoint."""
    print("\n🔧 Testing General WebSocket...")
    
    try:
        websocket_url = "ws://localhost:8000/ws/test_client_123"
        print(f"🔗 Connecting to: {websocket_url}")
        
        async with websockets.connect(websocket_url) as websocket:
            print("✅ Connected!")
            
            # Should receive welcome message
            message = await asyncio.wait_for(websocket.recv(), timeout=5.0)
            data = json.loads(message)
            
            print(f"📨 Received: {data['type']}")
            assert data['type'] == 'welcome'
            assert data['client_id'] == 'test_client_123'
            
            # Test echo
            await websocket.send("Hello WebSocket!")
            echo_message = await asyncio.wait_for(websocket.recv(), timeout=5.0)
            echo_data = json.loads(echo_message)
            
            print(f"📨 Echo: {echo_data['payload']['received']}")
            assert echo_data['type'] == 'echo'
            assert echo_data['payload']['received'] == "Hello WebSocket!"
            
        print("✅ General WebSocket test passed!")
        return True
        
    except Exception as e:
        print(f"❌ General WebSocket test failed: {e}")
        return False


async def main():
    """Main test function."""
    print("🧪 CSPBench WebSocket Monitor Test Suite")
    print("=" * 60)
    
    # Check if web interface is running
    try:
        import requests
        response = requests.get("http://localhost:8000/health", timeout=5)
        if response.status_code != 200:
            print("❌ Web interface not responding")
            print("   Please start the web interface first:")
            print("   python main.py web")
            return False
    except Exception as e:
        print("❌ Cannot connect to web interface")
        print("   Please start the web interface first:")
        print("   python main.py web")
        return False
        
    print("✅ Web interface is running")
    
    # Run tests
    tests_passed = 0
    total_tests = 2
    
    # Test 1: General WebSocket
    if await test_general_websocket():
        tests_passed += 1
    
    # Test 2: Work monitoring WebSocket
    print("\n" + "=" * 60)
    if await test_websocket_monitor():
        tests_passed += 1
    
    # Summary
    print("\n" + "=" * 60)
    print(f"🏁 Test Results: {tests_passed}/{total_tests} passed")
    
    if tests_passed == total_tests:
        print("🎉 All tests passed! WebSocket monitor is working correctly.")
        return True
    else:
        print("❌ Some tests failed. Check the logs above.")
        return False


if __name__ == "__main__":
    # Setup logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Run tests
    success = asyncio.run(main())
    sys.exit(0 if success else 1)
