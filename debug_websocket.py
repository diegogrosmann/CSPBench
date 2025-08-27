#!/usr/bin/env python3
"""
Script para debug de WebSocket - simula mensagens que o frontend recebe
"""
import json
import asyncio
import websockets
import sys
import time

async def test_websocket_client(work_id):
    uri = f"ws://localhost:8000/ws/work/{work_id}"
    print(f"Conectando ao WebSocket: {uri}")
    
    try:
        async with websockets.connect(uri) as websocket:
            print("✅ Conectado ao WebSocket")
            
            # Aguardar mensagens
            timeout_count = 0
            while timeout_count < 20:  # 20 tentativas de 2 segundos cada
                try:
                    message = await asyncio.wait_for(websocket.recv(), timeout=2.0)
                    data = json.loads(message)
                    
                    print(f"\n📨 Mensagem recebida:")
                    print(f"   Tipo: {data.get('type', 'unknown')}")
                    
                    if data.get('type') == 'error':
                        print(f"   ❌ Erro: {data.get('message', 'Erro desconhecido')}")
                        print(f"   📄 Dados completos: {json.dumps(data, indent=2)}")
                    
                    elif data.get('type') == 'event':
                        print(f"   📅 Evento: {data.get('event', {}).get('event_type', 'Desconhecido')}")
                        if data.get('event', {}).get('event_type') == 'WORK_STATUS_CHANGED':
                            old_status = data['event'].get('old_status')
                            new_status = data['event'].get('new_status')
                            print(f"      🔄 Status mudou: {old_status} → {new_status}")
                        print(f"   📄 Dados do evento: {json.dumps(data, indent=2)}")
                    
                    if data.get('payload'):
                        payload = data['payload']
                        if 'progress' in payload:
                            print(f"   Progress: ✅")
                        if 'executions' in payload:
                            executions = payload['executions']
                            print(f"   Executions: {len(executions)} execuções")
                            
                            # Agrupar por combination_id
                            by_combo = {}
                            for exec in executions:
                                combo_id = exec.get('combination_id', 'null')
                                if combo_id not in by_combo:
                                    by_combo[combo_id] = []
                                by_combo[combo_id].append(exec)
                            
                            for combo_id, execs in by_combo.items():
                                print(f"     Combo {combo_id}: {len(execs)} execuções")
                                statuses = {}
                                for exec in execs:
                                    status = exec.get('status', 'unknown')
                                    statuses[status] = statuses.get(status, 0) + 1
                                print(f"       Status: {statuses}")
                        
                        if 'combinations' in payload:
                            combinations = payload['combinations']
                            print(f"   Combinations: {len(combinations)} combinações")
                    
                    # Para mensagens de update, verificar estruturas diferentes
                    if data.get('type') == 'update':
                        if 'executions' in data:
                            executions = data['executions']
                            print(f"   Update Executions (root): {len(executions)} execuções")
                        if data.get('update') and 'executions' in data['update']:
                            executions = data['update']['executions']
                            print(f"   Update Executions (nested): {len(executions)} execuções")
                    
                    timeout_count = 0  # Reset timeout counter on successful message
                    
                except asyncio.TimeoutError:
                    timeout_count += 1
                    print(f"⏱️  Timeout {timeout_count}/20 - aguardando mensagens...")
                    
            print("🔚 Encerrando teste (timeout)")
            
    except Exception as e:
        print(f"❌ Erro: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uso: python debug_websocket.py <work_id>")
        sys.exit(1)
    
    work_id = sys.argv[1]
    print(f"🔍 Testando WebSocket para work_id: {work_id}")
    
    asyncio.run(test_websocket_client(work_id))
