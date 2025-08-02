#!/usr/bin/env python3
"""
Teste das funcionalidades de infraestrutura em .env
"""
import os
import sys
sys.path.insert(0, '.')

def test_infrastructure():
    print("=== Teste de Infraestrutura em .env ===")
    
    # Teste 1: Variáveis de ambiente básicas
    print("\n1. Variáveis de ambiente definidas:")
    env_vars = ['EXPORT_FORMAT', 'EXECUTOR_IMPL', 'DATASET_PATH', 'OUTPUT_BASE_DIRECTORY']
    for var in env_vars:
        value = os.getenv(var, 'NOT_SET')
        print(f"   {var}: {value}")
    
    # Teste 2: Modificar temporariamente variáveis
    print("\n2. Testando override de variáveis:")
    os.environ['EXPORT_FORMAT'] = 'csv'
    os.environ['EXECUTOR_IMPL'] = 'TestExecutor'
    
    print(f"   EXPORT_FORMAT modificado para: {os.getenv('EXPORT_FORMAT')}")
    print(f"   EXECUTOR_IMPL modificado para: {os.getenv('EXECUTOR_IMPL')}")
    
    # Teste 3: Importar e testar funções
    print("\n3. Testando funções de criação:")
    try:
        from main import load_config
        config = load_config()
        print("   ✅ Configuração carregada com sucesso")
        
        from main import create_exporter, create_executor
        
        # Testar exportador
        exporter = create_exporter(config)
        print(f"   ✅ Exportador criado: {type(exporter).__name__}")
        
        # Testar executor  
        executor = create_executor(config)
        print(f"   ✅ Executor criado: {type(executor).__name__}")
        
    except Exception as e:
        print(f"   ❌ Erro: {e}")
        import traceback
        traceback.print_exc()
    
    print("\n✅ Teste concluído!")

if __name__ == "__main__":
    test_infrastructure()
