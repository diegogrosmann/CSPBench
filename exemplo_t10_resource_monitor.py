#!/usr/bin/env python3
"""
Exemplo de teste do sistema de monitoramento de recursos T10.

Este script demonstra as melhorias:
- T10-1: Leitura de defaults de resource_limits_config
- T10-2: Coleta GC em cada check
- T10-3: Limpeza do gc_counter
"""

import sys
import time
from pathlib import Path

# Adicionar diretório raiz ao path
sys.path.insert(0, str(Path(__file__).parent.parent))

from csp_blfga.utils.resource_limits_config import (
    get_merged_gc_config,
    get_merged_memory_safety,
    get_merged_resource_limits,
)
from csp_blfga.utils.resource_monitor import (
    ResourceLimits,
    ResourceMonitor,
    force_garbage_collection,
)


def test_resource_config():
    """Testa T10-1: Leitura de configurações."""
    print("🔧 T10-1: Testando leitura de configurações")
    print("=" * 50)

    # Mostrar configurações básicas
    limits = get_merged_resource_limits()
    print("Limites de recursos:")
    for key, value in limits.items():
        print(f"  {key}: {value}")

    # Mostrar configurações de memória
    memory_config = get_merged_memory_safety()
    print("\nConfiguração de segurança de memória:")
    for key, value in memory_config.items():
        print(f"  {key}: {value}")

    # Mostrar configurações de GC
    gc_config = get_merged_gc_config()
    print("\nConfiguração de garbage collection:")
    for key, value in gc_config.items():
        print(f"  {key}: {value}")

    print("\n✅ T10-1: Configurações carregadas com sucesso!")


def test_gc_improvements():
    """Testa T10-2 e T10-3: Melhorias no GC."""
    print("\n🗑️ T10-2 e T10-3: Testando melhorias no GC")
    print("=" * 50)

    # Criar monitor com configurações
    limits = ResourceLimits.from_config()
    monitor = ResourceMonitor(limits)

    print(f"GC automático: {limits.gc_auto_collect}")
    print(f"Frequência do GC: {limits.gc_frequency}")
    print(f"GC forçado no limite: {limits.gc_force_on_limit}")
    print(f"Ratio do GC: {limits.gc_threshold_ratio}")

    # Teste de GC manual
    print("\nTestando força garbage collection...")
    force_garbage_collection()

    # Criar alguns objetos para GC
    print("\nCriando objetos temporários...")
    temp_objects = []
    for i in range(10000):
        temp_objects.append([i] * 100)

    print("Objetos criados, executando GC...")
    force_garbage_collection()

    # Limpar referências
    del temp_objects
    force_garbage_collection()

    # Testar estatísticas de memória
    stats = monitor.get_memory_stats()
    print("\nEstatísticas de memória:")
    for key, value in stats.items():
        if isinstance(value, float):
            print(f"  {key}: {value:.2f}")
        else:
            print(f"  {key}: {value}")

    print("\n✅ T10-2 e T10-3: Melhorias no GC testadas!")


def test_monitoring_loop():
    """Testa o loop de monitoramento melhorado."""
    print("\n📊 Testando loop de monitoramento")
    print("=" * 50)

    # Configurar limites baixos para teste
    limits = ResourceLimits.from_config()
    limits.check_interval = 1.0  # Verificar a cada 1 segundo
    limits.gc_frequency = 3  # GC a cada 3 checks

    monitor = ResourceMonitor(limits)

    # Callback para violações
    def on_violation(msg):
        print(f"⚠️ Violação detectada: {msg}")

    monitor.set_violation_callback(on_violation)

    print("Iniciando monitoramento por 5 segundos...")
    monitor.start_monitoring()

    try:
        # Aguardar e mostrar estatísticas
        for i in range(5):
            time.sleep(1)
            stats = monitor.get_memory_stats()
            print(
                f"Segundo {i+1}: Memória={stats['current_memory_mb']:.1f}MB, GC_counter={stats['gc_counter']}"
            )

    finally:
        monitor.stop_monitoring()

    print("✅ Loop de monitoramento testado!")


def main():
    """Executa todos os testes T10."""
    print("🚀 Testes T10 - Melhorias no ResourceMonitor")
    print("=" * 60)

    try:
        # T10-1: Configurações
        test_resource_config()

        # T10-2 e T10-3: GC
        test_gc_improvements()

        # Teste completo
        test_monitoring_loop()

        print("\n" + "=" * 60)
        print("🎯 Todos os testes T10 executados com sucesso!")
        print("✅ T10-1: Leitura de configurações aprimorada")
        print("✅ T10-2: Coleta GC em cada check")
        print("✅ T10-3: Limpeza do gc_counter")

    except Exception as e:
        print(f"\n❌ Erro durante os testes: {e}")
        import traceback

        traceback.print_exc()


if __name__ == "__main__":
    main()
