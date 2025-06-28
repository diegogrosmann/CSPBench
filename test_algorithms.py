#!/usr/bin/env python3
"""
Teste simples da estrutura de algoritmos registrados.

Funções:
    test_algorithms(): Executa todos os algoritmos do registry em um dataset de teste.
"""
from algorithms.base import global_registry
import algorithms

def test_algorithms():
    print("=== Teste da Estrutura de Algoritmos ===")
    print(f"Algoritmos registrados: {sorted(global_registry.keys())}")
    
    # Criar dataset pequeno para teste
    test_strings = ["ACGT", "ACTT", "AGGT"]
    alphabet = "ACGT"
    
    print(f"\nDataset de teste: {test_strings}")
    print(f"Alphabet: {alphabet}")
    
    # Testar cada algoritmo
    for alg_name, AlgClass in global_registry.items():
        print(f"\n--- Testando {alg_name} ---")
        try:
            alg = AlgClass(test_strings, alphabet)
            center, dist = alg.run()
            print(f"Centro: {center}")
            print(f"Distância: {dist}")
            print(f"✓ {alg_name} funcionou!")
        except Exception as e:
            print(f"✗ {alg_name} falhou: {e}")

if __name__ == "__main__":
    test_algorithms()
