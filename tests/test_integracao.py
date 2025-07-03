import subprocess
import sys
import os
import pytest

# Teste de integração: executa o programa principal simulando entradas do usuário

def test_main_integration(monkeypatch, tmp_path):
    """
    Executa o programa principal (main.py) em modo silencioso
    para testar a integração completa sem necessidade de inputs.
    """
    # Executa o main.py em modo silencioso
    env = os.environ.copy()
    env["CSP_AUTOMATED_TEST"] = "1"
    result = subprocess.run(
        [sys.executable, 'main.py', '--silent'],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        cwd=os.getcwd(),
        timeout=60,
        env=env
    )
    # Checa se rodou sem erro fatal
    assert result.returncode == 0
    # Checa se houve saída indicando execução dos algoritmos
    saida = result.stdout.decode(errors='ignore')
    # No modo silencioso, procura por indicadores de execução do algoritmo
    assert (
        'BLF-GA' in saida or
        'Criando população inicial' in saida or
        'Geração' in saida or
        'melhor=' in saida or
        'RESUMO RÁPIDO' in saida or
        'Algoritmo' in saida or
        'Melhor Dist' in saida or
        'tempo=' in saida or
        'dist=' in saida
    )
