import subprocess
import sys
import os
import pytest

# Teste de integração: executa o programa principal simulando entradas do usuário

def test_main_integration(monkeypatch, tmp_path):
    """
    Executa o programa principal (main.py) simulando entradas para:
    1. Habilitar modo debug? [n]
    2. Gerar dataset sintético
    3. Selecionar o primeiro algoritmo
    4. Não salvar dataset
    """
    # Entradas simuladas: [n] debug, [1] gerar sintético, [enter] para defaults, [1] primeiro algoritmo, [n] não salvar
    entradas = '\n'.join([
        'n',    # modo debug
        '1',    # menu: gerar dataset sintético
        '', '', '', '', '',  # defaults para n, L, alfabeto, noise, fully_random
        '',     # semente aleatória
        '1',    # selecionar primeiro algoritmo
        'n'     # não salvar dataset
    ]) + '\n'

    # Executa o main.py como subprocesso
    result = subprocess.run(
        [sys.executable, 'main.py'],
        input=entradas.encode(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        cwd=os.getcwd(),
        timeout=60
    )
    # Checa se rodou sem erro fatal
    assert result.returncode == 0
    # Checa se houve saída indicando execução dos algoritmos
    saida = result.stdout.decode(errors='ignore')
    assert 'Processo em execução' in saida or '[PID]' in saida
    assert (
        'Execução concluída' in saida or
        'centro' in saida or
        'distância' in saida or
        'resultado' in saida or
        'RESUMO RÁPIDO' in saida or
        'Gerando relatório detalhado' in saida
    )
