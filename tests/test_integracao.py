import os
import subprocess
import sys
from pathlib import Path

# Teste de integração: executa o programa principal simulando entradas do usuário


def test_main_integration(monkeypatch, tmp_path):
    """
    Executa o programa principal (main.py) em modo silencioso
    para testar a integração completa sem necessidade de inputs.
    """
    # Obter diretório do projeto
    project_root = Path(__file__).parent.parent

    # Executa o main.py em modo silencioso
    env = os.environ.copy()
    env["CSP_AUTOMATED_TEST"] = "1"
    result = subprocess.run(
        [sys.executable, "main.py", "--silent"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        cwd=str(project_root),
        timeout=60,
        env=env,
    )

    # Checa se rodou sem erro fatal
    # Note: main.py já faz sys.exit(0) ao final quando bem-sucedido
    saida = result.stdout.decode(errors="ignore")
    stderr_saida = result.stderr.decode(errors="ignore")

    # Em caso de erro, imprimir para debug
    if result.returncode != 0:
        print(f"STDOUT:\n{saida}")
        print(f"STDERR:\n{stderr_saida}")

    assert result.returncode == 0
    # Checa se houve saída indicando execução dos algoritmos
    saida = result.stdout.decode(errors="ignore")
    stderr_saida = result.stderr.decode(errors="ignore")
    # No modo silencioso, apenas verifica se não houve erro crítico
    assert "ERRO FATAL" not in saida
    assert "ERRO FATAL" not in stderr_saida
    # Se houver alguma saída, verifica se não há erros críticos
    if saida.strip():
        assert "Traceback" not in saida
    if stderr_saida.strip():
        assert "Traceback" not in stderr_saida
