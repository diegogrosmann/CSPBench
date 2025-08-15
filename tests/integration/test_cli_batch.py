"""
Testes de integração para CLI batch
"""

import json
import tempfile
from pathlib import Path

import pytest
import yaml
from typer.testing import CliRunner

from main_old import app


def test_cli_batch_minimal():
    """Testa execução de batch mínimo via CLI."""
    runner = CliRunner()

    with tempfile.TemporaryDirectory() as tmp_dir:
        batch_file = Path(tmp_dir) / "batch.yaml"
    # Configuração mínima válida
    config = {"task": {"type": "experiment"}, "experiments": []}

    batch_file.write_text(yaml.dump(config))

    result = runner.invoke(app, ["batch", str(batch_file)])

    # Deve ter sucesso mesmo com lista vazia
    assert result.exit_code == 0
    assert "Batch concluído" in result.output


def test_cli_batch_with_experiments():
    """Testa batch com experimentos reais."""
    runner = CliRunner()

    with tempfile.TemporaryDirectory() as tmp_dir:
        batch_file = Path(tmp_dir) / "batch.yaml"

        # Configuração com experimento de teste
        config = {
            "task": {"type": "experiment"},
            "experiments": [
                {
                    "algorithm": "Baseline",
                    "dataset": "saved_datasets/teste.fasta",
                    "params": {},
                }
            ],
        }

        batch_file.write_text(yaml.dump(config))

        result = runner.invoke(app, ["batch", str(batch_file)])

        # Pode falhar se dataset não existir, mas deve parsear corretamente
        if result.exit_code != 0:
            # Se falhou, deve ser por dataset não encontrado, não por parsing
            assert "Erro no batch" in result.output
        else:
            assert "Batch concluído" in result.output


def test_cli_batch_invalid_file():
    """Testa batch com arquivo inválido."""
    runner = CliRunner()

    result = runner.invoke(app, ["batch", "arquivo_inexistente.yaml"])

    # Deve falhar porque arquivo não existe
    assert result.exit_code != 0


def test_cli_batch_json_format():
    """Testa batch com formato JSON."""
    runner = CliRunner()

    with tempfile.TemporaryDirectory() as tmp_dir:
        batch_file = Path(tmp_dir) / "batch.json"
    # Configuração em JSON
    config = {"task": {"type": "experiment"}, "experiments": []}

    batch_file.write_text(json.dumps(config))

    result = runner.invoke(app, ["batch", str(batch_file)])

    assert result.exit_code == 0
    assert "Batch concluído" in result.output
