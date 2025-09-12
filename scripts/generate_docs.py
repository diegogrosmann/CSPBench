#!/usr/bin/env python3
"""
Gera documentação HTML a partir dos docstrings usando pdoc.

Saída padrão: docs/api
Gera documentação para os pacotes em `src/` e `algorithms/`.

Uso:
  .venv/bin/python scripts/generate_docs.py [--serve] [--port 8081]

Requisitos:
  - pdoc instalado no ambiente virtual
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT / "src"
ALG_DIR = ROOT / "algorithms"
OUTPUT_DIR = ROOT / "docs" / "api"


def run_pdoc(serve: bool = False, port: int = 8081) -> int:
    env = os.environ.copy()
    env.setdefault("PYTHONPATH", str(ROOT))
    cmd = [
        str(ROOT / ".venv" / "bin" / "python"),
        "-m",
        "pdoc",
        "--docformat",
        "google",
        "--footer-text",
        "CSPBench API Docs",
        "--logo",
        "https://raw.githubusercontent.com/diegogrosmann/csp-blfga/main/docs/logo.png",
    ]

    # Targets: todos os módulos principais em src e plugins em algorithms
    targets = []
    if SRC_DIR.exists():
        targets.append(str(SRC_DIR))
    if ALG_DIR.exists():
        targets.append(str(ALG_DIR))

    if not targets:
        print("Nada para documentar: diretórios 'src' e 'algorithms' não encontrados.")
        return 2

    if serve:
        cmd += ["--http", f"0.0.0.0:{port}"]
        cmd += targets
    else:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        cmd += ["-o", str(OUTPUT_DIR)] + targets

    print("Executando:", " ".join(cmd))
    proc = subprocess.run(cmd, env=env)
    return proc.returncode


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Gerar documentação com pdoc")
    parser.add_argument("--serve", action="store_true", help="Servir docs via HTTP")
    parser.add_argument("--port", type=int, default=8081, help="Porta para --serve")
    args = parser.parse_args(argv)

    return run_pdoc(serve=args.serve, port=args.port)


if __name__ == "__main__":
    raise SystemExit(main())
