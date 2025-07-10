#!/usr/bin/env python3
"""
Ponto de Entrada Principal do CSPBench

Este é o ponto de entrada principal do CSPBench (Closest String Problem Benchmark),
um framework experimental robusto para teste, comparação e análise de algoritmos
do Closest String Problem (CSP).

Funcionalidade:
    - Delegação para o sistema principal na pasta src/
    - Compatibilidade com execução direta via python main.py
    - Suporte a argumentos de linha de comando
    - Configuração de ambiente básica

Uso:
    ```bash
    # Execução interativa
    python main.py

    # Execução silenciosa
    python main.py --silent --dataset synthetic --algorithms BLF-GA

    # Execução em lote
    python main.py --batch config.yaml

    # Ajuda
    python main.py --help
    ```

Arquitetura:
    Este arquivo mantém a compatibilidade com execução existente,
    delegando toda a funcionalidade para src.ui.cli.app.main()
    que contém a implementação completa do sistema.

Estrutura do Sistema:
    - main.py: Ponto de entrada (este arquivo)
    - src/: Implementação principal do framework
    - algorithms/: Biblioteca de algoritmos CSP
    - batch_configs/: Configurações de execução em lote
    - tests/: Testes automatizados
    - docs/: Documentação completa

Autor: CSPBench Development Team
Data: 2024
Versão: 1.0.0
"""

from src.ui.cli.app import main

if __name__ == "__main__":
    main()
