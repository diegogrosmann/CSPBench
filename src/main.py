#!/usr/bin/env python3
"""
Ponto de Entrada Principal para a Aplicação CSP-BLFGA

Este script serve como o ponto de entrada central para a execução do projeto
a partir da linha de comando.

Função Principal:
- Inicia a interface de linha de comando (CLI) da aplicação.
- Invoca a função `main()` do módulo `src.ui.cli.app`, que é responsável
  por processar os argumentos, configurar e executar os algoritmos.

Como Usar:
- Para executar a aplicação, rode este script diretamente:
  ```bash
  python -m src.main
  ```
- Ou, se o pacote estiver instalado:
  ```bash
  csp-blfga --help
  ```

O design modular garante que este arquivo permaneça simples, delegando toda
a lógica complexa para os módulos apropriados dentro do diretório `src`.
"""

from src.ui.cli.app import main

if __name__ == "__main__":
    # Chama a função principal da aplicação CLI
    main()
