# Closest String Problem (CSP) - Plataforma Experimental

Este projeto implementa uma arquitetura modular e extensível para experimentação com algoritmos de solução do Closest String Problem (CSP), incluindo baseline, heurísticas avançadas e métodos exatos.

## Visão Geral

A aplicação permite:
- Geração, leitura e download de datasets (sintéticos, arquivos, NCBI)
- Execução interativa ou em lote de múltiplos algoritmos
- Relatórios detalhados e comparativos automáticos
- Adição de novos algoritmos sem modificar o código principal

## Fluxo Principal da Aplicação

O fluxo principal do sistema é:
1. **Seleção e leitura/geração do dataset**: O usuário pode escolher entre gerar um dataset sintético, carregar de arquivo, baixar do NCBI ou executar em lote.
2. **Seleção dos algoritmos**: O usuário seleciona quais algoritmos deseja executar sobre o dataset.
3. **Execução dos algoritmos**: Cada algoritmo é executado múltiplas vezes, com controle de timeout e monitoramento de recursos.
4. **Exibição e salvamento dos resultados**: Os resultados são exibidos no console, salvos em relatórios detalhados e exportados em CSV.

## Estrutura do Projeto

```
csp_blfga/                    # 📦 Pacote principal
├── main.py                   # Ponto de entrada do pacote
├── ui/                       # 🖥️ Interface de usuário
│   ├── cli/                  # Interface de linha de comando
│   │   ├── app.py           # Aplicação CLI principal
│   │   ├── console_manager.py # Gerenciamento thread-safe do console
│   │   └── menu.py          # Menus interativos
│   └── widgets/             # Placeholder para futuras interfaces gráficas
├── core/                     # ⚙️ Lógica principal do sistema
│   ├── exec/                # Execução de algoritmos
│   │   ├── algorithm_executor.py # Executor com controle de recursos
│   │   ├── batch_executor.py     # Execução em lote
│   │   └── runner.py            # Controle de execução e progresso
│   ├── io/                  # Entrada/saída de dados
│   │   ├── export_csv.py        # [DEPRECIADO] Proxy para CSPExporter
│   │   ├── export_csv_batch.py  # [DEPRECIADO] Proxy para CSPExporter
│   │   ├── exporter.py          # Sistema de exportação centralizado
│   │   └── results_formatter.py # Formatação de relatórios
│   └── report/              # Geração de relatórios
│       └── report_utils.py      # Utilitários de relatórios
└── utils/                    # 🔧 Utilitários gerais
    ├── config.py            # Configurações globais
    ├── distance.py          # Funções de distância
    ├── logging.py           # Sistema de logging padronizado
    ├── resource_monitor.py  # Monitoramento de recursos
    └── resource_limits_config.py # Configuração de limites

algorithms/                   # 🧬 Implementações dos algoritmos
├── baseline/                # Algoritmo de consenso ganancioso
├── blf_ga/                  # BLF-GA: Blockwise Learning Fusion + GA
├── csc/                     # CSC: Consensus String Clustering
├── dp_csp/                  # DP-CSP: Programação dinâmica exata
├── h3_csp/                  # H³-CSP: Híbrido hierárquico
└── README.md               # Guia para adicionar novos algoritmos

datasets/                     # 📊 Gerenciamento de datasets
├── dataset_file.py          # Leitura de arquivos
├── dataset_entrez.py        # Download do NCBI
├── dataset_synthetic.py     # Geração sintética
└── dataset_utils.py         # Utilitários

tests/                        # 🧪 Testes automatizados
main.py                      # 🚀 Ponto de entrada principal
results/                     # 📈 Relatórios gerados
logs/                        # 📝 Logs de execução
saved_datasets/              # 💾 Datasets salvos
batch_configs/               # ⚙️ Configurações de lote
```

## Como Executar

### 1. Instale o Python 3.10+ (recomendado: 3.10 ou 3.11)

No Ubuntu/Debian:
```bash
sudo apt update
sudo apt install python3.10 python3.10-venv python3.10-dev
```
No Windows:
Baixe em https://www.python.org/downloads/

### 2. Crie um ambiente virtual

No terminal (Linux/macOS):
```bash
python3.10 -m venv .venv
source .venv/bin/activate
```
No Windows (cmd):
```cmd
python -m venv .venv
.venv\Scripts\activate
```

### 3. Instale as dependências

```bash
pip install --upgrade pip
pip install -r requirements.txt
```

### 4. Execute a aplicação

```bash
python main.py
```

- Siga o menu interativo para escolher dataset, algoritmos e parâmetros.
- Para execução em lote, escolha a opção correspondente e selecione um arquivo YAML/XML em `batch_configs/`.

## Adicionando Novos Algoritmos

1. Crie uma nova pasta em `algorithms/` seguindo o padrão:
    ```
    algorithms/
      meu_algoritmo/
        __init__.py
        algorithm.py
        config.py
        implementation.py
    ```
2. Implemente a interface `Algorithm` e use o decorador `@register_algorithm`.
3. O algoritmo aparecerá automaticamente no menu, sem necessidade de alterar o main.

Veja `algorithms/README.md` para exemplos detalhados.

## Requisitos

- Python 3.10+
- Biopython
- scikit-learn
- tabulate

Instale dependências com:
```bash
pip install -r requirements.txt
```

## Relatórios e Resultados

- Relatórios detalhados são salvos em `results/` após cada execução.
- Resumos rápidos são exibidos no console.
- Execuções em lote geram relatórios consolidados.

## Suporte e Documentação

- Cada algoritmo possui README próprio explicando heurísticas, funcionamento e parâmetros.
- Consulte `REESTRUTURACAO.md` para detalhes sobre a arquitetura e modularização.

---

### Observações sobre o main.py

O arquivo `main.py` está totalmente documentado com docstrings no estilo Google, detalhando o fluxo, parâmetros e retornos de cada função. Consulte o código para detalhes de uso programático e integração.
