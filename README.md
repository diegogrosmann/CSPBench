# CSPBench - Closest Str### Cloud Deployment
```bash
# Deploy seguro no Google Cloud Run (recomendado)
make deploy-cloud-run-secure PROJECT_ID=meu-projeto

# Ou com variáveis de ambiente
export NCBI_EMAIL="seu-email@exemplo.com"
export NCBI_API_KEY="sua-chave-api"
make deploy-cloud-run-env PROJECT_ID=meu-projeto

# Ver todas as opções de cloud deployment
make help-cloud
```blem Benchmark

Framework para benchmarking e análise de algoritmos do Problema da String Mais Próxima (Closest String Problem).

## 🚀 Quick Start

### Execução Local
```bash
# 1. Instalar dependências
make install

# 2. Executar interface web
make run-web
```

### Docker (Desenvolvimento)
```bash
# Executar com Docker Compose
cd deployment/docker && docker-compose up --build

# Ou usando o Makefile
make deploy-local
```

### Cloud Deployment
```bash
# Deploy no Google Cloud Run
export NCBI_EMAIL="seu-email@exemplo.com"
cd deployment/cloud-run
./deploy.sh SEU_PROJECT_ID us-central1
```

## 📁 Estrutura do Projeto

```
CSPBench/
├── src/                    # Código fonte (Arquitetura Hexagonal)
│   ├── application/        # Casos de uso e serviços
│   ├── domain/            # Entidades e regras de negócio
│   ├── infrastructure/    # Adaptadores externos
│   └── presentation/      # CLI e Web Interface
├── algorithms/            # Implementações de algoritmos CSP
├── config/               # Configurações da aplicação
├── examples/             # Exemplos de datasets e batches
├── deployment/           # 🆕 Configurações de deployment
│   ├── cloud-run/       # Google Cloud Run
│   └── docker/          # Desenvolvimento local
└── tests/               # Testes automatizados
```

## 🛠️ Comandos Disponíveis

### Desenvolvimento
```bash
make install        # Instalar dependências
make test          # Executar testes
make lint          # Verificar qualidade do código
make run-local     # Executar CLI
make run-web       # Executar interface web
```

### Docker & Deployment
```bash
make docker-dev             # Build imagem de desenvolvimento
make docker-cloud           # Build imagem de produção
make deploy-local           # Docker Compose local
make deploy-cloud-run-secure PROJECT_ID=meu-projeto  # Deploy seguro no Cloud Run
make help-cloud             # Ver todas as opções de cloud deployment
```

## 📚 Documentação

- **[Deployment Guide](deployment/README.md)** - Configurações de deployment
- **[Cloud Run](deployment/cloud-run/README.md)** - Deploy específico para Cloud Run
- **[Web Interface](src/presentation/web/README.md)** - Interface web
- **[Algorithms](algorithms/README.md)** - Algoritmos disponíveis

## 🧬 Exemplos

### CLI
```bash
# Executar algoritmo específico
python main.py run baseline examples/datasets/example.fasta

# Executar batch de experimentos
python main.py batch examples/batches/experiment_example.yaml

# Iniciar menu interativo
python main.py
```

### Web Interface
Acesse `http://localhost:8000` após executar:
```bash
make run-web
# ou
python main.py web
```

## 🔧 Configuração

### Variáveis de Ambiente
Copie `.env.example` para `.env` e configure:
```bash
# APIs externas
NCBI_EMAIL=seu-email@exemplo.com
NCBI_API_KEY=sua-chave-ncbi

# Diretórios
DATASET_DIRECTORY=./datasets
OUTPUT_BASE_DIRECTORY=./outputs
LOG_DIRECTORY=./logs

# Web Interface
WEB_HOST=0.0.0.0
PORT=8000
```

## 🏗️ Arquitetura

O CSPBench usa **Arquitetura Hexagonal** (Ports & Adapters):

- **Domain**: Entidades core e regras de negócio
- **Application**: Casos de uso e orquestração
- **Infrastructure**: Adaptadores para I/O, persistência
- **Presentation**: CLI, Web API, interfaces externas

## 📦 Deploy

### Para desenvolvimento:
```bash
make deploy-local
```

### Para produção (Cloud Run):
```bash
# Opção 1: Deploy seguro (recomendado)
make deploy-cloud-run-secure PROJECT_ID=meu-projeto

# Opção 2: Com variáveis de ambiente
export NCBI_EMAIL="seu-email@exemplo.com"
export NCBI_API_KEY="sua-chave"
make deploy-cloud-run-env PROJECT_ID=meu-projeto

# Ver todas as opções
make help-cloud
```

Veja mais detalhes em [deployment/README.md](deployment/README.md).

## 🤝 Contribuindo

Veja [CONTRIBUTING.md](CONTRIBUTING.md) para diretrizes de contribuição.

## 📄 Licença

Este projeto está sob a licença especificada em [LICENSE](LICENSE).
