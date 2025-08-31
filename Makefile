# ===================================================================
# CSPBench Makefile
# ===================================================================
# Makefile para gerenciar build, test, deploy e operações do CSPBench

# Configurações padrão
APP_NAME = cspbench
VERSION ?= latest
REGISTRY ?= gcr.io
PROJECT_ID ?= your-gcp-project-id
IMAGE_NAME = $(REGISTRY)/$(PROJECT_ID)/$(APP_NAME)
IMAGE_TAG = $(IMAGE_NAME):$(VERSION)
CONTAINER_NAME = $(APP_NAME)-container
PORT ?= 8080
DATA_DIR ?= $(PWD)/data
CPU ?= 4
MEMORY ?= 2Gi

# Cores para output
RED = \033[0;31m
GREEN = \033[0;32m
YELLOW = \033[1;33m
BLUE = \033[0;34m
NC = \033[0m # No Color

# ===================================================================
# Help / Info
# ===================================================================

.DEFAULT_GOAL := help

help: ## Mostra esta mensagem de ajuda
	@echo "$(BLUE)===================================================================$(NC)"
	@echo "$(BLUE)                        CSPBench Makefile                        $(NC)"
	@echo "$(BLUE)===================================================================$(NC)"
	@echo ""
	@echo "$(YELLOW)Comandos disponíveis:$(NC)"
	@echo ""
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "  $(GREEN)%-20s$(NC) %s\n", $$1, $$2}' $(MAKEFILE_LIST)
	@echo ""
	@echo "$(YELLOW)Variáveis de ambiente:$(NC)"
	@echo "  $(GREEN)VERSION$(NC)     = $(VERSION)"
	@echo "  $(GREEN)PROJECT_ID$(NC)  = $(PROJECT_ID)"
	@echo "  $(GREEN)REGISTRY$(NC)    = $(REGISTRY)"
	@echo "  $(GREEN)PORT$(NC)        = $(PORT)"
	@echo "  $(GREEN)DATA_DIR$(NC)    = $(DATA_DIR)"
	@echo "  $(GREEN)CPU$(NC)         = $(CPU)"
	@echo "  $(GREEN)MEMORY$(NC)      = $(MEMORY)"

info: ## Mostra informações do ambiente
	@echo "$(BLUE)=== Informações do Sistema ===$(NC)"
	@echo "Docker version:"
	@docker --version
	@echo "Python version:"
	@python3 --version
	@echo "Make version:"
	@make --version | head -1

# ===================================================================
# Environment Setup
# ===================================================================

.env: ## Cria arquivo .env a partir do template
	@if [ ! -f .env ]; then \
		echo "$(YELLOW)Criando .env a partir de .env.example...$(NC)"; \
		cp .env.example .env; \
		echo "$(GREEN)✓ Arquivo .env criado!$(NC)"; \
		echo "$(YELLOW)⚠️  Configure as variáveis necessárias em .env$(NC)"; \
	else \
		echo "$(GREEN)✓ Arquivo .env já existe$(NC)"; \
	fi

setup: .env ## Configura ambiente de desenvolvimento
	@echo "$(BLUE)=== Configurando ambiente de desenvolvimento ===$(NC)"
	@if [ ! -d ".venv" ]; then \
		echo "$(YELLOW)Criando ambiente virtual...$(NC)"; \
		python3 -m venv .venv; \
	fi
	@echo "$(YELLOW)Instalando dependências...$(NC)"
	@.venv/bin/pip install --upgrade pip
	@.venv/bin/pip install -r requirements.txt
	@echo "$(GREEN)✓ Ambiente configurado!$(NC)"

clean-setup: ## Remove ambiente virtual e arquivos de setup
	@echo "$(YELLOW)Removendo ambiente virtual...$(NC)"
	@rm -rf .venv
	@echo "$(GREEN)✓ Ambiente limpo!$(NC)"

# ===================================================================
# Development
# ===================================================================

dev: setup ## Inicia ambiente de desenvolvimento
	@echo "$(BLUE)=== Iniciando CSPBench (desenvolvimento) ===$(NC)"
	@.venv/bin/python main.py

dev-web: setup ## Inicia interface web de desenvolvimento
	@echo "$(BLUE)=== Iniciando interface web (desenvolvimento) ===$(NC)"
	@.venv/bin/python main.py web

test: setup ## Executa testes
	@echo "$(BLUE)=== Executando testes ===$(NC)"
	@.venv/bin/python -m pytest tests/ -v

test-cov: setup ## Executa testes com cobertura
	@echo "$(BLUE)=== Executando testes com cobertura ===$(NC)"
	@.venv/bin/python -m pytest tests/ -v --cov=src --cov-report=html

format: setup ## Formata código
	@echo "$(BLUE)=== Formatando código ===$(NC)"
	@.venv/bin/python -m black .
	@.venv/bin/python -m isort .

lint: setup ## Executa linting
	@echo "$(BLUE)=== Executando linting ===$(NC)"
	@.venv/bin/python -m ruff check .

clean-dev: ## Limpa arquivos de desenvolvimento
	@echo "$(YELLOW)Limpando cache e arquivos temporários...$(NC)"
	@rm -rf __pycache__ *.pyc *.pyo */*.pyc */*.pyo */__pycache__
	@rm -rf .pytest_cache .coverage htmlcov
	@rm -rf logs outputs cache data
	@echo "$(GREEN)✓ Arquivos de desenvolvimento limpos!$(NC)"

# ===================================================================
# Docker Build
# ===================================================================

build: ## Constrói imagem Docker
	@echo "$(BLUE)=== Construindo imagem Docker ===$(NC)"
	@echo "$(YELLOW)Tag: $(IMAGE_TAG)$(NC)"
	@docker build -t $(IMAGE_TAG) .
	@docker tag $(IMAGE_TAG) $(IMAGE_NAME):latest
	@echo "$(GREEN)✓ Imagem construída: $(IMAGE_TAG)$(NC)"

build-no-cache: ## Constrói imagem Docker sem cache
	@echo "$(BLUE)=== Construindo imagem Docker (sem cache) ===$(NC)"
	@echo "$(YELLOW)Tag: $(IMAGE_TAG)$(NC)"
	@docker build --no-cache -t $(IMAGE_TAG) .
	@docker tag $(IMAGE_TAG) $(IMAGE_NAME):latest
	@echo "$(GREEN)✓ Imagem construída: $(IMAGE_TAG)$(NC)"

build-multi: ## Constrói imagem para múltiplas plataformas
	@echo "$(BLUE)=== Construindo imagem multi-plataforma ===$(NC)"
	@docker buildx build --platform linux/amd64,linux/arm64 -t $(IMAGE_TAG) --push .
	@echo "$(GREEN)✓ Imagem multi-plataforma construída e enviada!$(NC)"

# ===================================================================
# Docker Run
# ===================================================================

run: build ## Executa container localmente
	@echo "$(BLUE)=== Executando container ===$(NC)"
	@mkdir -p $(DATA_DIR)
	@docker run --rm -it \
		--name $(CONTAINER_NAME) \
		-p $(PORT):8080 \
		-v $(DATA_DIR):/data \
		-e PORT=8080 \
		$(IMAGE_TAG)

run-detached: build ## Executa container em background
	@echo "$(BLUE)=== Executando container em background ===$(NC)"
	@mkdir -p $(DATA_DIR)
	@docker run -d \
		--name $(CONTAINER_NAME) \
		-p $(PORT):8080 \
		-v $(DATA_DIR):/data \
		-e PORT=8080 \
		$(IMAGE_TAG)
	@echo "$(GREEN)✓ Container executando em background$(NC)"
	@echo "$(YELLOW)URL: http://localhost:$(PORT)$(NC)"

run-shell: build ## Executa container com shell interativo
	@echo "$(BLUE)=== Executando container com shell ===$(NC)"
	@mkdir -p $(DATA_DIR)
	@docker run --rm -it \
		--name $(CONTAINER_NAME)-shell \
		-v $(DATA_DIR):/data \
		--entrypoint /bin/bash \
		$(IMAGE_TAG)

stop: ## Para container em execução
	@echo "$(YELLOW)Parando container...$(NC)"
	@docker stop $(CONTAINER_NAME) || true
	@docker rm $(CONTAINER_NAME) || true
	@echo "$(GREEN)✓ Container parado$(NC)"

logs: ## Mostra logs do container
	@docker logs -f $(CONTAINER_NAME)

# ===================================================================
# Docker Registry
# ===================================================================

login: ## Faz login no registry (GCR)
	@echo "$(BLUE)=== Fazendo login no registry ===$(NC)"
	@gcloud auth configure-docker

push: build ## Envia imagem para registry
	@echo "$(BLUE)=== Enviando imagem para registry ===$(NC)"
	@docker push $(IMAGE_TAG)
	@docker push $(IMAGE_NAME):latest
	@echo "$(GREEN)✓ Imagem enviada: $(IMAGE_TAG)$(NC)"

pull: ## Baixa imagem do registry
	@echo "$(BLUE)=== Baixando imagem do registry ===$(NC)"
	@docker pull $(IMAGE_TAG)
	@echo "$(GREEN)✓ Imagem baixada: $(IMAGE_TAG)$(NC)"

# ===================================================================
# Cloud Deployment
# ===================================================================

deploy-cloud-run: ## Deploy no Google Cloud Run
	@echo "$(BLUE)=== Deploy no Google Cloud Run ===$(NC)"
	@gcloud run deploy $(APP_NAME) \
		--image $(IMAGE_TAG) \
		--platform managed \
		--region us-central1 \
		--allow-unauthenticated \
		--memory $(MEMORY) \
		--cpu $(CPU) \
		--max-instances 2 \
		--add-volume name=shared,type=cloud-storage,bucket=csp-bench \
		--add-volume-mount volume=shared,mount-path=/data
	@echo "$(GREEN)✓ Deploy realizado!$(NC)"

deploy-update: push ## Atualiza deployment no Cloud Run
	@echo "$(BLUE)=== Atualizando deployment ===$(NC)"
	@gcloud run services update $(APP_NAME) \
		--image $(IMAGE_TAG) \
		--region us-central1
	@echo "$(GREEN)✓ Deployment atualizado!$(NC)"

# ===================================================================
# Monitoring & Maintenance
# ===================================================================

health: ## Verifica saúde da aplicação
	@echo "$(BLUE)=== Verificando saúde da aplicação ===$(NC)"
	@curl -f http://localhost:$(PORT)/health || echo "$(RED)❌ Aplicação não está saudável$(NC)"

ps: ## Lista containers relacionados
	@echo "$(BLUE)=== Containers CSPBench ===$(NC)"
	@docker ps -a --filter name=$(APP_NAME)

images: ## Lista imagens relacionadas
	@echo "$(BLUE)=== Imagens CSPBench ===$(NC)"
	@docker images | grep $(APP_NAME) || echo "Nenhuma imagem encontrada"

stats: ## Mostra estatísticas do container
	@docker stats $(CONTAINER_NAME)

# ===================================================================
# Cleanup
# ===================================================================

clean: stop ## Limpa containers e imagens locais
	@echo "$(YELLOW)Limpando containers e imagens...$(NC)"
	@docker container prune -f
	@docker rmi $(IMAGE_TAG) $(IMAGE_NAME):latest || true
	@echo "$(GREEN)✓ Cleanup realizado!$(NC)"

clean-all: clean clean-dev ## Limpa tudo (containers, imagens, cache)
	@echo "$(YELLOW)Limpando sistema Docker...$(NC)"
	@docker system prune -f
	@echo "$(GREEN)✓ Cleanup completo realizado!$(NC)"

# ===================================================================
# Utilities
# ===================================================================

size: ## Mostra tamanho da imagem
	@echo "$(BLUE)=== Tamanho da imagem ===$(NC)"
	@docker images $(IMAGE_NAME) --format "table {{.Repository}}\t{{.Tag}}\t{{.Size}}"

inspect: ## Inspeciona imagem Docker
	@echo "$(BLUE)=== Inspecionando imagem ===$(NC)"
	@docker inspect $(IMAGE_TAG)

dive: ## Analisa camadas da imagem (requer dive)
	@echo "$(BLUE)=== Analisando camadas da imagem ===$(NC)"
	@dive $(IMAGE_TAG)

# ===================================================================
# Special Targets
# ===================================================================

.PHONY: help info setup clean-setup dev dev-web test test-cov format lint clean-dev \
	build build-no-cache build-multi run run-detached run-shell stop logs \
	login push pull deploy-cloud-run deploy-update health ps images stats \
	clean clean-all size inspect dive .env

# Previne a execução se variáveis críticas não estiverem definidas
check-project-id:
	@if [ "$(PROJECT_ID)" = "your-gcp-project-id" ]; then \
		echo "$(RED)❌ Configure PROJECT_ID antes de usar comandos de deploy!$(NC)"; \
		echo "$(YELLOW)Exemplo: make deploy-cloud-run PROJECT_ID=meu-projeto$(NC)"; \
		exit 1; \
	fi

deploy-cloud-run deploy-update push: check-project-id
