# CSPBench - Unified Makefile
# Development, testing, and deployment automation

.PHONY: help install clean test docker deploy

# Default help target
help:
	@echo "CSPBench - Available targets:"
	@echo "  install     - Setup Python environment and dependencies"
	@echo "  clean       - Clean temporary files and caches"
	@echo "  test        - Run tests with coverage"
	@echo "  lint        - Code quality checks (black, ruff, mypy)"
	@echo "  docker      - Build Docker image"
	@echo "  run-local   - Run locally with development settings"
	@echo "  run-web     - Start web interface"
	@echo "  deploy      - Deploy to cloud"

# Python Environment Setup
install:
	python3 -m venv .venv
	.venv/bin/pip install --upgrade pip
	.venv/bin/pip install -r requirements.txt

# Development
clean:
	find . -type d -name "__pycache__" -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
	rm -rf .pytest_cache htmlcov .coverage .mypy_cache

test:
	.venv/bin/pytest tests/ -v --cov=src --cov-report=html

lint:
	.venv/bin/black src/ tests/ algorithms/
	.venv/bin/ruff check src/ tests/ algorithms/
	.venv/bin/mypy src/

# Local execution
run-local:
	.venv/bin/python main.py

run-web:
	.venv/bin/uvicorn src.presentation.web.app:app --reload --host 0.0.0.0 --port 8000

# Docker operations
docker:
	docker build -t cspbench:latest .

docker-run:
	docker run --rm -it --env-file .env -p 8000:8000 \
		-v "$(PWD)/datasets:/app/datasets:ro" \
		-v "$(PWD)/config:/app/config:ro" \
		-v "$(PWD)/outputs:/app/outputs" \
		cspbench:latest

# Cloud deployment
deploy:
	./deploy/deploy-cloud.sh

# Docker Compose operations
compose-up:
	docker-compose -f deploy/docker-compose.yml up --build

compose-down:
	docker-compose -f deploy/docker-compose.yml down

web-up-detached:
	docker-compose up --build -d

web-down:
	docker-compose down

web-logs:
	docker-compose logs -f

web-restart:
	docker-compose restart

# Web Interface - Local Development
web-dev:
	.venv/bin/python main.py --web

web-dev-reload:
	.venv/bin/uvicorn src.presentation.web.app:app --reload --host 0.0.0.0 --port 8000

# Web Interface - Health Check
web-health:
	curl -f http://localhost:8000/api/health || echo "Web interface not running"

# Docker cleanup
docker-clean:
	docker system prune -f

docker-clean-all:
	docker system prune -a -f
	docker volume prune -f
