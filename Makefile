venv:
	python3 -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt

# Docker - Build and execution
docker-build:
	docker build -t cspbench:latest .

docker-run:
	docker run --rm -it --env-file .env cspbench:latest

docker-run-interactive:
	docker run --rm -it --env-file .env -v "$(PWD)/datasets:/app/datasets:ro" -v "$(PWD)/outputs:/app/outputs" cspbench:latest

# Docker Compose
compose-up:
	docker-compose up --build

compose-up-detached:
	docker-compose up --build -d

compose-down:
	docker-compose down

compose-logs:
	docker-compose logs -f

# Docker cleanup
docker-clean:
	docker system prune -f

docker-clean-all:
	docker system prune -a -f
	docker volume prune -f
