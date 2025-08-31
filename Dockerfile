# ===================================================================
# CSPBench - Dockerfile
# ===================================================================

# ==========================
# Stage 1: Builder
# ==========================
FROM python:3.11-slim AS builder

# Flags seguras p/ Python/pip no build
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

# Dependências de compilação (para wheels nativos)
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
 && rm -rf /var/lib/apt/lists/*

# Virtualenv isolado
RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Instala dependências Python (cache-friendly)
COPY requirements.txt /tmp/requirements.txt
RUN pip install --upgrade pip && \
    pip install -r /tmp/requirements.txt

# ==========================
# Stage 2: Runtime
# ==========================
FROM python:3.11-slim AS runtime

# Flags genéricas de runtime
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PATH="/opt/venv/bin:$PATH"

# Dependências básicas de runtime
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
 && rm -rf /var/lib/apt/lists/*

# Usuário não-root
RUN groupadd -r app && useradd -r -g app app

# Copia o virtualenv do builder
COPY --from=builder /opt/venv /opt/venv

# Diretório da aplicação
WORKDIR /app

# Copia o código-fonte
COPY --chown=app:app . /app/

# Diretório de dados (ex.: SQLite, logs); persista via volume
RUN mkdir -p /app/data && chown -R app:app /app

# Troca para usuário não-root
USER app

# Exponha a porta que sua app usa por padrão (PORT vem do .env)
EXPOSE 8000

# A aplicação deve carregar .env via python-dotenv internamente.
# Ex.: load_dotenv("/app/.env"), se .env for montado no contêiner.
#
# Exemplos de execução:
#   docker run --rm -it --env-file ./.env -p 8000:8000 -v cspbench_data:/app/data imagem
#   docker run --rm -it -p 8000:8000 -v "$(pwd)/.env:/app/.env:ro" -v cspbench_data:/app/data imagem

# Comando padrão
CMD ["python", "main.py", "web"]
