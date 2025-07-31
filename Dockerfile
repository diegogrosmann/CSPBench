# Dockerfile - CSPBench Web Interface (Cloud-Optimized)
# Multi-stage build for production deployment

# Build stage
FROM python:3.11-slim AS builder

# Install build dependencies
RUN apt-get update && apt-get install -y \
    --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

# Copy and install dependencies
COPY requirements.txt .
RUN pip install --user --no-cache-dir --upgrade pip && \
    pip install --user --no-cache-dir -r requirements.txt

# Production stage
FROM python:3.11-slim AS production

# Install runtime dependencies only
RUN apt-get update && apt-get install -y \
    --no-install-recommends \
    curl \
    && rm -rf /var/lib/apt/lists/* \
    && groupadd -r cspbench \
    && useradd -r -g cspbench -m cspbench

# Copy Python packages from builder
COPY --from=builder /root/.local /home/cspbench/.local

# Set application directory
WORKDIR /app

# Copy application code with proper ownership
COPY --chown=cspbench:cspbench src/ ./src/
COPY --chown=cspbench:cspbench algorithms/ ./algorithms/
COPY --chown=cspbench:cspbench config/ ./config/
COPY --chown=cspbench:cspbench datasets/ ./datasets/
COPY --chown=cspbench:cspbench main.py pyproject.toml ./

# Create runtime directories
RUN mkdir -p outputs tmp logs batches \
    && chown -R cspbench:cspbench /app

# Switch to non-root user
USER cspbench

# Set environment variables (cloud-agnostic)
ENV PYTHONPATH=/app \
    PATH=/home/cspbench/.local/bin:$PATH \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PORT=8000 \
    HOST=0.0.0.0 \
    WORKERS=1

# Expose port (configurable via ENV)
EXPOSE $PORT

# Health check (works with most cloud providers)
HEALTHCHECK --interval=30s --timeout=10s --start-period=40s --retries=3 \
    CMD curl -f http://localhost:${PORT}/api/health || exit 1

# Default command (respects cloud provider PORT env var)
CMD ["python", "-m", "uvicorn", "src.presentation.web.app:app", "--host", "0.0.0.0", "--port", "8000"]
