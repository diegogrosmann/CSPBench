# Dockerfile - CSPBench v0.1.0
FROM python:3.11-slim

# Container metadata
LABEL maintainer="Diego Grosmann <diego.grosmann@example.com>"
LABEL version="0.1.0"
LABEL description="CSPBench - Framework for Closest String Problem"

# Create non-root user
RUN groupadd -r cspbench && useradd -r -g cspbench cspbench

# Install system dependencies
RUN apt-get update && apt-get install -y \
    --no-install-recommends \
    git \
    && rm -rf /var/lib/apt/lists/*

# Create application directory
WORKDIR /app

# Copy and install Python dependencies first (for caching)
COPY requirements.txt .
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY src/ ./src/
COPY algorithms/ ./algorithms/
COPY config/ ./config/
COPY batches/ ./batches/
COPY main.py ./
COPY pyproject.toml ./

# Create necessary directories
RUN mkdir -p datasets outputs && \
    chown -R cspbench:cspbench /app

# Switch to non-root user
USER cspbench

# Healthcheck
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD python -c "import sys; sys.exit(0)"

# Default command
CMD ["python", "main.py"]
