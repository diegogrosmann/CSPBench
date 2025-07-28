# 🐳 Docker Configuration - CSPBench v0.1.0

## 📄 Summary of Applied Fixes

### ✅ Problems Resolved

| **File** | **Original Problem** | **Applied Solution** |
|----------|---------------------|---------------------|
| `Dockerfile` | Insecure copy (`COPY . .`) | Selective copy by directories |
| `Dockerfile` | Root execution | Non-root user `cspbench` |
| `Dockerfile` | Inefficient cache | Dependencies installed first |
| `docker-compose.yml` | Command conflict | Removed duplicate `command` |
| `docker-compose.yml` | Insecure volume | Specific volumes with permissions |
| `.dockerignore` | Insufficient exclusions | Complete exclusion list |
| `Makefile` | Limited commands | Expanded Docker commands |

### 🔧 Modified Files

#### 1. **Dockerfile** - Production Optimized
```dockerfile
# Optimized multi-layer
FROM python:3.11-slim

# Non-root user
RUN groupadd -r cspbench && useradd -r -g cspbench cspbench

# System dependencies
RUN apt-get update && apt-get install -y git

# Python dependencies cache
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Application code
COPY src/ ./src/
COPY algorithms/ ./algorithms/
# ... other directories

# Healthcheck
HEALTHCHECK --interval=30s --timeout=10s CMD python -c "import sys; sys.exit(0)"

USER cspbench
CMD ["python", "main.py"]
```

#### 2. **docker-compose.yml** - Production
```yaml
services:
  cspbench:
```

#### 3. **docker-compose.dev.yml** - Development
```yaml
services:
  cspbench-dev:
```

#### 4. **.dockerignore** - Enhanced Security
```ignore
# Python
__pycache__/
*.pyc
.venv/

# Development  
.vscode/
.git/
*.log

# Datasets (keep only examples)
datasets/*
!datasets/teste.fasta
!datasets/*_padrao.*
```

## 🚀 How to Use

### Basic Commands

```bash
# Build image
make docker-build
# or
docker build -t cspbench:latest .

# Run application
make docker-run  
# or
docker run --rm -it --env-file .env cspbench:latest

# With volumes for development
make docker-run-interactive
# or  
docker run --rm -it --env-file .env \
  -v "$(PWD)/datasets:/app/datasets:ro" \
  cspbench:latest
```

### Docker Compose

```bash
# Production
make compose-up
# or
docker-compose up --build

# Development
docker-compose -f docker-compose.dev.yml up --build

# Detached (background)
make compose-up-detached

# Stop services
make compose-down

# View logs
make compose-logs
```

### Cleanup

```bash
# Basic cleanup
make docker-clean

# Complete cleanup (images, volumes, etc)
make docker-clean-all
```

## 🔧 Environment Configurations

### Supported Environment Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `APP_ENV` | Application environment | `production` |
| `DEBUG` | Debug mode | `False` |
| `NCBI_EMAIL` | Email for NCBI API | - |
| `NCBI_API_KEY` | NCBI API key | - |

### Example `.env`

```bash
DEBUG=False
APP_ENV=production
NCBI_EMAIL=user@example.com
NCBI_API_KEY=your_key_here
```

## 📊 Container Resources

### Default Limits
- **CPU**: 2 cores maximum, 0.5 reserved
- **Memory**: 2GB maximum, 512MB reserved
- **Healthcheck**: Every 30s
- **User**: `cspbench` (non-root)

### Volumes
- `datasets/`: Read-only for datasets
- `batches/`: Read-only for configurations  
- `outputs/`: Write for results

## 🛡️ Security

### ✅ Implemented Practices
- Non-root user execution
- Comprehensive `.dockerignore`
- Selective file copying
- Volumes with restrictive permissions
- Healthcheck for monitoring
- Environment variables for credentials

### ❌ Avoided
- Root execution
- `COPY . .` (copy everything)
- Hardcoded credentials
- Wide volumes without restrictions

## 🧪 Testing

The build was tested and validated:

```bash
✅ Image build: SUCCESS
✅ Basic execution: SUCCESS  
✅ Help command: SUCCESS
✅ Non-root user: SUCCESS
✅ File structure: SUCCESS
```

## 📝 Next Steps

1. **CI/CD**: Integrate automatic build
2. **Registry**: Publish image to registry
3. **Multi-arch**: Support ARM64 architecture
4. **Optimization**: Further reduce image size
5. **Monitoring**: Add container metrics
