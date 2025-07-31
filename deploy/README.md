# 🚀 Deploy Directory

This directory contains all deployment-related files and configurations for CSPBench.

## 📁 Directory Structure

```
deploy/
├── README.md                   # This file
├── app.yaml                    # Google Cloud App Engine config
├── cloudbuild.yaml            # Google Cloud Build config
├── deploy-cloud.sh            # Automated deployment script
├── docker-compose.yml         # Docker Compose configuration
├── skaffold.yaml              # Skaffold configuration for K8s
└── k8s/                       # Kubernetes manifests
    └── cloudrun-service.yaml  # Cloud Run service definition
```

## 🔧 Deployment Methods

### 1. Automated Script (Recommended)
```bash
./deploy/deploy-cloud.sh
```

### 2. Docker Compose (Local/Development)
```bash
docker-compose -f deploy/docker-compose.yml up --build
```

### 3. Makefile Shortcuts
```bash
make deploy          # Automated cloud deployment
make compose-up      # Docker Compose up
make compose-down    # Docker Compose down
```

### 4. Manual Cloud Build
```bash
gcloud builds submit --config=deploy/cloudbuild.yaml
```

### 5. Skaffold (K8s Development)
```bash
skaffold dev -f deploy/skaffold.yaml
```

## 🌍 Supported Platforms

- ✅ Google Cloud Run
- ✅ Google Cloud App Engine
- ✅ Google Kubernetes Engine (GKE)
- ✅ Docker Compose (local)
- ✅ Any Docker-compatible platform

## 📋 Prerequisites

1. **Google Cloud SDK** (for cloud deployments)
2. **Docker** (for containerization)
3. **kubectl** (for Kubernetes deployments)
4. **Skaffold** (optional, for K8s development)

## 🔐 Environment Variables

Set these before deployment:

```bash
export PROJECT_ID="your-gcp-project-id"
export REGION="us-central1"
export SERVICE_NAME="cspbench-web"
```

## 📖 For More Information

See the main project documentation:
- `/docs/DEPLOY_GOOGLE_CLOUD.md` - Detailed Google Cloud deployment guide
- `/docs/QUICKSTART_CLOUD.md` - Quick start guide
- `README.md` - Main project documentation
