# 🚀 Guia de Deploy CSPBench na Google Cloud com VS Code

Este guia detalha como fazer deploy da aplicação CSPBench na Google Cloud usando VS Code com Cloud Code.

## 📋 Pré-requisitos

### 1. Extensões do VS Code
Instale as seguintes extensões:
- **Cloud Code** (`GoogleCloudTools.cloudcode`)
- **Docker** (`ms-azuretools.vscode-docker`)
- **Google Cloud Code** (inclui Kubernetes, Cloud Run support)

### 2. Ferramentas CLI
```bash
# Instalar Google Cloud CLI
curl https://sdk.cloud.google.com | bash
exec -l $SHELL

# Verificar instalação
gcloud --version

# Instalar Docker (se não estiver instalado)
sudo apt-get update
sudo apt-get install docker.io
sudo systemctl start docker
sudo usermod -aG docker $USER
```

## 🛠️ Configuração Inicial

### 1. Autenticação no Google Cloud
```bash
# Fazer login
gcloud auth login

# Configurar projeto (substitua PROJECT_ID pelo seu ID)
gcloud config set project YOUR_PROJECT_ID

# Habilitar APIs necessárias
gcloud services enable cloudbuild.googleapis.com
gcloud services enable run.googleapis.com
gcloud services enable containerregistry.googleapis.com
```

### 2. Configurar VS Code
1. Abra o VS Code no diretório do projeto
2. Abra a Command Palette (`Ctrl+Shift+P`)
3. Execute `Cloud Code: Sign In`
4. Selecione sua conta Google
5. Configure o projeto nas configurações:
   - Abra `.vscode/settings.json`
   - Substitua `""` em `"cloudcode.gcp.project"` pelo seu Project ID

### 3. Configurar Docker Registry
```bash
# Configurar autenticação do Docker para GCR
gcloud auth configure-docker
```

## 🔧 Deploy Usando Cloud Code (Recomendado)

### Método 1: Deploy Direto via Cloud Code

1. **Abrir Command Palette** (`Ctrl+Shift+P`)
2. **Executar** `Cloud Code: Deploy to Cloud Run`
3. **Selecionar configurações:**
   - Source: Current workspace
   - Dockerfile: `Dockerfile`
   - Service name: `cspbench-web`
   - Region: `us-central1`
   - Allow unauthenticated: `Yes`

4. **Monitorar o deploy** na aba "Cloud Code"

### Método 2: Deploy via Tasks do VS Code

1. **Abrir Terminal** (`Ctrl+Shift+~`)
2. **Executar Task** (`Ctrl+Shift+P` → `Tasks: Run Task`)
3. **Selecionar:** `Deploy to Cloud Run (gcloud)`

### Método 3: Deploy via Skaffold

1. **Abrir Command Palette** (`Ctrl+Shift+P`)
2. **Executar** `Cloud Code: Run on Cloud Run`
3. **Skaffold** usará o arquivo `skaffold.yaml` automaticamente

## 🧪 Teste Local Antes do Deploy

### 1. Build da imagem
```bash
# Via task do VS Code
Ctrl+Shift+P → Tasks: Run Task → "Build Cloud Image"

# Ou via terminal
docker build -t cspbench-web:latest .
```

### 2. Teste local
```bash
# Via task do VS Code
Ctrl+Shift+P → Tasks: Run Task → "Test Cloud Image"

# Ou via terminal
docker run --rm -p 8000:8000 cspbench-web:latest
```

### 3. Verificar aplicação
- Abra http://localhost:8000
- Teste o endpoint: http://localhost:8000/api/health

## 🔍 Deploy Manual via CLI

Se preferir deploy manual completo:

### 1. Build e Push
```bash
# Definir variáveis
export PROJECT_ID=your-project-id
export IMAGE_NAME=gcr.io/$PROJECT_ID/cspbench-web

# Build
docker build -t $IMAGE_NAME .

# Push
docker push $IMAGE_NAME
```

### 2. Deploy no Cloud Run
```bash
gcloud run deploy cspbench-web \\
  --image $IMAGE_NAME \\
  --region us-central1 \\
  --platform managed \\
  --allow-unauthenticated \\
  --memory 2Gi \\
  --cpu 1 \\
  --max-instances 10 \\
  --timeout 900 \\
  --port 8000
```

## 📊 Monitoramento e Debug

### 1. Logs via VS Code
1. **Abrir Command Palette** (`Ctrl+Shift+P`)
2. **Executar** `Cloud Code: View Logs`
3. **Selecionar** service `cspbench-web`

### 2. Logs via CLI
```bash
# Logs em tempo real
gcloud run services logs tail cspbench-web --region us-central1

# Logs específicos
gcloud logging read "resource.type=cloud_run_revision AND resource.labels.service_name=cspbench-web" --limit 50
```

### 3. Debugging
```bash
# Status do serviço
gcloud run services describe cspbench-web --region us-central1

# Listar revisions
gcloud run revisions list --service cspbench-web --region us-central1
```

## 🌐 Configurações Multi-Cloud

A aplicação foi configurada para ser cloud-agnostic:

### Variáveis de Ambiente Suportadas
- `PORT`: Porta da aplicação (padrão: 8000)
- `HOST`: Host binding (padrão: 0.0.0.0)
- `WORKERS`: Número de workers (padrão: 1)

### Deploy em Outras Clouds

#### AWS ECS/Fargate
```bash
# Use o mesmo Dockerfile
aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin 123456789012.dkr.ecr.us-east-1.amazonaws.com
docker tag cspbench-web:latest 123456789012.dkr.ecr.us-east-1.amazonaws.com/cspbench-web:latest
docker push 123456789012.dkr.ecr.us-east-1.amazonaws.com/cspbench-web:latest
```

#### Azure Container Instances
```bash
# Use o mesmo Dockerfile
az acr build --registry myregistry --image cspbench-web:latest .
az container create --resource-group myResourceGroup --name cspbench-web --image myregistry.azurecr.io/cspbench-web:latest
```

## 🔧 Solução de Problemas

### Erro: "Permission denied"
```bash
# Verificar autenticação
gcloud auth list
gcloud auth application-default login
```

### Erro: "Image not found"
```bash
# Verificar push da imagem
gcloud container images list --repository=gcr.io/YOUR_PROJECT_ID
```

### Erro: "Service timeout"
```bash
# Verificar logs para startup issues
gcloud run services logs tail cspbench-web --region us-central1
```

### Erro: "Port binding"
```bash
# Verificar se aplicação está escutando em 0.0.0.0:$PORT
# Verificar variável PORT no Cloud Run
```

## 🎯 URLs Importantes

Após deploy bem-sucedido:
- **Aplicação:** https://cspbench-web-xxxxx-uc.a.run.app
- **Health Check:** https://cspbench-web-xxxxx-uc.a.run.app/api/health
- **API Docs:** https://cspbench-web-xxxxx-uc.a.run.app/docs

## 📈 Otimizações de Produção

### 1. Configurar domínio customizado
```bash
gcloud run domain-mappings create --service cspbench-web --domain your-domain.com --region us-central1
```

### 2. Configurar SSL/TLS
Cloud Run gerencia SSL automaticamente para domínios customizados.

### 3. Configurar CI/CD
Use o arquivo `cloudbuild.yaml` para integração com Cloud Build:
```bash
gcloud builds submit --config cloudbuild.yaml .
```

## 🛡️ Segurança

- ✅ Container roda como usuário não-root
- ✅ Imagem multi-stage para menor superfície de ataque
- ✅ Health checks configurados
- ✅ Resource limits aplicados
- ✅ HTTPS enforced by Cloud Run
