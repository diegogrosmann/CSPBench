# 🚀 Quick Start: Deploy CSPBench na Google Cloud

## ⚡ Passo a Passo Rápido (VS Code + Cloud Code)

### 1. **Pré-requisitos** (5 min)
```bash
# Instalar Google Cloud CLI
curl https://sdk.cloud.google.com | bash && exec -l $SHELL

# Login e configurar projeto
gcloud auth login
gcloud config set project SEU_PROJECT_ID

# Configurar Docker
gcloud auth configure-docker
```

### 2. **Configurar VS Code** (2 min)
1. Instalar extensão **Cloud Code** no VS Code
2. `Ctrl+Shift+P` → `Cloud Code: Sign In`
3. Editar `.vscode/settings.json` → substituir `"cloudcode.gcp.project": "SEU_PROJECT_ID"`

### 3. **Deploy Automático** (1 comando)
```bash
# Opção 1: Script automatizado
./deploy-cloud.sh

# Opção 2: Via VS Code
# Ctrl+Shift+P → "Cloud Code: Deploy to Cloud Run"
# Selecionar: Dockerfile, região us-central1

# Opção 3: Via task do VS Code
# Ctrl+Shift+P → "Tasks: Run Task" → "Deploy to Cloud Run (gcloud)"
```

### 4. **Verificar Deploy** (30 seg)
```bash
# Sua aplicação estará disponível em:
# https://cspbench-web-xxxxx-uc.a.run.app

# Health check:
# https://cspbench-web-xxxxx-uc.a.run.app/api/health
```

## 🔧 Arquivos Criados/Modificados

### ✅ Novos Arquivos (Otimizados para Multi-Cloud)
- `Dockerfile` - Docker otimizado para múltiplas nuvens (único)
- `skaffold.yaml` - Configuração Skaffold para Cloud Code
- `k8s/cloudrun-service.yaml` - Configuração Kubernetes/Cloud Run
- `app.yaml` - Configuração App Engine (alternativa)
- `cloudbuild.yaml` - Pipeline CI/CD Google Cloud Build
- `deploy-cloud.sh` - Script automatizado de deploy
- `docs/DEPLOY_GOOGLE_CLOUD.md` - Documentação completa

### ✅ Arquivos Modificados
- `.vscode/settings.json` - Configurações Cloud Code
- `.vscode/tasks.json` - Tasks de build e deploy

## 🎯 Características do Deploy

### ✅ Cloud-Agnostic
- Funciona em Google Cloud, AWS, Azure
- Configuração via variáveis de ambiente
- Dockerfile multi-stage otimizado

### ✅ Produção-Ready
- Container não-root para segurança
- Health checks configurados
- Resource limits aplicados
- Auto-scaling habilitado
- HTTPS automático

### ✅ Desenvolvimento-Friendly
- Hot reload no desenvolvimento
- Debug remoto via Cloud Code
- Logs integrados no VS Code
- Deploy com 1 comando

## 🆘 Solução Rápida de Problemas

```bash
# Erro de autenticação
gcloud auth login && gcloud auth application-default login

# Erro de projeto
gcloud config set project SEU_PROJECT_ID

# Verificar APIs habilitadas
gcloud services enable cloudbuild.googleapis.com run.googleapis.com

# Verificar logs
gcloud run services logs tail cspbench-web --region us-central1

# Status do serviço
gcloud run services describe cspbench-web --region us-central1
```

## 📱 URLs Importantes (Após Deploy)

- **App:** https://cspbench-web-xxxxx-uc.a.run.app
- **Health:** https://cspbench-web-xxxxx-uc.a.run.app/api/health  
- **Docs:** https://cspbench-web-xxxxx-uc.a.run.app/docs
- **Console:** https://console.cloud.google.com/run

---

🎉 **Deploy concluído!** Sua aplicação CSPBench está rodando na Google Cloud!
