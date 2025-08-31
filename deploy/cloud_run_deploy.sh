#!/usr/bin/env bash

# ===================================================================
# CSPBench - Script de Deploy para Google Cloud Run
# ===================================================================
# Objetivo: Automatizar build, push e deploy da imagem Docker no Cloud Run
# Todas as configurações são definidas via arquivo de ambiente.
# -------------------------------------------------------------------
# Uso:
#   ./deploy/cloud_run_deploy.sh [--env-file .env.deploy] [--prepare-bucket] [--help]
#
# Exemplos:
#   ./deploy/cloud_run_deploy.sh
#   ./deploy/cloud_run_deploy.sh --env-file .env.production
#   ./deploy/cloud_run_deploy.sh --prepare-bucket
#
# Requisitos:
#   - gcloud autenticado (gcloud auth login)
#   - Permissão para Cloud Run & Artifact Registry / GCR
#   - Docker local funcionando para build
#   - Arquivo de ambiente configurado
# ===================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Variáveis padrão
ENV_FILE=".env.deploy"
PREPARE_BUCKET=false

usage() {
  cat << EOF
CSPBench - Deploy para Google Cloud Run

DESCRIÇÃO:
  Script automatizado para build, push e deploy no Cloud Run.
  Todas as configurações são definidas via arquivo de ambiente.

USO:
  $0 [FLAGS]

FLAGS:
  --env-file <arquivo>    Arquivo de ambiente (default: .env.deploy)
  --prepare-bucket        Prepara bucket GCS (cria estrutura + exemplos)
  --help                  Mostra esta ajuda

EXEMPLOS DE USO:

  # 1. Deploy básico (primeira vez)
  cp .env.deploy.example .env.deploy
  # Edite .env.deploy com seu PROJECT_ID
  $0 --prepare-bucket

  # 2. Deploy simples (configuração já pronta)
  $0

  # 3. Deploy com arquivo de ambiente customizado
  $0 --env-file .env.production

  # 4. Deploy apenas preparando bucket (sem deploy)
  $0 --prepare-bucket --env-file .env.staging

  # 5. Deploy para diferentes ambientes
  $0 --env-file .env.dev
  $0 --env-file .env.staging  
  $0 --env-file .env.production

CONFIGURAÇÃO INICIAL:
  1. Copie o arquivo exemplo:
     cp .env.deploy.example .env.deploy

  2. Configure as variáveis obrigatórias em .env.deploy:
     PROJECT_ID=seu-projeto-gcp
     BUCKET=seu-bucket-unico

  3. Execute o deploy:
     $0 --prepare-bucket

CONFIGURAÇÃO AVANÇADA:
  Para ambientes múltiplos, crie arquivos separados:
  - .env.deploy.dev
  - .env.deploy.staging  
  - .env.deploy.production

  E execute com:
  $0 --env-file .env.deploy.production

VARIÁVEIS PRINCIPAIS (definidas no arquivo .env):
  PROJECT_ID     - ID do projeto GCP (obrigatório)
  APP_NAME       - Nome do serviço (default: cspbench)
  REGION         - Região GCP (default: us-central1)
  BUCKET         - Bucket GCS (default: csp-bench)
  CPU/MEMORY     - Recursos (default: 4 CPUs, 2Gi RAM)
  DO_BUILD       - Fazer build local (default: true)
  DO_PUSH        - Push para registry (default: true)

EOF
}

# Parse de argumentos simplificado
while [[ $# -gt 0 ]]; do
  case "$1" in
    --env-file) ENV_FILE="$2"; shift 2;;
    --prepare-bucket) PREPARE_BUCKET=true; shift;;
    --help|-h) usage; exit 0;;
    *) echo "[ERRO] Flag desconhecida: $1"; usage; exit 1;;
  esac
done

# Carrega arquivo de ambiente
if [[ "${ENV_FILE}" != /* ]]; then
  ENV_FILE_ABS="${REPO_ROOT}/${ENV_FILE}"
else
  ENV_FILE_ABS="${ENV_FILE}"
fi

if [[ ! -f "${ENV_FILE_ABS}" ]]; then
  ENV_EXAMPLE_FILE="${ENV_FILE_ABS}.example"
  if [[ -f "${ENV_EXAMPLE_FILE}" ]]; then
    echo "[INFO] Arquivo ${ENV_FILE_ABS} não encontrado. Criando a partir de ${ENV_EXAMPLE_FILE}" >&2
    cp "${ENV_EXAMPLE_FILE}" "${ENV_FILE_ABS}"
    echo "[WARN] Arquivo ${ENV_FILE_ABS} criado com valores padrão. Revise as configurações antes de usar em produção." >&2
  else
    echo "[ERRO] Arquivo de ambiente não encontrado: ${ENV_FILE_ABS}" >&2
    echo "[ERRO] Arquivo exemplo também não encontrado: ${ENV_EXAMPLE_FILE}" >&2
    exit 1
  fi
fi

echo "[INFO] Carregando variáveis de ${ENV_FILE_ABS}" >&2
set -a
# shellcheck disable=SC1091
source "${ENV_FILE_ABS}"
set +a

# Validação de variáveis obrigatórias
required_vars=(PROJECT_ID APP_NAME REGION REGISTRY)
for var in "${required_vars[@]}"; do
  if [[ -z "${!var:-}" ]]; then
    echo "[ERRO] Variável obrigatória não definida: ${var}" >&2
    echo "[ERRO] Configure no arquivo: ${ENV_FILE_ABS}" >&2
    exit 1
  fi
done

# Defaults para variáveis opcionais
CPU="${CPU:-4}"
MEMORY="${MEMORY:-2Gi}"
BUCKET="${BUCKET:-csp-bench}"
MAX_INSTANCES="${MAX_INSTANCES:-2}"
MIN_INSTANCES="${MIN_INSTANCES:-0}"
CONCURRENCY="${CONCURRENCY:-80}"
ALLOW_UNAUTH="${ALLOW_UNAUTH:-true}"
DO_BUILD="${DO_BUILD:-true}"
DO_PUSH="${DO_PUSH:-true}"
MOUNT_PATH="${DATA_MOUNT_PATH:-/data}"
TAG="${TAG:-$(date +%Y%m%d-%H%M%S)}"

IMAGE_BASE="${REGISTRY}/${PROJECT_ID}/${APP_NAME}"
IMAGE_TAG="${IMAGE_BASE}:${TAG}"

echo "==================================================================="
echo " Deploy Cloud Run - CSPBench"
echo "==================================================================="
echo "APP_NAME       : ${APP_NAME}"
echo "PROJECT_ID     : ${PROJECT_ID}"
echo "REGION         : ${REGION}"
echo "ENV_FILE       : ${ENV_FILE_ABS}"
echo "IMAGE_TAG      : ${IMAGE_TAG}"
echo "CPU / MEMORY   : ${CPU} / ${MEMORY}"
echo "BUCKET         : ${BUCKET} (montado em ${MOUNT_PATH})"
echo "MAX_INSTANCES  : ${MAX_INSTANCES}"
echo "MIN_INSTANCES  : ${MIN_INSTANCES}"
echo "CONCURRENCY    : ${CONCURRENCY}"
echo "ALLOW_UNAUTH   : ${ALLOW_UNAUTH}"
echo "BUILD          : ${DO_BUILD}"
echo "PUSH           : ${DO_PUSH}"
echo "PREPARE_BUCKET : ${PREPARE_BUCKET}"
echo "==================================================================="

# Valida dependências
command -v gcloud >/dev/null 2>&1 || { echo "[ERRO] gcloud não encontrado no PATH"; exit 1; }
command -v docker >/dev/null 2>&1 || { echo "[ERRO] docker não encontrado no PATH"; exit 1; }

echo "[INFO] Verificando autenticação gcloud..."
if ! gcloud auth print-access-token >/dev/null 2>&1; then
  echo "[ERRO] Não autenticado. Execute: gcloud auth login" >&2
  exit 1
fi

echo "[INFO] Configurando projeto: ${PROJECT_ID}"
gcloud config set project "${PROJECT_ID}" >/dev/null

if [[ "${DO_BUILD}" == "true" ]]; then
  echo "[INFO] Build da imagem Docker..."
  docker build -t "${IMAGE_TAG}" -t "${IMAGE_BASE}:latest" "${REPO_ROOT}" || { echo "[ERRO] Build falhou"; exit 1; }
fi

if [[ "${DO_PUSH}" == "true" ]]; then
  echo "[INFO] Configurando Docker para autenticar no registry..."
  gcloud auth configure-docker --quiet || true
  echo "[INFO] Enviando imagem: ${IMAGE_TAG}"
  docker push "${IMAGE_TAG}"
  docker push "${IMAGE_BASE}:latest" || true
fi

if [[ "${PREPARE_BUCKET}" == "true" ]]; then
  if command -v gsutil >/dev/null 2>&1; then
    echo "[INFO] Preparando bucket ${BUCKET}..."
    "${SCRIPT_DIR}/prepare_bucket.sh" -p "${PROJECT_ID}" -b "${BUCKET}" -r "${REGION}" || { echo "[ERRO] Falha ao preparar bucket"; exit 1; }
  else
    echo "[WARN] gsutil não encontrado; pulando preparação de bucket. Instale via gcloud components install gsutil" >&2
  fi
fi

# Preparar variáveis de ambiente da aplicação
APP_ENV_VARS="PYTHONUNBUFFERED=1,LOG_TO_STDOUT=true,DATA_MOUNT_PATH=${MOUNT_PATH}"

# Adicionar variáveis de aplicação críticas se definidas no arquivo .env
if [[ -n "${SETTINGS_PATH:-}" ]]; then
  APP_ENV_VARS="${APP_ENV_VARS},SETTINGS_PATH=${SETTINGS_PATH}"
fi
if [[ -n "${NCBI_EMAIL:-}" ]]; then
  APP_ENV_VARS="${APP_ENV_VARS},NCBI_EMAIL=${NCBI_EMAIL}"
fi
if [[ -n "${NCBI_API_KEY:-}" ]]; then
  APP_ENV_VARS="${APP_ENV_VARS},NCBI_API_KEY=${NCBI_API_KEY}"
fi
if [[ -n "${WEB_HOST:-}" ]]; then
  APP_ENV_VARS="${APP_ENV_VARS},WEB_HOST=${WEB_HOST}"
fi
if [[ -n "${WEB_DEBUG:-}" ]]; then
  APP_ENV_VARS="${APP_ENV_VARS},WEB_DEBUG=${WEB_DEBUG}"
fi
if [[ -n "${WEB_LOG_LEVEL:-}" ]]; then
  APP_ENV_VARS="${APP_ENV_VARS},WEB_LOG_LEVEL=${WEB_LOG_LEVEL}"
fi
if [[ -n "${WEB_ACCESS_LOG:-}" ]]; then
  APP_ENV_VARS="${APP_ENV_VARS},WEB_ACCESS_LOG=${WEB_ACCESS_LOG}"
fi
if [[ -n "${LOG_LEVEL:-}" ]]; then
  APP_ENV_VARS="${APP_ENV_VARS},LOG_LEVEL=${LOG_LEVEL}"
fi
if [[ -n "${FORCE_CLEANUP:-}" ]]; then
  APP_ENV_VARS="${APP_ENV_VARS},FORCE_CLEANUP=${FORCE_CLEANUP}"
fi
if [[ -n "${DATASET_DIRECTORY:-}" ]]; then
  APP_ENV_VARS="${APP_ENV_VARS},DATASET_DIRECTORY=${DATASET_DIRECTORY}"
fi
if [[ -n "${BATCH_DIRECTORY:-}" ]]; then
  APP_ENV_VARS="${APP_ENV_VARS},BATCH_DIRECTORY=${BATCH_DIRECTORY}"
fi
if [[ -n "${OUTPUT_BASE_DIRECTORY:-}" ]]; then
  APP_ENV_VARS="${APP_ENV_VARS},OUTPUT_BASE_DIRECTORY=${OUTPUT_BASE_DIRECTORY}"
fi
if [[ -n "${WORK_DB_PATH:-}" ]]; then
  APP_ENV_VARS="${APP_ENV_VARS},WORK_DB_PATH=${WORK_DB_PATH}"
fi

echo "[INFO] Variáveis de ambiente da aplicação: ${APP_ENV_VARS}"

DEPLOY_CMD=(
  gcloud run deploy "${APP_NAME}" \
    --image "${IMAGE_TAG}" \
    --platform managed \
    --region "${REGION}" \
    --memory "${MEMORY}" \
    --cpu "${CPU}" \
    --max-instances "${MAX_INSTANCES}" \
    --min-instances "${MIN_INSTANCES}" \
    --concurrency "${CONCURRENCY}" \
    --set-env-vars "${APP_ENV_VARS}" \
    --add-volume "name=shared,type=cloud-storage,bucket=${BUCKET}" \
    --add-volume-mount "volume=shared,mount-path=${MOUNT_PATH}"
)

if [[ "${ALLOW_UNAUTH}" == "true" ]]; then
  DEPLOY_CMD+=(--allow-unauthenticated)
else
  DEPLOY_CMD+=(--no-allow-unauthenticated)
fi

echo "[INFO] Realizando deploy no Cloud Run..."
"${DEPLOY_CMD[@]}"

echo "[INFO] Obtendo URL do serviço..."
SERVICE_URL=$(gcloud run services describe "${APP_NAME}" --region "${REGION}" --format 'value(status.url)')
echo "==================================================================="
echo "✓ Deploy concluído com sucesso!"
echo "Service URL: ${SERVICE_URL}"
echo "Imagem      : ${IMAGE_TAG}"
echo "==================================================================="

exit 0
echo "[INFO] Realizando deploy no Cloud Run..."
"${DEPLOY_CMD[@]}"

echo "[INFO] Obtendo URL do serviço..."
SERVICE_URL=$(gcloud run services describe "${APP_NAME}" --region "${REGION}" --format 'value(status.url)')
echo "==================================================================="
echo "✓ Deploy concluído com sucesso!"
echo "Service URL: ${SERVICE_URL}"
echo "Imagem      : ${IMAGE_TAG}"
echo "==================================================================="

exit 0
