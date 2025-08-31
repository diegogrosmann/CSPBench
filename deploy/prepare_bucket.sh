#!/usr/bin/env bash

# ===================================================================
# CSPBench - Preparação de Bucket GCS
# ===================================================================
# Cria o bucket (se não existir), estrutura de diretórios esperada e
# copia exemplos (datasets, batches) e cria estrutura de pastas usada pelas variáveis de ambiente.
# -------------------------------------------------------------------
# Uso:
#   ./deploy/prepare_bucket.sh -p <PROJECT_ID> -b <bucket-name> [-r us-central1]
# Exemplos:
#   ./deploy/prepare_bucket.sh -p meu-projeto -b csp-bench
#   ./deploy/prepare_bucket.sh -p meu-projeto -b csp-bench --no-examples
# Flags:
#   -p|--project    ID do projeto (obrigatório)
#   -b|--bucket     Nome do bucket (obrigatório, sem gs://)
#   -r|--region     Região/Location (default: us-central1)
#      --no-examples  Não copia exemplos
#      --overwrite    Sobrescreve arquivos existentes
#      --help         Ajuda
# ===================================================================
set -euo pipefail

REGION="us-central1"
PROJECT_ID=""
BUCKET=""
COPY_EXAMPLES=true
OVERWRITE=false

usage(){ grep '^#' "$0" | sed -E 's/^# ?//' | sed -n '1,/^$/p'; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    -p|--project) PROJECT_ID="$2"; shift 2;;
    -b|--bucket) BUCKET="$2"; shift 2;;
    -r|--region) REGION="$2"; shift 2;;
    --no-examples) COPY_EXAMPLES=false; shift;;
    --overwrite) OVERWRITE=true; shift;;
    --cpu|--memory|--prepare-bucket)
      echo "[ERRO] A flag '$1' pertence ao script cloud_run_deploy.sh (não use em prepare_bucket.sh)." >&2
      echo "       Exemplo correto: ./deploy/cloud_run_deploy.sh -p <PROJECT_ID> --bucket <bucket> --cpu 4 --memory 2Gi" >&2
      usage; exit 1;;
    --help|-h) usage; exit 0;;
    *) echo "[ERRO] Flag desconhecida: $1"; usage; exit 1;;
  esac
done

[[ -z "$PROJECT_ID" ]] && { echo "[ERRO] --project é obrigatório"; exit 1; }
[[ -z "$BUCKET" ]] && { echo "[ERRO] --bucket é obrigatório"; exit 1; }

command -v gcloud >/dev/null || { echo "[ERRO] gcloud não encontrado"; exit 1; }
command -v gsutil >/dev/null || { echo "[ERRO] gsutil não encontrado (instale com gcloud components)"; exit 1; }

echo "[INFO] Projeto: $PROJECT_ID | Bucket: $BUCKET | Região: $REGION"
gcloud config set project "$PROJECT_ID" >/dev/null

if gsutil ls -b "gs://$BUCKET" >/dev/null 2>&1; then
  echo "[INFO] Bucket já existe: gs://$BUCKET"
else
  echo "[INFO] Criando bucket gs://$BUCKET ..."
  gsutil mb -l "$REGION" -p "$PROJECT_ID" "gs://$BUCKET"
fi

# Estrutura esperada
DIRS=(batches datasets outputs logs)
for d in "${DIRS[@]}"; do
  echo "[INFO] Garantindo diretório: $d/"
  # Criar marcador vazio se não existir
  if ! gsutil ls "gs://$BUCKET/$d" >/dev/null 2>&1; then
    echo "[INFO] Criando gs://$BUCKET/$d/"
    echo "CSPBench placeholder" | gsutil cp - "gs://$BUCKET/$d/.keep"
  fi
done

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if $COPY_EXAMPLES; then
  echo "[INFO] Copiando exemplos..."
  if [[ -d "$ROOT_DIR/examples/datasets" ]]; then
    echo "[INFO] Datasets de exemplo -> gs://$BUCKET/datasets/"
    gsutil -m cp -r $( $OVERWRITE && echo "-n" ) "$ROOT_DIR/examples/datasets/*" "gs://$BUCKET/datasets/" 2>/dev/null || true
  fi
  if [[ -d "$ROOT_DIR/examples/batches" ]]; then
    echo "[INFO] Batches de exemplo -> gs://$BUCKET/batches/"
    gsutil -m cp -r $( $OVERWRITE && echo "-n" ) "$ROOT_DIR/examples/batches/*" "gs://$BUCKET/batches/" 2>/dev/null || true
  fi
fi

echo "[OK] Bucket preparado."
exit 0
