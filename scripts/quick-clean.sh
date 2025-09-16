#!/usr/bin/env bash

# ===================================================================
# CSPBench - Limpeza Rápida
# ===================================================================
# Script simples para limpeza rápida quando o HD estiver cheio
# -------------------------------------------------------------------

set -e

echo "🧹 CSPBench - Limpeza Rápida"
echo "=============================="

# Função para mostrar espaço
show_space() {
    echo "📊 Espaço em disco:"
    df -h | grep -E "(Filesystem|/dev/loop4|overlay)" | head -2
    echo ""
}

echo "📊 Espaço antes da limpeza:"
show_space

# Parar e remover containers Docker
echo "🛑 Parando containers Docker..."
docker stop $(docker ps -aq) 2>/dev/null || true
docker rm $(docker ps -aq) 2>/dev/null || true

# Remover imagens Docker não utilizadas
echo "🗑️  Removendo imagens Docker..."
docker image prune -af >/dev/null 2>&1 || true

# Remover volumes Docker
echo "📦 Removendo volumes Docker..."
docker volume prune -f >/dev/null 2>&1 || true

# Remover build cache
echo "⚡ Removendo build cache..."
docker builder prune -af >/dev/null 2>&1 || true

# Limpeza geral Docker
echo "🔧 Limpeza geral Docker..."
docker system prune -af --volumes >/dev/null 2>&1 || true

# Limpar cache Python
echo "🐍 Limpando cache Python..."
find . -name "*.pyc" -delete 2>/dev/null || true
find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true

# Limpar arquivos temporários
echo "🗂️  Removendo arquivos temporários..."
rm -rf __pycache__ *.pyc *.pyo .pytest_cache logs outputs cache data *.db 2>/dev/null || true

# Limpar /tmp
echo "🧽 Limpando /tmp..."
sudo rm -rf /tmp/* 2>/dev/null || true

echo ""
echo "✅ Limpeza concluída!"
echo ""
show_space

echo "💡 Para limpeza mais detalhada, use:"
echo "   ./deploy/cleanup.sh"
echo "   make clean-full"
