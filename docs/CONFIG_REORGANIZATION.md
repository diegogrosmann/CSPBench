# ===================================================================
# CSPBENCH CONFIGURATION REORGANIZATION SUMMARY
# ===================================================================
# Data: 2025-08-01
# Autor: Sistema CSPBench
# 
# Este documento resume as mudanças realizadas na organização das 
# configurações e variáveis de ambiente do CSPBench
# ===================================================================

## 📋 OBJETIVO

Reorganizar as configurações para:
1. Seguir a estrutura do TEMPLATE.yaml 
2. Eliminar duplicatas entre settings.yaml e .env
3. Garantir uso correto de variáveis de ambiente
4. Centralizar configurações globais no .env
5. Manter configurações padrão de algoritmos nos próprios algoritmos

## 🔧 MUDANÇAS REALIZADAS

### 1. ARQUIVO .env REORGANIZADO

**Antes:**
- Configurações desordenadas
- Duplicatas (WEB_HOST, WEB_PORT)
- Nomes inconsistentes (EXPORT_FMT vs EXPORT_FORMAT)

**Depois:**
- Organizado em 7 seções lógicas
- Eliminadas duplicatas
- Nomes padronizados
- Documentação clara de cada seção

**Novas variáveis adicionadas:**
- OUTPUT_PATH=outputs (controla diretório base de saída)
- DEBUG=false (para modo debug da web interface)

**Variáveis renomeadas:**
- EXPORT_FMT → EXPORT_FORMAT

### 2. ARQUIVO settings.yaml REESTRUTURADO

**Antes:**
- Estrutura antiga sem seguir TEMPLATE.yaml
- Configurações espalhadas (logging, result, etc.)
- Mistura de configurações de infraestrutura com output

**Depois:**
- Segue exatamente a estrutura do TEMPLATE.yaml
- Seções 1, 2, 7-12 implementadas
- Seções 3-6 excluídas (específicas de batch)
- Configuração unificada de output
- Suporte a override via variáveis de ambiente

### 3. CÓDIGO ATUALIZADO

**Arquivos modificados:**
- main.py: Uso correto de LOG_LEVEL env var
- src/application/services/experiment_service.py: 
  - Uso de LOG_LEVEL
  - Uso de OUTPUT_PATH para diretórios
- src/presentation/cli/commands.py: 
  - EXPORT_FMT → EXPORT_FORMAT

**Melhorias implementadas:**
- Variáveis de ambiente têm precedência sobre configurações
- Caminhos hardcoded substituídos por env vars
- Configurações de logging unificadas

### 4. ARQUIVOS CRIADOS

- .env.template: Template limpo com comentários
- batches/test_config.yaml: Teste das configurações
- docs/CONFIG_REORGANIZATION.md: Esta documentação

## 📁 ESTRUTURA FINAL

```
.env                    # Configurações globais (produção)
.env.example           # Exemplo completo
.env.template          # Template limpo
config/settings.yaml   # Configurações padrão (seções 1,2,7-12)
batches/*.yaml         # Configurações específicas (seções 3-6)
```

## 🌍 VARIÁVEIS DE AMBIENTE

### Obrigatórias:
- NCBI_EMAIL (para datasets Entrez)

### Opcionais com defaults sensatos:
- LOG_LEVEL=INFO
- EXPORT_FORMAT=json
- DATASET_PATH=./datasets
- OUTPUT_PATH=outputs
- WEB_HOST=0.0.0.0
- WEB_PORT=8000
- EXECUTOR_IMPL=Executor
- INTERNAL_WORKERS=1

## ✅ TESTES REALIZADOS

1. ✅ Comando config-info funciona
2. ✅ Execução de batch funciona 
3. ✅ Logs são salvos corretamente
4. ✅ Resultados são exportados
5. ✅ Estrutura de diretórios está correta
6. ✅ Variáveis de ambiente são respeitadas

## 🔄 COMPATIBILIDADE

- ✅ Mantém compatibilidade com batches existentes
- ✅ Configurações antigas ainda funcionam
- ✅ Novos batches podem usar estrutura atualizada
- ✅ Variáveis de ambiente sobrescrevem configurações

## 🎯 BENEFÍCIOS

1. **Organização:** Configurações organizadas logicamente
2. **Flexibilidade:** Env vars permitem customização fácil
3. **Consistência:** Segue estrutura do TEMPLATE.yaml
4. **Manutenibilidade:** Código mais limpo e organizado
5. **Deployment:** Facilita deployment em diferentes ambientes
6. **Documentação:** Configurações bem documentadas

## 🔧 PRÓXIMOS PASSOS

1. Verificar se todos os algoritmos usam env vars corretamente
2. Revisar deploy configs para usar novas variáveis
3. Atualizar documentação do usuário
4. Considerar adicionar validação de configurações
