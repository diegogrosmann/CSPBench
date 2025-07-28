# Refatoração dos Parâmetros do Batch - CSPBench

## Resumo das Mudanças

Este documento resume a refatoração completa dos parâmetros de batch do CSPBench para usar os nomes padronizados do `TEMPLATE_PADRONIZADO.yaml`.

## 1. Parâmetros Mapeados e Implementados

### ✅ Parâmetros já implementados (mantidos):
- **metadados**: Informações básicas do batch
- **datasets**: Configuração de datasets
- **algorithms**: Configuração de algoritmos
- **task**: Tipo de tarefa (execution/optimization/sensitivity)
- **execution**: Configurações de execução
- **optimization**: Configurações de otimização
- **sensitivity**: Configurações de análise de sensibilidade
- **resources**: Configurações de recursos (parcialmente)

### 🆕 Parâmetros recém-implementados:
- **infrastructure**: Configurações de infraestrutura (histórico e resultados)
- **export**: Configurações de exportação (formatos, destino, filtros)
- **plots**: Configurações de visualização (tipos de gráficos, estilos)
- **monitoring**: Configurações de monitoramento (interface, intervalos)
- **logging**: Configurações de logging (níveis, saídas, formatos)
- **system**: Configurações de sistema (reprodutibilidade, checkpoints, erro handling)

## 2. Arquivos Modificados

### 2.1 `src/application/services/config_parser.py`
**Mudanças principais:**
- Adicionadas novas dataclasses para cada seção de configuração:
  - `InfrastructureConfig`
  - `ExportConfig`
  - `PlotsConfig`
  - `MonitoringConfig`
  - `LoggingConfig`
  - `SystemConfig`
- Adicionados novos métodos de parsing:
  - `parse_infrastructure_config()`
  - `parse_export_config()`
  - `parse_plots_config()`
  - `parse_monitoring_config()`
  - `parse_logging_config()`
  - `parse_system_config()`

### 2.2 `src/application/services/experiment_service.py`
**Mudanças principais:**
- Atualizado `run_batch()` para usar todos os novos parsers
- Adicionados novos métodos auxiliares:
  - `_update_batch_logging_with_config()`
  - `_configure_monitoring()`
  - `_export_batch_results_with_config()`
  - `_generate_plots_with_config()`
- Assinaturas dos métodos de processamento atualizadas para aceitar `resources_config`
- Importações atualizadas para incluir as novas classes

## 3. Estrutura de Parsing

### 3.1 Fluxo de Parsing Atualizado
```python
# No ExperimentService.run_batch()
batch_config = self._parse_batch_config(batch_cfg)

# Parse de todas as seções
metadata = ConfigurationParser.parse_metadata(batch_config)
infrastructure_config = ConfigurationParser.parse_infrastructure_config(batch_config)
export_config = ConfigurationParser.parse_export_config(batch_config)
plots_config = ConfigurationParser.parse_plots_config(batch_config)
monitoring_config = ConfigurationParser.parse_monitoring_config(batch_config)
logging_config = ConfigurationParser.parse_logging_config(batch_config)
system_config = ConfigurationParser.parse_system_config(batch_config)
resources_config = ConfigurationParser.parse_resources_config(batch_config)
```

### 3.2 Configuração Aplicada
```python
# Logging
self._update_batch_logging_with_config(logging_config)

# Monitoring
if monitoring_config.enabled and self._monitoring_service:
    self._configure_monitoring(monitoring_config)

# Processamento com recursos
results = self._process_execution_batch(batch_config, resources_config)

# Export
if export_config.enabled:
    self._export_batch_results_with_config(export_data, export_config, task_type)

# Plots
if plots_config.enabled:
    self._generate_plots_with_config(results, plots_config, task_type)
```

## 4. Exemplo de Uso

Um arquivo de exemplo completo foi criado em `batches/exemplo_refatorado.yaml` demonstrando o uso de todos os novos parâmetros.

## 5. Retrocompatibilidade

A implementação mantém retrocompatibilidade através de:
- Valores padrão para todos os novos parâmetros
- Fallbacks para configurações legadas
- Graceful degradation quando seções não estão presentes

## 6. Implementação Futura

### 6.1 TODOs Identificados
- **Plots**: Implementação completa de geração de gráficos
- **Infrastructure**: Uso efetivo das configurações de histórico
- **Logging**: Configurações avançadas de arquivo e formatação
- **System**: Implementação de checkpointing e early stopping
- **Export**: Suporte para formatos Parquet e Pickle

### 6.2 Pontos de Extensão
Todas as implementações incluem pontos de extensão para funcionalidades futuras, permitindo evolução gradual sem breaking changes.

## 7. Teste e Validação

- ✅ Syntax check: Todos os arquivos validados
- ✅ Exemplo funcional: `exemplo_refatorado.yaml` criado
- ✅ Importações: Todas as novas classes importáveis
- ✅ Estrutura: Mantém padrões do hexagonal architecture

## 8. Benefícios Alcançados

1. **Padronização**: Alinhamento com template oficial
2. **Extensibilidade**: Estrutura preparada para funcionalidades futuras
3. **Modularidade**: Cada seção tem parser e validação próprios
4. **Flexibilidade**: Configurações granulares para cada aspecto
5. **Manutenibilidade**: Código organizado e bem documentado
