# Changelog

Todas as mudanças notáveis neste projeto serão documentadas neste arquivo.

O formato é baseado em [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
e este projeto adere ao [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Documentação técnica completa
- Guia de contribuição detalhado
- Estrutura de changelog
- Exemplos de uso no README

### Changed
- README.md completamente reestruturado com arquitetura visual
- Documentação do src/ expandida com detalhes técnicos
- Melhorias na estrutura de documentação

### Improved
- Clareza na documentação de APIs
- Exemplos práticos de uso
- Guias de instalação e configuração

## [1.0.0] - 2025-01-08

### Added
- 🎯 **Sistema de Agendamento Inteligente**
  - ExecutionScheduler com fila FIFO absoluta
  - Controle automático de recursos (CPU/memória)
  - Monitoramento de processos filhos
  - Timeout configurável por tarefa
  - Delay inteligente entre execuções

- 🧬 **Algoritmos CSP Implementados**
  - Baseline: Algoritmo de consenso ganancioso
  - BLF-GA: Algoritmo genético híbrido com aprendizado por blocos
  - CSC: Constraint Satisfaction with Clustering
  - DP-CSP: Programação Dinâmica para CSP
  - H3-CSP: Heurística H3 para CSP

- 📊 **Sistema de Datasets**
  - Geração sintética com controle de ruído
  - Carregamento de arquivos FASTA
  - Download automático do NCBI/Entrez
  - Validação e normalização automática

- 🖥️ **Interface de Usuário**
  - CLI interativa com menus guiados
  - Interface curses para monitoramento em tempo real
  - Modo silencioso para automação
  - Sistema de progresso visual

- 📈 **Sistema de Relatórios**
  - Geração automática de relatórios JSON/CSV
  - Análise estatística comparativa
  - Visualizações de performance
  - Exportação flexível de dados

- 🔧 **Utilitários e Ferramentas**
  - Sistema de logging estruturado
  - Monitoramento de recursos do sistema
  - Configuração centralizada
  - Funções de distância otimizadas

- 🧪 **Otimização e Análise**
  - Otimização de hiperparâmetros com Optuna
  - Análise de sensibilidade de parâmetros
  - Execução em lote com YAML
  - Processamento paralelo controlado

### Technical Features
- **Arquitetura Modular**: Separação clara de responsabilidades
- **Extensibilidade**: Sistema de registro automático de algoritmos
- **Thread Safety**: Operações seguras em ambiente concorrente
- **Error Handling**: Tratamento robusto de erros
- **Type Safety**: Type hints completos
- **Testing**: Framework de testes abrangente

### Performance
- **Scheduler Otimizado**: Balanceamento automático de carga
- **Memory Management**: Limpeza automática de memória
- **Resource Monitoring**: Prevenção de sobrecarga do sistema
- **Parallel Execution**: Execução paralela eficiente

### Documentation
- **README Completo**: Guia abrangente de uso
- **Documentação Técnica**: Detalhes de implementação
- **API Documentation**: Documentação completa das APIs
- **Exemplos**: Casos de uso práticos
- **Guias**: Instalação, configuração e contribuição

### Quality Assurance
- **Code Standards**: Seguindo PEP 8 e boas práticas
- **Type Checking**: Verificação estática de tipos
- **Linting**: Análise de código com Ruff
- **Formatting**: Formatação automática com Black
- **Testing**: Testes unitários e de integração

## [0.9.0] - 2024-12-15

### Added
- Implementação inicial do BLF-GA
- Sistema básico de execução
- Interface CLI rudimentar
- Suporte básico a datasets sintéticos

### Changed
- Refatoração da estrutura de algoritmos
- Melhoria na organização de módulos

### Fixed
- Correções em bugs de execução
- Melhorias de estabilidade

## [0.8.0] - 2024-11-20

### Added
- Algoritmo Baseline implementado
- Sistema de logging básico
- Estrutura inicial do projeto

### Technical Debt
- Código legado removido
- Refatoração de interfaces
- Padronização de nomenclatura

## [0.7.0] - 2024-10-10

### Added
- Protótipo inicial
- Estrutura básica de classes
- Primeiros testes de conceito

## Tipos de Mudanças

- **Added**: para novas funcionalidades
- **Changed**: para mudanças em funcionalidades existentes
- **Deprecated**: para funcionalidades que serão removidas
- **Removed**: para funcionalidades removidas
- **Fixed**: para correções de bugs
- **Security**: para correções de segurança
- **Improved**: para melhorias gerais
- **Performance**: para otimizações de performance
- **Documentation**: para mudanças na documentação
- **Technical**: para mudanças técnicas internas

## Convenções de Versionamento

Este projeto usa [Semantic Versioning](https://semver.org/):

- **MAJOR**: Mudanças incompatíveis na API
- **MINOR**: Funcionalidades adicionadas de forma compatível
- **PATCH**: Correções de bugs compatíveis

### Exemplos:
- `1.0.0`: Primeira versão estável
- `1.1.0`: Nova funcionalidade compatível
- `1.1.1`: Correção de bug
- `2.0.0`: Mudança incompatível na API

## Links Úteis

- [Releases](https://github.com/seu-usuario/csp-blfga/releases)
- [Issues](https://github.com/seu-usuario/csp-blfga/issues)
- [Pull Requests](https://github.com/seu-usuario/csp-blfga/pulls)
- [Discussões](https://github.com/seu-usuario/csp-blfga/discussions)
