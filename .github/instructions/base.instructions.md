---
applyTo: '**'
---

# 📋 Diretrizes de Desenvolvimento – CSPBench v0.1.0

> **🔒 IMUTÁVEL:** Estas diretrizes devem ser mantidas em **todas** as interações futuras com IAs. Elas consolidam a arquitetura, padrões e práticas operacionais implementadas. Qualquer proposta de mudança **deve** ser avaliada e aprovada antes da implementação.

---

## 📑 Índice

1. [🏗️ Arquitetura e Estrutura](#1-arquitetura-e-estrutura)
2. [🧩 Sistema de Plugins (Algoritmos)](#2-sistema-de-plugins-algoritmos)
3. [⚙️ Configuração e Ambiente](#3-configuração-e-ambiente)
4. [🔧 Comandos e Fluxos Críticos](#4-comandos-e-fluxos-críticos)
5. [📝 Convenções de Código](#5-convenções-de-código)
6. [🧪 Testes e Qualidade](#6-testes-e-qualidade)
7. [🔐 Segurança e Credenciais](#7-segurança-e-credenciais)
8. [📄 Documentação](#8-documentação)
9. [🔄 Processo de Desenvolvimento](#9-processo-de-desenvolvimento)
10. [🌐 Internacionalização](#10-internacionalização)
11. [🤖 Diretrizes para IAs](#11-diretrizes-para-ias)

---

## 1. 🏗️ Arquitetura e Estrutura

### 1.1 Arquitetura Hexagonal (Clean Architecture)

| Camada | Diretório | Responsabilidade | Dependências Permitidas |
|--------|-----------|------------------|------------------------|
| **🔷 Domain** | `src/domain/` | Regras de negócio, entidades, interfaces | **APENAS** Python StdLib |
| **🔶 Application** | `src/application/` | Casos de uso, orquestração, portas | Domain + interfaces |
| **🔸 Infrastructure** | `src/infrastructure/` | Adaptadores externos, I/O, APIs | Application Ports + StdLib |
| **🔹 Presentation** | `src/presentation/` | CLI, TUI, interfaces de usuário | Application Services |
| **🧩 Plugins** | `algorithms/` | Implementações de algoritmos | Domain (apenas interfaces) |

### 1.2 Princípios Fundamentais

- ✅ **Inversão de Dependências**: Application Services recebem implementações via DI
- ✅ **Isolamento do Domain**: Camada Domain é **pura** - sem I/O, frameworks ou APIs externas
- ✅ **Plugin System**: Algoritmos são descobertos dinamicamente via `global_registry`
- ❌ **Imports Diretos**: NUNCA `from algorithms.xyz import ...` na aplicação

### 1.3 Estrutura de Diretórios Obrigatória

```
project_root/
├── src/                    # Core da aplicação
│   ├── domain/            # Regras de negócio
│   ├── application/       # Casos de uso
│   ├── infrastructure/    # Adaptadores
│   └── presentation/      # Interfaces
├── algorithms/            # Plugins de algoritmos
├── config/               # Configurações
├── batches/              # Configurações de batch
├── tests/                # Testes organizados por tipo
└── datasets/             # Datasets para experimentos
```

---

## 2. 🧩 Sistema de Plugins (Algoritmos)

### 2.1 Estrutura Padrão de Plugin

```
algorithms/algorithm_name/
├── __init__.py          # Importa e reexporta a classe principal
├── algorithm.py         # Classe wrapper com @register_algorithm
├── implementation.py    # Lógica core isolada
├── config.py            # Parâmetros padrão
└── README.md            # Documentação específica
```

### 2.2 Template de Implementação

```python
# algorithm.py
from src.domain.algorithms import CSPAlgorithm, register_algorithm
from .implementation import AlgorithmCore
from .config import DEFAULT_PARAMS

@register_algorithm
class MyAlgorithm(CSPAlgorithm):
    """Descrição clara do algoritmo."""
    
    name = "MyAlgorithm"
    default_params = DEFAULT_PARAMS
    
    def run(self, dataset, **kwargs):
        """
        Executa o algoritmo.
        
        Args:
            dataset: Dataset a processar
            **kwargs: Parâmetros adicionais
            
        Returns:
            tuple[str, int, dict]: (best_string, max_distance, metadata)
        """
        # Usar self._report_progress(msg) para feedback
        return AlgorithmCore().solve(dataset, **kwargs)
```

### 2.3 Regras para Plugins

- ✅ **Determinismo**: Mesmo `seed` = mesmo resultado
- ✅ **Isolamento**: Um plugin não pode depender de outro
- ✅ **Interface Única**: Implementar `CSPAlgorithm.run()`
- ✅ **Registro Automático**: Usar `@register_algorithm`
- ✅ **Feedback**: Usar `self._report_progress()` para progresso

---

## 3. ⚙️ Configuração e Ambiente

### 3.1 Ambiente Virtual Python

O projeto utiliza um ambiente virtual Python localizado em `.venv/`:

| Componente | Localização | Descrição |
|------------|-------------|-----------|
| **Executável Python** | `.venv/bin/python` | Interpretador Python do ambiente |
| **Pip** | `.venv/bin/pip` | Gerenciador de pacotes |
| **Scripts** | `.venv/bin/` | Executáveis e scripts do ambiente |

#### Ativação do Ambiente

```bash
# Ativar ambiente virtual (Linux/macOS)
source .venv/bin/activate

# Verificar ambiente ativo
which python  # Deve retornar .venv/bin/python
```

#### Execução Direta

```bash
# Executar sem ativar o ambiente
.venv/bin/python main.py

# Instalar dependências
.venv/bin/pip install -r requirements.txt
```

### 3.2 Hierarquia de Configuração

| Prioridade | Fonte | Localização | Propósito |
|------------|-------|-------------|-----------|
| 1 (Maior) | Variáveis de Ambiente | Sistema | Override runtime |
| 2 | Arquivo Batch | `batches/*.yaml` | Configuração específica |
| 3 (Menor) | Settings Global | `config/settings.yaml` | Configuração base |

### 3.3 Arquivos de Configuração

```yaml
# config/settings.yaml - Configuração base do sistema
infrastructure:
  dataset_repository:
    implementation: "FileDatasetRepository"
    config:
      base_path: "datasets/"
  
# batches/exemplo.yaml - Configuração de experimento
execution:
  algorithm: "BLF_GA"
  dataset: "test_dataset.fasta"
  params:
    population_size: 100
```

### 3.4 Variáveis de Ambiente

```bash
# .env (não versionado)
NCBI_EMAIL=user@example.com
NCBI_API_KEY=your_key_here
```

### 3.5 Credenciais e Segurança

- ✅ **Versionamento**: Commit apenas `.env.example`
- ❌ **Dados Sensíveis**: NUNCA commit credenciais reais
- ✅ **Carregamento**: Use `python-dotenv` para carregar `.env`

---

## 4. 🔧 Comandos e Fluxos Críticos

### 4.1 CLI Principal

| Comando | Descrição | Exemplo |
|---------|-----------|---------|
| `python main.py` | Menu interativo | - |
| `python main.py --help` | Ajuda geral | - |
| `python main.py --algorithms` | Lista algoritmos disponíveis | - |
| `python main.py --datasetsave` | Gera/salva datasets | - |
| `python main.py <file.yaml>` | Executa batch | `python main.py batches/teste.yaml` |

#### Comandos Específicos

| Comando | Descrição | Exemplo |
|---------|-----------|---------|
| `test` | Teste básico do sistema | `python main.py test` |
| `run <algorithm> <dataset>` | Executa algoritmo em dataset | `python main.py run Baseline test.fasta` |
| `batch <file.yaml>` | Executa batch | `python main.py batch batches/example.yaml` |
| `algorithms` | Lista algoritmos | `python main.py algorithms` |
| `config-info` | Mostra configuração | `python main.py config-info` |
| `sessions` | Lista sessões | `python main.py sessions` |
| `cleanup` | Remove sessões antigas | `python main.py cleanup` |
| `show-session <name>` | Mostra detalhes da sessão | `python main.py show-session session_name` |
| `view-report <name>` | Abre relatório no browser | `python main.py view-report session_name` |

#### Uso com Ambiente Virtual

```bash
# Forma recomendada - usando ambiente virtual diretamente
.venv/bin/python main.py
.venv/bin/python main.py test
.venv/bin/python main.py run BLF_GA test_dataset.fasta
.venv/bin/python main.py batches/example.yaml

# Alternativa - ativando ambiente primeiro
source .venv/bin/activate
python main.py
python main.py test
deactivate  # Para desativar o ambiente
```

### 4.2 Tasks do VS Code

| Task | Propósito | Comando |
|------|-----------|---------|
| Run CSPBench | Executar aplicação | `python main.py` |
| Run Tests | Executar testes | `pytest tests/ -v` |
| Format Code | Formatar código | `black .` |
| Lint Code | Análise estática | `ruff check .` |
| Type Check | Verificação de tipos | `mypy src/` |
| Run Coverage | Cobertura de testes | `pytest --cov=src --cov-report=html tests/` |

**Nota**: As tasks do VS Code utilizam o ambiente virtual automaticamente quando configurado corretamente.

### 4.3 Fluxo de Batch

```yaml
# Tipos de batch suportados
execution:     # Execução simples
optimization:  # Otimização de parâmetros  
sensitivity:   # Análise de sensibilidade
```

---

## 5. 📝 Convenções de Código

### 5.1 Ferramentas de Qualidade

| Ferramenta | Propósito | Configuração | Obrigatório |
|------------|-----------|--------------|-------------|
| `black` | Formatação | `pyproject.toml` | ✅ |
| `ruff` | Linting | `pyproject.toml` | ✅ |
| `mypy` | Type checking | `--strict` | ✅ |
| `isort` | Ordenação imports | `pyproject.toml` | ✅ |

### 5.2 Contratos de Interface

```python
# Interface CSPAlgorithm
def run(self, dataset, **kwargs) -> tuple[str, int, dict]:
    """
    Returns:
        tuple: (best_string, max_distance, metadata_dict)
    """
```

### 5.3 Tratamento de Erros

- ✅ **Exceções Customizadas**: Use `src.domain.errors`
- ✅ **Propagação**: Deixe erros subirem com contexto
- ✅ **Logging**: Log erros com nível apropriado

---

## 6. 🧪 Testes e Qualidade

### 6.1 Estrutura de Testes

```
tests/
├── unit/              # Testes unitários (mocks/fakes)
├── integration/       # Testes de integração (filesystem real)
├── algorithms/        # Testes específicos de algoritmos
└── fixtures/          # Dados de teste compartilhados
```

### 6.2 Diretrizes por Tipo

| Tipo | Localização | Características | Ferramentas |
|------|-------------|----------------|-------------|
| **Unit** | `tests/unit/` | Isolados, rápidos, mocks | pytest, unittest.mock |
| **Integration** | `tests/integration/` | I/O real, tmpdir | pytest, tempfile |
| **Algorithm** | `tests/algorithms/` | Determinismo, contratos | pytest, fixtures |

### 6.3 Checklist de Teste

- ✅ **Determinismo**: Algoritmos com mesmo seed = mesmo resultado
- ✅ **Contratos**: Validar retorno `(str, int, dict)`
- ✅ **Edge Cases**: Testar casos limite
- ✅ **Cobertura**: Manter cobertura > 80%

---

## 7. 🔐 Segurança e Credenciais

### 7.1 Gestão de Credenciais

- ✅ **Arquivo .env**: Para desenvolvimento local
- ✅ **Variáveis Sistema**: Para produção
- ❌ **Hardcode**: NUNCA credenciais no código
- ❌ **Versionamento**: NUNCA commit `.env` real

### 7.2 Template de Credenciais

```bash
# .env.example (versionado)
NCBI_EMAIL=your.email@example.com
NCBI_API_KEY=your_ncbi_api_key_here
```

### 7.3 Validação de Credenciais

```python
def load_credentials():
    """Carrega e valida credenciais necessárias."""
    required_vars = ["NCBI_EMAIL", "NCBI_API_KEY"]
    missing = [var for var in required_vars if not os.getenv(var)]
    
    if missing:
        raise EnvironmentError(f"Variáveis obrigatórias: {missing}")
```

---

## 8. 📄 Documentação

### 8.1 Versionamento

- ✅ **Ponto Único**: `src/__init__.py` contém `__version__`
- ✅ **Sincronização**: Atualizar `pyproject.toml`, docs, CLI
- ✅ **Changelog**: Documentar mudanças por versão

### 8.2 Documentação de Algoritmo

```markdown
# Algorithm Name

## Descrição
Breve descrição do algoritmo e seu propósito.

## Parâmetros
- `param1`: Descrição do parâmetro
- `param2`: Descrição do parâmetro

## Exemplo de Uso
```yaml
algorithm: "AlgorithmName"
params:
  param1: value1
  param2: value2
```

## Referências
- Paper original ou fonte
```

---

## 9. 🔄 Processo de Desenvolvimento

### 9.1 Fluxo de Mudanças

1. **🔍 Análise**
   - Identificar claramente o que fazer
   - Verificar alinhamento com diretrizes
   - Avaliar impacto no sistema

2. **📖 Compreensão**
   - Ler e analisar código existente
   - Consultar documentação e testes
   - Entender contexto e dependências

3. **📋 Planejamento**
   - Definir arquivos a alterar
   - Planejar implementação
   - Considerar necessidade de testes/docs

4. **✋ Aprovação**
   - Documentar alterações planejadas
   - Solicitar aprovação antes de implementar
   - Aguardar confirmação para prosseguir

5. **⚡ Implementação**
   - Aplicar mudanças conforme aprovado
   - Seguir diretrizes de estilo e arquitetura
   - Criar/atualizar testes necessários

6. **🧹 Limpeza**
   - Remover código obsoleto
   - Verificar que arquivos removidos não são necessários
   - Garantir que não há impactos não intencionais

7. **✅ Validação**
   - Executar todos os testes
   - Verificar que nada foi quebrado
   - Validar funcionamento da funcionalidade

### 9.2 Branches e Versionamento

```bash
# Padrão de branches
feature/<descricao>     # Novas funcionalidades
bugfix/<descricao>      # Correções
refactor/<descricao>    # Refatorações
```

### 9.3 Checklist de PR

- ✅ Testes passando
- ✅ Cobertura mantida/melhorada
- ✅ Documentação atualizada
- ✅ Diretrizes seguidas
- ✅ Código limpo e bem comentado

---

## 10. 🌐 Internacionalização

### 10.1 Política Obrigatória

**REGRA FUNDAMENTAL**: Todo código deve estar em inglês - sem exceções.

| Elemento | Idioma | Exemplo |
|----------|--------|---------|
| **Variáveis/Funções** | Inglês | `calculate_distance()` |
| **Classes/Métodos** | Inglês | `class Algorithm:` |
| **Docstrings** | Inglês | `"""Execute algorithm and return..."""` |
| **Comentários** | Inglês | `# Save final state to history` |
| **Mensagens** | Inglês | `"Starting CSP algorithm..."` |
| **Metadados** | Inglês | `"center_found": center` |

### 10.2 Mapeamento de Termos

| Português | Inglês | Contexto |
|-----------|--------|----------|
| `algoritmo` | `algorithm` | Geral |
| `execução` | `execution` | Logs/código |
| `configuração` | `configuration` | Config |
| `parâmetros` | `parameters` | Documentação |
| `centro_encontrado` | `center_found` | Metadados |
| `melhor_distancia` | `best_distance` | Metadados |
| `iteracoes` | `iterations` | Metadados |
| `Iniciando` | `Starting` | Mensagens |
| `finalizado com sucesso` | `completed successfully` | Mensagens |

### 10.3 Auditoria

```bash
# Verificar termos em português
grep -r "algoritmo\|execução\|configuração" src/
grep -r "Iniciando\|Executando\|Finalizando" src/
grep -r "[àáâãäçéêëíîïóôõöúûüÀÁÂÃÄÇÉÊËÍÎÏÓÔÕÖÚÛÜ]" src/
```

---

## 11. 🤖 Diretrizes para IAs

### 11.1 Comunicação
- ✅ **Resposta**: Sempre em português brasileiro
- ✅ **Clareza**: Manter precisão técnica
- ✅ **Contexto**: Considerar arquitetura existente

### 11.2 Implementação
- ✅ **Ambiente Virtual**: Sempre usar `.venv/bin/python` para execução
- ✅ **Ferramentas Internas**: Priorizar APIs/funções internas vs comandos externos
- ✅ **Incrementalidade**: Mudanças pequenas e testáveis
- ✅ **Validação**: Testar após mudanças significativas
- ✅ **Cleanup**: Remover código obsoleto

### 11.3 Restrições Críticas
- ❌ **Dependências Cruzadas**: Plugins NÃO podem depender uns dos outros
- ❌ **Imports Diretos**: Aplicação NÃO pode importar plugins diretamente
- ❌ **Domain Purity**: Domain layer NÃO pode ter I/O ou dependências externas
- ❌ **Credenciais**: NUNCA hardcodar ou versionar credenciais

### 11.4 Internacionalização Obrigatória

**AÇÃO AUTOMÁTICA**: Sempre traduzir português → inglês ao modificar código.

#### Checklist Obrigatório:
- [ ] **Variáveis/Funções**: snake_case em inglês
- [ ] **Classes**: PascalCase em inglês  
- [ ] **Docstrings**: Google Style em inglês
- [ ] **Comentários**: Inglês explicativo
- [ ] **Mensagens**: Logs/erros em inglês
- [ ] **Metadados**: Chaves em inglês

#### Traduções Automáticas:
```python
TRANSLATIONS = {
    "algoritmo": "algorithm", "execução": "execution",
    "configuração": "configuration", "parâmetros": "parameters",
    "centro_encontrado": "center_found", "iteracoes": "iterations",
    "Iniciando": "Starting", "finalizado com sucesso": "completed successfully"
}
```

### 11.5 Qualidade
- ✅ **Code Style**: black, ruff, mypy
- ✅ **Testes**: Criar/atualizar para mudanças
- ✅ **Type Hints**: Tipagem estática apropriada

---

## 🔒 Cláusula de Imutabilidade

> **ATENÇÃO:** Estas diretrizes são **IMUTÁVEIS** até haver consenso explícito para nova versão. Qualquer proposta de alteração deve:
> 1. Ser aprovada pela equipe
> 2. Resultar em atualização deste documento
> 3. Manter compatibilidade com código existente
> 4. Incluir plano de migração se necessário

**Versão:** 0.1.0  
**Data:** Julho 2025  
**Status:** Ativo