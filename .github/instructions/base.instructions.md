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
11. [🌐 Interface Web](#11-interface-web)
12. [📄 Diretrizes para Publicação Acadêmica (JOSS/JORS)](#12-diretrizes-para-publicação-acadêmica-jossjors)
13. [🤖 Diretrizes para IAs](#13-diretrizes-para-ias)

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

### 4.2 Interface Web - Comandos Específicos

#### Inicialização da Interface Web

| Método | Comando | Uso |
|--------|---------|-----|
| **Desenvolvimento** | `.venv/bin/python src/presentation/web/run_web.py` | Local com reload |
| **Uvicorn Direto** | `.venv/bin/uvicorn src.presentation.web.app:app --host 0.0.0.0 --port 8000 --reload` | Debug/desenvolvimento |
| **Script Shell** | `./start-web.sh` | Docker/produção |
| **Docker Compose** | `docker-compose -f docker-compose.web.yml up -d` | Containerizado |

#### URLs e Acesso

| URL | Descrição | Uso |
|-----|-----------|-----|
| `http://localhost:8000` | Interface principal | Acesso browser |
| `http://localhost:8000/api/algorithms` | Lista algoritmos | API REST |
| `http://localhost:8000/api/datasets` | Lista datasets | API REST |
| `http://localhost:8000/docs` | Documentação FastAPI | API docs |

#### Gerenciamento Docker

```bash
# Iniciar interface web
docker-compose -f docker-compose.web.yml up -d

# Parar interface web
docker-compose -f docker-compose.web.yml down

# Ver logs
docker-compose -f docker-compose.web.yml logs -f

# Rebuild e restart
docker-compose -f docker-compose.web.yml up --build -d
```

### 4.3 Tasks do VS Code

| Task | Propósito | Comando |
|------|-----------|---------|
| Run CSPBench | Executar aplicação | `python main.py` |
| **Start Web Interface** | Iniciar interface web | `uvicorn src.presentation.web.app:app --reload --host 0.0.0.0 --port 8000` |
| Run Tests | Executar testes | `pytest tests/ -v` |
| Format Code | Formatar código | `black .` |
| Lint Code | Análise estática | `ruff check .` |
| Type Check | Verificação de tipos | `mypy src/` |
| Run Coverage | Cobertura de testes | `pytest --cov=src --cov-report=html tests/` |
| Build Docker Web | Build imagem web | `docker build -f Dockerfile.web -t cspbench-web:latest .` |

**Nota**: As tasks do VS Code utilizam o ambiente virtual automaticamente quando configurado corretamente.

### 4.4 Fluxo de Batch

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

## 11. 🌐 Interface Web

### 11.1 Arquitetura da Interface Web

A interface web segue os princípios da arquitetura hexagonal e está localizada em `src/presentation/web/`:

| Componente | Localização | Responsabilidade | Tecnologias |
|------------|-------------|------------------|-------------|
| **FastAPI App** | `app.py` | API REST e servir templates | FastAPI, Uvicorn |
| **Templates** | `templates/` | Interface HTML | Jinja2, HTML5 |
| **Static Assets** | `static/` | CSS, JS, imagens | CSS3, ES6+ |
| **Launcher** | `run_web.py` | Script de inicialização | Python |

### 11.2 Configuração e Inicialização

#### Métodos de Inicialização

| Método | Comando | Uso Recomendado |
|--------|---------|-----------------|
| **Desenvolvimento** | `.venv/bin/python src/presentation/web/run_web.py` | Desenvolvimento local |
| **Uvicorn Direto** | `.venv/bin/uvicorn src.presentation.web.app:app --host 0.0.0.0 --port 8000 --reload` | Debug/teste |
| **VS Code Task** | `Start Web Interface` | Desenvolvimento integrado |
| **Docker** | `./start-web.sh` ou `docker-compose -f docker-compose.web.yml up` | Produção/containerizado |

#### Configuração Obrigatória

```yaml
# config/settings.yaml - Seção web (se necessária)
web:
  host: "0.0.0.0"  # NUNCA hardcodar IPs específicos
  port: 8000       # Usar variável de ambiente PORT se disponível
  debug: false     # NUNCA true em produção
  
infrastructure:
  # Configurações existentes devem ser respeitadas
  dataset_repository:
    implementation: "FileDatasetRepository"
    config:
      base_path: "datasets/"
```

### 11.3 Estrutura de Componentes

#### Frontend Modular - Estrutura Real

```
src/presentation/web/
├── app.py                          # Aplicação FastAPI principal
├── run_web.py                      # Script de inicialização
├── README.md                       # Documentação da interface
├── templates/
│   ├── base.html                   # Layout base com CSS/JS imports
│   ├── index.html                  # Seleção de tipo de execução
│   └── single_execution.html       # Workflow de execução simples
└── static/
    ├── css/
    │   ├── style.css               # Estilos principais globais
    │   └── components/             # CSS por componente
    │       ├── algorithm-selector.css
    │       ├── dataset-selector.css
    │       ├── algorithm-config.css
    │       └── execution-monitor.css
    └── js/
        ├── app.js                  # Aplicação principal (legacy)
        └── components/             # Componentes modulares
            ├── api-client.js       # Cliente API centralizado
            ├── algorithm-selector.js # Seleção de algoritmos
            ├── dataset-selector.js   # Upload/seleção de datasets
            ├── algorithm-config.js   # Configuração de algoritmos
            ├── parameter-config.js   # Configuração de parâmetros
            ├── execution-monitor.js  # Monitor de execução
            └── results-viewer.js     # Visualização de resultados
```

#### Princípios de Organização

- ✅ **Separação CSS/JS**: Cada componente tem CSS e JS dedicados
- ✅ **Modularidade**: Componentes independentes e reutilizáveis
- ✅ **Carregamento Base**: `base.html` carrega todos os componentes
- ✅ **API Centralizada**: `api-client.js` centraliza comunicação
- ✅ **Responsivo**: CSS mobile-first em todos os componentes

#### Convenções de Naming

| Tipo | Padrão | Exemplo |
|------|--------|---------|
| **Componentes JS** | `kebab-case.js` | `algorithm-selector.js` |
| **Classes CSS** | `component-name__element` | `algorithm-selector__dropdown` |
| **CSS Files** | `component-name.css` | `algorithm-selector.css` |
| **Templates** | `snake_case.html` | `single_execution.html` |

#### Migração de app.js (Legacy)

**AÇÃO OBRIGATÓRIA**: Refatorar `app.js` monolítico para componentes modulares:

```javascript
// ❌ DEPRECADO - app.js monolítico
class CSPBenchApp {
    constructor() {
        // Lógica misturada
    }
}

// ✅ CORRETO - Componentes especializados
class AlgorithmSelector extends Component {
    constructor(container, apiClient) {
        super(container);
        this.apiClient = apiClient;
    }
}
```

### 11.4 Arquitetura de Componentes JavaScript

#### Estrutura Base de Componente

```javascript
/**
 * Base Component Class
 * Padrão para todos os componentes da interface
 */
class Component {
    constructor(container, options = {}) {
        this.container = container;
        this.options = { ...this.defaultOptions, ...options };
        this.state = {};
        this.init();
    }
    
    get defaultOptions() {
        return {};
    }
    
    init() {
        this.render();
        this.bindEvents();
    }
    
    render() {
        // Implementar em subclasses
        throw new Error('render() must be implemented');
    }
    
    bindEvents() {
        // Implementar em subclasses se necessário
    }
    
    setState(newState) {
        this.state = { ...this.state, ...newState };
        this.render();
    }
    
    destroy() {
        // Cleanup resources
    }
}
```

#### Comunicação Entre Componentes

```javascript
// ✅ CORRETO - Event-driven communication
class ComponentManager {
    constructor() {
        this.components = new Map();
        this.eventBus = new EventBus();
    }
    
    register(name, component) {
        this.components.set(name, component);
        component.eventBus = this.eventBus;
    }
    
    emit(event, data) {
        this.eventBus.emit(event, data);
    }
}

// Uso nos componentes
class AlgorithmSelector extends Component {
    onAlgorithmSelected(algorithm) {
        this.eventBus.emit('algorithm:selected', { algorithm });
    }
}
```

### 11.5 Integração com Core

#### Uso do ExperimentService

```python
# app.py - SEMPRE usar dependency injection
from src.application.services.experiment_service import ExperimentService

# ❌ PROIBIDO - Hardcoded paths
dataset_repository = FileDatasetRepository("/hardcoded/path")

# ✅ CORRETO - Usar configuração
def load_config():
    config_path = Path("config/settings.yaml")
    with open(config_path, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)
    
    dataset_repository = FileDatasetRepository(
        config["infrastructure"]["dataset_repository"]["config"]["base_path"]
    )
```

#### Descoberta de Algoritmos

```python
# ❌ PROIBIDO - Import direto de algoritmos
from algorithms.blf_ga import BLF_GA

# ✅ CORRETO - Usar global_registry
from src.domain.algorithms import global_registry

@app.get("/api/algorithms")
async def get_algorithms():
    """Lista algoritmos disponíveis via registry."""
    algorithms = []
    for name, algo_class in global_registry.get_all().items():
        algorithms.append({
            "name": name,
            "description": algo_class.__doc__ or "No description",
            "default_params": algo_class.default_params
        })
    return algorithms
```

### 11.6 APIs e Endpoints

#### Estrutura Padrão de Endpoints

| Endpoint | Método | Propósito | Resposta |
|----------|--------|-----------|----------|
| `/` | GET | Página principal | HTML |
| `/execution/{type}` | GET | Páginas de execução | HTML |
| `/api/algorithms` | GET | Lista algoritmos | JSON |
| `/api/datasets` | GET | Lista datasets | JSON |
| `/api/execute` | POST | Executa algoritmo | JSON |
| `/api/status/{session_id}` | GET | Status execução | JSON |
| `/api/results/{session_id}` | GET | Resultados | JSON |
| `/api/download/{session_id}` | GET | Download ZIP | Binary |

#### Modelos Pydantic Obrigatórios

```python
# ✅ CORRETO - Modelos tipados
class ExecutionRequest(BaseModel):
    algorithm: str
    dataset_content: Optional[str] = None
    dataset_name: Optional[str] = "uploaded_dataset.fasta"
    parameters: Dict = {}  # NUNCA hardcodar parâmetros
    save_history: bool = False
    timeout: int = 300  # Timeout configurável, não fixo

class ExecutionResult(BaseModel):
    session_id: str
    status: str
    result: Optional[Dict] = None
    error: Optional[str] = None
    download_url: Optional[str] = None

class AlgorithmInfo(BaseModel):
    name: str
    description: str
    default_params: Dict
    is_deterministic: bool = True
    supports_parallel: bool = False
    category: Optional[str] = None
```

#### Validação e Sanitização

```python
# ✅ CORRETO - Validação completa
@app.post("/api/execute")
async def execute_algorithm(request: ExecutionRequest):
    # 1. Validar algoritmo existe
    if request.algorithm not in global_registry.get_all():
        raise HTTPException(400, f"Unknown algorithm: {request.algorithm}")
    
    # 2. Validar parâmetros contra schema do algoritmo
    algo_class = global_registry.get(request.algorithm)
    validated_params = validate_algorithm_params(
        request.parameters, 
        algo_class.default_params
    )
    
    # 3. Sanitizar nome do dataset
    safe_filename = sanitize_filename(request.dataset_name)
    
    # 4. Validar tamanho do dataset
    if request.dataset_content and len(request.dataset_content) > MAX_DATASET_SIZE:
        raise HTTPException(413, "Dataset too large")
    
    # 5. Executar com timeout
    return await execute_with_timeout(request, validated_params, safe_filename)
```

### 11.7 Gerenciamento de Estado e Sessões

#### Sessões e Temporários

```python
# ✅ CORRETO - Usar SessionManager
from src.infrastructure.orchestrators.session_manager import SessionManager

class WebSessionManager:
    def __init__(self, config):
        self.session_manager = SessionManager(config)
        # NUNCA hardcodar paths temporários
        self.temp_dir = Path(config.get("temp_dir", "/tmp"))
    
    def create_session(self) -> str:
        """Cria sessão web com ID único."""
        session_id = str(uuid.uuid4())
        # Usar infraestrutura existente
        return self.session_manager.create_session(session_id)
```

### 11.7 Gerenciamento de Estado e Sessões

#### WebSocket para Execução Real-time

```python
# ✅ NOVO - WebSocket para progresso em tempo real
from fastapi import WebSocket

@app.websocket("/ws/{session_id}")
async def websocket_endpoint(websocket: WebSocket, session_id: str):
    await websocket.accept()
    
    # Configurar callback de progresso
    def progress_callback(message: str, progress: float):
        asyncio.create_task(
            websocket.send_json({
                "type": "progress",
                "message": message,
                "progress": progress,
                "timestamp": datetime.now().isoformat()
            })
        )
    
    # Executar algoritmo com callback
    await execute_algorithm_with_progress(session_id, progress_callback)
```

#### Sessões Web com Cleanup Automático

```python
# ✅ CORRETO - Gerenciamento de sessão robusto
class WebSessionManager:
    def __init__(self, config):
        self.session_manager = SessionManager(config)
        self.temp_dir = Path(config.get("temp_dir", "/tmp/cspbench"))
        self.max_session_age = config.get("max_session_age_minutes", 60)
        self.cleanup_interval = config.get("cleanup_interval_minutes", 10)
        
        # Agendar limpeza automática
        self.schedule_cleanup()
    
    async def create_session(self) -> str:
        session_id = str(uuid.uuid4())
        session_dir = self.temp_dir / session_id
        session_dir.mkdir(parents=True, exist_ok=True)
        
        # Registrar sessão com timestamp
        await self.session_manager.create_session(session_id, {
            "created_at": datetime.now(),
            "work_dir": str(session_dir),
            "status": "created"
        })
        
        return session_id
    
    async def cleanup_expired_sessions(self):
        """Remove sessões expiradas automaticamente."""
        cutoff = datetime.now() - timedelta(minutes=self.max_session_age)
        expired_sessions = await self.session_manager.get_expired_sessions(cutoff)
        
        for session_id in expired_sessions:
            await self.cleanup_session(session_id)
```

#### Download Streaming e Compressão

```python
# ✅ MELHORADO - Download otimizado
from fastapi.responses import StreamingResponse
import aiofiles

@app.get("/api/download/{session_id}")
async def download_results_stream(session_id: str):
    """Download ZIP com streaming para arquivos grandes."""
    session_info = await session_manager.get_session(session_id)
    if not session_info:
        raise HTTPException(404, "Session not found")
    
    zip_path = await create_results_zip_async(session_id, session_info)
    
    async def stream_file():
        async with aiofiles.open(zip_path, 'rb') as f:
            chunk_size = 8192
            while chunk := await f.read(chunk_size):
                yield chunk
    
    return StreamingResponse(
        stream_file(),
        media_type="application/zip",
        headers={
            "Content-Disposition": f"attachment; filename=cspbench_results_{session_id}.zip",
            "Content-Length": str(zip_path.stat().st_size)
        }
    )

async def create_results_zip_async(session_id: str, session_info: Dict) -> Path:
    """Criar ZIP assíncrono com estrutura inteligente."""
    work_dir = Path(session_info["work_dir"])
    zip_path = work_dir / f"results_{session_id}.zip"
    
    async with aiofiles.open(zip_path, 'wb') as zip_file:
        with zipfile.ZipFile(zip_file, 'w', zipfile.ZIP_DEFLATED, compresslevel=6) as zipf:
            # Estrutura baseada no conteúdo real da sessão
            results_file = work_dir / "results.json"
            if results_file.exists():
                zipf.write(results_file, "results.json")
            
            # Adicionar logs se existirem
            log_file = work_dir / "execution.log"
            if log_file.exists():
                zipf.write(log_file, "execution.log")
            
            # Adicionar outputs específicos do algoritmo
            output_dir = work_dir / "algorithm_output"
            if output_dir.exists():
                for file_path in output_dir.rglob("*"):
                    if file_path.is_file():
                        arc_name = f"algorithm_output/{file_path.relative_to(output_dir)}"
                        zipf.write(file_path, arc_name)
    
    return zip_path
```

### 11.8 Containerização e Deploy

### 11.8 Containerização e Deploy

#### Dockerfile Multi-Stage Otimizado

```dockerfile
# Dockerfile.web - Produção otimizada
FROM python:3.11-slim as builder

# Build dependencies
RUN apt-get update && apt-get install -y \
    --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build
COPY requirements.web.txt .
RUN pip install --user --no-cache-dir -r requirements.web.txt

# Production stage
FROM python:3.11-slim as production

# Runtime dependencies apenas
RUN apt-get update && apt-get install -y \
    --no-install-recommends \
    curl \
    && rm -rf /var/lib/apt/lists/* \
    && groupadd -r cspbench \
    && useradd -r -g cspbench cspbench

# Copy Python packages from builder
COPY --from=builder /root/.local /home/cspbench/.local
ENV PATH=/home/cspbench/.local/bin:$PATH

WORKDIR /app
COPY --chown=cspbench:cspbench . .

# ✅ CORRETO - Configuração via ENV com fallbacks
ENV PYTHONPATH=/app
ENV WEB_HOST=0.0.0.0
ENV WEB_PORT=8000
ENV WORKERS=1
ENV MAX_SESSION_AGE_MINUTES=60

USER cspbench
EXPOSE 8000

# Health check robusto
HEALTHCHECK --interval=30s --timeout=10s --start-period=40s --retries=3 \
  CMD curl -f http://localhost:${WEB_PORT}/api/health || exit 1

CMD ["uvicorn", "src.presentation.web.app:app", "--host", "0.0.0.0", "--port", "8000"]
```

#### Docker Compose para Diferentes Ambientes

```yaml
# docker-compose.web.yml - Produção
version: '3.8'

services:
  cspbench-web:
    build:
      context: .
      dockerfile: Dockerfile.web
      target: production
    ports:
      - "${WEB_PORT:-8000}:8000"
    environment:
      # ✅ CORRETO - ENV vars com fallbacks seguros
      - PYTHONPATH=/app
      - WEB_HOST=0.0.0.0
      - WEB_PORT=8000
      - WORKERS=${WEB_WORKERS:-1}
      - MAX_SESSION_AGE_MINUTES=${SESSION_AGE:-60}
      - CLEANUP_INTERVAL_MINUTES=${CLEANUP_INTERVAL:-10}
      # Credenciais via env vars
      - NCBI_EMAIL=${NCBI_EMAIL:-}
      - NCBI_API_KEY=${NCBI_API_KEY:-}
    volumes:
      - ./datasets:/app/datasets:ro
      - ./config:/app/config:ro
      - ./outputs:/app/outputs
      - web_sessions:/app/tmp
    restart: unless-stopped
    networks:
      - cspbench-network
    labels:
      - "traefik.enable=true"
      - "traefik.http.routers.cspbench.rule=Host(`${DOMAIN:-cspbench.local}`)"

  # Opcional: Redis para sessões distribuídas
  redis:
    image: redis:7-alpine
    volumes:
      - redis_data:/data
    networks:
      - cspbench-network
    command: redis-server --appendonly yes

volumes:
  web_sessions:
  redis_data:

networks:
  cspbench-network:
    driver: bridge
```

### 11.9 Segurança e Performance

### 11.9 Segurança e Performance

#### Rate Limiting e Throttling

```python
# ✅ NOVO - Rate limiting por IP e endpoint
from slowapi import Limiter, _rate_limit_exceeded_handler
from slowapi.util import get_remote_address
from slowapi.errors import RateLimitExceeded

limiter = Limiter(key_func=get_remote_address)
app.state.limiter = limiter
app.add_exception_handler(RateLimitExceeded, _rate_limit_exceeded_handler)

@app.post("/api/execute")
@limiter.limit("3/minute")  # Max 3 execuções por minuto por IP
async def execute_algorithm(request: Request, execution_request: ExecutionRequest):
    """Execute algorithm with rate limiting."""
    pass

@app.post("/api/upload")
@limiter.limit("10/minute")  # Max 10 uploads por minuto por IP
async def upload_dataset(request: Request, file: UploadFile):
    """Upload dataset with rate limiting."""
    pass
```

#### Validação e Sanitização Avançada

```python
# ✅ CORRETO - Validação multicamadas
import bleach
from pathlib import Path

class SecurityConfig:
    MAX_DATASET_SIZE = 50 * 1024 * 1024  # 50MB
    MAX_FILENAME_LENGTH = 255
    ALLOWED_EXTENSIONS = {'.fasta', '.fa', '.txt'}
    BLOCKED_PATTERNS = ['..', '/', '\\', '<', '>', '|', '*', '?']

def sanitize_filename(filename: str) -> str:
    """Sanitizar nome de arquivo de forma rigorosa."""
    # Remover caracteres perigosos
    safe_chars = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-"
    sanitized = "".join(c for c in filename if c in safe_chars)
    
    # Limitar tamanho
    if len(sanitized) > SecurityConfig.MAX_FILENAME_LENGTH:
        name, ext = sanitized.rsplit('.', 1) if '.' in sanitized else (sanitized, '')
        sanitized = name[:SecurityConfig.MAX_FILENAME_LENGTH-len(ext)-1] + '.' + ext
    
    return sanitized or "uploaded_file.fasta"

def validate_dataset_content(content: str) -> bool:
    """Validar conteúdo do dataset."""
    if not content or len(content) > SecurityConfig.MAX_DATASET_SIZE:
        return False
    
    # Verificar se é FASTA válido básico
    lines = content.strip().split('\n')
    has_header = any(line.startswith('>') for line in lines)
    has_sequence = any(line and not line.startswith('>') for line in lines)
    
    return has_header and has_sequence

@app.post("/api/execute")
async def execute_algorithm_secure(request: ExecutionRequest):
    # 1. Validar algoritmo
    if request.algorithm not in global_registry.get_all():
        raise HTTPException(400, "Invalid algorithm")
    
    # 2. Sanitizar e validar dataset
    if request.dataset_content:
        if not validate_dataset_content(request.dataset_content):
            raise HTTPException(400, "Invalid dataset format")
    
    # 3. Sanitizar filename
    safe_filename = sanitize_filename(request.dataset_name)
    
    # 4. Validar parâmetros contra schema
    algo_class = global_registry.get(request.algorithm)
    validated_params = validate_algorithm_parameters(
        request.parameters, 
        algo_class.default_params
    )
    
    return await execute_with_security_context(request, validated_params, safe_filename)
```

#### Monitoramento e Observabilidade

```python
# ✅ NOVO - Métricas e monitoring
import time
from prometheus_client import Counter, Histogram, generate_latest

# Métricas
REQUEST_COUNT = Counter('web_requests_total', 'Total requests', ['method', 'endpoint'])
REQUEST_DURATION = Histogram('web_request_duration_seconds', 'Request duration')
EXECUTION_COUNT = Counter('algorithm_executions_total', 'Total executions', ['algorithm'])
EXECUTION_DURATION = Histogram('algorithm_execution_seconds', 'Execution duration', ['algorithm'])

@app.middleware("http")
async def add_process_time_header(request: Request, call_next):
    start_time = time.time()
    response = await call_next(request)
    process_time = time.time() - start_time
    
    # Registrar métricas
    REQUEST_COUNT.labels(method=request.method, endpoint=request.url.path).inc()
    REQUEST_DURATION.observe(process_time)
    
    response.headers["X-Process-Time"] = str(process_time)
    return response

@app.get("/api/metrics")
async def get_metrics():
    """Endpoint de métricas para Prometheus."""
    return Response(generate_latest(), media_type="text/plain")

@app.get("/api/health")
async def health_check():
    """Health check detalhado."""
    try:
        # Verificar componentes essenciais
        algorithms_available = len(global_registry.get_all()) > 0
        config_loaded = config is not None
        
        status = "healthy" if algorithms_available and config_loaded else "degraded"
        
        return {
            "status": status,
            "timestamp": datetime.now().isoformat(),
            "components": {
                "algorithms": algorithms_available,
                "config": config_loaded,
                "session_manager": session_manager is not None
            }
        }
    except Exception as e:
        raise HTTPException(503, f"Health check failed: {str(e)}")
```

### 11.10 Testing e Quality Assurance

### 11.10 Testing e Quality Assurance

#### Testes Frontend

```javascript
// ✅ CORRETO - Testes unitários de componentes
// tests/web/components/test-algorithm-selector.js
import { AlgorithmSelector } from '../../../src/presentation/web/static/js/components/algorithm-selector.js';

describe('AlgorithmSelector', () => {
    let container, apiClient, selector;
    
    beforeEach(() => {
        container = document.createElement('div');
        apiClient = { getAlgorithms: jest.fn() };
        selector = new AlgorithmSelector(container, apiClient);
    });
    
    test('should load algorithms on init', async () => {
        const mockAlgorithms = [
            { name: 'BLF_GA', description: 'Test algorithm' }
        ];
        apiClient.getAlgorithms.mockResolvedValue(mockAlgorithms);
        
        await selector.loadAlgorithms();
        
        expect(apiClient.getAlgorithms).toHaveBeenCalled();
        expect(selector.algorithms).toEqual(mockAlgorithms);
    });
});
```

#### Testes E2E

```python
# ✅ CORRETO - Testes end-to-end
# tests/web/test_web_integration.py
import pytest
from fastapi.testclient import TestClient
from src.presentation.web.app import app

@pytest.fixture
def client():
    return TestClient(app)

def test_web_interface_loads(client):
    """Test that main page loads correctly."""
    response = client.get("/")
    assert response.status_code == 200
    assert "CSPBench Web Interface" in response.text

def test_algorithms_api(client):
    """Test algorithms API endpoint."""
    response = client.get("/api/algorithms")
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, list)
    assert len(data) > 0

def test_execute_algorithm_validation(client):
    """Test algorithm execution with validation."""
    # Invalid algorithm
    response = client.post("/api/execute", json={
        "algorithm": "NonExistentAlgorithm",
        "dataset_content": ">test\nATCG"
    })
    assert response.status_code == 400
    
    # Valid execution
    response = client.post("/api/execute", json={
        "algorithm": "Baseline",
        "dataset_content": ">seq1\nATCG\n>seq2\nATGC",
        "parameters": {}
    })
    assert response.status_code == 200
```

#### Performance Testing

```python
# ✅ CORRETO - Testes de performance
# tests/web/test_performance.py
import asyncio
import time
from fastapi.testclient import TestClient

async def test_concurrent_executions():
    """Test multiple concurrent algorithm executions."""
    client = TestClient(app)
    
    async def execute_algorithm():
        start_time = time.time()
        response = client.post("/api/execute", json={
            "algorithm": "Baseline",
            "dataset_content": ">seq1\nATCG\n>seq2\nATGC"
        })
        duration = time.time() - start_time
        return response.status_code, duration
    
    # Execute 5 concurrent requests
    tasks = [execute_algorithm() for _ in range(5)]
    results = await asyncio.gather(*tasks)
    
    # All should succeed
    for status_code, duration in results:
        assert status_code == 200
        assert duration < 30  # Should complete within 30 seconds
```

#### Browser Testing

```javascript
// ✅ CORRETO - Testes de browser
// tests/web/browser/test-user-workflow.js
const { chromium } = require('playwright');

describe('User Workflow', () => {
    let browser, page;
    
    beforeAll(async () => {
        browser = await chromium.launch();
    });
    
    afterAll(async () => {
        await browser.close();
    });
    
    beforeEach(async () => {
        page = await browser.newPage();
        await page.goto('http://localhost:8000');
    });
    
    test('complete algorithm execution workflow', async () => {
        // 1. Navigate to single execution
        await page.click('a[href="/execution/single"]');
        
        // 2. Select algorithm
        await page.selectOption('#algorithm-select', 'Baseline');
        
        // 3. Upload dataset
        const fileInput = page.locator('input[type="file"]');
        await fileInput.setInputFiles('./tests/fixtures/test.fasta');
        
        // 4. Execute
        await page.click('#execute-button');
        
        // 5. Wait for results
        await page.waitForSelector('#results-section', { timeout: 30000 });
        
        // 6. Verify results displayed
        const results = await page.textContent('#results-section');
        expect(results).toContain('Best String');
        expect(results).toContain('Maximum Distance');
    });
});
```

## 12. 📄 Diretrizes para Publicação Acadêmica (JOSS/JORS)

### 12.1 Padrões de Qualidade para Publicação

#### Checklist JOSS (Journal of Open Source Software)

| Critério | Requisito | Status CSPBench |
|----------|-----------|-----------------|
| **Software Functionality** | Software deve ter funcionalidade clara e bem definida | ✅ Benchmarking CSP |
| **Open Source License** | Licença OSI aprovada | ✅ MIT License |
| **Documentation** | README claro, documentação de API, tutoriais | ✅ Documentação completa |
| **Tests** | Testes automatizados com cobertura adequada | ✅ pytest + coverage |
| **Installation Instructions** | Instruções claras de instalação | ✅ Requirements + Docker |
| **Example Usage** | Exemplos funcionais de uso | ✅ Datasets exemplo + CLI |
| **Community Guidelines** | CONTRIBUTING.md, CODE_OF_CONDUCT.md | ✅ Implementado |
| **Sustainable Development** | Desenvolvimento ativo e sustentável | ✅ Arquitetura modular |

#### Checklist JORS (Journal of Open Research Software)

| Critério | Requisito | Implementação |
|----------|-----------|---------------|
| **Research Purpose** | Software deve ter propósito de pesquisa científica | ✅ Algoritmos CSP para pesquisa |
| **Technical Documentation** | Documentação técnica detalhada | ✅ Documentação por algoritmo |
| **Reproducibility** | Resultados reproduzíveis | ✅ Seeds determinísticos |
| **Citation Information** | CITATION.cff, DOI, metadados | ✅ CITATION.cff presente |
| **Performance Evaluation** | Benchmarks e métricas | ✅ Sistema de métricas |
| **Extensibility** | Facilidade para extensão | ✅ Plugin system |

### 12.2 Estrutura de Documentação Acadêmica

#### Arquivo paper.md (JOSS)

```markdown
# CSPBench: A Comprehensive Benchmarking Framework for Closest String Problem Algorithms

## Summary

CSPBench is a Python framework for benchmarking algorithms that solve the Closest String Problem (CSP). 
It provides a plugin-based architecture that allows researchers to easily implement, compare, and evaluate 
different CSP algorithms using standardized datasets and metrics.

## Statement of Need

The Closest String Problem is fundamental in bioinformatics applications such as primer design, 
motif finding, and consensus sequence generation. However, comparing algorithm performance across 
different implementations has been challenging due to lack of standardized benchmarking tools.

## Key Features

- Plugin-based algorithm architecture
- Standardized benchmarking metrics
- Web interface for algorithm execution
- Reproducible experimental setup
- Docker containerization for portability

## Research Applications

CSPBench has been designed to facilitate research in:
- Comparative analysis of CSP algorithms
- Development of new heuristic approaches
- Performance evaluation on biological datasets
- Educational use in bioinformatics courses
```

#### Arquivo CITATION.cff

```yaml
# CITATION.cff - Citation File Format
cff-version: 1.2.0
title: "CSPBench: A Comprehensive Benchmarking Framework for Closest String Problem Algorithms"
message: "If you use this software, please cite it using the metadata from this file."
type: software
authors:
  - given-names: "Diego"
    family-names: "Grosmann"
    email: "diego.grosmann@example.com"
    orcid: "https://orcid.org/0000-0000-0000-0000"
repository-code: "https://github.com/diegogrosmann/CSPBench"
url: "https://cspbench.readthedocs.io"
abstract: >
  CSPBench is a Python framework for benchmarking algorithms that solve 
  the Closest String Problem (CSP). It provides a plugin-based architecture 
  that allows researchers to easily implement, compare, and evaluate 
  different CSP algorithms.
keywords:
  - bioinformatics
  - closest-string-problem
  - benchmarking
  - algorithms
  - computational-biology
license: MIT
version: "0.1.0"
date-released: "2025-07-30"
```

### 12.3 Padrões de Código para Publicação

#### Documentação de Algoritmos

```python
# ✅ PADRÃO ACADÊMICO - Documentação completa
@register_algorithm
class MyAlgorithm(CSPAlgorithm):
    """
    Implementation of [Algorithm Name] for the Closest String Problem.
    
    This algorithm implements the approach described in [Reference Paper].
    
    References:
        Author, A. et al. (Year). "Paper Title". Journal Name, Volume(Issue), pages.
        DOI: https://doi.org/10.xxxx/xxxxx
    
    Time Complexity:
        O(n*m*k) where n=number of strings, m=string length, k=parameter
        
    Space Complexity:
        O(n*m) for storing input and intermediate results
        
    Parameters:
        population_size (int): Size of the population (default: 100)
            Valid range: [10, 1000]
            
        max_generations (int): Maximum number of generations (default: 500)
            Valid range: [1, 10000]
            
        mutation_rate (float): Probability of mutation (default: 0.1)
            Valid range: [0.0, 1.0]
    
    Returns:
        tuple: (best_string, max_distance, metadata)
            best_string (str): The closest string found
            max_distance (int): Maximum Hamming distance to input strings
            metadata (dict): Algorithm-specific information including:
                - 'iterations': Number of iterations performed
                - 'convergence_generation': Generation where best solution was found
                - 'execution_time': Runtime in seconds
                - 'memory_usage': Peak memory usage in MB
    
    Examples:
        >>> algorithm = MyAlgorithm()
        >>> dataset = load_dataset("example.fasta")
        >>> result = algorithm.run(dataset, population_size=50)
        >>> best_string, distance, metadata = result
        >>> print(f"Best string: {best_string}, Distance: {distance}")
    
    Note:
        This implementation is deterministic when seed is set.
        For reproducible results, use the same seed value.
    """
    
    name = "MyAlgorithm"
    default_params = {
        "population_size": 100,
        "max_generations": 500,
        "mutation_rate": 0.1,
        "seed": 42
    }
```

#### Testes para Reprodutibilidade

```python
# ✅ PADRÃO ACADÊMICO - Testes de reprodutibilidade
class TestAlgorithmReproducibility:
    """Test suite ensuring algorithmic reproducibility for academic publication."""
    
    @pytest.mark.parametrize("algorithm_name", ["BLF_GA", "CSC", "DP_CSP"])
    def test_deterministic_behavior(self, algorithm_name):
        """Test that algorithm produces identical results with same seed."""
        algorithm_class = global_registry.get(algorithm_name)
        algorithm1 = algorithm_class()
        algorithm2 = algorithm_class()
        
        dataset = create_test_dataset(num_strings=5, string_length=10)
        params = {"seed": 42}
        
        result1 = algorithm1.run(dataset, **params)
        result2 = algorithm2.run(dataset, **params)
        
        assert result1[0] == result2[0]  # Same best string
        assert result1[1] == result2[1]  # Same distance
        assert result1[2]["iterations"] == result2[2]["iterations"]  # Same iterations
    
    def test_algorithm_metadata_completeness(self):
        """Ensure all algorithms return complete metadata for academic analysis."""
        required_metadata_keys = [
            "iterations", "execution_time", "convergence_info", 
            "algorithm_version", "parameters_used"
        ]
        
        for name, algo_class in global_registry.get_all().items():
            algorithm = algo_class()
            dataset = create_test_dataset(num_strings=3, string_length=8)
            
            _, _, metadata = algorithm.run(dataset)
            
            for key in required_metadata_keys:
                assert key in metadata, f"Algorithm {name} missing metadata key: {key}"
    
    def test_performance_benchmarking(self):
        """Test performance characteristics for academic comparison."""
        algorithm_class = global_registry.get("BLF_GA")
        algorithm = algorithm_class()
        
        # Test scalability
        for size in [5, 10, 20]:
            dataset = create_test_dataset(num_strings=size, string_length=15)
            start_time = time.time()
            
            _, _, metadata = algorithm.run(dataset, max_generations=100)
            
            execution_time = time.time() - start_time
            
            # Ensure reasonable performance
            assert execution_time < 60, f"Algorithm too slow for size {size}"
            assert "execution_time" in metadata
            assert abs(metadata["execution_time"] - execution_time) < 1.0
```

### 12.4 Métricas e Avaliação Científica

#### Sistema de Métricas Padronizado

```python
# ✅ PADRÃO ACADÊMICO - Métricas científicas
class ScientificMetrics:
    """
    Standardized metrics for academic evaluation of CSP algorithms.
    
    Implements metrics commonly used in CSP literature for fair comparison.
    """
    
    @staticmethod
    def solution_quality_score(best_string: str, input_strings: List[str]) -> Dict[str, float]:
        """
        Calculate comprehensive solution quality metrics.
        
        Returns:
            dict: Metrics including hamming distances, consensus score, entropy
        """
        distances = [hamming_distance(best_string, s) for s in input_strings]
        
        return {
            "max_distance": max(distances),
            "avg_distance": np.mean(distances),
            "std_distance": np.std(distances),
            "consensus_score": calculate_consensus_score(best_string, input_strings),
            "solution_entropy": calculate_sequence_entropy(best_string)
        }
    
    @staticmethod
    def algorithm_efficiency_metrics(metadata: Dict) -> Dict[str, float]:
        """
        Calculate algorithm efficiency metrics for academic comparison.
        
        Returns:
            dict: Efficiency metrics including convergence rate, resource usage
        """
        return {
            "convergence_rate": metadata.get("convergence_generation", 0) / metadata.get("iterations", 1),
            "iterations_per_second": metadata.get("iterations", 0) / metadata.get("execution_time", 1),
            "memory_efficiency": metadata.get("peak_memory_mb", 0) / metadata.get("dataset_size_mb", 1),
            "solution_stability": calculate_solution_stability(metadata)
        }
    
    @staticmethod
    def generate_comparison_report(results: Dict[str, Dict]) -> str:
        """
        Generate standardized comparison report for academic publication.
        
        Args:
            results: Dict mapping algorithm names to their results
            
        Returns:
            str: Formatted comparison report with statistical analysis
        """
        report = []
        report.append("# Algorithm Comparison Report")
        report.append("## Methodology")
        report.append("All algorithms tested on identical datasets with fixed seeds for reproducibility.")
        
        # Statistical analysis
        quality_scores = {}
        efficiency_scores = {}
        
        for algo_name, result in results.items():
            best_string, distance, metadata = result
            quality_scores[algo_name] = ScientificMetrics.solution_quality_score(
                best_string, metadata["input_strings"]
            )
            efficiency_scores[algo_name] = ScientificMetrics.algorithm_efficiency_metrics(metadata)
        
        # Generate comparative tables
        report.append("## Solution Quality Comparison")
        report.append(generate_quality_table(quality_scores))
        
        report.append("## Efficiency Comparison")
        report.append(generate_efficiency_table(efficiency_scores))
        
        # Statistical significance tests
        report.append("## Statistical Analysis")
        report.append(perform_statistical_tests(quality_scores, efficiency_scores))
        
        return "\n".join(report)
```

### 12.5 Reprodutibilidade e Versionamento

#### Controle de Versão Científico

```python
# ✅ PADRÃO ACADÊMICO - Versionamento para reprodutibilidade
class ReproducibilityManager:
    """
    Manage reproducibility requirements for academic publication.
    """
    
    @staticmethod
    def create_experiment_manifest(algorithm_name: str, dataset_name: str, 
                                  parameters: Dict) -> Dict:
        """
        Create detailed manifest for experiment reproduction.
        
        Returns:
            dict: Complete experiment specification including versions
        """
        return {
            "experiment_id": generate_experiment_id(),
            "timestamp": datetime.now().isoformat(),
            "cspbench_version": get_cspbench_version(),
            "python_version": sys.version,
            "algorithm": {
                "name": algorithm_name,
                "version": get_algorithm_version(algorithm_name),
                "parameters": parameters,
                "implementation_hash": get_algorithm_hash(algorithm_name)
            },
            "dataset": {
                "name": dataset_name,
                "hash": calculate_dataset_hash(dataset_name),
                "size": get_dataset_size(dataset_name)
            },
            "environment": {
                "os": platform.platform(),
                "cpu_info": get_cpu_info(),
                "memory_gb": get_total_memory_gb(),
                "dependencies": get_dependency_versions()
            }
        }
    
    @staticmethod
    def validate_reproduction(original_manifest: Dict, current_result: Dict) -> bool:
        """
        Validate that current execution matches original experiment.
        
        Returns:
            bool: True if reproduction is valid within tolerances
        """
        # Verify versions match
        if original_manifest["cspbench_version"] != get_cspbench_version():
            logger.warning("CSPBench version mismatch")
        
        # Verify algorithm implementation hasn't changed
        current_hash = get_algorithm_hash(original_manifest["algorithm"]["name"])
        if original_manifest["algorithm"]["implementation_hash"] != current_hash:
            logger.error("Algorithm implementation has changed")
            return False
        
        # Verify dataset integrity
        current_dataset_hash = calculate_dataset_hash(original_manifest["dataset"]["name"])
        if original_manifest["dataset"]["hash"] != current_dataset_hash:
            logger.error("Dataset has been modified")
            return False
        
        return True
```

### 12.6 Estrutura de Arquivos para Publicação

#### Organização de Projeto Acadêmico

```
project_root/
├── paper/                          # Artigos e documentação acadêmica
│   ├── paper.md                   # Artigo principal (JOSS)
│   ├── paper.bib                  # Bibliografia
│   ├── figures/                   # Figuras para o artigo
│   └── supplementary/             # Material suplementar
├── CITATION.cff                   # Metadados de citação
├── CODE_OF_CONDUCT.md            # Código de conduta
├── CONTRIBUTING.md               # Guidelines para contribuição
├── LICENSE                       # Licença MIT
├── README.md                     # Documentação principal
├── CHANGELOG.md                  # Histórico de mudanças
├── benchmarks/                   # Resultados de benchmarks
│   ├── datasets/                 # Datasets de teste padronizados
│   ├── results/                  # Resultados reproduzíveis
│   └── analysis/                 # Análises estatísticas
├── docs/                         # Documentação completa
│   ├── api/                      # Documentação da API
│   ├── tutorials/                # Tutoriais passo a passo
│   ├── algorithms/               # Documentação por algoritmo
│   └── reproducibility/          # Guias de reprodução
└── tests/                        # Testes para reprodutibilidade
    ├── reproducibility/          # Testes de reprodução
    ├── benchmarks/               # Testes de performance
    └── academic/                 # Testes específicos para publicação
```

#### Arquivo de Configuração Acadêmica

```yaml
# config/academic.yaml - Configurações para publicação acadêmica
academic:
  reproducibility:
    enforce_deterministic: true
    require_seed: true
    validate_versions: true
    
  metrics:
    standard_datasets: 
      - "bio_test_set_1.fasta"
      - "synthetic_benchmark.fasta"
      - "real_world_dataset.fasta"
    
    required_metrics:
      - "solution_quality"
      - "execution_time"
      - "memory_usage"
      - "convergence_analysis"
    
  publication:
    joss:
      require_tests: true
      coverage_threshold: 0.85
      documentation_completeness: true
      
    jors:
      require_benchmarks: true
      performance_analysis: true
      comparative_evaluation: true
      
  export:
    formats: ["csv", "json", "latex"]
    include_metadata: true
    statistical_analysis: true
```

## 13. 🤖 Diretrizes para IAs

### 13.1 Comunicação
- ✅ **Resposta**: Sempre em português brasileiro
- ✅ **Clareza**: Manter precisão técnica
- ✅ **Contexto**: Considerar arquitetura existente

### 13.2 Implementação
- ✅ **Ambiente Virtual**: Sempre usar `.venv/bin/python` para execução
- ✅ **Ferramentas Internas**: Priorizar APIs/funções internas vs comandos externos
- ✅ **Incrementalidade**: Mudanças pequenas e testáveis
- ✅ **Validação**: Testar após mudanças significativas
- ✅ **Cleanup**: Remover código obsoleto

### 13.3 Restrições Críticas
- ❌ **Dependências Cruzadas**: Plugins NÃO podem depender uns dos outros
- ❌ **Imports Diretos**: Aplicação NÃO pode importar plugins diretamente
- ❌ **Domain Purity**: Domain layer NÃO pode ter I/O ou dependências externas
- ❌ **Credenciais**: NUNCA hardcodar ou versionar credenciais
- ❌ **Paths Hardcoded**: NUNCA usar paths absolutos ou específicos de máquina
- ❌ **IPs/Hosts Fixos**: NUNCA hardcodar endereços IP ou hostnames
- ❌ **Dados Arbitrários**: NUNCA inserir dados de exemplo fixos sem configuração

### 13.4 Interface Web - Restrições Específicas

#### PROIBIÇÕES ABSOLATAS:
- ❌ **Hardcoded Paths**: `/home/user/data`, `C:\Users\John`, etc.
- ❌ **Hardcoded IPs**: `192.168.1.100`, `10.0.0.5`, etc.
- ❌ **Credenciais Fixas**: API keys, senhas, tokens reais
- ❌ **Dados Mock Fixos**: Resultados fake, datasets exemplo hardcoded
- ❌ **Configuração Inline**: Parâmetros direto no código vs arquivo config

#### PRÁTICAS OBRIGATÓRIAS:
- ✅ **Configuração Externa**: Tudo via `config/settings.yaml` ou ENV vars
- ✅ **Paths Relativos**: Sempre relativo ao project root
- ✅ **Fallbacks Seguros**: Valores padrão não sensíveis
- ✅ **Validação**: Input validation em todos os endpoints
- ✅ **Tipo Safety**: Pydantic models para todas as APIs

### 13.5 Padrões Acadêmicos (JOSS/JORS)

#### REQUISITOS OBRIGATÓRIOS:
- ✅ **Documentação Científica**: Docstrings com referências, complexidade, exemplos
- ✅ **Reprodutibilidade**: Seeds determinísticos, controle de versão
- ✅ **Métricas Padronizadas**: Sistema de métricas comparáveis
- ✅ **Testes de Reprodução**: Validação de resultados idênticos
- ✅ **Metadados Completos**: CITATION.cff, paper.md, benchmarks

#### DIRETRIZES DE IMPLEMENTAÇÃO:
- ✅ **Algoritmos**: Documentar referências, complexidade, parâmetros
- ✅ **Experimentos**: Manifest completo para reprodução
- ✅ **Análises**: Comparações estatísticas e métricas científicas
- ✅ **Versionamento**: Hash de implementações, controle de dataset
- ✅ **Relatórios**: Formato acadêmico com análise estatística

### 13.6 Internacionalização Obrigatória

**AÇÃO AUTOMÁTICA**: Sempre traduzir português → inglês ao modificar código.

#### Checklist Obrigatório:
- [ ] **Variáveis/Funções**: snake_case em inglês
- [ ] **Classes**: PascalCase em inglês  
- [ ] **Docstrings**: Google Style em inglês
- [ ] **Comentários**: Inglês explicativo
- [ ] **Mensagens**: Logs/erros em inglês
- [ ] **Metadados**: Chaves em inglês
- [ ] **API Responses**: Todas em inglês
- [ ] **HTML/JS**: Variáveis e funções em inglês

#### Traduções Automáticas:
```python
TRANSLATIONS = {
    "algoritmo": "algorithm", "execução": "execution",
    "configuração": "configuration", "parâmetros": "parameters",
    "centro_encontrado": "center_found", "iteracoes": "iterations",
    "Iniciando": "Starting", "finalizado com sucesso": "completed successfully",
    "interface_web": "web_interface", "servidor": "server",
    "requisição": "request", "resposta": "response"
}
```

### 13.7 Qualidade
- ✅ **Code Style**: black, ruff, mypy
- ✅ **Testes**: Criar/atualizar para mudanças
- ✅ **Type Hints**: Tipagem estática apropriada
- ✅ **Web Standards**: HTML5, CSS3, ES6+ compatível
- ✅ **Academic Standards**: JOSS/JORS compliance

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