---
applyTo: '**'
---

# 📋 Diretrizes CSPBench v0.1.0

> **🔒 IMUTÁVEL:** Diretrizes fundamentais para todas as interações com IAs. Alterações devem ser aprovadas.

## 🏗️ Arquitetura

**Arquitetura Hexagonal (Clean Architecture)**
- **Domain** (`src/domain/`): Regras de negócio - APENAS Python StdLib
- **Application** (`src/application/`): Casos de uso
- **Infrastructure** (`src/infrastructure/`): Adaptadores, I/O
- **Presentation** (`src/presentation/`): CLI, TUI, Web
- **Plugins** (`algorithms/`): Algoritmos via `@register_algorithm`

**Ver detalhes:** `src/domain/algorithms.py`, `README.md`

## 🧩 Sistema de Plugins

**Estrutura obrigatória:**
```
algorithms/nome_algoritmo/
├── __init__.py
├── algorithm.py         # @register_algorithm
├── implementation.py    # Lógica core
├── config.py           # Parâmetros
└── README.md
```

**Template mínimo:**
```python
@register_algorithm
class MeuAlgoritmo(CSPAlgorithm):
    name = "MeuAlgoritmo"
    def run(self, dataset, **kwargs):
        return implementacao.solve(dataset, **kwargs)
```

**Ver exemplos:** `algorithms/baseline/`, `algorithms/blf_ga/`

## ⚙️ Ambiente

**Python Virtual Environment:**
- Executável: `.venv/bin/python`
- Pip: `.venv/bin/pip`
- **SEMPRE** usar `.venv/bin/python` (não apenas `python`)

**Configuração:**
- Global: `config/settings.yaml`
- Experimentos: `batches/*.yaml`
- Credenciais: `.env` (NÃO versionar)

**Ver detalhes:** `config/settings.yaml`, `.env.example`

## 🔧 Comandos Principais

**CLI:**
```bash
.venv/bin/python main.py              # Menu interativo
.venv/bin/python main.py run Baseline example.fasta
.venv/bin/python main.py batch batches/experimento.yaml
```

**Web Interface:**
```bash
.venv/bin/python -m uvicorn src.presentation.web.app:app --reload --host 0.0.0.0 --port 8000
```

**Ver tasks:** VS Code Tasks (Ctrl+Shift+P → "Tasks: Run Task")

## 📝 Convenções

**Código:**
- **Inglês:** Variáveis, funções, classes, docstrings
- **Formatação:** black, ruff, mypy
- **Testes:** pytest com cobertura >85%

**Nomenclatura:**
```python
class AlgorithmName:        # PascalCase
def function_name():        # snake_case
CONSTANT_VALUE = 42         # UPPER_SNAKE_CASE
```

**Ver padrões:** `pyproject.toml`, `tests/`

## 🚫 Restrições Críticas

**PROIBIDO:**
- ❌ Imports diretos de plugins na aplicação
- ❌ I/O ou dependências externas no Domain
- ❌ Hardcoded paths, IPs, credenciais
- ❌ Dados mock fixos em código
- ❌ Usar `python` direto (sempre `.venv/bin/python`)

**OBRIGATÓRIO:**
- ✅ Plugins via `@register_algorithm`
- ✅ Configuração via arquivos YAML/ENV
- ✅ Internacionalização (código em inglês)
- ✅ Testes para mudanças
- ✅ Tipagem estática (type hints)

## 🌐 Interface Web

**Estrutura:** `src/presentation/web/`
- `app.py`: FastAPI app principal
- `templates/`: HTML templates
- `static/`: CSS, JS, assets

**APIs:** RESTful com Pydantic models
**Frontend:** HTML5, CSS3, vanilla JS

**Ver implementação:** `src/presentation/web/`

## 🧪 Qualidade

**Testes:**
```bash
.venv/bin/python -m pytest tests/ -v
.venv/bin/python -m pytest --cov=src tests/
```

**Formatação:**
```bash
.venv/bin/python -m black .
.venv/bin/python -m ruff check .
.venv/bin/python -m mypy src/
```

**Ver configuração:** `pyproject.toml`

## 📄 Acadêmico (JOSS/JORS)

**Requisitos:**
- Documentação científica completa
- Reprodutibilidade (seeds determinísticos)
- Métricas padronizadas
- Testes de reprodução
- Metadados (`CITATION.cff`)

**Ver standards:** `CITATION.cff`, `docs/`

---

**Versão:** 0.1.0 | **Data:** Julho 2025 | **Status:** Ativo

> Para detalhes completos, consulte os arquivos específicos mencionados em cada seção.
