# Guia de Contribuição - CSP-BLFGA

## 🎯 Visão Geral

Obrigado por considerar contribuir para o CSP-BLFGA! Este projeto é uma plataforma experimental para o Closest String Problem que visa ser robusta, extensível e útil para a comunidade científica.

## 🚀 Como Contribuir

### Tipos de Contribuições

1. **🐛 Correção de Bugs**: Encontrou um problema? Ajude-nos a corrigi-lo!
2. **✨ Novas Funcionalidades**: Ideias para melhorar o projeto
3. **📊 Novos Algoritmos**: Implementações de algoritmos CSP
4. **📚 Documentação**: Melhorias na documentação
5. **🧪 Testes**: Ampliação da cobertura de testes
6. **🔧 Otimizações**: Melhorias de performance

### Processo de Contribuição

1. **Fork** o repositório
2. **Crie uma branch** para sua feature: `git checkout -b feature/nome-da-feature`
3. **Faça as alterações** seguindo os padrões do projeto
4. **Adicione testes** para suas mudanças
5. **Execute os testes** para garantir que tudo funciona
6. **Commit** suas mudanças: `git commit -m 'Adiciona nova feature'`
7. **Push** para sua branch: `git push origin feature/nome-da-feature`
8. **Abra um Pull Request** com descrição detalhada

## 📋 Preparando o Ambiente

### Configuração Inicial

```bash
# 1. Clone o repositório
git clone https://github.com/seu-usuario/csp-blfga.git
cd csp-blfga

# 2. Crie um ambiente virtual
python -m venv .venv
source .venv/bin/activate  # Linux/macOS
# ou
.venv\Scripts\activate     # Windows

# 3. Instale dependências de desenvolvimento
pip install -r requirements.txt
pip install -e .[dev]

# 4. Configure pre-commit hooks
pre-commit install
```

### Verificação do Ambiente

```bash
# Executar testes
python -m pytest tests/ -v

# Verificar formatação
python -m black --check .

# Verificar linting
python -m ruff check .

# Verificar tipos
python -m mypy src/

# Executar aplicação
python main.py --help
```

## 🧪 Executando Testes

### Testes Unitários

```bash
# Executar todos os testes
pytest tests/

# Executar testes específicos
pytest tests/unit/test_scheduler.py

# Executar com cobertura
pytest tests/ --cov=src --cov-report=html --cov-report=term-missing

# Executar testes paralelos
pytest tests/ -n auto
```

### Testes de Integração

```bash
# Testes de integração completos
pytest tests/integration/

# Teste de execução completa
pytest tests/integration/test_full_flow.py -v
```

### Testes de Performance

```bash
# Benchmarks
python -m pytest tests/benchmarks/ --benchmark-only

# Profiling
python -m pytest tests/performance/ --profile
```

## 📝 Padrões de Código

### Estilo de Código

O projeto segue **PEP 8** com algumas extensões:

```python
# ✅ Bom
def calculate_distance(center: str, sequences: List[str]) -> int:
    """
    Calcula distância máxima entre centro e sequências.
    
    Args:
        center: String centro
        sequences: Lista de sequências
        
    Returns:
        Distância máxima
    """
    return max(hamming_distance(center, seq) for seq in sequences)

# ❌ Ruim
def calc_dist(c, s):
    return max([hamming_distance(c, seq) for seq in s])
```

### Naming Conventions

- **Classes**: `PascalCase` (ex: `ExecutionScheduler`)
- **Funções/Métodos**: `snake_case` (ex: `submit_task`)
- **Constantes**: `UPPER_SNAKE_CASE` (ex: `MAX_WORKERS`)
- **Variáveis**: `snake_case` (ex: `task_id`)
- **Arquivos**: `snake_case` (ex: `execution_scheduler.py`)

### Type Hints

Sempre use type hints:

```python
from typing import List, Dict, Optional, Tuple, Any, Union

def process_algorithm_results(
    results: List[Dict[str, Any]], 
    algorithm_name: str,
    timeout: Optional[float] = None
) -> Tuple[float, Dict[str, Any]]:
    """Processa resultados de algoritmo."""
    pass
```

### Documentação

Use docstrings no formato Google:

```python
def complex_function(param1: str, param2: int, param3: Optional[bool] = None) -> Dict[str, Any]:
    """
    Função complexa com múltiplos parâmetros.
    
    Args:
        param1: Descrição do primeiro parâmetro
        param2: Descrição do segundo parâmetro
        param3: Parâmetro opcional com valor padrão
        
    Returns:
        Dicionário com resultados processados
        
    Raises:
        ValueError: Se param1 estiver vazio
        TypeError: Se param2 não for inteiro
        
    Examples:
        >>> result = complex_function("test", 42)
        >>> assert "processed" in result
    """
    if not param1:
        raise ValueError("param1 não pode estar vazio")
    
    # Implementação...
    return {"processed": True}
```

## 🧬 Adicionando Novos Algoritmos

### Estrutura Básica

```python
# algorithms/meu_algoritmo/algorithm.py
from algorithms.base import CSPAlgorithm, register_algorithm

@register_algorithm
class MeuAlgoritmo(CSPAlgorithm):
    """
    Descrição do algoritmo.
    
    Referências:
        - Paper original: [link]
        - Implementação baseada em: [fonte]
    """
    
    name = "MeuAlgoritmo"
    default_params = {
        'param1': 10,
        'param2': 0.5,
        'param3': True
    }
    is_deterministic = False  # ou True se for determinístico
    supports_internal_parallel = False  # ou True se suportar paralelismo
    
    def __init__(self, strings: List[str], alphabet: str, **params):
        super().__init__(strings, alphabet, **params)
        
        # Validar parâmetros específicos
        self._validate_params()
        
        # Inicializar estruturas do algoritmo
        self._initialize_algorithm()
    
    def _validate_params(self) -> None:
        """Valida parâmetros específicos do algoritmo."""
        if self.params['param1'] <= 0:
            raise ValueError("param1 deve ser positivo")
    
    def _initialize_algorithm(self) -> None:
        """Inicializa estruturas específicas do algoritmo."""
        self.state = {}
        self.iteration_count = 0
    
    def run(self) -> Tuple[str, int, Dict[str, Any]]:
        """
        Executa o algoritmo.
        
        Returns:
            Tuple contendo (centro, distância, metadados)
        """
        self._report_progress("Iniciando algoritmo...")
        
        try:
            # Algoritmo principal
            center = self._find_center()
            distance = self._calculate_distance(center)
            
            # Coletar metadados
            metadata = self._collect_metadata()
            
            return center, distance, metadata
            
        except Exception as e:
            self._report_warning(f"Erro durante execução: {e}")
            raise
    
    def _find_center(self) -> str:
        """Implementa a lógica principal do algoritmo."""
        # Sua implementação aqui
        pass
    
    def _calculate_distance(self, center: str) -> int:
        """Calcula distância do centro."""
        from src.utils.distance import hamming_distance
        return max(hamming_distance(center, seq) for seq in self.strings)
    
    def _collect_metadata(self) -> Dict[str, Any]:
        """Coleta metadados da execução."""
        return {
            'iterations': self.iteration_count,
            'algorithm_specific_metric': self.state.get('metric', 0),
            'convergence_info': self.state.get('convergence', {})
        }
```

### Estrutura de Diretórios

```
algorithms/meu_algoritmo/
├── __init__.py
├── algorithm.py          # Classe principal
├── implementation.py     # Implementação do algoritmo
├── config.py            # Configurações padrão
├── README.md            # Documentação do algoritmo
└── tests/
    ├── __init__.py
    ├── test_algorithm.py
    └── test_implementation.py
```

### Testando Novos Algoritmos

```python
# algorithms/meu_algoritmo/tests/test_algorithm.py
import pytest
from algorithms.meu_algoritmo.algorithm import MeuAlgoritmo

class TestMeuAlgoritmo:
    def setup_method(self):
        self.sequences = [
            "ACGTACGT",
            "ACGTACGA",
            "ACGTACGC"
        ]
        self.alphabet = "ACGT"
    
    def test_initialization(self):
        """Testa inicialização do algoritmo."""
        algo = MeuAlgoritmo(self.sequences, self.alphabet)
        assert algo.strings == self.sequences
        assert algo.alphabet == self.alphabet
    
    def test_run_basic(self):
        """Testa execução básica."""
        algo = MeuAlgoritmo(self.sequences, self.alphabet)
        center, distance, metadata = algo.run()
        
        assert isinstance(center, str)
        assert isinstance(distance, int)
        assert isinstance(metadata, dict)
        assert distance >= 0
    
    def test_progress_callback(self):
        """Testa callback de progresso."""
        progress_messages = []
        
        def progress_callback(message):
            progress_messages.append(message)
        
        algo = MeuAlgoritmo(self.sequences, self.alphabet)
        algo.set_progress_callback(progress_callback)
        algo.run()
        
        assert len(progress_messages) > 0
    
    def test_invalid_params(self):
        """Testa validação de parâmetros inválidos."""
        with pytest.raises(ValueError):
            MeuAlgoritmo(self.sequences, self.alphabet, param1=-1)
```

## 📊 Adicionando Novos Datasets

### Estrutura de Dataset Loader

```python
# src/datasets/meu_dataset.py
from typing import List, Dict, Any, Tuple
from src.datasets.base import DatasetLoader

class MeuDatasetLoader(DatasetLoader):
    """
    Carregador para meu tipo de dataset.
    
    Suporta:
        - Formato específico
        - Validação customizada
        - Metadados específicos
    """
    
    def __init__(self, **params):
        super().__init__(**params)
        self.required_params = ['source_path', 'format_type']
        self._validate_params()
    
    def load(self) -> Tuple[List[str], Dict[str, Any]]:
        """
        Carrega dataset do arquivo.
        
        Returns:
            Tuple com (sequências, metadados)
        """
        sequences = self._load_sequences()
        metadata = self._extract_metadata()
        
        # Validar sequências
        self._validate_sequences(sequences)
        
        return sequences, metadata
    
    def _load_sequences(self) -> List[str]:
        """Carrega sequências do arquivo."""
        # Implementação específica
        pass
    
    def _extract_metadata(self) -> Dict[str, Any]:
        """Extrai metadados do dataset."""
        return {
            'source': self.params['source_path'],
            'format': self.params['format_type'],
            'loader': 'MeuDatasetLoader'
        }
    
    def _validate_sequences(self, sequences: List[str]) -> None:
        """Valida sequências carregadas."""
        if not sequences:
            raise ValueError("Nenhuma sequência encontrada")
        
        # Validação específica do formato
        # ...
```

## 🎨 Melhorando a Interface

### Adicionando Novo Menu

```python
# src/ui/cli/menu.py
def meu_novo_menu():
    """Menu para nova funcionalidade."""
    questions = [
        {
            'type': 'list',
            'name': 'option',
            'message': 'Selecione uma opção:',
            'choices': [
                {'name': 'Opção 1', 'value': 'option1'},
                {'name': 'Opção 2', 'value': 'option2'},
                {'name': 'Voltar', 'value': 'back'}
            ]
        }
    ]
    
    answers = prompt(questions)
    
    if answers['option'] == 'option1':
        # Implementar opção 1
        pass
    elif answers['option'] == 'option2':
        # Implementar opção 2
        pass
    else:
        return  # Voltar
```

### Melhorando Interface Curses

```python
# src/ui/curses_integration.py
class NovaTelaCurses:
    """Nova tela para funcionalidade específica."""
    
    def __init__(self, stdscr):
        self.stdscr = stdscr
        self.setup_colors()
    
    def setup_colors(self):
        """Configura cores para a tela."""
        curses.start_color()
        curses.init_pair(1, curses.COLOR_GREEN, curses.COLOR_BLACK)
        curses.init_pair(2, curses.COLOR_RED, curses.COLOR_BLACK)
    
    def render(self, data):
        """Renderiza a tela com dados."""
        self.stdscr.clear()
        
        # Renderizar conteúdo
        self.stdscr.addstr(0, 0, "Título da Tela", curses.color_pair(1))
        
        # Atualizar tela
        self.stdscr.refresh()
```

## 📈 Adicionando Métricas

### Novas Métricas de Performance

```python
# src/utils/metrics.py
class NovaMetrica:
    """Nova métrica para avaliar algoritmos."""
    
    def __init__(self, name: str):
        self.name = name
        self.values = []
    
    def record(self, value: float):
        """Registra novo valor."""
        self.values.append(value)
    
    def calculate(self) -> Dict[str, float]:
        """Calcula estatísticas da métrica."""
        if not self.values:
            return {}
        
        return {
            'mean': sum(self.values) / len(self.values),
            'min': min(self.values),
            'max': max(self.values),
            'std': self._calculate_std(),
            'count': len(self.values)
        }
    
    def _calculate_std(self) -> float:
        """Calcula desvio padrão."""
        if len(self.values) < 2:
            return 0.0
        
        mean = sum(self.values) / len(self.values)
        variance = sum((x - mean) ** 2 for x in self.values) / len(self.values)
        return variance ** 0.5
```

## 🔧 Debugging e Profiling

### Ferramentas de Debug

```python
# src/utils/debug.py
import functools
import time
from typing import Callable, Any

def debug_execution(func: Callable) -> Callable:
    """Decorador para debug de execução."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        print(f"🐛 Executando {func.__name__}")
        print(f"   Args: {args}")
        print(f"   Kwargs: {kwargs}")
        
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        
        print(f"   Tempo: {end_time - start_time:.4f}s")
        print(f"   Resultado: {type(result).__name__}")
        
        return result
    return wrapper

def profile_memory(func: Callable) -> Callable:
    """Decorador para profiling de memória."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        import tracemalloc
        tracemalloc.start()
        
        result = func(*args, **kwargs)
        
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        
        print(f"📊 {func.__name__} - Memória: {current / 1024 / 1024:.2f} MB (pico: {peak / 1024 / 1024:.2f} MB)")
        
        return result
    return wrapper
```

## 📋 Checklist de Pull Request

Antes de submeter um PR, verifique:

### ✅ Código
- [ ] Código segue padrões do projeto (PEP 8)
- [ ] Type hints adicionados
- [ ] Documentação (docstrings) atualizada
- [ ] Comentários explicativos onde necessário
- [ ] Nomes de variáveis/funções são descritivos

### ✅ Testes
- [ ] Testes unitários adicionados
- [ ] Testes de integração (se aplicável)
- [ ] Todos os testes passam
- [ ] Cobertura de testes mantida/melhorada

### ✅ Documentação
- [ ] README atualizado (se necessário)
- [ ] Documentação técnica atualizada
- [ ] Exemplos de uso fornecidos
- [ ] Changelog atualizado

### ✅ Funcionalidade
- [ ] Feature funciona como esperado
- [ ] Não quebra funcionalidades existentes
- [ ] Performance é aceitável
- [ ] Tratamento de erros implementado

### ✅ Git
- [ ] Commit messages são descritivos
- [ ] Branch tem nome apropriado
- [ ] Sem arquivos desnecessários no commit
- [ ] Conflitos resolvidos

## 🏷️ Convenções de Commit

Use conventional commits:

```
feat: adiciona novo algoritmo BFS-CSP
fix: corrige memory leak no scheduler
docs: atualiza documentação do BLF-GA
test: adiciona testes para dataset loader
refactor: melhora estrutura do executor
perf: otimiza cálculo de distância
style: formata código com black
chore: atualiza dependências
```

## 📞 Comunicação

### Canais de Comunicação

- **Issues**: Para bugs e feature requests
- **Discussions**: Para discussões gerais
- **Pull Requests**: Para contribuições de código
- **Email**: Para questões privadas

### Relatando Bugs

Use o template:

```markdown
**Descrição do Bug**
Descrição clara do problema.

**Passos para Reproduzir**
1. Vá para '...'
2. Clique em '...'
3. Execute '...'
4. Veja o erro

**Comportamento Esperado**
O que deveria acontecer.

**Comportamento Atual**
O que acontece na realidade.

**Ambiente**
- OS: [Ubuntu 20.04]
- Python: [3.11.0]
- Versão: [1.0.0]

**Informações Adicionais**
Logs, screenshots, etc.
```

### Solicitando Features

Use o template:

```markdown
**Problema a Resolver**
Descreva o problema que a feature resolve.

**Solução Proposta**
Descreva a solução que você gostaria.

**Alternativas Consideradas**
Outras soluções que você considerou.

**Contexto Adicional**
Informações extras relevantes.
```

## 🎓 Recursos de Aprendizado

### Documentação Interna
- `README.md`: Visão geral do projeto
- `TECHNICAL_DOCUMENTATION.md`: Documentação técnica
- `src/README.md`: Arquitetura do código
- Docstrings: Documentação inline

### Recursos Externos
- [PEP 8](https://pep8.org/): Style guide para Python
- [Type Hints](https://docs.python.org/3/library/typing.html): Documentação oficial
- [Pytest](https://docs.pytest.org/): Framework de testes
- [Black](https://black.readthedocs.io/): Formatador de código

### Algoritmos CSP
- Livros sobre bioinformática
- Papers sobre Closest String Problem
- Implementações de referência

---

Obrigado por contribuir para o CSP-BLFGA! Sua participação ajuda a tornar a plataforma melhor para toda a comunidade científica. 🚀
