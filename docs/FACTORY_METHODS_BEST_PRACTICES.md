# Factory Methods nos Persistence Wrappers

## Visão Geral

Os factory methods foram implementados nos persistence wrappers para facilitar a navegação e criação de instâncias entre os diferentes níveis hierárquicos. Esta é considerada uma **excelente prática** por várias razões.

## Factory Methods Implementados

### WorkScopedPersistence

```python
def for_combination(self, combination_id: int) -> CombinationScopedPersistence:
    """Cria CombinationScopedPersistence para uma combination específica deste work."""

def for_execution(self, unit_id: str) -> ExecutionScopedPersistence:
    """Cria ExecutionScopedPersistence para uma execution específica deste work."""

def get_all_combinations(self) -> list[CombinationScopedPersistence]:
    """Retorna todas as combinations deste work como instâncias."""

def get_all_executions(self) -> list[ExecutionScopedPersistence]:
    """Retorna todas as executions deste work como instâncias."""
```

### CombinationScopedPersistence

```python
def for_execution(self, unit_id: str) -> ExecutionScopedPersistence:
    """Cria ExecutionScopedPersistence para uma execution específica desta combination."""

def get_all_executions(self) -> list[ExecutionScopedPersistence]:
    """Retorna todas as executions desta combination como instâncias."""
```

## Por que é uma Boa Prática?

### 1. **Fluent Interface** 🔗
```python
# Navegação natural e expressiva
work = WorkScopedPersistence("work_123")
combination = work.for_combination(456)
execution = combination.for_execution("unit_789")
```

### 2. **Encapsulamento** 🛡️
```python
# Esconde a complexidade de criação
# Em vez de:
combination = CombinationScopedPersistence(456, work.store)

# Use:
combination = work.for_combination(456)
```

### 3. **Validação Automática** ✅
```python
# Valida automaticamente se a combination pertence ao work
try:
    combo = work.for_combination(999)  # ID inválido
except ValueError as e:
    print(f"Erro: {e}")  # Combination 999 not found in work work_123
```

### 4. **Type Safety** 🎯
```python
# IDEs podem inferir tipos corretamente
work: WorkScopedPersistence = get_work()
combo: CombinationScopedPersistence = work.for_combination(456)  # Tipo garantido
```

### 5. **Composição sobre Herança** 🧩
```python
# Permite composição flexível
def process_work_data(work: WorkScopedPersistence):
    for combo in work.get_all_combinations():
        combo.log_warning("Processing...")
        for exec in combo.get_all_executions():
            exec.update_execution_status("running")
```

## Padrões de Uso Recomendados

### Navegação Hierárquica
```python
def process_work_hierarchy(work_id: str):
    work = WorkScopedPersistence(work_id)
    
    for combo_data in work.get_combinations():
        combo = work.for_combination(combo_data["id"])
        combo.log_warning("Starting combination processing")
        
        for exec_data in combo.get_executions():
            exec = combo.for_execution(exec_data["unit_id"])
            exec.update_execution_status("running")
```

### Processamento em Lote
```python
def process_all_executions(work_id: str):
    work = WorkScopedPersistence(work_id)
    
    # Forma concisa de obter todas as executions
    executions = work.get_all_executions()
    
    for execution in executions:
        execution.add_progress(0.0, "Starting...")
        # ... processamento
        execution.add_progress(1.0, "Completed!")
```

### Pipeline de Dados
```python
def create_execution_pipeline(work_id: str):
    work = WorkScopedPersistence(work_id)
    
    return (
        work.for_combination(combo_id)
        .for_execution(unit_id)  # Chain de factory methods
        .add_progress(0.5, "Pipeline created")
    )
```

## Benefícios Técnicos

### Performance
- ✅ Reutiliza a mesma instância de `WorkPersistence`
- ✅ Evita re-conexões desnecessárias ao banco
- ✅ Cache de dados compartilhado

### Manutenibilidade  
- ✅ Mudanças centralizadas nos factory methods
- ✅ Evolução da API sem quebrar código existente
- ✅ Facilita refatoração

### Testabilidade
- ✅ Fácil criação de mocks/stubs
- ✅ Testes isolados por nível
- ✅ Validação de relacionamentos

## Comparação: Antes vs Depois

### Antes (Manual)
```python
from src.infrastructure.persistence.work_state.core import WorkPersistence
from src.infrastructure.persistence.work_state.wrappers.work_scoped import WorkScopedPersistence
from src.infrastructure.persistence.work_state.wrappers.combination_scoped import CombinationScopedPersistence
from src.infrastructure.persistence.work_state.wrappers.execution_scoped import ExecutionScopedPersistence

# Criação manual - verbosa e propensa a erros
store = WorkPersistence()
work = WorkScopedPersistence("work_123", store)
combination = CombinationScopedPersistence(456, store)  # Sem validação!
execution = ExecutionScopedPersistence("unit_789", "work_123", store)  # IDs podem não corresponder!
```

### Depois (Factory Methods)
```python
from src.infrastructure.persistence.work_state.wrappers.work_scoped import WorkScopedPersistence

# Criação elegante com validação automática
work = WorkScopedPersistence("work_123")
combination = work.for_combination(456)  # Valida se pertence ao work
execution = combination.for_execution("unit_789")  # Valida se pertence à combination
```

## Alinhamento com Padrões de Design

### Factory Method Pattern
- ✅ Encapsula lógica de criação de objetos
- ✅ Permite extensibilidade futura
- ✅ Segue princípios SOLID

### Fluent Interface Pattern
- ✅ APIs expressivas e legíveis
- ✅ Método chaining natural
- ✅ Reduz verbosidade

### Builder Pattern Elements
- ✅ Construção incremental de objetos complexos
- ✅ Validação durante a construção
- ✅ Interface clara e intuitiva

## Conclusão

A implementação de factory methods nos persistence wrappers é uma **excelente prática** que:

1. **Melhora a experiência do desenvolvedor** com APIs mais intuitivas
2. **Reduz bugs** através de validação automática
3. **Facilita manutenção** centralizando lógica de criação
4. **Aumenta produtividade** com código mais conciso
5. **Segue padrões estabelecidos** da engenharia de software

Esta implementação estabelece uma base sólida para expansão futura do sistema de persistence, mantendo alta qualidade de código e usabilidade.
