# Reestruturação da Aplicação CSP - Sumário

## O que foi implementado

### 1. Interface Padronizada de Algoritmos
- **Classe base abstrata** `Algorithm` em `algorithms/base.py`
- **Sistema de registro automático** com decorator `@register_algorithm`
- **Registry global** que mantém todos os algoritmos disponíveis

### 2. Estrutura Modular por Algoritmo
Cada algoritmo agora possui sua própria pasta com:
```
algorithms/
├── {nome_algoritmo}/
│   ├── __init__.py          # Expõe o algoritmo
│   ├── algorithm.py         # Wrapper com interface padronizada
│   ├── config.py            # Configurações específicas
│   └── implementation.py    # Implementação original
```

### 3. Algoritmos Migrados
- ✅ **Baseline** - `algorithms/baseline/`
- ✅ **BLF-GA** - `algorithms/blf_ga/`
- ✅ **CSC** - `algorithms/csc/`
- ✅ **DP-CSP** - `algorithms/dp_csp/`
- ✅ **H³-CSP** - `algorithms/h3_csp/`

### 4. Separação de Configurações
- **Configurações da aplicação** permanecem em `utils/config.py`
- **Configurações específicas** de cada algoritmo estão em suas respectivas pastas
- Removidas dependências cruzadas entre algoritmos e config global

### 5. Main.py Dinâmico
- **Menu automático** lista algoritmos do registry
- **Execução genérica** funciona com qualquer algoritmo registrado
- **Zero modificações necessárias** para adicionar novos algoritmos

## Como Adicionar Novos Algoritmos

### Passos Simples:
1. Criar pasta `algorithms/meu_algoritmo/`
2. Adicionar os 4 arquivos (config.py, implementation.py, algorithm.py, __init__.py)
3. Usar `@register_algorithm` no wrapper
4. Pronto! Aparece automaticamente no menu

### Exemplo Mínimo:
```python
# algorithms/meu_algoritmo/algorithm.py
from algorithms.base import Algorithm, register_algorithm

@register_algorithm
class MeuAlgoritmo(Algorithm):
    name = "Meu Algoritmo"
    default_params = {}
    
    def __init__(self, strings, alphabet, **params):
        self.strings = strings
        self.alphabet = alphabet
    
    def run(self) -> tuple[str, int]:
        # Sua implementação aqui
        return center_string, max_distance
```

## Benefícios Alcançados

### ✅ **Modularidade Total**
- Cada algoritmo é independente
- Configurações isoladas
- Implementações preservadas

### ✅ **Facilidade de Expansão**
- Novos algoritmos não requerem modificação do main
- Interface padronizada garantida
- Auto-registro automático

### ✅ **Manutenibilidade**
- Código organizado por responsabilidade
- Configurações centralizadas por algoritmo
- Dependências claras

### ✅ **Compatibilidade**
- Todas as implementações existentes preservadas
- Funcionalidade inalterada
- Performance mantida

## Arquivos Modificados

### Criados:
- `algorithms/base.py` - Interface e registry
- `algorithms/*/algorithm.py` - Wrappers padronizados
- `algorithms/*/config.py` - Configurações isoladas
- `algorithms/README.md` - Documentação

### Modificados:
- `main.py` - Lógica dinâmica de execução
- `utils/config.py` - Removidas configs de algoritmos
- `algorithms/__init__.py` - Auto-import de subpacotes

### Preservados:
- `algorithms/*/implementation.py` - Código original dos algoritmos
- Toda funcionalidade existente
- Estrutura de datasets inalterada

## Teste de Validação

✅ Todos os 5 algoritmos funcionando corretamente
✅ Registry automático operacional  
✅ Novo algoritmo pode ser adicionado em minutos
✅ Main.py completamente genérico

A aplicação agora possui uma arquitetura extensível e modular que facilita muito a adição de novos algoritmos!
