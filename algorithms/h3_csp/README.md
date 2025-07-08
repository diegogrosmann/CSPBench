# H³-CSP: Hybrid Hierarchical Hamming Search

O **H³-CSP** (Hybrid Hierarchical Hamming Search) é um algoritmo híbrido avançado que combina decomposição hierárquica, técnicas especializadas por bloco e refinamento global para resolver eficientemente o **Closest String Problem** (CSP).

## 🏗️ Arquitetura do Algoritmo

O H³-CSP opera em três fases principais:

### 1. **B-Splitter** (Divisão Hierárquica)
- Divide as strings em **~√L blocos contíguos** usando a "regra √L"
- Cada bloco tem tamanho aproximadamente uniforme
- Permite processamento paralelo conceitual dos blocos

### 2. **Smart-Core** (Núcleo Inteligente)
Seleciona a técnica ótima para cada bloco baseada na **dificuldade do bloco** (d_b):

| Dificuldade | Técnica | Descrição |
|-------------|---------|-----------|
| d_b ≤ 2 | **Busca Exaustiva** | Solução ótima para blocos pequenos |
| 2 < d_b ≤ 4 | **Beam Search Reduzido** | Beam width reduzido para eficiência |
| d_b > 4 | **Beam Search Completo** | Beam width completo para blocos difíceis |

### 3. **Global Refine** (Refinamento Global)
- **Fusão**: Combina os melhores candidatos de cada bloco
- **Hill-Climbing**: Refinamento iterativo por busca local
- Convergência garantida para ótimo local

## 🧠 Heurísticas Avançadas

### Decomposição Adaptativa
- **Regra √L**: Número ótimo de blocos baseado no comprimento
- **Blocos Contíguos**: Preserva localidade espacial dos dados
- **Tamanho Uniforme**: Distribui carga computacional

### Seleção Inteligente de Técnicas
- **Análise de Dificuldade**: Calcula d_b = max_distance(consenso, blocos)
- **Adaptação Dinâmica**: Técnica escolhida por bloco individual
- **Eficiência Otimizada**: Técnicas caras apenas quando necessário

### Geração de Candidatos Diversos
- **Múltiplos Candidatos**: k candidatos por bloco (padrão: 5)
- **Estratégias Híbridas**: Combina exaustiva, beam search e consenso
- **Fallback Inteligente**: Usa dataset original para blocos muito grandes

## ⚙️ Parâmetros Principais

### Divisão de Blocos
```python
"auto_blocks": True,        # Usa divisão automática √L
"min_block_size": 2,       # Tamanho mínimo de bloco
"max_blocks": None,        # Máximo de blocos (None = automático)
```

### Limiares de Dificuldade
```python
"block_small": 2,          # Limite para busca exaustiva
"block_medium": 4,         # Limite para beam search reduzido
"block_large": 8,          # Limite para beam search completo
```

### Técnicas de Busca
```python
"exhaustive_limit": 10000, # Limite |Σ|^m para busca exaustiva
"beam_width": 32,          # Largura do beam search
"k_candidates": 5,         # Candidatos por bloco
```

### Refinamento e Controle
```python
"local_iters": 3,          # Iterações de hill-climbing
"max_time": 300,           # Timeout em segundos
"seed": None,              # Semente para reprodutibilidade
```

## 🎯 Casos de Uso Ideais

### ✅ Quando Usar H³-CSP
- **Instâncias médias**: L entre 50-500 caracteres
- **Dados estruturados**: Presença de padrões locais
- **Qualidade vs Eficiência**: Necessidade de equilíbrio
- **Recursos limitados**: Quando algoritmos exaustivos são inviáveis

### ❌ Limitações
- **Overhead para instâncias pequenas** (L < 20)
- **Múltiplos parâmetros** requerem ajuste fino
- **Complexidade adicional** comparado a algoritmos simples
- **Sem garantia de ótimo global** (apenas ótimo local)

## 📊 Complexidade Computacional

- **Tempo**: O(√L × |Σ|^(√L) + L × |Σ|_local + k × refinamento)
- **Espaço**: O(L + k × √L)
- **Escalabilidade**: Sublinear em L para a maioria dos casos

## 💻 Exemplo de Uso

```python
from algorithms.h3_csp.algorithm import H3CSPAlgorithm

# Instanciar o algoritmo
strings = ["ACGTTAGC", "AGGTTAGC", "ACTTTAGC"]
alphabet = "ACGT"

alg = H3CSPAlgorithm(
    strings, 
    alphabet,
    beam_width=16,          # Beam search mais focado
    k_candidates=3,         # Menos candidatos por bloco
    local_iters=5,          # Mais refinamento
    max_time=120            # Timeout de 2 minutos
)

# Executar
center, distance, metadata = alg.run()

print(f"Solução encontrada: {center}")
print(f"Distância máxima: {distance}")
print(f"Iterações: {metadata['iteracoes']}")
print(f"Parâmetros: {metadata['parametros_usados']}")
```

## 🔄 Callback de Progresso

```python
def progress_handler(message):
    print(f"[H³-CSP] {message}")

alg.set_progress_callback(progress_handler)
center, distance, metadata = alg.run()
```

## 📚 Documentação Técnica

### Arquivos do Módulo
- `algorithm.py`: Wrapper de integração ao framework
- `implementation.py`: Implementação completa do algoritmo
- `config.py`: Configurações e parâmetros padrão
- `README.md`: Documentação do usuário

### Estilo de Documentação
- **Docstrings**: Formato Google Style
- **Type Hints**: Tipagem completa
- **Logging**: Registro detalhado de execução
- **Testes**: Integração com framework de testes

## 🔗 Integração com Framework

O H³-CSP está totalmente integrado ao framework CSP-BLFGA através do decorador `@register_algorithm`, permitindo:

- **Execução automatizada** em batches
- **Comparação** com outros algoritmos
- **Análise de performance** e métricas
- **Visualização** de resultados

## 🎨 Visualização de Resultados

O algoritmo gera metadados ricos que podem ser usados para:
- **Análise de convergência** do refinamento local
- **Comparação de eficácia** por bloco
- **Profiling de performance** por fase
- **Visualização** da divisão hierárquica

---

*H³-CSP: Quando você precisa de qualidade e eficiência em harmonia.*
