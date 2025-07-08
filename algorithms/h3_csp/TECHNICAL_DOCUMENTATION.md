# H³-CSP: Documentação Técnica Detalhada

## Visão Geral

O **H³-CSP** (Hybrid Hierarchical Hamming Search) é um algoritmo híbrido inovador que resolve o Closest String Problem através de uma abordagem de três camadas. Este documento fornece uma análise técnica aprofundada do algoritmo, suas heurísticas e implementação.

## Fundamentação Teórica

### O Problema do Closest String (CSP)

**Definição**: Dado um conjunto de strings S = {s₁, s₂, ..., sₙ} sobre um alfabeto Σ, encontrar uma string c que minimize:
```
d(c, S) = max{d(c, sᵢ) : i = 1, ..., n}
```

onde d(c, sᵢ) é a distância de Hamming entre c e sᵢ.

### Motivação para Abordagem Híbrida

1. **Complexidade Exponencial**: O CSP é NP-hard, exigindo heurísticas eficientes
2. **Localidade Espacial**: Sequências biológicas apresentam padrões locais
3. **Escalabilidade**: Necessidade de lidar com sequências de tamanho variável
4. **Qualidade vs Eficiência**: Equilíbrio entre qualidade da solução e tempo computacional

## Arquitetura do Algoritmo

### Fase 1: B-Splitter (Divisão Hierárquica)

#### Regra √L
```python
B = ⌈√L⌉
block_size = ⌈L/B⌉
```

**Justificativa Matemática**:
- Minimiza o produto B × block_size
- Otimiza o trade-off entre número de blocos e tamanho por bloco
- Reduz complexidade de O(|Σ|^L) para O(√L × |Σ|^(√L))

#### Algoritmo de Divisão
```python
def split_in_blocks(L: int) -> list[Block]:
    B = math.ceil(math.sqrt(L))
    base_size = math.ceil(L / B)
    blocks = []
    cur = 0
    while cur < L:
        blocks.append((cur, min(cur + base_size, L)))
        cur += base_size
    return blocks
```

### Fase 2: Smart-Core (Núcleo Inteligente)

#### Métrica de Dificuldade
Para cada bloco b = (l, r), calcula-se:
```
d_b = max_distance(consensus_b, {s[l:r] : s ∈ S})
```

onde `consensus_b` é a string consenso por maioria para o bloco.

#### Seleção Adaptativa de Técnicas

| Condição | Técnica | Complexidade | Justificativa |
|----------|---------|-------------|---------------|
| d_b ≤ block_small | Busca Exaustiva | O(\|Σ\|^m) | Solução ótima para blocos fáceis |
| block_small < d_b ≤ block_medium | Beam Search Reduzido | O(beam_width/2 × m × \|Σ\|) | Eficiência para blocos médios |
| d_b > block_medium | Beam Search Completo | O(beam_width × m × \|Σ\|) | Exploração ampla para blocos difíceis |

### Fase 3: Global Refine (Refinamento Global)

#### Fusão de Blocos
```python
def _fuse_blocks(self, chosen: list[String]) -> String:
    return "".join(chosen)
```

#### Hill-Climbing Local
```python
def _local_search(candidate: String, strings: Sequence[String]) -> String:
    # Algoritmo guloso que para em ótimo local
    # Explora vizinhança de Hamming distância 1
    # Usa apenas caracteres presentes nas strings originais
```

**Propriedades**:
- **Convergência**: Garantida em número finito de passos
- **Ótimo Local**: Sem garantia de ótimo global
- **Complexidade**: O(L × |Σ|_local × iterações)

## Análise de Complexidade

### Complexidade Temporal

#### Fase 1: B-Splitter
- **Tempo**: O(L)
- **Operação**: Divisão aritmética simples

#### Fase 2: Smart-Core
Para cada bloco b de tamanho m:

1. **Busca Exaustiva**: O(|Σ|^m)
2. **Beam Search**: O(beam_width × m × |Σ|)
3. **Consenso**: O(n × m)

**Total**: O(√L × max(|Σ|^(√L), beam_width × √L × |Σ|))

#### Fase 3: Global Refine
- **Fusão**: O(L)
- **Hill-Climbing**: O(local_iters × L × |Σ|_local)

**Complexidade Total**: O(√L × |Σ|^(√L) + L × |Σ|_local)

### Complexidade Espacial

- **Armazenamento de strings**: O(n × L)
- **Candidatos por bloco**: O(k × √L × √L) = O(k × L)
- **Estruturas auxiliares**: O(L)

**Total**: O(n × L + k × L) = O((n + k) × L)

## Heurísticas Avançadas

### 1. Decomposição Hierárquica

**Princípio**: Divide et impera aplicado ao CSP
- Reduz dimensionalidade do problema
- Permite paralelização conceitual
- Preserva localidade espacial

### 2. Seleção Adaptativa

**Princípio**: Economia de recursos computacionais
- Técnicas caras apenas quando necessário
- Adaptação baseada em características dos dados
- Balanceamento dinâmico qualidade/eficiência

### 3. Diversidade de Candidatos

**Princípio**: Exploração múltipla do espaço de soluções
- k candidatos por bloco aumentam robustez
- Combinação de estratégias (exaustiva, beam, consenso)
- Fallback inteligente para casos extremos

### 4. Refinamento Global

**Princípio**: Correção de inconsistências entre blocos
- Hill-climbing corrige junções subótimas
- Busca local melhora solução global
- Vizinhança restrita acelera convergência

## Implementação Técnica

### Estruturas de Dados

```python
String = str                    # Alias para strings
Block = tuple[int, int]        # (início, fim) - fim exclusivo
```

### Padrões de Design

1. **Strategy Pattern**: Seleção de técnicas por bloco
2. **Template Method**: Estrutura geral do algoritmo
3. **Observer Pattern**: Callbacks de progresso
4. **Factory Pattern**: Criação de candidatos

### Tratamento de Erros

```python
try:
    # Execução do algoritmo
    return center, best_val
except Exception as e:
    logger.error("Erro no H3CSP: %s", str(e))
    if self.progress_callback:
        self.progress_callback(f"Erro: {str(e)}")
    raise e  # Re-propaga para wrapper
```

## Otimizações e Melhorias

### Otimizações Implementadas

1. **Lazy Evaluation**: Consenso calculado apenas quando necessário
2. **Early Stopping**: Hill-climbing para sem melhorias
3. **Timeout Protection**: Controle de tempo máximo
4. **Memory Efficiency**: Reutilização de estruturas

### Melhorias Futuras

1. **Paralelização**: Processamento simultâneo de blocos
2. **Algoritmos Genéticos**: Evolução de população de candidatos
3. **Machine Learning**: Aprendizado de padrões para seleção de técnicas
4. **Aproximação Estocástica**: Técnicas probabilísticas para blocos grandes

## Casos de Teste e Validação

### Datasets de Teste

1. **Sintéticos**: Controle total sobre características
2. **Biológicos**: DNA, RNA, proteínas reais
3. **Adversários**: Casos extremos e corner cases

### Métricas de Avaliação

1. **Qualidade**: Distância da solução ótima
2. **Eficiência**: Tempo de execução
3. **Escalabilidade**: Comportamento com L crescente
4. **Robustez**: Estabilidade em cenários diversos

## Comparação com Outros Algoritmos

### Vantagens do H³-CSP

1. **Escalabilidade**: Sublinear em L
2. **Adaptabilidade**: Técnicas baseadas em características
3. **Qualidade**: Bom equilíbrio com eficiência
4. **Flexibilidade**: Múltiplos parâmetros ajustáveis

### Desvantagens

1. **Complexidade**: Mais complexo que algoritmos simples
2. **Parâmetros**: Requer ajuste fino
3. **Overhead**: Ineficiente para instâncias muito pequenas
4. **Não-ótimo**: Sem garantia de solução ótima

## Conclusão

O H³-CSP representa uma abordagem inovadora para o Closest String Problem, combinando teoria algorítmica sólida com heurísticas práticas eficientes. Sua arquitetura híbrida permite adaptação dinâmica às características dos dados, oferecendo um excelente equilíbrio entre qualidade da solução e eficiência computacional.

A implementação técnica segue padrões de design estabelecidos, garantindo manutenibilidade e extensibilidade. As otimizações implementadas e as melhorias futuras sugeridas demonstram o potencial de evolução contínua do algoritmo.

---

*Para mais informações técnicas, consulte o código-fonte e os docstrings detalhados.*
