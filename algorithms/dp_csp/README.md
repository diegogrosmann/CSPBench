# DP-CSP (Dynamic Programming Closest String Problem)

O algoritmo **DP-CSP** é uma solução **exata** para o Closest String Problem baseada em **programação dinâmica**. Diferentemente de heurísticas que buscam soluções aproximadas, o DP-CSP **garante encontrar a solução ótima** - a string central com o menor raio possível - usando uma busca sistemática sobre todos os possíveis valores de distância.

## 📋 Índice

- [Estratégia Algorítmica](#estratégia-algorítmica)
- [Funcionamento Detalhado](#funcionamento-detalhado)
- [Parâmetros e Configuração](#parâmetros-e-configuração)
- [Casos de Uso](#casos-de-uso)
- [Análise Algorítmica](#análise-algorítmica)
- [Exemplos de Uso](#exemplos-de-uso)
- [Limitações](#limitações)
- [Integração com CSPBench](#integração-com-cspbench)

## 🎯 Estratégia Algorítmica

### Abordagem Principal
O DP-CSP utiliza uma estratégia de **busca binária incremental** combinada com **programação dinâmica decisória**:

1. **Busca Incremental**: Testa valores crescentes de d (0, 1, 2, ...) até encontrar uma solução
2. **DP Decisório**: Para cada d, verifica se existe uma string central com raio ≤ d
3. **Estados DP**: Mantém vetores de "erros restantes" para cada string do dataset
4. **Construção**: Reconstrói a string central ótima quando encontrada

### Vantagens
- **Solução Exata**: Garante encontrar o raio mínimo possível (solução ótima)
- **Determinístico**: Sempre produz o mesmo resultado para a mesma entrada
- **Matematicamente Rigoroso**: Baseado em teoria sólida de programação dinâmica
- **Verificável**: Resultados podem ser validados facilmente

### Filosofia
O DP-CSP sacrifica tempo de execução em favor de **precisão absoluta**. É a referência para verificar a qualidade de algoritmos heurísticos.

## ⚙️ Funcionamento Detalhado

### Algoritmo Principal: Busca Incremental
```
Para d = 0, 1, 2, ..., max_d:
    Se existe_centro_com_raio(d):
        Retorna centro encontrado
    Senão:
        Tenta próximo d
Se nenhum d funciona:
    Falha (não deveria acontecer)
```

### Subproblema: DP Decisório
**Entrada**: Conjunto de strings S, alfabeto Σ, raio d  
**Saída**: String central c tal que max(H(c,s)) ≤ d, ou NULL se não existir  

### Estados da Programação Dinâmica
```
Estado: (posição, vetor_erros_restantes)
onde:
- posição: índice atual na string sendo construída (0 a L-1)
- vetor_erros_restantes: [r₁, r₂, ..., rₙ]
  rᵢ = número máximo de erros que ainda podemos cometer com string i
```

### Transições de Estado
```
Estado atual: (pos, [r₁, r₂, ..., rₙ])
Para cada caractere σ ∈ Σ:
    Calcula desconto dᵢ = 1 se strings[i][pos] ≠ σ, senão 0
    Novo estado: (pos+1, [r₁-d₁, r₂-d₂, ..., rₙ-dₙ])
    Se min(rᵢ-dᵢ) ≥ 0: estado é viável
```

### Exemplo Detalhado
```
Strings: ["ACG", "ATG", "AAG"]
Alfabeto: "ACGT"
Testando d = 1:

Estado inicial: (0, [1,1,1])  # posição 0, 1 erro permitido para cada

Posição 0:
├── Testa 'A': strings[0][0]='A', strings[1][0]='A', strings[2][0]='A'
│   └── Desconto: [0,0,0] → novo estado: (1, [1,1,1])
├── Testa 'C': descontos [1,1,1] → estado: (1, [0,0,0])
├── Testa 'G': descontos [1,1,1] → estado: (1, [0,0,0])
└── Testa 'T': descontos [1,1,1] → estado: (1, [0,0,0])

Posição 1 (a partir do melhor estado anterior):
├── Testa 'A': descontos [1,1,1] → inviável (estados negativos)
├── Testa 'C': descontos [0,1,1] → estado: (2, [1,0,0])
├── Testa 'G': descontos [1,0,0] → estado: (2, [0,1,1])
└── Testa 'T': descontos [1,0,1] → estado: (2, [0,1,0])

Posição 2:
└── Escolhe 'G' para chegar a estado final viável

Centro encontrado: "ACG" ou "ATG" (ambos têm raio 1)
```

### Reconstrução da Solução
```
A partir do estado final, segue backtrack:
- Estado final: (3, [0,0,0])
- Posição 2: caractere escolhido = 'G'
- Posição 1: caractere escolhido = 'T'  
- Posição 0: caractere escolhido = 'A'
Centro: "ATG"
```

## 🔧 Parâmetros e Configuração

### Parâmetros Principais

| Parâmetro | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| `max_d` | int | Auto | Limite superior para busca de d (usa baseline se None) |
| `max_time` | int | 300 | Timeout em segundos para evitar execução infinita |
| `warn_threshold` | int | 9 | Alerta se (d+1)^n > 10^9 estados |

### Cálculo Automático do max_d
```python
# Se max_d não fornecido, usa baseline como upper bound
baseline = max_distance(strings[0], strings)  # distância da primeira string
max_d = baseline  # garante que pelo menos uma solução existe
```

### Limites de Segurança
```python
# Memória: Estima estados como (d+1)^n
# Para n=5, d=10: (10+1)^5 = 161M estados (~1GB RAM)
# Para n=6, d=10: (10+1)^6 = 1.77B estados (~14GB RAM)

# Tempo: Monitora elapsed_time < max_time
# Cancelamento: Permite interrupção via SIGTERM
```

## 📊 Casos de Uso

### 🟢 Ideal Para:
- **Verificação de Referência**: Validar qualidade de algoritmos heurísticos
- **Instâncias Pequenas**: n ≤ 5-8 strings, comprimentos ≤ 20-50
- **Análise Teórica**: Estudar propriedades exatas do CSP
- **Benchmarking**: Estabelecer lower bounds para comparação

### 🟡 Adequado Para:
- **Protótipos**: Desenvolvimento e teste de novos algoritmos
- **Casos Críticos**: Quando precisão absoluta é essencial
- **Datasets Específicos**: Strings curtas com alfabetos pequenos
- **Pesquisa Acadêmica**: Experimentos controlados

### 🔴 Limitado Para:
- **Instâncias Grandes**: n > 10 strings (explosão exponencial)
- **Strings Longas**: L > 100 caracteres (muitos estados)
- **Aplicações em Tempo Real**: Execução pode ser muito lenta
- **Datasets Reais**: Geralmente muito grandes para DP

## 📈 Análise Algorítmica

### Complexidade Temporal
- **Busca Externa**: O(d_ótimo × F(n,d,L))
- **DP Interno**: O(L × |Σ| × estados_únicos)
- **Estados Únicos**: O((d+1)<sup>n</sup>) no pior caso
- **Total**: O(d × L × |Σ| × (d+1)<sup>n</sup>)

### Complexidade Espacial
- **Armazenamento de Estados**: O((d+1)<sup>n</sup>)
- **Tabela de Transições**: O(L × |Σ| × (d+1)^n)
- **Backtracking**: O(L)
- **Total**: O(L × |Σ| × (d+1)<sup>n</sup>)

### Explosão Exponencial
```
n=3: (d+1)³ estados máximos
n=4: (d+1)⁴ estados máximos  
n=5: (d+1)⁵ estados máximos
n=6: (d+1)⁶ estados máximos → 1B+ estados para d≥10
n=7: (d+1)⁷ estados máximos → inviável para d>5
```

### Performance Estimada
```
n=3, L=10, d≤5:    < 1 segundo
n=4, L=20, d≤8:    1-10 segundos  
n=5, L=30, d≤6:    10 segundos - 2 minutos
n=6, L=40, d≤4:    1-10 minutos
n=7, L≥50:         Provavelmente inviável
```

## 💡 Exemplos de Uso

### Exemplo 1: Instância Pequena
```python
from algorithms.dp_csp import DPCSPAlgorithm

strings = ["ACG", "ATG", "AAG"]
algorithm = DPCSPAlgorithm(strings, alphabet="ACGT")
center, distance, metadata = algorithm.run()

print(f"Centro ótimo: {center}")
print(f"Raio mínimo: {distance}")
print(f"Solução exata: {metadata['solucao_exata']}")
```

### Exemplo 2: Com Limite de max_d
```python
# Limitar busca para evitar timeout
algorithm = DPCSPAlgorithm(
    strings, 
    alphabet="ACGT",
    max_d=5,  # não testa d > 5
    max_time=60  # timeout em 1 minuto
)

try:
    center, distance, metadata = algorithm.run()
    print(f"Solução encontrada: {center} com d={distance}")
except RuntimeError as e:
    print(f"DP-CSP falhou: {e}")
```

### Exemplo 3: Verificação de Benchmark
```python
# Comparar DP-CSP (exato) vs Baseline (heurística)
from algorithms.baseline import BaselineAlgorithm

# Solução exata
dp_center, dp_dist, _ = DPCSPAlgorithm(strings, "ACGT").run()

# Solução heurística  
baseline_center, baseline_dist, _ = BaselineAlgorithm(strings, "ACGT").run()

print(f"DP-CSP (ótimo): d={dp_dist}")
print(f"Baseline: d={baseline_dist}")
print(f"Gap: {baseline_dist - dp_dist} ({100*(baseline_dist-dp_dist)/dp_dist:.1f}%)")
```

### Exemplo 4: Análise de Limites
```python
import time

def test_limits():
    for n in range(3, 8):
        strings = ["A"*10, "T"*10, "G"*10, "C"*10][:n]
        
        try:
            start = time.time()
            center, dist, meta = DPCSPAlgorithm(strings, "ACGT", max_time=30).run()
            elapsed = time.time() - start
            print(f"n={n}: SUCESSO d={dist} em {elapsed:.2f}s")
        except RuntimeError as e:
            print(f"n={n}: FALHOU - {e}")

test_limits()
```

## ⚠️ Limitações

### Limitações Fundamentais
1. **Explosão Exponencial**: (d+1)<sup>n</sup> estados crescem exponencialmente
2. **Limite de Memória**: Pode consumir gigabytes de RAM rapidamente
3. **Timeout**: Execução pode demorar horas/dias para instâncias grandes
4. **Escalabilidade**: Impraticável para n > 8-10 strings

### Limitações Práticas
1. **Datasets Reais**: A maioria é grande demais para DP exato
2. **Aplicações Online**: Latência inaceitável para uso interativo
3. **Recursos Computacionais**: Requer máquinas potentes para n > 6
4. **Implementação**: Complexidade de código vs algoritmos simples

### Cenários Problemáticos
```python
# Caso 1: Muitas strings (explosão exponencial)
strings = ["ACGT"] * 10  # n=10, inviável

# Caso 2: Raio alto necessário (muitos estados)
strings = ["AAAA", "TTTT", "GGGG", "CCCC"]  # d=4, ainda ok para n=4

# Caso 3: Strings longas (muitas posições)
strings = ["A"*1000, "T"*1000, "G"*1000]  # L=1000, muitas iterações

# Caso 4: Alfabeto grande (mais transições)
strings = ["ABC", "DEF", "GHI"]  # alfabeto com 9 letras
```

### Workarounds
```python
# Limite conservador de max_d
algorithm = DPCSPAlgorithm(strings, alphabet, max_d=min(5, len(strings[0])//4))

# Timeout agressivo
algorithm = DPCSPAlgorithm(strings, alphabet, max_time=30)

# Pré-filtro por tamanho
if len(strings) > 8 or len(strings[0]) > 50:
    print("Instância muito grande para DP-CSP, use heurística")
else:
    result = algorithm.run()
```

## 🔗 Integração com CSPBench

### Registro Automático
```python
@register_algorithm
class DPCSPAlgorithm(CSPAlgorithm):
    name = "DP-CSP"
    supports_internal_parallel = False
    is_deterministic = True
```

### Configuração via YAML
```yaml
algorithm:
  name: "DP-CSP"
  params:
    max_d: 5
    max_time: 120
```

### Execução via CLI
```bash
python main.py --algorithm DP-CSP --dataset small_synthetic --max_d 3
```

### Suporte a Paralelização
- **Paralelismo Interno**: ❌ Não suportado (algoritmo sequencial)
- **Paralelismo de Runs**: ✅ Múltiplas execuções independentes
- **Compatibilidade**: ⚠️ Cuidado com consumo de memória em paralelo

### Metadados Retornados
```python
metadata = {
    "iteracoes": 1,
    "max_d_usado": 5,
    "solucao_exata": True,
    "centro_encontrado": "ACGT"
}
```

### Handling de Erros
```python
try:
    center, distance, metadata = algorithm.run()
    # Sucesso: solução ótima encontrada
except RuntimeError as e:
    # Falha: timeout, limite de memória, ou max_d insuficiente
    print(f"DP-CSP não conseguiu resolver: {e}")
    # Fallback para algoritmo heurístico
```

### Troubleshooting

**Problema**: "uso de memória excedeu limite seguro"
```
Solução: Reduzir max_d ou usar máquina com mais RAM
```

**Problema**: "tempo de execução excedeu Xs"
```
Solução: Aumentar max_time ou reduzir tamanho da instância
```

**Problema**: "(d+1)^n excede limite prático"
```
Solução: Usar menos strings ou algoritmo heurístico
```

---

**Desenvolvido para CSPBench** - Framework de Experimentação para o Closest String Problem  
📚 Para mais informações, consulte a [documentação principal](../../README.md) do framework.
