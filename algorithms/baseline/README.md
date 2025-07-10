# Baseline: Algoritmo de Consenso Ganancioso

O **Baseline** é um algoritmo determinístico simples e eficiente que implementa uma estratégia de consenso ganancioso para resolver o Closest String Problem. Serve como referência fundamental para comparação com métodos mais sofisticados.

## 📊 Visão Geral

### **Estratégia Principal**
- **Consenso por Posição**: Para cada posição, escolhe o símbolo mais frequente
- **Decisão Gananciosa**: Toma decisões localmente ótimas sem considerar impacto global
- **Determinístico**: Sempre produz o mesmo resultado para a mesma entrada
- **Eficiência**: Execução linear em O(n × L × |Σ|)

### **Funcionamento**
1. Para cada posição i ∈ [0, L-1]:
   - Conta a frequência de cada símbolo do alfabeto
   - Seleciona o símbolo com maior frequência
   - Em caso de empate, escolhe o primeiro símbolo alfabeticamente
2. Constrói a string consenso concatenando os símbolos escolhidos
3. Calcula a distância máxima de Hamming para todas as strings de entrada

## 🔧 Características Técnicas

### **Complexidade**
- **Temporal**: O(n × L × |Σ|)
  - n: número de strings
  - L: comprimento das strings  
  - |Σ|: tamanho do alfabeto
- **Espacial**: O(|Σ|) para contadores + O(L) para resultado

### **Propriedades**
- ✅ **Determinístico**: Sempre produz o mesmo resultado
- ✅ **Rápido**: Execução quase instantânea
- ✅ **Simples**: Implementação direta e compreensível
- ✅ **Estável**: Não há parâmetros para ajustar
- ❌ **Qualidade**: Pode não encontrar o ótimo global
- ❌ **Independência**: Não considera dependências entre posições

## 🎯 Casos de Uso

### **✅ Quando Usar**
- **Baseline de Comparação**: Estabelecer linha de base para outros algoritmos
- **Execução Rápida**: Quando tempo é extremamente limitado
- **Dados com Consenso Forte**: Sequências com posições bem conservadas
- **Pré-processamento**: Solução inicial para algoritmos iterativos
- **Validação**: Verificar funcionamento básico do framework

### **❌ Limitações**
- **Ótimos Locais**: Pode ficar preso em soluções subótimas
- **Empates**: Resolução arbitrária pode impactar qualidade
- **Dados Ruidosos**: Performance degradada com muito ruído
- **Dependências**: Ignora correlações entre posições

## 🧮 Parâmetros

O algoritmo Baseline **não possui parâmetros configuráveis**, garantindo:
- Reprodutibilidade total
- Simplicidade de uso
- Ausência de tuning necessário
- Comportamento consistente

## 💻 Exemplo de Uso

### **Uso Básico**
```python
from algorithms.baseline.algorithm import BaselineAlg

# Dataset de exemplo
strings = ["ACGTACGT", "AGGTACGT", "ACGTAAGT"]
alphabet = "ACGT"

# Criar e executar algoritmo
algorithm = BaselineAlg(strings, alphabet)
center, distance, metadata = algorithm.run()

print(f"Centro encontrado: {center}")
print(f"Distância máxima: {distance}")
print(f"Metadados: {metadata}")
```

### **Via Framework**
```bash
# Execução via CLI
python main.py --algorithms Baseline --dataset synthetic

# Execução silenciosa
python main.py --silent --algorithms Baseline --dataset synthetic --num-execs 1
```

### **Em Lote (YAML)**
```yaml
algorithms: ["Baseline"]
task:
  type: "execution"
  execution:
    executions:
      - nome: "Teste Baseline"
        dataset: dataset_1
        runs_per_algorithm_per_base: 1  # Determinístico
        timeout: 30
```

## 🔬 Análise Algorítmica

### **Pseudocódigo**
```
function baseline_consensus(strings, alphabet):
    L = length(strings[0])
    consensus = ""
    
    for position in range(L):
        # Contar frequências
        counts = {}
        for symbol in alphabet:
            counts[symbol] = 0
        
        for string in strings:
            symbol = string[position]
            counts[symbol] += 1
        
        # Encontrar símbolo mais frequente
        max_count = 0
        best_symbol = alphabet[0]  # Tie-breaking
        
        for symbol in alphabet:
            if counts[symbol] > max_count:
                max_count = counts[symbol]
                best_symbol = symbol
        
        consensus += best_symbol
    
    return consensus
```

### **Análise Matemática**
Para uma posição i, seja f(s,i) a frequência do símbolo s na posição i:
- Escolha: argmax_s f(s,i)
- Distância esperada por posição: ≈ n × (1 - max_s(f(s,i)/n))
- Distância total esperada: Σ_i n × (1 - max_s(f(s,i)/n))

## 🎨 Visualizações

### **Análise de Consenso**
```python
# Gerar heatmap de consenso por posição
import matplotlib.pyplot as plt
import numpy as np

def visualize_consensus(strings, alphabet):
    L = len(strings[0])
    n = len(strings)
    
    # Matriz de frequências
    freq_matrix = np.zeros((len(alphabet), L))
    
    for i, symbol in enumerate(alphabet):
        for pos in range(L):
            count = sum(1 for s in strings if s[pos] == symbol)
            freq_matrix[i, pos] = count / n
    
    # Plotar heatmap
    plt.imshow(freq_matrix, aspect='auto', cmap='viridis')
    plt.xlabel('Posição')
    plt.ylabel('Símbolo')
    plt.colorbar(label='Frequência')
    plt.title('Consenso por Posição')
    plt.show()
```

## 🔗 Integração com CSPBench

O Baseline está totalmente integrado ao framework através de:

- **Registro Automático**: Detectado via `@register_algorithm`
- **Interface Padronizada**: Implementa `CSPAlgorithm`
- **Execução Paralela**: Compatível com sistema de execução
- **Relatórios**: Gera metadados estruturados
- **Monitoramento**: Suporte a callbacks de progresso

---

*Baseline: A base sólida para comparação de algoritmos CSP - simples, rápido e confiável.*
