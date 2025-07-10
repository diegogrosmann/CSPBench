"""
Pacote de Utilitários do CSPBench

Este pacote contém módulos utilitários essenciais para o funcionamento
do framework CSPBench. Oferece funcionalidades de apoio para configuração,
cálculos, logging, monitoramento e gestão de recursos.

MÓDULOS DISPONÍVEIS:
===================
- config: Gestão de configurações e parâmetros
- distance: Cálculos de distância entre strings
- logging: Sistema de logging estruturado
- curses_console: Interface curses para monitoramento
- signal_manager: Gestão de sinais do sistema
- worker_calculator: Cálculos de paralelização

FUNCIONALIDADES:
===============
- Configuração centralizada e flexível
- Cálculos de distância otimizados
- Logging estruturado com níveis configuráveis
- Interface visual para monitoramento em tempo real
- Gestão robusta de sinais e interrupções
- Cálculo inteligente de workers para paralelização

EXEMPLO DE USO:
==============
```python
from src.utils.distance import hamming_distance
from src.utils.config import load_config
from src.utils.logging import setup_logging

# Configurar logging
setup_logging(level="INFO")

# Calcular distância
dist = hamming_distance("ACGT", "AGCT")

# Carregar configuração
config = load_config("config.yaml")
```

DESIGN PATTERNS:
===============
- Singleton para configurações globais
- Factory para criação de loggers
- Observer para monitoramento de recursos
- Strategy para diferentes tipos de cálculos
"""
