"""
Pacote de Gestão de Datasets - CSPBench

Este pacote fornece um sistema completo para gestão, geração e processamento
de datasets para problemas de Closest String Problem (CSP). Oferece suporte
a múltiplas fontes de dados e formatos, com validação automática e persistência.

Módulos Incluídos:
    - dataset_file: Carregamento de arquivos FASTA e texto
    - dataset_synthetic: Geração de dados sintéticos controlados
    - dataset_entrez: Download de dados biológicos do NCBI
    - dataset_utils: Utilitários para persistência e gerenciamento

Funcionalidades:
    - Carregamento de arquivos em múltiplos formatos
    - Geração de datasets sintéticos parametrizáveis
    - Download automático de dados biológicos reais
    - Validação e normalização automática
    - Persistência em formatos padrão
    - Cache inteligente para otimização

Arquitetura:
    O pacote implementa uma arquitetura modular com:
    - Interfaces padronizadas para diferentes fontes
    - Sistema de validação unificado
    - Tratamento robusto de erros
    - Logging detalhado de operações
    - Configuração flexível via parâmetros

Exemplo de Uso:
    ```python
    # Importar módulos necessários
    from src.datasets import dataset_file, dataset_synthetic, dataset_entrez

    # Carregar de arquivo
    sequences, params = dataset_file.load_dataset_with_params({
        'filepath': 'data/sequences.fasta'
    })

    # Gerar dados sintéticos
    sequences, params = dataset_synthetic.generate_dataset_with_params({
        'n': 100, 'L': 50, 'alphabet': 'ACGT', 'noise': 0.15
    })

    # Download do NCBI
    sequences, params = dataset_entrez.fetch_dataset_silent({
        'email': 'user@example.com',
        'db': 'nucleotide',
        'term': 'COVID-19 spike',
        'n': 50
    })
    ```

Tipos de Dados Suportados:
    - Sequências de DNA/RNA
    - Sequências de proteínas
    - Dados sintéticos parametrizáveis
    - Arquivos FASTA e texto
    - Dados do NCBI via API Entrez

Validação e Qualidade:
    - Verificação de formato e integridade
    - Normalização automática (maiúsculas)
    - Filtragem por comprimento uniforme
    - Validação de alfabeto
    - Relatórios de qualidade dos dados

Autor: CSPBench Development Team
Data: 2024
"""

# Pacote de gestão de datasets do CSPBench
