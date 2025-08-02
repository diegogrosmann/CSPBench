# Interface de Download de Resultados - CSPBench

## Visão Geral

Foi implementada uma interface web completa para gerenciar e baixar resultados de execuções do CSPBench. A funcionalidade inclui:

## Funcionalidades Implementadas

### 1. Página de Resultados Principal (`/results`)
- **Listagem de sessões**: Mostra todas as execuções com informações como algoritmo, dataset, timestamp, status, etc.
- **Filtros avançados**: Pesquisa por ID de sessão, algoritmo, status, e intervalo de datas
- **Paginação**: Suporte para grandes quantidades de resultados
- **Seleção múltipla**: Permite selecionar várias sessões para operações em lote
- **Resumo estatístico**: Mostra totais de sessões, tamanho total, sessões completadas/falhadas

### 2. Downloads Individuais
- **Download completo**: ZIP com todos os arquivos da sessão
- **Downloads específicos**: 
  - Resultados JSON
  - Resultados CSV
  - Arquivo de log
  - Plots e gráficos
  - Relatório HTML

### 3. Downloads em Lote
- **Múltiplas sessões**: Baixar várias sessões selecionadas em um único ZIP
- **Limite de segurança**: Máximo de 20 sessões por lote

### 4. Página de Detalhes da Sessão (`/results/{session_id}`)
- **Informações detalhadas**: Dados completos da sessão
- **Lista de arquivos**: Todos os arquivos gerados na execução
- **Downloads rápidos**: Botões para diferentes tipos de arquivos
- **Resultados detalhados**: Exibição dos dados JSON da execução

### 5. Gerenciamento de Sessões
- **Exclusão**: Deletar sessões individuais ou múltiplas
- **Confirmação**: Modal de confirmação para operações destrutivas

## APIs Implementadas

### Endpoints de Resultados (`/api/results`)

#### `GET /api/results/`
Lista todas as sessões com filtros opcionais:
- `limit`, `offset`: Paginação
- `algorithm`: Filtrar por algoritmo
- `status`: Filtrar por status
- `date_from`, `date_to`: Filtrar por período

#### `GET /api/results/{session_id}`
Detalhes completos de uma sessão específica.

#### `GET /api/results/{session_id}/download`
Download do ZIP completo da sessão.

#### `GET /api/results/{session_id}/files/{file_type}`
Download de tipo específico de arquivo:
- `json`: Arquivos JSON de resultados
- `csv`: Arquivos CSV de resultados
- `log`: Arquivo de log da execução
- `plots`: Gráficos e visualizações
- `report`: Relatório HTML

#### `POST /api/results/download/batch`
Download em lote de múltiplas sessões.
Body: `["session_id1", "session_id2", ...]`

#### `DELETE /api/results/{session_id}`
Exclui uma sessão específica.

## Componentes Frontend

### ResultsManager
Classe JavaScript que gerencia a interface principal de resultados:
- Carregamento e filtros de dados
- Paginação e seleção
- Operações em lote
- Interface reativa

### ResultDetailsManager
Classe JavaScript para a página de detalhes:
- Carregamento de informações detalhadas
- Downloads individuais
- Ações de gerenciamento

## Estrutura de Arquivos

```
src/presentation/web/
├── routes/
│   └── results.py              # APIs de gerenciamento de resultados
├── templates/
│   ├── results.html            # Página principal de resultados
│   └── result_details.html     # Página de detalhes da sessão
└── static/js/components/
    └── results-manager.js      # Componente JavaScript principal
```

## Recursos de Segurança

- **Validação de entrada**: Todos os parâmetros são validados
- **Limitação de lote**: Máximo de 20 sessões por download em lote
- **Sanitização de paths**: Prevenção de ataques de path traversal
- **Gerenciamento de arquivos temporários**: Limpeza automática de ZIPs temporários

## Recursos de Usabilidade

- **Interface responsiva**: Funciona em desktop e mobile
- **Feedback visual**: Loading states, toasts de sucesso/erro
- **Navegação intuitiva**: Breadcrumbs e navegação contextual
- **Filtros persistentes**: Mantém filtros aplicados durante navegação
- **Paginação inteligente**: Navegação eficiente entre páginas

## Como Usar

### Acessar a Interface
1. Inicie a aplicação web: `python src/presentation/web/run_web.py`
2. Acesse: http://localhost:8000/results

### Operações Básicas
1. **Visualizar resultados**: A página carrega automaticamente todos os resultados
2. **Filtrar**: Use os filtros no topo da página
3. **Baixar sessão**: Clique no botão de download na linha da sessão
4. **Ver detalhes**: Clique no ID da sessão ou no botão de visualização
5. **Download em lote**: Selecione múltiplas sessões e clique em "Download Selected"

### Gerenciamento
- **Deletar**: Use o botão de lixeira para deletar sessões
- **Atualizar**: Use o botão "Refresh" para recarregar dados
- **Navegação**: Use os links de navegação e breadcrumbs

## Benefícios

1. **Centralização**: Um local único para gerenciar todos os resultados
2. **Eficiência**: Downloads em lote e filtros avançados
3. **Flexibilidade**: Múltiplos formatos de download
4. **Usabilidade**: Interface intuitiva e responsiva
5. **Escalabilidade**: Suporte para grandes quantidades de dados
6. **Manutenibilidade**: Código modular e bem organizado

A interface fornece uma solução completa para gerenciamento de resultados do CSPBench, permitindo que usuários encontrem, visualizem e baixem resultados de forma eficiente e intuitiva.
