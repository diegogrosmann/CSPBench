# CSPBench Web Interface - Component Architecture

## Overview
The web interface has been completely redesigned with a modular, component-based architecture that provides an improved user experience.

## Key Features

### 🎯 Execution Type Selection
- **Batch Manager**: Manage batch configuration files  
- **Algorithm Comparison**: Compare multiple algorithms on the same dataset
- **Benchmark Suite**: Run comprehensive benchmarks
- **Custom Workflow**: Create custom execution workflows

### 🧩 Modular Components

#### JavaScript Components
- **APIClient** (`api-client.js`): Centralized API communication
- **DatasetSelector** (`dataset-selector.js`): File upload, samples, generation
- **AlgorithmSelector** (`algorithm-selector.js`): Algorithm selection with filtering
- **ParameterConfig** (`parameter-config.js`): Parameter configuration with validation
- **ExecutionMonitor** (`execution-monitor.js`): Real-time execution tracking
- **ResultsViewer** (`results-viewer.js`): Results display and analysis

#### Templates
Templates are organized by functionality:

- **Base Layout** (`base.html`): Common structure and navigation
- **Dashboard** (`index.html`): Main dashboard with metrics and controls  
- **Results** (`results.html`): Algorithm execution results
- **Batch Manager** (`batch_manager.html`): Batch file management

### 🎨 Modern Design
- CSS custom properties for consistent theming
- Responsive grid layouts
- Card-based UI components
- Loading states and transitions
- Mobile-friendly responsive design

## Usage

### Starting the Web Interface

```bash
# Method 1: Using the launcher script
python src/presentation/web/run_web.py

# Method 2: Direct uvicorn
uvicorn src.presentation.web.app:app --host 0.0.0.0 --port 8000 --reload
```

### Accessing the Interface
Open your browser to: http://localhost:8000

### Workflow

1. **Select Execution Type**: Choose from 6 different execution workflows
2. **Configure Dataset**: Upload file, select samples, or generate synthetic data
3. **Choose Algorithm**: Browse and filter available algorithms with detailed info
4. **Set Parameters**: Configure algorithm parameters with validation and presets
5. **Monitor Execution**: Real-time progress tracking with logs and metrics
6. **View Results**: Comprehensive results display with download options

## API Endpoints

### Core Endpoints
- `GET /` - Main page with execution type selection
- `GET /execution/{type}` - Specific execution workflow pages
- `GET /api/algorithms` - List available algorithms with metadata
- `POST /api/execute` - Execute algorithm with parameters

### Dataset Endpoints  
- `GET /api/sample-datasets` - Get predefined sample datasets
- `POST /api/upload-dataset` - Upload and parse FASTA files
- `POST /api/generate-dataset` - Generate synthetic datasets

### Results Endpoints
- `GET /api/download/{session_id}` - Download results as ZIP file

## Component Integration

### Template Inheritance
```html
<!-- base.html provides common layout -->
{% extends "base.html" %}

<!-- Individual pages extend base -->
{% block content %}
<!-- Page-specific content -->
{% endblock %}
```

### JavaScript Component Usage
```javascript
// Initialize components
window.apiClient = new APIClient('/api');
const datasetSelector = new DatasetSelector('#dataset-container');
const algorithmSelector = new AlgorithmSelector('#algorithm-container');
const parameterConfig = new ParameterConfig('#parameter-container');
const executionMonitor = new ExecutionMonitor('#execution-container');
const resultsViewer = new ResultsViewer('#results-container');
```

### Component Communication
Components communicate through:
- Event callbacks (`onSelect`, `onUpdate`, etc.)
- Shared API client for data consistency
- DOM events for loose coupling
- Local storage for persistence

## Styling System

### CSS Custom Properties
```css
:root {
    --primary-color: #007acc;
    --spacing-md: 1rem;
    --radius-lg: 0.5rem;
    --shadow-md: 0 3px 6px rgba(0, 0, 0, 0.15);
}
```

### Component Classes
```css
.card { /* Base card component */ }
.btn { /* Button component with variants */ }
.execution-card { /* Execution type selection cards */ }
.algorithm-card { /* Algorithm selection cards */ }
```

## Docker Support

The interface includes Docker configuration for containerized deployment:

```dockerfile
# Dockerfile.web
FROM python:3.11-slim
COPY . /app
WORKDIR /app
RUN pip install -r requirements.txt
EXPOSE 8000
CMD ["python", "src/presentation/web/run_web.py"]
```

## Extension Points

### Adding New Execution Types
1. Create new template in `templates/`
2. Add route handler in `app.py`
3. Update execution type list in index template

### Adding New Components
1. Create JavaScript class in `components/`
2. Include in base template
3. Initialize in specific pages

### Customizing Styling
1. Modify CSS custom properties in `style.css`
2. Add component-specific styles
3. Update responsive breakpoints as needed

## Future Enhancements

- Real-time WebSocket communication for live updates
- Advanced visualization with Chart.js/D3.js
- User authentication and session management
- Algorithm performance profiling
- Result comparison matrices
- Export to additional formats (PDF, CSV, etc.)

The modular architecture makes these enhancements straightforward to implement.
