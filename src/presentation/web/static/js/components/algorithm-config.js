/**
 * AlgorithmConfig - Component for algorithm parameter configuration
 */
class AlgorithmConfig {
    constructor(containerId, apiClient) {
        this.container = document.getElementById(containerId);
        this.apiClient = apiClient;
        this.algorithms = [];
        this.selectedAlgorithm = null;
        this.parameters = {};
        
        this.init();
    }
    
    async init() {
        try {
            this.algorithms = await this.apiClient.getAlgorithms();
            this.render();
            this.setupEventListeners();
        } catch (error) {
            console.error('Error loading algorithms:', error);
            this.showError('Failed to load algorithms');
        }
    }
    
    render() {
        this.container.innerHTML = `
            <div class="algorithm-config">
                <div class="form-group">
                    <label class="form-label">Select Algorithm</label>
                    <select class="form-select" id="algorithm-select">
                        <option value="">Choose an algorithm...</option>
                        ${this.algorithms.map(alg => 
                            `<option value="${alg.name}">${alg.name}</option>`
                        ).join('')}
                    </select>
                </div>
                
                <div id="algorithm-info" class="algorithm-info" style="display: none;">
                    <div class="info-section">
                        <h4 class="info-title">Algorithm Information</h4>
                        <div class="algorithm-details">
                            <div class="detail-item">
                                <span class="detail-label">Type:</span>
                                <span class="detail-value" id="algorithm-type"></span>
                            </div>
                            <div class="detail-item">
                                <span class="detail-label">Deterministic:</span>
                                <span class="detail-value" id="algorithm-deterministic"></span>
                            </div>
                            <div class="detail-item">
                                <span class="detail-label">Parallel Support:</span>
                                <span class="detail-value" id="algorithm-parallel"></span>
                            </div>
                        </div>
                        <div class="algorithm-description" id="algorithm-description"></div>
                    </div>
                </div>
                
                <div id="parameters-section" class="parameters-section" style="display: none;">
                    <h4 class="section-title">Algorithm Parameters</h4>
                    <div id="parameters-container" class="parameters-container"></div>
                    <div class="parameter-actions">
                        <button type="button" class="btn btn-secondary" id="reset-params">Reset to Defaults</button>
                        <button type="button" class="btn btn-secondary" id="export-params">Export Config</button>
                        <button type="button" class="btn btn-secondary" id="import-params">Import Config</button>
                    </div>
                </div>
            </div>
        `;
    }
    
    setupEventListeners() {
        const algorithmSelect = document.getElementById('algorithm-select');
        const resetBtn = document.getElementById('reset-params');
        const exportBtn = document.getElementById('export-params');
        const importBtn = document.getElementById('import-params');
        
        algorithmSelect.addEventListener('change', (e) => {
            this.selectAlgorithm(e.target.value);
        });
        
        resetBtn?.addEventListener('click', () => {
            this.resetParameters();
        });
        
        exportBtn?.addEventListener('click', () => {
            this.exportConfiguration();
        });
        
        importBtn?.addEventListener('click', () => {
            this.importConfiguration();
        });
    }
    
    selectAlgorithm(algorithmName) {
        if (!algorithmName) {
            this.selectedAlgorithm = null;
            this.hideAlgorithmInfo();
            this.hideParameters();
            return;
        }
        
        this.selectedAlgorithm = this.algorithms.find(alg => alg.name === algorithmName);
        if (this.selectedAlgorithm) {
            this.showAlgorithmInfo();
            this.renderParameters();
            this.notifySelectionChange();
        }
    }
    
    showAlgorithmInfo() {
        const infoSection = document.getElementById('algorithm-info');
        const typeElement = document.getElementById('algorithm-type');
        const deterministicElement = document.getElementById('algorithm-deterministic');
        const parallelElement = document.getElementById('algorithm-parallel');
        const descriptionElement = document.getElementById('algorithm-description');
        
        if (this.selectedAlgorithm) {
            typeElement.textContent = this.selectedAlgorithm.is_deterministic ? 'Deterministic' : 'Stochastic';
            deterministicElement.textContent = this.selectedAlgorithm.is_deterministic ? 'Yes' : 'No';
            parallelElement.textContent = this.selectedAlgorithm.supports_internal_parallel ? 'Yes' : 'No';
            descriptionElement.textContent = this.selectedAlgorithm.description;
            
            infoSection.style.display = 'block';
        }
    }
    
    hideAlgorithmInfo() {
        const infoSection = document.getElementById('algorithm-info');
        infoSection.style.display = 'none';
    }
    
    renderParameters() {
        const parametersSection = document.getElementById('parameters-section');
        const container = document.getElementById('parameters-container');
        
        if (!this.selectedAlgorithm || !this.selectedAlgorithm.default_params) {
            this.hideParameters();
            return;
        }
        
        const params = this.selectedAlgorithm.default_params;
        this.parameters = { ...params }; // Clone default parameters
        
        container.innerHTML = '';
        
        Object.entries(params).forEach(([key, value]) => {
            const paramElement = this.createParameterElement(key, value);
            container.appendChild(paramElement);
        });
        
        parametersSection.style.display = 'block';
    }
    
    createParameterElement(key, value) {
        const wrapper = document.createElement('div');
        wrapper.className = 'parameter-item';
        
        const type = this.getParameterType(value);
        const inputId = `param-${key}`;
        
        wrapper.innerHTML = `
            <label class="parameter-label" for="${inputId}">
                ${this.formatParameterName(key)}
            </label>
            <div class="parameter-input-wrapper">
                ${this.createParameterInput(inputId, key, value, type)}
                <span class="parameter-type">${type}</span>
            </div>
            <div class="parameter-description">
                ${this.getParameterDescription(key)}
            </div>
        `;
        
        // Setup event listener for parameter changes
        const input = wrapper.querySelector('input, select');
        input.addEventListener('change', (e) => {
            this.updateParameter(key, e.target.value, type);
        });
        
        return wrapper;
    }
    
    createParameterInput(id, key, value, type) {
        switch (type) {
            case 'boolean':
                return `
                    <select id="${id}" class="parameter-input">
                        <option value="true" ${value === true ? 'selected' : ''}>True</option>
                        <option value="false" ${value === false ? 'selected' : ''}>False</option>
                    </select>
                `;
            case 'string':
                return `<input type="text" id="${id}" class="parameter-input" value="${value || ''}" placeholder="Enter text...">`;
            case 'integer':
                return `<input type="number" id="${id}" class="parameter-input" value="${value || 0}" step="1">`;
            case 'float':
                return `<input type="number" id="${id}" class="parameter-input" value="${value || 0}" step="0.001">`;
            case 'null':
                return `<input type="text" id="${id}" class="parameter-input" value="" placeholder="null (auto)">`;
            default:
                return `<input type="text" id="${id}" class="parameter-input" value="${value || ''}" placeholder="Enter value...">`;
        }
    }
    
    getParameterType(value) {
        if (value === null) return 'null';
        if (typeof value === 'boolean') return 'boolean';
        if (typeof value === 'string') return 'string';
        if (Number.isInteger(value)) return 'integer';
        if (typeof value === 'number') return 'float';
        return 'string';
    }
    
    formatParameterName(key) {
        return key.replace(/_/g, ' ')
                 .replace(/\b\w/g, l => l.toUpperCase());
    }
    
    getParameterDescription(key) {
        const descriptions = {
            pop_size: 'Population size multiplier',
            min_pop_size: 'Minimum population size',
            seed: 'Random seed for reproducibility (null for random)',
            max_gens: 'Maximum number of generations',
            max_time: 'Maximum execution time in seconds',
            cross_prob: 'Crossover probability (0.0 to 1.0)',
            mut_prob: 'Mutation probability (0.0 to 1.0)',
            elite_rate: 'Elite selection rate (0.0 to 1.0)',
            beam_width: 'Beam search width',
            max_d: 'Maximum distance threshold',
            warn_threshold: 'Warning threshold for large problems',
            min_d: 'Minimum distance parameter',
            d_factor: 'Distance factor',
            local_iters: 'Number of local search iterations'
        };
        
        return descriptions[key] || 'Algorithm parameter';
    }
    
    updateParameter(key, value, type) {
        let parsedValue = value;
        
        switch (type) {
            case 'boolean':
                parsedValue = value === 'true';
                break;
            case 'integer':
                parsedValue = value === '' ? null : parseInt(value);
                break;
            case 'float':
                parsedValue = value === '' ? null : parseFloat(value);
                break;
            case 'null':
                parsedValue = value === '' ? null : value;
                break;
        }
        
        this.parameters[key] = parsedValue;
        this.notifyParameterChange(key, parsedValue);
    }
    
    resetParameters() {
        if (this.selectedAlgorithm) {
            this.parameters = { ...this.selectedAlgorithm.default_params };
            this.renderParameters();
            this.notifyParametersReset();
        }
    }
    
    exportConfiguration() {
        if (!this.selectedAlgorithm) {
            alert('Please select an algorithm first');
            return;
        }
        
        const config = {
            algorithm: this.selectedAlgorithm.name,
            parameters: this.parameters,
            timestamp: new Date().toISOString()
        };
        
        const blob = new Blob([JSON.stringify(config, null, 2)], { type: 'application/json' });
        const url = URL.createObjectURL(blob);
        
        const a = document.createElement('a');
        a.href = url;
        a.download = `${this.selectedAlgorithm.name.toLowerCase()}_config.json`;
        a.click();
        
        URL.revokeObjectURL(url);
    }
    
    importConfiguration() {
        const input = document.createElement('input');
        input.type = 'file';
        input.accept = '.json';
        
        input.onchange = (e) => {
            const file = e.target.files[0];
            if (!file) return;
            
            const reader = new FileReader();
            reader.onload = (e) => {
                try {
                    const config = JSON.parse(e.target.result);
                    
                    if (config.algorithm && config.parameters) {
                        // Select algorithm
                        const algorithmSelect = document.getElementById('algorithm-select');
                        algorithmSelect.value = config.algorithm;
                        this.selectAlgorithm(config.algorithm);
                        
                        // Update parameters
                        this.parameters = { ...config.parameters };
                        this.renderParameters();
                        
                        alert('Configuration imported successfully');
                    } else {
                        alert('Invalid configuration file format');
                    }
                } catch (error) {
                    alert('Error reading configuration file');
                }
            };
            reader.readAsText(file);
        };
        
        input.click();
    }
    
    hideParameters() {
        const parametersSection = document.getElementById('parameters-section');
        parametersSection.style.display = 'none';
    }
    
    showError(message) {
        this.container.innerHTML = `
            <div class="error-message">
                <div class="error-icon">⚠️</div>
                <div class="error-text">${message}</div>
            </div>
        `;
    }
    
    // Event notification methods
    notifySelectionChange() {
        const event = new CustomEvent('algorithmSelected', {
            detail: {
                algorithm: this.selectedAlgorithm,
                parameters: this.parameters
            }
        });
        this.container.dispatchEvent(event);
    }
    
    notifyParameterChange(key, value) {
        const event = new CustomEvent('parameterChanged', {
            detail: {
                algorithm: this.selectedAlgorithm,
                parameter: key,
                value: value,
                parameters: this.parameters
            }
        });
        this.container.dispatchEvent(event);
    }
    
    notifyParametersReset() {
        const event = new CustomEvent('parametersReset', {
            detail: {
                algorithm: this.selectedAlgorithm,
                parameters: this.parameters
            }
        });
        this.container.dispatchEvent(event);
    }
    
    // Public methods for external access
    getSelectedAlgorithm() {
        return this.selectedAlgorithm;
    }
    
    getParameters() {
        return { ...this.parameters };
    }
    
    setParameters(params) {
        this.parameters = { ...params };
        this.renderParameters();
    }
    
    getConfiguration() {
        return {
            algorithm: this.selectedAlgorithm?.name || null,
            parameters: this.parameters
        };
    }
    
    isValid() {
        return this.selectedAlgorithm !== null;
    }
}

// Export para uso global
window.AlgorithmConfig = AlgorithmConfig;
