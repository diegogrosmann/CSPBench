// CSPBench Web Interface JavaScript

class CSPBenchApp {
    constructor() {
        this.algorithms = [];
        this.currentExecution = null;
        this.init();
    }

    async init() {
        await this.loadAlgorithms();
        this.setupEventListeners();
        this.setupFileUpload();
    }

    async loadAlgorithms() {
        try {
            const response = await fetch('/api/algorithms');
            this.algorithms = await response.json();
            this.populateAlgorithmSelect();
        } catch (error) {
            console.error('Failed to load algorithms:', error);
            this.showError('Failed to load algorithms. Please refresh the page.');
        }
    }

    populateAlgorithmSelect() {
        const select = document.getElementById('algorithm-select');
        select.innerHTML = '<option value="">Select an algorithm...</option>';
        
        this.algorithms.forEach(algorithm => {
            const option = document.createElement('option');
            option.value = algorithm.name;
            option.textContent = algorithm.name;
            select.appendChild(option);
        });
    }

    setupEventListeners() {
        // Algorithm selection change
        document.getElementById('algorithm-select').addEventListener('change', (e) => {
            this.onAlgorithmChange(e.target.value);
        });

        // Form submission
        document.getElementById('execution-form').addEventListener('submit', (e) => {
            e.preventDefault();
            this.executeAlgorithm();
        });
    }

    setupFileUpload() {
        const fileInput = document.getElementById('dataset-file');
        const textArea = document.getElementById('dataset-input');

        fileInput.addEventListener('change', (e) => {
            const file = e.target.files[0];
            if (file) {
                const reader = new FileReader();
                reader.onload = (e) => {
                    textArea.value = e.target.result;
                };
                reader.readAsText(file);
            }
        });
    }

    onAlgorithmChange(algorithmName) {
        const algorithm = this.algorithms.find(a => a.name === algorithmName);
        if (algorithm) {
            this.showAlgorithmInfo(algorithm);
            this.generateParametersForm(algorithm);
        } else {
            this.clearParametersForm();
        }
    }

    showAlgorithmInfo(algorithm) {
        // You could show algorithm description somewhere
        console.log('Selected algorithm:', algorithm);
    }

    generateParametersForm(algorithm) {
        const container = document.getElementById('parameters-container');
        container.innerHTML = '';

        if (Object.keys(algorithm.default_params).length === 0) {
            container.innerHTML = '<p>This algorithm has no configurable parameters.</p>';
            return;
        }

        Object.entries(algorithm.default_params).forEach(([key, value]) => {
            const paramDiv = document.createElement('div');
            paramDiv.className = 'parameter-item';

            const label = document.createElement('label');
            label.textContent = this.formatParameterName(key);
            label.setAttribute('for', `param-${key}`);

            const input = document.createElement('input');
            input.type = this.getInputType(value);
            input.id = `param-${key}`;
            input.name = key;
            input.value = value;
            input.className = 'parameter-input';

            if (input.type === 'number') {
                input.step = this.isInteger(value) ? '1' : '0.01';
                if (key.includes('rate') || key.includes('probability')) {
                    input.min = '0';
                    input.max = '1';
                    input.step = '0.01';
                }
            }

            paramDiv.appendChild(label);
            paramDiv.appendChild(input);
            container.appendChild(paramDiv);
        });
    }

    clearParametersForm() {
        const container = document.getElementById('parameters-container');
        container.innerHTML = '<p>Select an algorithm to see available parameters</p>';
    }

    formatParameterName(name) {
        return name.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase());
    }

    getInputType(value) {
        if (typeof value === 'number') return 'number';
        if (typeof value === 'boolean') return 'checkbox';
        return 'text';
    }

    isInteger(value) {
        return Number.isInteger(value);
    }

    async executeAlgorithm() {
        const form = document.getElementById('execution-form');
        const formData = new FormData(form);
        
        // Validate form
        const algorithm = formData.get('algorithm');
        const dataset = formData.get('dataset');
        
        if (!algorithm) {
            this.showError('Please select an algorithm.');
            return;
        }
        
        if (!dataset.trim()) {
            this.showError('Please provide a dataset.');
            return;
        }

        // Collect parameters
        const parameters = {};
        const paramInputs = document.querySelectorAll('.parameter-input');
        paramInputs.forEach(input => {
            const value = input.type === 'checkbox' ? input.checked : 
                         input.type === 'number' ? parseFloat(input.value) : 
                         input.value;
            parameters[input.name] = value;
        });

        // Prepare request
        const request = {
            algorithm: algorithm,
            dataset_content: dataset,
            dataset_name: 'uploaded_dataset.fasta',
            parameters: parameters,
            save_history: document.getElementById('save-history').checked,
            timeout: parseInt(document.getElementById('timeout').value)
        };

        // Show progress
        this.showProgress();
        this.disableForm();

        try {
            const response = await fetch('/api/execute', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify(request)
            });

            const result = await response.json();
            
            if (result.status === 'completed') {
                this.showResults(result);
            } else {
                this.showError(result.error || 'Execution failed');
            }
        } catch (error) {
            console.error('Execution failed:', error);
            this.showError('Network error. Please try again.');
        } finally {
            this.hideProgress();
            this.enableForm();
        }
    }

    showProgress() {
        const resultsSection = document.getElementById('results-section');
        const progressIndicator = document.getElementById('progress-indicator');
        const resultsContent = document.getElementById('results-content');
        
        resultsSection.style.display = 'block';
        progressIndicator.style.display = 'block';
        resultsContent.style.display = 'none';
        
        resultsSection.scrollIntoView({ behavior: 'smooth' });
    }

    hideProgress() {
        const progressIndicator = document.getElementById('progress-indicator');
        progressIndicator.style.display = 'none';
    }

    showResults(result) {
        const resultsContent = document.getElementById('results-content');
        
        const html = `
            <div class="success-message fade-in">
                <strong>✅ Algorithm executed successfully!</strong>
            </div>
            
            <div class="result-card fade-in">
                <h3>🎯 Best Solution</h3>
                <div class="result-value">${result.result.best_string || 'N/A'}</div>
            </div>
            
            <div class="result-card fade-in">
                <h3>📏 Maximum Distance</h3>
                <div class="result-value">${result.result.max_distance || 'N/A'}</div>
            </div>
            
            <div class="result-card fade-in">
                <h3>⏱️ Execution Time</h3>
                <div class="result-value">${result.result.execution_time ? result.result.execution_time.toFixed(4) + ' seconds' : 'N/A'}</div>
            </div>
            
            <div class="result-card fade-in">
                <h3>📊 Metadata</h3>
                <div class="result-value"><pre>${JSON.stringify(result.result.metadata, null, 2)}</pre></div>
            </div>
            
            <div class="result-card fade-in">
                <h3>💾 Download Results</h3>
                <p>Download a ZIP file containing all results and metadata:</p>
                <a href="${result.download_url}" class="download-button" download>
                    📥 Download Results ZIP
                </a>
            </div>
        `;
        
        resultsContent.innerHTML = html;
        resultsContent.style.display = 'block';
    }

    showError(message) {
        const resultsSection = document.getElementById('results-section');
        const resultsContent = document.getElementById('results-content');
        
        resultsContent.innerHTML = `
            <div class="error-message fade-in">
                <strong>❌ Error:</strong> ${message}
            </div>
        `;
        
        resultsSection.style.display = 'block';
        resultsContent.style.display = 'block';
        
        resultsSection.scrollIntoView({ behavior: 'smooth' });
    }

    disableForm() {
        const form = document.getElementById('execution-form');
        const inputs = form.querySelectorAll('input, select, textarea, button');
        inputs.forEach(input => {
            input.disabled = true;
        });
        form.classList.add('loading');
    }

    enableForm() {
        const form = document.getElementById('execution-form');
        const inputs = form.querySelectorAll('input, select, textarea, button');
        inputs.forEach(input => {
            input.disabled = false;
        });
        form.classList.remove('loading');
    }
}

// Initialize app when DOM is loaded
document.addEventListener('DOMContentLoaded', () => {
    new CSPBenchApp();
});
