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
