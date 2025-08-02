/**
 * Modern Algorithm Selector Component
 * Bootstrap-styled algorithm selection with enhanced UX
 */

class AlgorithmSelector {
    constructor(container, options = {}) {
        this.container = typeof container === 'string' ? document.querySelector(container) : container;
        this.options = {
            apiClient: window.apiClient,
            showDescription: true,
            showParameters: true,
            autoSelect: false,
            onSelectionChange: null,
            ...options
        };
        
        this.algorithms = [];
        this.selectedAlgorithm = null;
        this.isLoading = false;
        
        this.init();
    }

    async init() {
        this.render();
        await this.loadAlgorithms();
        this.bindEvents();
    }

    render() {
        this.container.innerHTML = `
            <div class="algorithm-selector">
                <div class="algorithm-selector__header">
                    <h5 class="mb-0">
                        <i class="bi bi-cpu me-2 text-primary"></i>
                        Select Algorithm
                    </h5>
                </div>
                
                <div class="algorithm-selector__content">
                    <div class="algorithm-selector__loading" style="display: none;">
                        <div class="d-flex align-items-center justify-content-center py-4">
                            <div class="spinner-border text-primary me-3" role="status">
                                <span class="visually-hidden">Loading algorithms...</span>
                            </div>
                            <span class="text-muted">Loading available algorithms...</span>
                        </div>
                    </div>
                    
                    <div class="algorithm-selector__dropdown-container" style="display: none;">
                        <select class="form-select algorithm-selector__dropdown" id="algorithm-select">
                            <option value="">Choose an algorithm...</option>
                        </select>
                    </div>
                    
                    <div class="algorithm-selector__error" style="display: none;">
                        <div class="alert alert-danger">
                            <i class="bi bi-exclamation-triangle me-2"></i>
                            <span class="error-message">Failed to load algorithms</span>
                        </div>
                    </div>
                    
                    <div class="algorithm-selector__details mt-3" style="display: none;">
                        <div class="algorithm-description">
                            <h6 class="text-primary">
                                <i class="bi bi-info-circle me-2"></i>
                                Description
                            </h6>
                            <p class="text-muted algorithm-description__content"></p>
                        </div>
                        
                        <div class="algorithm-parameters mt-3">
                            <h6 class="text-primary">
                                <i class="bi bi-sliders me-2"></i>
                                Default Parameters
                            </h6>
                            <div class="parameters-list"></div>
                        </div>
                    </div>
                </div>
            </div>
        `;
    }

    async loadAlgorithms() {
        try {
            this.setLoading(true);
            this.algorithms = await this.options.apiClient.getAlgorithms();
            this.populateDropdown();
            this.showDropdown();
        } catch (error) {
            console.error('Failed to load algorithms:', error);
            this.showError('Failed to load algorithms. Please try again.');
        } finally {
            this.setLoading(false);
        }
    }

    populateDropdown() {
        const dropdown = this.container.querySelector('.algorithm-selector__dropdown');
        
        // Clear existing options except the first one
        while (dropdown.children.length > 1) {
            dropdown.removeChild(dropdown.lastChild);
        }
        
        // Add algorithm options
        this.algorithms.forEach(algorithm => {
            const option = document.createElement('option');
            option.value = algorithm.name;
            option.textContent = `${algorithm.name} - ${algorithm.description || 'No description'}`;
            dropdown.appendChild(option);
        });

        // Auto-select first algorithm if requested
        if (this.options.autoSelect && this.algorithms.length > 0) {
            dropdown.value = this.algorithms[0].name;
            this.selectAlgorithm(this.algorithms[0].name);
        }
    }

    bindEvents() {
        const dropdown = this.container.querySelector('.algorithm-selector__dropdown');
        
        dropdown.addEventListener('change', (event) => {
            const algorithmName = event.target.value;
            this.selectAlgorithm(algorithmName);
        });
    }

    selectAlgorithm(algorithmName) {
        if (!algorithmName) {
            this.selectedAlgorithm = null;
            this.hideDetails();
            this.notifySelectionChange();
            return;
        }

        const algorithm = this.algorithms.find(alg => alg.name === algorithmName);
        if (!algorithm) {
            console.warn(`Algorithm not found: ${algorithmName}`);
            return;
        }

        this.selectedAlgorithm = algorithm;
        this.showDetails(algorithm);
        this.notifySelectionChange();
    }

    showDetails(algorithm) {
        if (!this.options.showDescription && !this.options.showParameters) {
            return;
        }

        const detailsContainer = this.container.querySelector('.algorithm-selector__details');
        
        // Update description
        if (this.options.showDescription) {
            const descriptionContent = this.container.querySelector('.algorithm-description__content');
            descriptionContent.textContent = algorithm.description || 'No description available.';
        }

        // Update parameters
        if (this.options.showParameters) {
            const parametersList = this.container.querySelector('.parameters-list');
            this.renderParameters(parametersList, algorithm.default_params || {});
        }

        detailsContainer.style.display = 'block';
    }

    hideDetails() {
        const detailsContainer = this.container.querySelector('.algorithm-selector__details');
        detailsContainer.style.display = 'none';
    }

    renderParameters(container, parameters) {
        if (!parameters || Object.keys(parameters).length === 0) {
            container.innerHTML = '<p class="text-muted small">No default parameters specified.</p>';
            return;
        }

        const parameterCards = Object.entries(parameters).map(([key, value]) => {
            return `
                <div class="parameter-item d-flex justify-content-between align-items-center py-2 border-bottom">
                    <span class="parameter-name fw-medium">${key}</span>
                    <span class="parameter-value badge bg-light text-dark">${this.formatParameterValue(value)}</span>
                </div>
            `;
        }).join('');

        container.innerHTML = parameterCards;
    }

    formatParameterValue(value) {
        if (typeof value === 'boolean') {
            return value ? 'true' : 'false';
        }
        if (typeof value === 'number') {
            return value.toString();
        }
        if (typeof value === 'string') {
            return value.length > 20 ? value.substring(0, 17) + '...' : value;
        }
        return JSON.stringify(value);
    }

    setLoading(loading) {
        this.isLoading = loading;
        const loadingElement = this.container.querySelector('.algorithm-selector__loading');
        const dropdownContainer = this.container.querySelector('.algorithm-selector__dropdown-container');
        
        if (loading) {
            loadingElement.style.display = 'block';
            dropdownContainer.style.display = 'none';
        } else {
            loadingElement.style.display = 'none';
        }
    }

    showDropdown() {
        const dropdownContainer = this.container.querySelector('.algorithm-selector__dropdown-container');
        dropdownContainer.style.display = 'block';
    }

    showError(message) {
        const errorContainer = this.container.querySelector('.algorithm-selector__error');
        const errorMessage = errorContainer.querySelector('.error-message');
        
        errorMessage.textContent = message;
        errorContainer.style.display = 'block';
    }

    hideError() {
        const errorContainer = this.container.querySelector('.algorithm-selector__error');
        errorContainer.style.display = 'none';
    }

    notifySelectionChange() {
        if (this.options.onSelectionChange) {
            this.options.onSelectionChange(this.selectedAlgorithm);
        }

        // Dispatch custom event
        const event = new CustomEvent('algorithm:selected', {
            detail: { algorithm: this.selectedAlgorithm }
        });
        this.container.dispatchEvent(event);
    }

    // Public API methods

    getSelectedAlgorithm() {
        return this.selectedAlgorithm;
    }

    setSelectedAlgorithm(algorithmName) {
        const dropdown = this.container.querySelector('.algorithm-selector__dropdown');
        dropdown.value = algorithmName;
        this.selectAlgorithm(algorithmName);
    }

    getAlgorithms() {
        return [...this.algorithms];
    }

    refresh() {
        return this.loadAlgorithms();
    }

    disable() {
        const dropdown = this.container.querySelector('.algorithm-selector__dropdown');
        dropdown.disabled = true;
    }

    enable() {
        const dropdown = this.container.querySelector('.algorithm-selector__dropdown');
        dropdown.disabled = false;
    }

    destroy() {
        this.container.innerHTML = '';
    }
}

// Auto-initialize algorithm selectors with data-algorithm-selector attribute
document.addEventListener('DOMContentLoaded', () => {
    const selectors = document.querySelectorAll('[data-algorithm-selector]');
    
    selectors.forEach(element => {
        const options = JSON.parse(element.dataset.algorithmSelectorOptions || '{}');
        new AlgorithmSelector(element, options);
    });
});

// Export for manual instantiation
if (typeof window !== 'undefined') {
    window.AlgorithmSelector = AlgorithmSelector;
}
