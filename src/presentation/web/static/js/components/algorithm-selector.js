/**
 * Algorithm Selector Component
 * 
 * Reusable component for algorithm selection with detailed information
 */

class AlgorithmSelector {
    constructor(container, options = {}) {
        this.container = typeof container === 'string' ? document.querySelector(container) : container;
        
        if (!this.container) {
            console.error('AlgorithmSelector: Container not found:', container);
            return;
        }
        
        this.options = {
            allowMultiple: false,
            showDetails: true,
            onSelect: () => {},
            ...options
        };
        
        this.algorithms = [];
        this.selectedAlgorithm = null;
        this.selectedAlgorithms = [];
        
        this.init();
    }

    async init() {
        await this.loadAlgorithms();
        this.render();
        this.setupEventListeners();
    }

    async loadAlgorithms() {
        try {
            this.algorithms = await window.apiClient.getAlgorithms();
        } catch (error) {
            console.error('Failed to load algorithms:', error);
            this.algorithms = [];
        }
    }

    render() {
        if (this.algorithms.length === 0) {
            this.container.innerHTML = `
                <div class="loading-state">
                    <div class="spinner"></div>
                    <p>Loading algorithms...</p>
                </div>
            `;
            return;
        }

        this.container.innerHTML = `
            <div class="algorithm-selector-modern">
                <div class="selector-header">
                    <h3 class="selector-title">
                        <span class="title-icon">🧠</span>
                        Choose Your Algorithm
                    </h3>
                    <p class="selector-subtitle">Select a CSP algorithm to solve your problem</p>
                </div>
                
                <div class="algorithm-search-bar">
                    <div class="search-input-wrapper">
                        <span class="search-icon">🔍</span>
                        <input type="text" id="algorithm-search" class="search-input" 
                               placeholder="Search algorithms by name or features..." />
                    </div>
                    <div class="algorithm-filters">
                        <button class="filter-chip" data-filter="all" data-active="true">
                            All Algorithms
                        </button>
                        <button class="filter-chip" data-filter="deterministic">
                            🎯 Deterministic
                        </button>
                        <button class="filter-chip" data-filter="parallel">
                            ⚡ Parallel
                        </button>
                        <button class="filter-chip" data-filter="fast">
                            🚀 Fast
                        </button>
                    </div>
                </div>
                
                <div class="algorithm-grid">
                    ${this.renderModernAlgorithmCards()}
                </div>
                
                ${this.options.showDetails ? this.renderModernAlgorithmDetails() : ''}
            </div>
        `;
    }

    renderModernAlgorithmCards() {
        return this.algorithms.map(algorithm => {
            const isSelected = this.isSelected(algorithm);
            const algorithmType = this.getAlgorithmType(algorithm);
            const complexity = this.getAlgorithmComplexity(algorithm);
            
            return `
                <div class="algorithm-card-modern ${isSelected ? 'selected' : ''}" 
                     data-algorithm="${algorithm.name}"
                     data-deterministic="${algorithm.is_deterministic}"
                     data-parallel="${algorithm.supports_internal_parallel}">
                    
                    <div class="card-header">
                        <div class="algorithm-icon">${this.getAlgorithmIcon(algorithm)}</div>
                        <div class="card-title-section">
                            <h4 class="algorithm-name">${algorithm.name}</h4>
                            <span class="algorithm-type">${algorithmType}</span>
                        </div>
                        <div class="selection-indicator">
                            <div class="checkmark">✓</div>
                        </div>
                    </div>
                    
                    <div class="card-body">
                        <div class="algorithm-summary">
                            ${this.getAlgorithmSummary(algorithm.description)}
                        </div>
                        
                        <div class="algorithm-features">
                            <div class="feature-tags">
                                <span class="feature-tag ${algorithm.is_deterministic ? 'active' : 'inactive'}">
                                    <span class="tag-icon">🎯</span>
                                    Deterministic
                                </span>
                                <span class="feature-tag ${algorithm.supports_internal_parallel ? 'active' : 'inactive'}">
                                    <span class="tag-icon">⚡</span>
                                    Parallel
                                </span>
                                <span class="feature-tag complexity-${complexity}">
                                    <span class="tag-icon">⚙️</span>
                                    ${complexity}
                                </span>
                            </div>
                        </div>
                        
                        <div class="algorithm-metrics">
                            <div class="metric">
                                <span class="metric-label">Parameters</span>
                                <span class="metric-value">${Object.keys(algorithm.default_params).length}</span>
                            </div>
                            <div class="metric">
                                <span class="metric-label">Type</span>
                                <span class="metric-value">${algorithm.is_deterministic ? 'Deterministic' : 'Stochastic'}</span>
                            </div>
                        </div>
                    </div>
                    
                    <div class="card-footer">
                        <button type="button" class="select-button" data-algorithm="${algorithm.name}">
                            <span class="button-icon">
                                ${isSelected ? '✓' : '+'}
                            </span>
                            <span class="button-text">
                                ${isSelected ? 'Selected' : (this.options.allowMultiple ? 'Add Algorithm' : 'Select Algorithm')}
                            </span>
                        </button>
                        ${this.options.showDetails ? `
                            <button type="button" class="details-button" data-algorithm="${algorithm.name}">
                                <span class="button-icon">ℹ️</span>
                                Details
                            </button>
                        ` : ''}
                    </div>
                </div>
            `;
        }).join('');
    }

    // Funções auxiliares para renderização moderna
    getAlgorithmIcon(algorithm) {
        const iconMap = {
            'Baseline': '📏',
            'BLF-GA': '🧬',
            'CSC': '🎯',
            'DP-CSP': '🏗️',
            'H³-CSP': '🔺'
        };
        return iconMap[algorithm.name] || '🔬';
    }

    getAlgorithmType(algorithm) {
        if (algorithm.name.includes('GA') || algorithm.name.includes('Genetic')) return 'Genetic Algorithm';
        if (algorithm.name.includes('DP')) return 'Dynamic Programming';
        if (algorithm.name.includes('CSC')) return 'Clustering';
        if (algorithm.name.includes('H³')) return 'Hierarchical Search';
        if (algorithm.name === 'Baseline') return 'Greedy Heuristic';
        return 'Optimization Algorithm';
    }

    getAlgorithmComplexity(algorithm) {
        const paramCount = Object.keys(algorithm.default_params).length;
        if (paramCount <= 3) return 'Simple';
        if (paramCount <= 8) return 'Moderate';
        return 'Advanced';
    }

    getAlgorithmSummary(description) {
        // Extrai a primeira frase ou até 120 caracteres
        const firstSentence = description.split('.')[0];
        if (firstSentence.length > 120) {
            return firstSentence.substring(0, 120) + '...';
        }
        return firstSentence + '.';
    }

    renderModernAlgorithmDetails() {
        return `
            <div class="algorithm-details-modern" id="algorithm-details" style="display: none;">
                <div class="details-overlay" onclick="window.algorithmSelector.hideDetails()"></div>
                <div class="details-modal">
                    <div class="details-header">
                        <div class="details-title-section">
                            <span class="details-icon" id="details-icon">🔬</span>
                            <h3 class="details-title" id="details-title">Algorithm Details</h3>
                        </div>
                        <button type="button" class="close-button" onclick="window.algorithmSelector.hideDetails()">
                            <span>✕</span>
                        </button>
                    </div>
                    
                    <div class="details-content">
                        <div class="details-section">
                            <h4 class="section-title">
                                <span class="section-icon">📋</span>
                                Description
                            </h4>
                            <p class="details-description" id="details-description"></p>
                        </div>
                        
                        <div class="details-section">
                            <h4 class="section-title">
                                <span class="section-icon">⚙️</span>
                                Properties
                            </h4>
                            <div class="properties-grid" id="details-properties"></div>
                        </div>
                        
                        <div class="details-section">
                            <h4 class="section-title">
                                <span class="section-icon">🎛️</span>
                                Parameters
                            </h4>
                            <div class="parameters-list" id="details-parameters"></div>
                        </div>
                    </div>
                </div>
            </div>
        `;
    }

    setupEventListeners() {
        // Search functionality
        const searchInput = this.container.querySelector('#algorithm-search');
        if (searchInput) {
            searchInput.addEventListener('input', (e) => {
                this.filterAlgorithms(e.target.value);
            });
        }

        // Filter chips functionality
        this.container.querySelectorAll('.filter-chip').forEach(chip => {
            chip.addEventListener('click', (e) => {
                // Remove active from all chips
                this.container.querySelectorAll('.filter-chip').forEach(c => c.removeAttribute('data-active'));
                // Add active to clicked chip
                e.target.setAttribute('data-active', 'true');
                
                const filter = e.target.dataset.filter;
                this.applyChipFilter(filter);
            });
        });

        // Card selection (click anywhere on card)
        this.container.querySelectorAll('.algorithm-card-modern').forEach(card => {
            card.addEventListener('click', (e) => {
                // Don't trigger if clicking on details button
                if (e.target.closest('.details-button')) return;
                
                const algorithmName = card.dataset.algorithm;
                this.selectAlgorithm(algorithmName);
            });
        });

        // Select button clicks
        this.container.querySelectorAll('.select-button').forEach(btn => {
            btn.addEventListener('click', (e) => {
                e.stopPropagation(); // Prevent card click
                const algorithmName = e.target.closest('[data-algorithm]').dataset.algorithm;
                this.selectAlgorithm(algorithmName);
            });
        });

        // Details button clicks
        this.container.querySelectorAll('.details-button').forEach(btn => {
            btn.addEventListener('click', (e) => {
                e.stopPropagation(); // Prevent card click
                const algorithmName = e.target.closest('[data-algorithm]').dataset.algorithm;
                this.showDetails(algorithmName);
            });
        });
    }

    selectAlgorithm(algorithmName) {
        const algorithm = this.algorithms.find(a => a.name === algorithmName);
        if (!algorithm) return;

        if (this.options.allowMultiple) {
            const index = this.selectedAlgorithms.findIndex(a => a.name === algorithmName);
            if (index === -1) {
                this.selectedAlgorithms.push(algorithm);
            } else {
                this.selectedAlgorithms.splice(index, 1);
            }
            this.updateMultipleSelection();
        } else {
            this.selectedAlgorithm = algorithm;
            this.updateSingleSelection();
        }

        this.options.onSelect(this.options.allowMultiple ? this.selectedAlgorithms : this.selectedAlgorithm);
    }

    updateSingleSelection() {
        this.container.querySelectorAll('.algorithm-card-modern').forEach(card => {
            const algorithmName = card.dataset.algorithm;
            const isSelected = this.selectedAlgorithm?.name === algorithmName;
            
            card.classList.toggle('selected', isSelected);
            
            const btn = card.querySelector('.select-button');
            if (btn) {
                const buttonIcon = btn.querySelector('.button-icon');
                const buttonText = btn.querySelector('.button-text');
                
                if (isSelected) {
                    buttonIcon.textContent = '✓';
                    buttonText.textContent = 'Selected';
                    btn.classList.add('selected');
                } else {
                    buttonIcon.textContent = '+';
                    buttonText.textContent = 'Select Algorithm';
                    btn.classList.remove('selected');
                }
            }
        });
    }

    updateMultipleSelection() {
        this.container.querySelectorAll('.algorithm-card-modern').forEach(card => {
            const algorithmName = card.dataset.algorithm;
            const isSelected = this.selectedAlgorithms.some(a => a.name === algorithmName);
            
            card.classList.toggle('selected', isSelected);
            
            const btn = card.querySelector('.select-button');
            if (btn) {
                const buttonIcon = btn.querySelector('.button-icon');
                const buttonText = btn.querySelector('.button-text');
                
                if (isSelected) {
                    buttonIcon.textContent = '✓';
                    buttonText.textContent = 'Selected';
                    btn.classList.add('selected');
                } else {
                    buttonIcon.textContent = '+';
                    buttonText.textContent = 'Add Algorithm';
                    btn.classList.remove('selected');
                }
            }
        });
    }

    applyChipFilter(filter) {
        const cards = this.container.querySelectorAll('.algorithm-card-modern');
        
        cards.forEach(card => {
            let show = true;
            
            switch (filter) {
                case 'all':
                    show = true;
                    break;
                case 'deterministic':
                    show = card.dataset.deterministic === 'true';
                    break;
                case 'parallel':
                    show = card.dataset.parallel === 'true';
                    break;
                case 'fast':
                    // Consider simple algorithms as fast
                    const algorithmName = card.dataset.algorithm;
                    show = ['Baseline', 'CSC'].includes(algorithmName);
                    break;
            }
            
            card.style.display = show ? 'block' : 'none';
        });
    }

    isSelected(algorithm) {
        if (this.options.allowMultiple) {
            return this.selectedAlgorithms.some(a => a.name === algorithm.name);
        } else {
            return this.selectedAlgorithm?.name === algorithm.name;
        }
    }

    showDetails(algorithmName) {
        const algorithm = this.algorithms.find(a => a.name === algorithmName);
        if (!algorithm) return;

        // Update modal content
        const detailsIcon = document.getElementById('details-icon');
        const detailsTitle = document.getElementById('details-title');
        const detailsDescription = document.getElementById('details-description');
        
        if (detailsIcon) detailsIcon.textContent = this.getAlgorithmIcon(algorithm);
        if (detailsTitle) detailsTitle.textContent = `${algorithm.name}`;
        if (detailsDescription) detailsDescription.textContent = algorithm.description;
        
        // Properties Grid
        const propertiesGrid = document.getElementById('details-properties');
        if (propertiesGrid) {
            propertiesGrid.innerHTML = `
                <div class="property-item">
                    <div class="property-label">Type</div>
                    <div class="property-value">${algorithm.is_deterministic ? 'Deterministic' : 'Stochastic'}</div>
                </div>
                <div class="property-item">
                    <div class="property-label">Parallel Support</div>
                    <div class="property-value">${algorithm.supports_internal_parallel ? 'Yes' : 'No'}</div>
                </div>
                <div class="property-item">
                    <div class="property-label">Algorithm Class</div>
                    <div class="property-value">${this.getAlgorithmType(algorithm)}</div>
                </div>
                <div class="property-item">
                    <div class="property-label">Complexity</div>
                    <div class="property-value">${this.getAlgorithmComplexity(algorithm)}</div>
                </div>
            `;
        }
        
        // Parameters List
        const parametersList = document.getElementById('details-parameters');
        if (parametersList) {
            const params = algorithm.default_params;
            
            if (Object.keys(params).length === 0) {
                parametersList.innerHTML = '<p style="color: var(--text-muted); text-align: center;">No configurable parameters</p>';
            } else {
                parametersList.innerHTML = Object.entries(params).map(([key, value]) => `
                    <div class="parameter-item">
                        <span class="parameter-name">${this.formatParameterName(key)}</span>
                        <span class="parameter-value">${this.formatParameterValue(value)}</span>
                    </div>
                `).join('');
            }
        }

        // Show modal
        const modal = document.getElementById('algorithm-details');
        if (modal) modal.style.display = 'flex';
    }

    hideDetails() {
        const modal = document.getElementById('algorithm-details');
        if (modal) modal.style.display = 'none';
    }

    renderParametersList(params) {
        return `
            <div class="parameters-list">
                ${Object.entries(params).map(([key, value]) => `
                    <div class="parameter-item">
                        <span class="param-name">${this.formatParameterName(key)}</span>
                        <span class="param-value">${this.formatParameterValue(value)}</span>
                        <span class="param-type">${this.getParameterType(value)}</span>
                    </div>
                `).join('')}
            </div>
        `;
    }

    formatParameterName(name) {
        return name.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase());
    }

    formatParameterValue(value) {
        if (value === null) return 'null';
        if (typeof value === 'boolean') return value ? 'true' : 'false';
        if (typeof value === 'string') return `"${value}"`;
        return String(value);
    }

    getParameterType(value) {
        if (value === null) return 'nullable';
        return typeof value;
    }

    truncateDescription(description) {
        const maxLength = 150;
        if (description.length <= maxLength) return description;
        return description.substring(0, maxLength) + '...';
    }

    filterAlgorithms(searchTerm) {
        const term = searchTerm.toLowerCase();
        this.container.querySelectorAll('.algorithm-card').forEach(card => {
            const algorithmName = card.dataset.algorithm.toLowerCase();
            const description = card.querySelector('.algorithm-description').textContent.toLowerCase();
            const matches = algorithmName.includes(term) || description.includes(term);
            card.style.display = matches ? 'block' : 'none';
        });
    }

    applyFilters() {
        const filters = Array.from(this.container.querySelectorAll('.algorithm-filters input:checked'))
            .map(input => input.value);

        this.container.querySelectorAll('.algorithm-card').forEach(card => {
            const algorithmName = card.dataset.algorithm;
            const algorithm = this.algorithms.find(a => a.name === algorithmName);
            
            let shouldShow = true;
            
            if (filters.includes('deterministic') && !algorithm.is_deterministic) {
                shouldShow = false;
            }
            
            if (filters.includes('parallel') && !algorithm.supports_internal_parallel) {
                shouldShow = false;
            }
            
            card.style.display = shouldShow ? 'block' : 'none';
        });
    }

    getSelectedAlgorithm() {
        return this.selectedAlgorithm;
    }

    getSelectedAlgorithms() {
        return this.selectedAlgorithms;
    }

    reset() {
        this.selectedAlgorithm = null;
        this.selectedAlgorithms = [];
        this.updateSingleSelection();
        this.updateMultipleSelection();
        this.hideDetails();
    }
}

// Export para uso global
window.AlgorithmSelector = AlgorithmSelector;
