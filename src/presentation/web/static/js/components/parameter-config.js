/**
 * Parameter Configuration Component
 * 
 * Reusable component for configuring algorithm parameters
 */

class ParameterConfig {
    constructor(container, options = {}) {
        this.container = typeof container === 'string' ? document.querySelector(container) : container;
        this.options = {
            showDefaults: true,
            showDescriptions: true,
            validateInput: true,
            onUpdate: () => {},
            ...options
        };
        
        this.parameters = {};
        this.values = {};
        this.validators = {};
        
        this.init();
    }

    init() {
        this.container.classList.add('parameter-config');
    }

    loadParameters(algorithm) {
        this.algorithm = algorithm;
        this.parameters = algorithm.default_params || {};
        this.values = { ...this.parameters };
        this.setupValidators();
        this.render();
    }

    setupValidators() {
        this.validators = {
            population_size: (value) => {
                const num = parseInt(value);
                return num > 0 && num <= 10000 ? null : 'Must be between 1 and 10000';
            },
            max_generations: (value) => {
                const num = parseInt(value);
                return num > 0 && num <= 10000 ? null : 'Must be between 1 and 10000';
            },
            mutation_rate: (value) => {
                const num = parseFloat(value);
                return num >= 0 && num <= 1 ? null : 'Must be between 0 and 1';
            },
            crossover_rate: (value) => {
                const num = parseFloat(value);
                return num >= 0 && num <= 1 ? null : 'Must be between 0 and 1';
            },
            elite_size: (value) => {
                const num = parseInt(value);
                return num >= 0 && num <= 100 ? null : 'Must be between 0 and 100';
            },
            tournament_size: (value) => {
                const num = parseInt(value);
                return num > 0 && num <= 20 ? null : 'Must be between 1 and 20';
            },
            seed: (value) => {
                if (value === '' || value === null) return null;
                const num = parseInt(value);
                return !isNaN(num) && num >= 0 ? null : 'Must be a non-negative integer or empty';
            }
        };
    }

    render() {
        if (!this.algorithm || Object.keys(this.parameters).length === 0) {
            this.container.innerHTML = `
                <div class="no-parameters">
                    <p>This algorithm has no configurable parameters.</p>
                </div>
            `;
            return;
        }

        this.container.innerHTML = `
            <div class="parameter-config-header">
                <h3>Configure Parameters for ${this.algorithm.name}</h3>
                <div class="parameter-actions">
                    <button class="btn btn-outline" id="reset-params">Reset to Defaults</button>
                    <button class="btn btn-outline" id="load-preset">Load Preset</button>
                    <button class="btn btn-outline" id="save-preset">Save Preset</button>
                </div>
            </div>
            
            <div class="parameter-groups">
                ${this.renderParameterGroups()}
            </div>
            
            <div class="parameter-summary">
                <h4>Current Configuration</h4>
                <pre id="config-preview">${JSON.stringify(this.values, null, 2)}</pre>
            </div>
        `;

        this.setupEventListeners();
    }

    renderParameterGroups() {
        const groups = this.groupParameters();
        
        return Object.entries(groups).map(([groupName, params]) => `
            <div class="parameter-group">
                <h4 class="group-title">${this.formatGroupName(groupName)}</h4>
                <div class="parameter-fields">
                    ${params.map(param => this.renderParameterField(param)).join('')}
                </div>
            </div>
        `).join('');
    }

    groupParameters() {
        const groups = {
            basic: [],
            genetic: [],
            performance: [],
            advanced: []
        };

        Object.keys(this.parameters).forEach(key => {
            if (['population_size', 'max_generations', 'seed'].includes(key)) {
                groups.basic.push(key);
            } else if (['mutation_rate', 'crossover_rate', 'elite_size', 'tournament_size'].includes(key)) {
                groups.genetic.push(key);
            } else if (['parallel_processes', 'max_time', 'convergence_threshold'].includes(key)) {
                groups.performance.push(key);
            } else {
                groups.advanced.push(key);
            }
        });

        // Remove empty groups
        Object.keys(groups).forEach(key => {
            if (groups[key].length === 0) delete groups[key];
        });

        return groups;
    }

    renderParameterField(paramKey) {
        const value = this.values[paramKey];
        const defaultValue = this.parameters[paramKey];
        const type = this.getParameterType(value);
        
        return `
            <div class="parameter-field" data-param="${paramKey}">
                <label class="parameter-label">
                    ${this.formatParameterName(paramKey)}
                    ${this.options.showDefaults ? `<span class="default-value">(default: ${this.formatValue(defaultValue)})</span>` : ''}
                </label>
                
                ${this.renderInput(paramKey, value, type)}
                
                ${this.options.showDescriptions ? this.renderParameterDescription(paramKey) : ''}
                
                <div class="validation-message" id="validation-${paramKey}"></div>
            </div>
        `;
    }

    renderInput(paramKey, value, type) {
        switch (type) {
            case 'boolean':
                return `
                    <div class="checkbox-group">
                        <label class="checkbox-label">
                            <input type="checkbox" id="param-${paramKey}" ${value ? 'checked' : ''} />
                            <span class="checkbox-text">Enabled</span>
                        </label>
                    </div>
                `;
            
            case 'number':
                const step = this.isFloatParameter(paramKey) ? '0.01' : '1';
                const min = this.getParameterMin(paramKey);
                const max = this.getParameterMax(paramKey);
                
                return `
                    <div class="number-input-group">
                        <input type="number" 
                               id="param-${paramKey}" 
                               value="${value}" 
                               step="${step}"
                               ${min !== null ? `min="${min}"` : ''}
                               ${max !== null ? `max="${max}"` : ''}
                               class="parameter-input" />
                        ${this.renderSlider(paramKey, value, min, max)}
                    </div>
                `;
            
            case 'string':
                return `
                    <input type="text" 
                           id="param-${paramKey}" 
                           value="${value || ''}" 
                           class="parameter-input" />
                `;
            
            default:
                return `
                    <input type="text" 
                           id="param-${paramKey}" 
                           value="${JSON.stringify(value)}" 
                           class="parameter-input" />
                `;
        }
    }

    renderSlider(paramKey, value, min, max) {
        if (min === null || max === null) return '';
        
        // Only show slider for parameters with reasonable ranges
        const range = max - min;
        if (range > 1000) return '';
        
        return `
            <input type="range" 
                   class="parameter-slider" 
                   id="slider-${paramKey}"
                   min="${min}" 
                   max="${max}" 
                   value="${value}"
                   step="${this.isFloatParameter(paramKey) ? '0.01' : '1'}" />
        `;
    }

    renderParameterDescription(paramKey) {
        const descriptions = {
            population_size: 'Number of individuals in each generation',
            max_generations: 'Maximum number of generations to evolve',
            mutation_rate: 'Probability of mutation for each gene',
            crossover_rate: 'Probability of crossover between parents',
            elite_size: 'Number of best individuals to preserve',
            tournament_size: 'Number of individuals in tournament selection',
            seed: 'Random seed for reproducible results (leave empty for random)'
        };
        
        const description = descriptions[paramKey];
        return description ? `<small class="parameter-description">${description}</small>` : '';
    }

    setupEventListeners() {
        // Parameter inputs
        this.container.querySelectorAll('.parameter-input').forEach(input => {
            input.addEventListener('input', (e) => {
                this.updateParameter(e.target);
            });
            
            input.addEventListener('blur', (e) => {
                this.validateParameter(e.target);
            });
        });

        // Sliders
        this.container.querySelectorAll('.parameter-slider').forEach(slider => {
            slider.addEventListener('input', (e) => {
                const paramKey = e.target.id.replace('slider-', '');
                const input = document.getElementById(`param-${paramKey}`);
                input.value = e.target.value;
                this.updateParameter(input);
            });
        });

        // Checkboxes
        this.container.querySelectorAll('input[type="checkbox"]').forEach(checkbox => {
            checkbox.addEventListener('change', (e) => {
                this.updateParameter(e.target);
            });
        });

        // Action buttons
        document.getElementById('reset-params')?.addEventListener('click', () => {
            this.resetToDefaults();
        });

        document.getElementById('load-preset')?.addEventListener('click', () => {
            this.showPresetDialog();
        });

        document.getElementById('save-preset')?.addEventListener('click', () => {
            this.savePreset();
        });
    }

    updateParameter(input) {
        const paramKey = input.id.replace('param-', '');
        const type = this.getParameterType(this.parameters[paramKey]);
        
        let newValue;
        switch (type) {
            case 'boolean':
                newValue = input.checked;
                break;
            case 'number':
                newValue = this.isFloatParameter(paramKey) ? 
                    parseFloat(input.value) : parseInt(input.value);
                break;
            case 'string':
                newValue = input.value || null;
                break;
            default:
                try {
                    newValue = JSON.parse(input.value);
                } catch {
                    newValue = input.value;
                }
        }
        
        this.values[paramKey] = newValue;
        
        // Update slider if exists
        const slider = document.getElementById(`slider-${paramKey}`);
        if (slider) {
            slider.value = newValue;
        }
        
        this.updateConfigPreview();
        this.options.onUpdate(this.values);
    }

    validateParameter(input) {
        const paramKey = input.id.replace('param-', '');
        const validator = this.validators[paramKey];
        const messageEl = document.getElementById(`validation-${paramKey}`);
        
        if (validator) {
            const error = validator(input.value);
            if (error) {
                messageEl.textContent = error;
                messageEl.className = 'validation-message error';
                input.classList.add('invalid');
                return false;
            }
        }
        
        messageEl.textContent = '';
        messageEl.className = 'validation-message';
        input.classList.remove('invalid');
        return true;
    }

    updateConfigPreview() {
        const preview = document.getElementById('config-preview');
        if (preview) {
            preview.textContent = JSON.stringify(this.values, null, 2);
        }
    }

    resetToDefaults() {
        this.values = { ...this.parameters };
        this.render();
    }

    savePreset() {
        const name = prompt('Enter preset name:');
        if (name) {
            const presets = JSON.parse(localStorage.getItem('csp-presets') || '{}');
            presets[`${this.algorithm.name}-${name}`] = { ...this.values };
            localStorage.setItem('csp-presets', JSON.stringify(presets));
            alert('Preset saved successfully!');
        }
    }

    showPresetDialog() {
        const presets = JSON.parse(localStorage.getItem('csp-presets') || '{}');
        const algorithmPresets = Object.keys(presets)
            .filter(key => key.startsWith(`${this.algorithm.name}-`))
            .map(key => key.replace(`${this.algorithm.name}-`, ''));
        
        if (algorithmPresets.length === 0) {
            alert('No saved presets found for this algorithm.');
            return;
        }
        
        const selected = prompt(`Select preset:\n${algorithmPresets.join('\n')}`);
        if (selected && algorithmPresets.includes(selected)) {
            this.values = { ...presets[`${this.algorithm.name}-${selected}`] };
            this.render();
        }
    }

    // Utility methods
    getParameterType(value) {
        if (typeof value === 'boolean') return 'boolean';
        if (typeof value === 'number') return 'number';
        if (typeof value === 'string') return 'string';
        return 'object';
    }

    isFloatParameter(paramKey) {
        return ['mutation_rate', 'crossover_rate', 'convergence_threshold'].includes(paramKey);
    }

    getParameterMin(paramKey) {
        const mins = {
            population_size: 1,
            max_generations: 1,
            mutation_rate: 0,
            crossover_rate: 0,
            elite_size: 0,
            tournament_size: 1,
            seed: 0
        };
        return mins[paramKey] ?? null;
    }

    getParameterMax(paramKey) {
        const maxs = {
            population_size: 10000,
            max_generations: 10000,
            mutation_rate: 1,
            crossover_rate: 1,
            elite_size: 100,
            tournament_size: 20
        };
        return maxs[paramKey] ?? null;
    }

    formatParameterName(name) {
        return name.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase());
    }

    formatGroupName(name) {
        return name.replace(/\b\w/g, l => l.toUpperCase());
    }

    formatValue(value) {
        if (value === null) return 'null';
        if (typeof value === 'string') return `"${value}"`;
        return String(value);
    }

    getValues() {
        return { ...this.values };
    }

    setValues(values) {
        this.values = { ...values };
        this.render();
    }

    validate() {
        let isValid = true;
        this.container.querySelectorAll('.parameter-input').forEach(input => {
            if (!this.validateParameter(input)) {
                isValid = false;
            }
        });
        return isValid;
    }
}
