/**
 * Dataset Generator Component
 * 
 * Handles synthetic dataset generation and NCBI dataset downloading
 * with real-time progress tracking and result management.
 */

class DatasetGenerator {
    constructor() {
        this.apiClient = window.apiClient || new APIClient();
        this.currentSession = null;
        this.currentType = 'synthetic';
        
        this.initializeComponents();
        this.bindEvents();
        this.loadSavedDatasets();
        
        console.log('Dataset Generator initialized');
    }
    
    initializeComponents() {
        // Get DOM elements
        this.elements = {
            // Type selector
            syntheticTab: document.getElementById('synthetic-tab'),
            ncbiTab: document.getElementById('ncbi-tab'),
            syntheticGenerator: document.getElementById('synthetic-generator'),
            ncbiGenerator: document.getElementById('ncbi-generator'),
            
            // Forms
            syntheticForm: document.getElementById('synthetic-form'),
            ncbiForm: document.getElementById('ncbi-form'),
            
            // Synthetic form elements
            syntheticAlphabet: document.getElementById('synthetic-alphabet'),
            customAlphabet: document.getElementById('custom-alphabet'),
            syntheticMethod: document.getElementById('synthetic-method'),
            methodParams: document.getElementById('method-params'),
            
            // Status and results
            statusSection: document.getElementById('generation-status'),
            resultsSection: document.getElementById('results-section'),
            progressFill: document.getElementById('progress-fill'),
            statusMessage: document.getElementById('status-message'),
            statusDetails: document.getElementById('status-details'),
            cancelButton: document.getElementById('cancel-generation'),
            
            // Results elements
            datasetDetails: document.getElementById('dataset-details'),
            sequencePreview: document.getElementById('sequence-preview'),
            downloadButton: document.getElementById('download-dataset'),
            useButton: document.getElementById('use-in-execution'),
            generateAnother: document.getElementById('generate-another'),
            
            // Saved datasets
            savedDatasetsList: document.getElementById('saved-datasets-list'),
            refreshButton: document.getElementById('refresh-datasets')
        };
        
        this.setupMethodParameters();
    }
    
    bindEvents() {
        // Type switching
        this.elements.syntheticTab?.addEventListener('click', () => this.switchType('synthetic'));
        this.elements.ncbiTab?.addEventListener('click', () => this.switchType('ncbi'));
        
        // Form submissions
        this.elements.syntheticForm?.addEventListener('submit', (e) => this.handleSyntheticSubmit(e));
        this.elements.ncbiForm?.addEventListener('submit', (e) => this.handleNCBISubmit(e));
        
        // Alphabet selection
        this.elements.syntheticAlphabet?.addEventListener('change', (e) => this.handleAlphabetChange(e));
        
        // Method selection
        this.elements.syntheticMethod?.addEventListener('change', (e) => this.updateMethodParameters(e));
        
        // Cancel generation
        this.elements.cancelButton?.addEventListener('click', () => this.cancelGeneration());
        
        // Results actions
        this.elements.downloadButton?.addEventListener('click', () => this.downloadDataset());
        this.elements.useButton?.addEventListener('click', () => this.useInExecution());
        this.elements.generateAnother?.addEventListener('click', () => this.resetForm());
        
        // Saved datasets
        this.elements.refreshButton?.addEventListener('click', () => this.loadSavedDatasets());
        
        // Real-time validation
        this.setupFormValidation();
    }
    
    switchType(type) {
        this.currentType = type;
        
        // Update tabs
        this.elements.syntheticTab?.classList.toggle('active', type === 'synthetic');
        this.elements.ncbiTab?.classList.toggle('active', type === 'ncbi');
        
        // Update generators
        this.elements.syntheticGenerator?.classList.toggle('active', type === 'synthetic');
        this.elements.ncbiGenerator?.classList.toggle('active', type === 'ncbi');
        
        // Hide status and results when switching
        this.hideStatus();
        this.hideResults();
    }
    
    handleAlphabetChange(event) {
        const value = event.target.value;
        const customInput = this.elements.customAlphabet;
        
        if (value === 'custom') {
            customInput.style.display = 'block';
            customInput.required = true;
        } else {
            customInput.style.display = 'none';
            customInput.required = false;
            customInput.value = '';
        }
    }
    
    setupMethodParameters() {
        // Initialize with default method
        this.updateMethodParameters({ target: { value: 'random' } });
    }
    
    updateMethodParameters(event) {
        const method = event.target.value;
        const container = this.elements.methodParams;
        
        if (!container) return;
        
        container.innerHTML = this.getMethodParametersHTML(method);
        this.bindMethodParameterEvents();
    }
    
    getMethodParametersHTML(method) {
        switch (method) {
            case 'random':
                return `
                    <h4>Random Generation Parameters</h4>
                    <p>Generates completely random sequences using the selected alphabet.</p>
                `;
                
            case 'noise':
                return `
                    <h4>Center-based Generation Parameters</h4>
                    <div class="form-grid">
                        <div class="form-group">
                            <label for="center-sequence">Center Sequence</label>
                            <input type="text" id="center-sequence" name="center_sequence" 
                                   placeholder="Enter center sequence or leave empty for auto-generation" 
                                   pattern="[A-Za-z]+">
                            <small>Base sequence around which variations will be generated</small>
                        </div>
                        <div class="form-group">
                            <label for="noise-rate">Noise Rate</label>
                            <input type="number" id="noise-rate" name="noise_rate" 
                                   min="0" max="1" step="0.01" value="0.1" required>
                            <small>Probability of mutation per position (0.0-1.0)</small>
                        </div>
                    </div>
                `;
                
            case 'clustered':
                return `
                    <h4>Clustered Generation Parameters</h4>
                    <div class="form-grid">
                        <div class="form-group">
                            <label for="num-clusters">Number of Clusters</label>
                            <input type="number" id="num-clusters" name="num_clusters" 
                                   min="2" max="10" value="3" required>
                            <small>Number of sequence clusters to generate</small>
                        </div>
                        <div class="form-group">
                            <label for="cluster-noise">Cluster Noise Rate</label>
                            <input type="number" id="cluster-noise" name="cluster_noise" 
                                   min="0" max="1" step="0.01" value="0.1" required>
                            <small>Noise within each cluster (0.0-1.0)</small>
                        </div>
                    </div>
                `;
                
            case 'mutations':
                return `
                    <h4>Mutation-based Generation Parameters</h4>
                    <div class="form-grid">
                        <div class="form-group">
                            <label for="base-sequence">Base Sequence</label>
                            <input type="text" id="base-sequence" name="base_sequence" 
                                   placeholder="Enter base sequence or leave empty for auto-generation" 
                                   pattern="[A-Za-z]+">
                            <small>Starting sequence for mutations</small>
                        </div>
                        <div class="form-group">
                            <label for="mutation-rate">Mutation Rate</label>
                            <input type="number" id="mutation-rate" name="mutation_rate" 
                                   min="0" max="1" step="0.01" value="0.05" required>
                            <small>Probability of mutation per sequence (0.0-1.0)</small>
                        </div>
                    </div>
                    <div class="form-group">
                        <label for="mutation-types">Mutation Types</label>
                        <div class="checkbox-group">
                            <label class="checkbox-label">
                                <input type="checkbox" name="mutation_types" value="substitution" checked>
                                Substitution
                            </label>
                            <label class="checkbox-label">
                                <input type="checkbox" name="mutation_types" value="insertion">
                                Insertion
                            </label>
                            <label class="checkbox-label">
                                <input type="checkbox" name="mutation_types" value="deletion">
                                Deletion
                            </label>
                        </div>
                        <small>Types of mutations to apply</small>
                    </div>
                `;
                
            default:
                return '<p>Select a generation method to see parameters.</p>';
        }
    }
    
    bindMethodParameterEvents() {
        // Add any specific event listeners for method parameters
        // This could include validation, dynamic updates, etc.
    }
    
    setupFormValidation() {
        // Add real-time validation for both forms
        const forms = [this.elements.syntheticForm, this.elements.ncbiForm];
        
        forms.forEach(form => {
            if (!form) return;
            
            const inputs = form.querySelectorAll('input, select');
            inputs.forEach(input => {
                input.addEventListener('blur', () => this.validateField(input));
                input.addEventListener('input', () => this.clearFieldError(input));
            });
        });
    }
    
    validateField(field) {
        const formGroup = field.closest('.form-group');
        if (!formGroup) return;
        
        let isValid = field.checkValidity();
        let errorMessage = '';
        
        // Custom validation rules
        if (field.name === 'sequences' && field.value) {
            const value = parseInt(field.value);
            if (value < 3 || value > 1000) {
                isValid = false;
                errorMessage = 'Number of sequences must be between 3 and 1000';
            }
        }
        
        if (field.name === 'length' && field.value) {
            const value = parseInt(field.value);
            if (value < 5 || value > 10000) {
                isValid = false;
                errorMessage = 'Sequence length must be between 5 and 10000';
            }
        }
        
        if (field.name === 'custom_alphabet' && field.style.display !== 'none' && field.value) {
            if (field.value.length < 2) {
                isValid = false;
                errorMessage = 'Custom alphabet must have at least 2 characters';
            }
        }
        
        // Update visual feedback
        formGroup.classList.remove('error', 'success');
        const existingError = formGroup.querySelector('.error-message');
        if (existingError) existingError.remove();
        
        if (!isValid) {
            formGroup.classList.add('error');
            if (errorMessage) {
                const errorDiv = document.createElement('div');
                errorDiv.className = 'error-message';
                errorDiv.innerHTML = `⚠️ ${errorMessage}`;
                formGroup.appendChild(errorDiv);
            }
        } else if (field.value) {
            formGroup.classList.add('success');
        }
        
        return isValid;
    }
    
    clearFieldError(field) {
        const formGroup = field.closest('.form-group');
        if (!formGroup) return;
        
        formGroup.classList.remove('error');
        const errorMessage = formGroup.querySelector('.error-message');
        if (errorMessage) errorMessage.remove();
    }
    
    async handleSyntheticSubmit(event) {
        event.preventDefault();
        
        if (!this.validateForm(this.elements.syntheticForm)) {
            return;
        }
        
        const formData = new FormData(this.elements.syntheticForm);
        const params = this.extractSyntheticParams(formData);
        
        try {
            this.showStatus('Generating synthetic dataset...');
            const result = await this.apiClient.generateSyntheticDataset(params);
            await this.handleGenerationResult(result, 'synthetic');
        } catch (error) {
            this.showError('Failed to generate synthetic dataset', error);
        }
    }
    
    async handleNCBISubmit(event) {
        event.preventDefault();
        
        if (!this.validateForm(this.elements.ncbiForm)) {
            return;
        }
        
        const formData = new FormData(this.elements.ncbiForm);
        const params = this.extractNCBIParams(formData);
        
        try {
            this.showStatus('Downloading dataset from NCBI...');
            const result = await this.apiClient.downloadNCBIDataset(params);
            await this.handleGenerationResult(result, 'ncbi');
        } catch (error) {
            this.showError('Failed to download NCBI dataset', error);
        }
    }
    
    validateForm(form) {
        if (!form) return false;
        
        const inputs = form.querySelectorAll('input[required], select[required]');
        let isValid = true;
        
        inputs.forEach(input => {
            if (!this.validateField(input)) {
                isValid = false;
            }
        });
        
        return isValid;
    }
    
    extractSyntheticParams(formData) {
        const params = {
            n: parseInt(formData.get('sequences')),
            length: parseInt(formData.get('length')),
            alphabet: formData.get('alphabet') === 'custom' ? 
                     formData.get('custom_alphabet') : formData.get('alphabet'),
            method: formData.get('method'),
            filename: formData.get('filename') || null
        };
        
        // Add seed if provided
        const seed = formData.get('seed');
        if (seed) {
            params.seed = parseInt(seed);
        }
        
        // Add method-specific parameters
        this.addMethodParams(params, formData);
        
        return params;
    }
    
    extractNCBIParams(formData) {
        return {
            query: formData.get('query'),
            database: formData.get('database'),
            max_sequences: parseInt(formData.get('max_sequences')),
            min_length: parseInt(formData.get('min_length')),
            max_length: parseInt(formData.get('max_length')),
            filename: formData.get('filename') || null
        };
    }
    
    addMethodParams(params, formData) {
        switch (params.method) {
            case 'noise':
                const centerSequence = formData.get('center_sequence');
                if (centerSequence) params.center_sequence = centerSequence;
                params.noise_rate = parseFloat(formData.get('noise_rate'));
                break;
                
            case 'clustered':
                params.num_clusters = parseInt(formData.get('num_clusters'));
                params.cluster_noise = parseFloat(formData.get('cluster_noise'));
                break;
                
            case 'mutations':
                const baseSequence = formData.get('base_sequence');
                if (baseSequence) params.base_sequence = baseSequence;
                params.mutation_rate = parseFloat(formData.get('mutation_rate'));
                
                // Get checked mutation types
                const mutationTypes = [];
                const checkboxes = document.querySelectorAll('input[name="mutation_types"]:checked');
                checkboxes.forEach(cb => mutationTypes.push(cb.value));
                params.mutation_types = mutationTypes;
                break;
        }
    }
    
    showStatus(message, details = '') {
        this.elements.statusSection.style.display = 'block';
        this.elements.resultsSection.style.display = 'none';
        this.elements.statusMessage.textContent = message;
        this.elements.statusDetails.textContent = details;
        this.elements.progressFill.style.width = '10%';
        this.elements.cancelButton.style.display = 'block';
        
        // Scroll to status section
        this.elements.statusSection.scrollIntoView({ behavior: 'smooth' });
    }
    
    updateProgress(progress, message = '', details = '') {
        this.elements.progressFill.style.width = `${progress}%`;
        if (message) this.elements.statusMessage.textContent = message;
        if (details) this.elements.statusDetails.textContent = details;
    }
    
    hideStatus() {
        this.elements.statusSection.style.display = 'none';
    }
    
    async handleGenerationResult(result, type) {
        if (result.session_id) {
            this.currentSession = result.session_id;
            
            // Monitor progress
            await this.monitorGeneration(result.session_id);
        } else {
            // Immediate result
            this.showResults(result.result);
        }
    }
    
    async monitorGeneration(sessionId) {
        let attempts = 0;
        const maxAttempts = 60; // 5 minutes at 5-second intervals
        
        const checkStatus = async () => {
            try {
                const status = await this.apiClient.getGenerationStatus(sessionId);
                
                this.updateProgress(
                    status.progress || 50,
                    status.message || 'Processing...',
                    status.details || ''
                );
                
                if (status.status === 'completed' && status.result) {
                    this.showResults(status.result);
                    return;
                }
                
                if (status.status === 'failed') {
                    throw new Error(status.error || 'Generation failed');
                }
                
                attempts++;
                if (attempts < maxAttempts) {
                    setTimeout(checkStatus, 5000);
                } else {
                    throw new Error('Generation timed out');
                }
                
            } catch (error) {
                this.showError('Generation monitoring failed', error);
            }
        };
        
        setTimeout(checkStatus, 2000);
    }
    
    showResults(result) {
        this.hideStatus();
        
        this.elements.resultsSection.style.display = 'block';
        
        // Populate dataset details
        this.populateDatasetDetails(result.dataset_info);
        
        // Show sequence preview
        this.populateSequencePreview(result.sequences);
        
        // Set up download
        if (result.download_url) {
            this.elements.downloadButton.onclick = () => {
                window.open(result.download_url, '_blank');
            };
        }
        
        // Scroll to results
        this.elements.resultsSection.scrollIntoView({ behavior: 'smooth' });
        
        // Refresh saved datasets list
        this.loadSavedDatasets();
    }
    
    populateDatasetDetails(datasetInfo) {
        if (!this.elements.datasetDetails || !datasetInfo) return;
        
        const details = [
            { label: 'Filename', value: datasetInfo.filename },
            { label: 'Number of Sequences', value: datasetInfo.num_sequences },
            { label: 'Sequence Length', value: datasetInfo.sequence_length },
            { label: 'Alphabet', value: datasetInfo.alphabet },
            { label: 'Type', value: datasetInfo.type },
            { label: 'File Size', value: datasetInfo.file_size },
            { label: 'Generated At', value: new Date(datasetInfo.created_at).toLocaleString() }
        ];
        
        this.elements.datasetDetails.innerHTML = details.map(detail => `
            <div class="dataset-detail-item">
                <span class="dataset-detail-label">${detail.label}:</span>
                <span class="dataset-detail-value">${detail.value}</span>
            </div>
        `).join('');
    }
    
    populateSequencePreview(sequences) {
        if (!this.elements.sequencePreview || !sequences) return;
        
        const preview = sequences.slice(0, 5).map((seq, index) => `
            <div class="sequence-item">
                <div class="sequence-header">&gt;seq_${index + 1}</div>
                <div class="sequence-content">${seq}</div>
            </div>
        `).join('');
        
        this.elements.sequencePreview.innerHTML = preview;
    }
    
    showError(message, error) {
        this.hideStatus();
        
        console.error('Dataset Generator Error:', error);
        
        alert(`${message}\n\nError: ${error.message || error.toString()}`);
    }
    
    cancelGeneration() {
        if (this.currentSession) {
            this.apiClient.cancelGeneration(this.currentSession);
        }
        this.hideStatus();
        this.currentSession = null;
    }
    
    downloadDataset() {
        if (this.currentSession) {
            window.open(`/api/dataset/download/${this.currentSession}`, '_blank');
        }
    }
    
    useInExecution() {
        if (this.currentSession) {
            // Store dataset info for use in execution
            localStorage.setItem('selectedDataset', JSON.stringify({
                type: 'generated',
                session_id: this.currentSession
            }));
            
            // Navigate to single execution page
            window.location.href = '/execution/single';
        }
    }
    
    resetForm() {
        this.hideResults();
        this.currentSession = null;
        
        // Reset form values to defaults
        if (this.currentType === 'synthetic') {
            this.elements.syntheticForm.reset();
            document.getElementById('synthetic-sequences').value = 20;
            document.getElementById('synthetic-length').value = 50;
            document.getElementById('synthetic-alphabet').value = 'ACGT';
            document.getElementById('synthetic-method').value = 'random';
        } else {
            this.elements.ncbiForm.reset();
            document.getElementById('ncbi-query').value = 'COI[Gene]';
            document.getElementById('ncbi-database').value = 'nucleotide';
            document.getElementById('ncbi-max-sequences').value = 50;
            document.getElementById('ncbi-min-length').value = 100;
            document.getElementById('ncbi-max-length').value = 2000;
        }
        
        this.setupMethodParameters();
        this.clearFormValidation();
    }
    
    hideResults() {
        this.elements.resultsSection.style.display = 'none';
    }
    
    clearFormValidation() {
        const formGroups = document.querySelectorAll('.form-group');
        formGroups.forEach(group => {
            group.classList.remove('error', 'success');
            const errorMessage = group.querySelector('.error-message');
            if (errorMessage) errorMessage.remove();
        });
    }
    
    async loadSavedDatasets() {
        if (!this.elements.savedDatasetsList) return;
        
        try {
            this.elements.savedDatasetsList.innerHTML = '<div class="loading-message">Loading saved datasets...</div>';
            
            const datasets = await this.apiClient.getSavedDatasets();
            
            if (datasets.length === 0) {
                this.elements.savedDatasetsList.innerHTML = '<div class="loading-message">No saved datasets found.</div>';
                return;
            }
            
            this.elements.savedDatasetsList.innerHTML = datasets.map(dataset => 
                this.createDatasetItemHTML(dataset)
            ).join('');
            
            // Bind dataset action events
            this.bindDatasetActions();
            
        } catch (error) {
            console.error('Failed to load saved datasets:', error);
            this.elements.savedDatasetsList.innerHTML = '<div class="loading-message">Failed to load saved datasets.</div>';
        }
    }
    
    createDatasetItemHTML(dataset) {
        return `
            <div class="dataset-item" data-id="${dataset.id}">
                <div class="dataset-item-info">
                    <div class="dataset-item-name">${dataset.filename}</div>
                    <div class="dataset-item-details">
                        <span>📊 ${dataset.num_sequences} sequences</span>
                        <span>📏 ${dataset.sequence_length} bp</span>
                        <span>🧬 ${dataset.alphabet}</span>
                        <span>📅 ${new Date(dataset.created_at).toLocaleDateString()}</span>
                        <span>💾 ${dataset.file_size}</span>
                    </div>
                </div>
                <div class="dataset-item-actions">
                    <button class="dataset-action-button download" onclick="window.datasetGenerator.downloadSavedDataset('${dataset.id}')">
                        💾 Download
                    </button>
                    <button class="dataset-action-button" onclick="window.datasetGenerator.useSavedDataset('${dataset.id}')">
                        🚀 Use
                    </button>
                    <button class="dataset-action-button delete" onclick="window.datasetGenerator.deleteSavedDataset('${dataset.id}')">
                        🗑️ Delete
                    </button>
                </div>
            </div>
        `;
    }
    
    bindDatasetActions() {
        // Actions are bound via onclick attributes in HTML
        // This could be refactored to use proper event delegation
    }
    
    async downloadSavedDataset(datasetId) {
        try {
            const downloadUrl = await this.apiClient.getDatasetDownloadUrl(datasetId);
            window.open(downloadUrl, '_blank');
        } catch (error) {
            alert('Failed to download dataset: ' + error.message);
        }
    }
    
    async useSavedDataset(datasetId) {
        try {
            // Store dataset info for use in execution
            localStorage.setItem('selectedDataset', JSON.stringify({
                type: 'saved',
                dataset_id: datasetId
            }));
            
            // Navigate to single execution page
            window.location.href = '/execution/single';
        } catch (error) {
            alert('Failed to select dataset: ' + error.message);
        }
    }
    
    async deleteSavedDataset(datasetId) {
        if (!confirm('Are you sure you want to delete this dataset? This action cannot be undone.')) {
            return;
        }
        
        try {
            await this.apiClient.deleteDataset(datasetId);
            this.loadSavedDatasets(); // Refresh list
        } catch (error) {
            alert('Failed to delete dataset: ' + error.message);
        }
    }
}

// Global function for template onclick handlers
function switchDatasetType(type) {
    if (window.datasetGenerator) {
        window.datasetGenerator.switchType(type);
    }
}

// Make DatasetGenerator globally available
window.DatasetGenerator = DatasetGenerator;

// Auto-initialize when DOM is ready
document.addEventListener('DOMContentLoaded', function() {
    if (document.querySelector('.dataset-generator-container')) {
        window.datasetGenerator = new DatasetGenerator();
    }
});
