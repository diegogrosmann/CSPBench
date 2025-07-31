/**
 * Dataset Selector Component
 * 
 * Reusable component for dataset selection and configuration
 */

class DatasetSelector {
    constructor(container, options = {}) {
        this.container = typeof container === 'string' ? document.querySelector(container) : container;
        this.options = {
            allowUpload: true,
            allowGenerate: true,
            allowSamples: true,
            onConfig: () => {},
            ...options
        };
        
        this.config = null;
        this.init();
    }

    init() {
        this.render();
        this.setupEventListeners();
    }

    render() {
        this.container.innerHTML = `
            <div class="dataset-selector">
                <div class="dataset-source-tabs">
                    ${this.options.allowUpload ? '<button class="tab-btn active" data-source="upload">📁 Upload File</button>' : ''}
                    ${this.options.allowSamples ? '<button class="tab-btn" data-source="samples">📊 Sample Datasets</button>' : ''}
                    ${this.options.allowGenerate ? '<button class="tab-btn" data-source="generate">🧬 Generate Dataset</button>' : ''}
                </div>
                
                <div class="dataset-source-content">
                    ${this.renderUploadTab()}
                    ${this.renderSamplesTab()}
                    ${this.renderGenerateTab()}
                </div>
                
                <div class="dataset-preview" style="display: none;">
                    <h4>Dataset Preview</h4>
                    <div class="preview-stats">
                        <span class="stat">Sequences: <span id="preview-sequences">-</span></span>
                        <span class="stat">Length: <span id="preview-length">-</span></span>
                        <span class="stat">Alphabet: <span id="preview-alphabet">-</span></span>
                    </div>
                    <div class="preview-content">
                        <textarea id="preview-text" readonly rows="6"></textarea>
                    </div>
                </div>
            </div>
        `;
    }

    renderUploadTab() {
        if (!this.options.allowUpload) return '';
        
        return `
            <div class="source-tab active" data-source="upload">
                <div class="upload-area">
                    <div class="upload-drop-zone" id="drop-zone">
                        <div class="upload-icon">📁</div>
                        <p>Drag and drop a FASTA file here, or click to browse</p>
                        <input type="file" id="file-input" accept=".fasta,.fa,.txt" style="display: none;">
                        <button type="button" class="btn btn-outline" onclick="document.getElementById('file-input').click()">
                            Choose File
                        </button>
                    </div>
                    
                    <div class="upload-text-option">
                        <h4>Or paste FASTA content:</h4>
                        <textarea id="dataset-text" placeholder=">sequence1&#10;ACGTACGTACGT&#10;>sequence2&#10;ACGTACGTACGG" rows="8"></textarea>
                    </div>
                </div>
            </div>
        `;
    }

    renderSamplesTab() {
        if (!this.options.allowSamples) return '';
        
        return `
            <div class="source-tab" data-source="samples">
                <div class="samples-grid" id="samples-grid">
                    <div class="loading-placeholder">Loading sample datasets...</div>
                </div>
            </div>
        `;
    }

    renderGenerateTab() {
        if (!this.options.allowGenerate) return '';
        
        return `
            <div class="source-tab" data-source="generate">
                <div class="generate-form">
                    <div class="form-group">
                        <label for="gen-sequences">Number of sequences:</label>
                        <input type="number" id="gen-sequences" value="10" min="2" max="100">
                    </div>
                    <div class="form-group">
                        <label for="gen-length">Sequence length:</label>
                        <input type="number" id="gen-length" value="20" min="5" max="500">
                    </div>
                    <div class="form-group">
                        <label for="gen-alphabet">Alphabet:</label>
                        <select id="gen-alphabet">
                            <option value="DNA">DNA (A, C, G, T)</option>
                            <option value="RNA">RNA (A, C, G, U)</option>
                            <option value="PROTEIN">Protein (20 amino acids)</option>
                            <option value="BINARY">Binary (0, 1)</option>
                        </select>
                    </div>
                    <div class="form-group">
                        <label for="gen-diversity">Diversity level:</label>
                        <select id="gen-diversity">
                            <option value="low">Low (similar sequences)</option>
                            <option value="medium">Medium (moderate variation)</option>
                            <option value="high">High (diverse sequences)</option>
                        </select>
                    </div>
                    <button type="button" class="btn btn-primary" onclick="this.generateDataset()" id="generate-btn">
                        Generate Dataset
                    </button>
                </div>
            </div>
        `;
    }

    setupEventListeners() {
        // Tab switching
        this.container.querySelectorAll('.tab-btn').forEach(btn => {
            btn.addEventListener('click', (e) => {
                this.switchTab(e.target.dataset.source);
            });
        });

        // File upload
        const fileInput = this.container.querySelector('#file-input');
        if (fileInput) {
            fileInput.addEventListener('change', (e) => {
                this.handleFileUpload(e.target.files[0]);
            });
        }

        // Drag and drop
        const dropZone = this.container.querySelector('#drop-zone');
        if (dropZone) {
            dropZone.addEventListener('dragover', (e) => {
                e.preventDefault();
                dropZone.classList.add('drag-over');
            });

            dropZone.addEventListener('dragleave', () => {
                dropZone.classList.remove('drag-over');
            });

            dropZone.addEventListener('drop', (e) => {
                e.preventDefault();
                dropZone.classList.remove('drag-over');
                this.handleFileUpload(e.dataTransfer.files[0]);
            });
        }

        // Text input
        const textArea = this.container.querySelector('#dataset-text');
        if (textArea) {
            textArea.addEventListener('input', () => {
                this.handleTextInput(textArea.value);
            });
        }

        // Load sample datasets
        if (this.options.allowSamples) {
            this.loadSampleDatasets();
        }
    }

    switchTab(source) {
        // Update tab buttons
        this.container.querySelectorAll('.tab-btn').forEach(btn => {
            btn.classList.toggle('active', btn.dataset.source === source);
        });

        // Update tab content
        this.container.querySelectorAll('.source-tab').forEach(tab => {
            tab.classList.toggle('active', tab.dataset.source === source);
        });
    }

    async handleFileUpload(file) {
        if (!file) return;

        try {
            const content = await this.readFile(file);
            this.processDatasetContent(content, file.name);
        } catch (error) {
            this.showError('Failed to read file: ' + error.message);
        }
    }

    handleTextInput(content) {
        if (content.trim()) {
            this.processDatasetContent(content, 'pasted_dataset.fasta');
        }
    }

    readFile(file) {
        return new Promise((resolve, reject) => {
            const reader = new FileReader();
            reader.onload = (e) => resolve(e.target.result);
            reader.onerror = (e) => reject(new Error('File reading failed'));
            reader.readAsText(file);
        });
    }

    processDatasetContent(content, name) {
        try {
            const sequences = this.parseFASTA(content);
            
            if (sequences.length === 0) {
                throw new Error('No valid sequences found');
            }

            this.config = {
                type: 'upload',
                content: content,
                name: name,
                sequences: sequences.length,
                length: sequences[0]?.length || 0,
                alphabet: this.inferAlphabet(sequences)
            };

            this.showPreview();
            this.options.onConfig(this.config);

        } catch (error) {
            this.showError('Invalid FASTA format: ' + error.message);
        }
    }

    parseFASTA(content) {
        const sequences = [];
        const lines = content.trim().split('\n');
        let currentSequence = '';

        for (const line of lines) {
            if (line.startsWith('>')) {
                if (currentSequence) {
                    sequences.push(currentSequence);
                    currentSequence = '';
                }
            } else {
                currentSequence += line.trim();
            }
        }

        if (currentSequence) {
            sequences.push(currentSequence);
        }

        return sequences;
    }

    inferAlphabet(sequences) {
        const chars = new Set();
        sequences.forEach(seq => {
            for (const char of seq.toUpperCase()) {
                chars.add(char);
            }
        });
        return Array.from(chars).sort().join('');
    }

    showPreview() {
        const preview = this.container.querySelector('.dataset-preview');
        preview.style.display = 'block';

        document.getElementById('preview-sequences').textContent = this.config.sequences;
        document.getElementById('preview-length').textContent = this.config.length;
        document.getElementById('preview-alphabet').textContent = this.config.alphabet;
        
        // Show first few lines
        const lines = this.config.content.split('\n').slice(0, 10);
        document.getElementById('preview-text').value = lines.join('\n') + 
            (this.config.content.split('\n').length > 10 ? '\n...' : '');
    }

    async loadSampleDatasets() {
        const samplesGrid = this.container.querySelector('#samples-grid');
        if (!samplesGrid) return;

        try {
            // Mock sample datasets for now
            const samples = [
                {
                    name: 'Small DNA Example',
                    description: '5 sequences, 10bp each',
                    sequences: 5,
                    length: 10,
                    type: 'DNA'
                },
                {
                    name: 'Medium Protein Set',
                    description: '20 sequences, 50aa each',
                    sequences: 20,
                    length: 50,
                    type: 'Protein'
                },
                {
                    name: 'Large DNA Dataset',
                    description: '100 sequences, 200bp each',
                    sequences: 100,
                    length: 200,
                    type: 'DNA'
                }
            ];

            samplesGrid.innerHTML = samples.map(sample => `
                <div class="sample-card" data-sample="${sample.name}">
                    <h4>${sample.name}</h4>
                    <p>${sample.description}</p>
                    <div class="sample-stats">
                        <span>Sequences: ${sample.sequences}</span>
                        <span>Length: ${sample.length}</span>
                        <span>Type: ${sample.type}</span>
                    </div>
                    <button class="btn btn-outline" onclick="datasetSelector.selectSample('${sample.name}')">
                        Use This Dataset
                    </button>
                </div>
            `).join('');

        } catch (error) {
            samplesGrid.innerHTML = '<div class="error">Failed to load sample datasets</div>';
        }
    }

    selectSample(sampleName) {
        // Mock sample data - in real implementation, this would fetch from API
        const sampleContent = `>seq1\nACGTACGTACGT\n>seq2\nACGTACGTACGG\n>seq3\nACGTACGTACGA`;
        
        this.config = {
            type: 'sample',
            content: sampleContent,
            name: `${sampleName}.fasta`,
            sequences: 3,
            length: 12,
            alphabet: 'ACGT'
        };

        this.showPreview();
        this.options.onConfig(this.config);
    }

    showError(message) {
        const errorDiv = document.createElement('div');
        errorDiv.className = 'error-message';
        errorDiv.textContent = message;
        
        // Remove any existing error messages
        this.container.querySelectorAll('.error-message').forEach(el => el.remove());
        
        // Add new error message
        this.container.appendChild(errorDiv);
        
        // Auto-remove after 5 seconds
        setTimeout(() => errorDiv.remove(), 5000);
    }

    getConfig() {
        return this.config;
    }

    reset() {
        this.config = null;
        this.container.querySelector('.dataset-preview').style.display = 'none';
        this.container.querySelectorAll('input, textarea').forEach(input => {
            input.value = '';
        });
    }
}

// Export para uso global
window.DatasetSelector = DatasetSelector;
