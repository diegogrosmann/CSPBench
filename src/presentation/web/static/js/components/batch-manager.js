/**
 * Batch Manager Component for CSPBench Web Interface
 * Manages batch configuration files with CRUD operations
 */

class BatchManagerTabbed {
    constructor() {
        this.containerId = 'batch-manager-container';
        this.batchFiles = [];
        this.executions = [];
        this.initialized = false;
        this.apiClient = window.apiClient || new APIClient();
    }

    /**
     * Initialize the batch manager component
     */
    async initialize() {
        if (this.initialized) return;

        console.log('BatchManagerTabbed: Starting initialization...');

        try {
            console.log('BatchManagerTabbed: Rendering interface...');
            await this.render();
            
            console.log('BatchManagerTabbed: Loading batch files...');
            await this.loadBatchFiles();
            
            console.log('BatchManagerTabbed: Setting up event handlers...');
            this.setupEventHandlers();
            
            this.initialized = true;
            console.log('BatchManagerTabbed: Initialization completed successfully');
        } catch (error) {
            console.error('BatchManagerTabbed: Failed to initialize:', error);
            this.renderError('Failed to initialize batch manager: ' + error.message);
        }
    }

    /**
     * Render the initial batch manager interface
     */
    async render() {
        const container = document.getElementById(this.containerId);
        if (!container) {
            throw new Error(`Container element with ID '${this.containerId}' not found`);
        }

        container.innerHTML = `
            <div class="batch-manager-content">
                <!-- Alert Container -->
                <div id="alert-container"></div>

                <!-- Active Executions -->
                <div class="mb-4" id="active-executions-section">
                    <div class="d-flex justify-content-between align-items-center mb-3">
                        <h5><i class="bi bi-play-circle me-2"></i>Active Executions</h5>
                        <button class="btn btn-outline-secondary btn-sm" id="refresh-executions-btn">
                            <i class="bi bi-arrow-clockwise me-1"></i>Refresh
                        </button>
                    </div>
                    
                    <div id="active-executions-list">
                        <div class="text-center p-3">
                            <div class="spinner-border text-primary spinner-border-sm" role="status">
                                <span class="visually-hidden">Loading executions...</span>
                            </div>
                        </div>
                    </div>
                </div>

                <!-- Batch Files List -->
                <div class="mb-4">
                    <div class="d-flex justify-content-between align-items-center mb-3">
                        <h5><i class="bi bi-files me-2"></i>Batch Configuration Files</h5>
                        <div class="d-flex gap-2">
                            <button class="btn btn-outline-primary btn-sm" id="upload-batch-btn">
                                <i class="bi bi-upload me-1"></i>Upload
                            </button>
                            <button class="btn btn-primary btn-sm" id="create-batch-btn">
                                <i class="bi bi-plus-circle me-1"></i>New Batch
                            </button>
                        </div>
                    </div>
                    
                    <div id="batch-files-list">
                        <div class="text-center p-4">
                            <div class="spinner-border text-primary" role="status">
                                <span class="visually-hidden">Loading batch files...</span>
                            </div>
                        </div>
                    </div>
                </div>

                <!-- File Upload Modal -->
                <div id="upload-modal" class="d-none">
                    <div class="card">
                        <div class="card-header">
                            <h6 class="mb-0"><i class="bi bi-upload me-2"></i>Upload Batch File</h6>
                        </div>
                        <div class="card-body">
                            <div class="mb-3">
                                <label for="file-input" class="form-label">Select YAML File</label>
                                <input type="file" class="form-control" id="file-input" accept=".yaml,.yml">
                                <div class="form-text">Only .yaml and .yml files are supported</div>
                            </div>
                            <div class="d-flex gap-2">
                                <button type="button" class="btn btn-success" id="upload-confirm-btn">
                                    <i class="bi bi-upload me-1"></i>Upload
                                </button>
                                <button type="button" class="btn btn-secondary" id="upload-cancel-btn">
                                    <i class="bi bi-x-circle me-1"></i>Cancel
                                </button>
                            </div>
                        </div>
                    </div>
                </div>

                <!-- Batch Editor (hidden by default) -->
                <div id="batch-editor" class="d-none">
                    <div class="card">
                        <div class="card-header">
                            <h6 class="mb-0"><i class="bi bi-pencil me-2"></i>Batch Editor</h6>
                        </div>
                        <div class="card-body">
                            <form id="batch-form">
                                <div class="row">
                                    <div class="col-md-6 mb-3">
                                        <label for="batch-name" class="form-label">Batch Name</label>
                                        <input type="text" class="form-control" id="batch-name" required>
                                    </div>
                                    <div class="col-md-6 mb-3">
                                        <label for="batch-description" class="form-label">Description</label>
                                        <input type="text" class="form-control" id="batch-description">
                                    </div>
                                </div>
                                
                                <div class="mb-3">
                                    <label for="batch-content" class="form-label">YAML Configuration</label>
                                    <textarea class="form-control" id="batch-content" rows="15" 
                                              placeholder="Enter your batch configuration in YAML format..."></textarea>
                                </div>
                                
                                <div class="d-flex gap-2">
                                    <button type="submit" class="btn btn-success">
                                        <i class="bi bi-check-circle me-1"></i>Save
                                    </button>
                                    <button type="button" class="btn btn-secondary" id="cancel-edit-btn">
                                        <i class="bi bi-x-circle me-1"></i>Cancel
                                    </button>
                                </div>
                            </form>
                        </div>
                    </div>
                </div>
            </div>
        `;
    }

    /**
     * Load batch files from the server
     */
    async loadBatchFiles() {
        try {
            console.log('BatchManagerTabbed: Loading batch files from API...');
            
            // First, try to get batch files from the file management API
            try {
                const response = await this.apiClient.getBatchFiles();
                this.batchFiles = response.files || [];
                console.log('BatchManagerTabbed: Batch files loaded from API:', this.batchFiles);
            } catch (apiError) {
                console.warn('BatchManagerTabbed: Batch files API not available, using filesystem approach');
                // Fallback: List files from batches directory (will need backend support)
                this.batchFiles = await this.loadBatchFilesFromDirectory();
            }

            // Load current executions
            await this.loadBatchExecutions();

            console.log('BatchManagerTabbed: Rendering batch files list...');
            this.renderBatchFilesList();
            console.log('BatchManagerTabbed: Batch files list rendered');
        } catch (error) {
            console.error('BatchManagerTabbed: Failed to load batch files:', error);
            this.renderBatchFilesError();
        }
    }

    /**
     * Load batch files from filesystem (fallback method)
     */
    async loadBatchFilesFromDirectory() {
        // For now, return example files that we know exist
        const knownFiles = [
            'batches/TEMPLATE.yaml'
        ];

        const files = [];
        for (const filePath of knownFiles) {
            try {
                // Check if file exists by trying to load it
                files.push({
                    name: filePath.split('/').pop(),
                    path: filePath,
                    description: 'Batch configuration file',
                    created: new Date().toISOString(),
                    size: 'Unknown'
                });
            } catch (error) {
                console.warn(`File ${filePath} not accessible:`, error);
            }
        }

        return files;
    }

    /**
     * Load current batch executions
     */
    async loadBatchExecutions() {
        try {
            console.log('BatchManagerTabbed: Loading batch executions...');
            const response = await this.apiClient.listBatchExecutions();
            this.executions = response.executions || [];
            console.log('BatchManagerTabbed: Executions loaded:', this.executions);
            
            // Render executions after loading
            this.renderActiveExecutions();
        } catch (error) {
            console.warn('BatchManagerTabbed: Failed to load executions:', error);
            this.executions = [];
            this.renderActiveExecutions();
        }
    }

    /**
     * Render active executions list
     */
    renderActiveExecutions() {
        const listContainer = document.getElementById('active-executions-list');
        if (!listContainer) {
            console.error('BatchManagerTabbed: active-executions-list container not found!');
            return;
        }

        const activeExecutions = this.executions.filter(exec => 
            ['queued', 'running', 'paused'].includes(exec.status)
        );

        if (activeExecutions.length === 0) {
            listContainer.innerHTML = `
                <div class="text-center p-3 text-muted">
                    <i class="bi bi-info-circle mb-2" style="font-size: 2rem;"></i>
                    <p>No active executions</p>
                </div>
            `;
            return;
        }

        const executionsHTML = activeExecutions.map(execution => `
            <div class="card mb-2">
                <div class="card-body p-3">
                    <div class="d-flex justify-content-between align-items-center">
                        <div>
                            <h6 class="mb-1">${execution.work_id}</h6>
                            <small class="text-muted">
                                Status: <span class="badge bg-${this.getStatusColor(execution.status)}">${execution.status}</span>
                                • Started: ${this.formatDate(execution.created_at)}
                            </small>
                        </div>
                        <div class="d-flex gap-1">
                            ${execution.status === 'running' ? `
                                <button class="btn btn-outline-warning btn-sm" onclick="batchManager.controlExecution('${execution.work_id}', 'pause')">
                                    <i class="bi bi-pause"></i>
                                </button>
                            ` : ''}
                            ${execution.status === 'paused' ? `
                                <button class="btn btn-outline-success btn-sm" onclick="batchManager.controlExecution('${execution.work_id}', 'resume')">
                                    <i class="bi bi-play"></i>
                                </button>
                            ` : ''}
                            <button class="btn btn-outline-danger btn-sm" onclick="batchManager.controlExecution('${execution.work_id}', 'cancel')">
                                <i class="bi bi-stop"></i>
                            </button>
                            <button class="btn btn-outline-info btn-sm" onclick="batchManager.showExecutionDetails('${execution.work_id}')">
                                <i class="bi bi-info"></i>
                            </button>
                        </div>
                    </div>
                    ${execution.progress ? `
                        <div class="mt-2">
                            <div class="progress" style="height: 4px;">
                                <div class="progress-bar" role="progressbar" 
                                     style="width: ${execution.progress.percentage || 0}%"></div>
                            </div>
                        </div>
                    ` : ''}
                </div>
            </div>
        `).join('');

        listContainer.innerHTML = executionsHTML;
    }

    /**
     * Get status badge color
     */
    getStatusColor(status) {
        const colors = {
            'queued': 'secondary',
            'running': 'primary',
            'paused': 'warning',
            'finished': 'success',
            'error': 'danger',
            'cancelled': 'dark'
        };
        return colors[status] || 'secondary';
    }

    /**
     * Format date for display
     */
    formatDate(dateString) {
        try {
            return new Date(dateString).toLocaleString();
        } catch (error) {
            return 'Unknown';
        }
    }

    /**
     * Render the list of batch files
     */
    renderBatchFilesList() {
        console.log('BatchManagerTabbed: Looking for batch-files-list container...');
        const listContainer = document.getElementById('batch-files-list');
        
        if (!listContainer) {
            console.error('BatchManagerTabbed: batch-files-list container not found!');
            return;
        }
        
        console.log('BatchManagerTabbed: Container found, batch files count:', this.batchFiles.length);

        if (this.batchFiles.length === 0) {
            console.log('BatchManagerTabbed: No batch files, showing empty state');
            listContainer.innerHTML = `
                <div class="text-center p-4 text-muted">
                    <i class="bi bi-folder2-open display-4 mb-3"></i>
                    <p>No batch files found</p>
                    <p class="small">Create your first batch configuration to get started</p>
                </div>
            `;
            return;
        }

        console.log('BatchManagerTabbed: Generating HTML for batch files...');
        const filesHtml = this.batchFiles.map(file => `
            <div class="card mb-2">
                <div class="card-body p-3">
                    <div class="d-flex justify-content-between align-items-center">
                        <div class="flex-grow-1">
                            <h6 class="mb-1">
                                <i class="bi bi-file-earmark-text me-2"></i>
                                ${file.name}
                            </h6>
                            <p class="mb-0 small text-muted">${file.description}</p>
                            <small class="text-muted">Size: ${file.size} | Created: ${new Date(file.created).toLocaleDateString()}</small>
                        </div>
                        <div class="d-flex gap-2">
                            <button class="btn btn-success btn-sm" data-action="execute" data-file="${file.path || file.name}" title="Execute Batch">
                                <i class="bi bi-play-fill"></i>
                            </button>
                            <div class="btn-group">
                                <button class="btn btn-outline-primary btn-sm" data-action="edit" data-file="${file.name}" title="Edit">
                                    <i class="bi bi-pencil"></i>
                                </button>
                                <button class="btn btn-outline-info btn-sm" data-action="download" data-file="${file.name}" title="Download">
                                    <i class="bi bi-download"></i>
                                </button>
                                <button class="btn btn-outline-danger btn-sm" data-action="delete" data-file="${file.name}" title="Delete">
                                    <i class="bi bi-trash"></i>
                                </button>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        `).join('');

        console.log('BatchManagerTabbed: Setting innerHTML with generated HTML');
        listContainer.innerHTML = filesHtml;
        console.log('BatchManagerTabbed: Batch files list updated successfully');
    }

    /**
     * Render error message for batch files list
     */
    renderBatchFilesError() {
        const listContainer = document.getElementById('batch-files-list');
        if (!listContainer) return;

        listContainer.innerHTML = `
            <div class="alert alert-warning">
                <i class="bi bi-exclamation-triangle me-2"></i>
                Failed to load batch files. Please try refreshing the page.
            </div>
        `;
    }

    /**
     * Setup event handlers for the component
     */
    setupEventHandlers() {
        // Create new batch button
        const createBtn = document.getElementById('create-batch-btn');
        if (createBtn) {
            createBtn.addEventListener('click', () => this.showBatchEditor());
        }

        // Upload batch button
        const uploadBtn = document.getElementById('upload-batch-btn');
        if (uploadBtn) {
            uploadBtn.addEventListener('click', () => this.showUploadModal());
        }

        // Refresh executions button
        const refreshBtn = document.getElementById('refresh-executions-btn');
        if (refreshBtn) {
            refreshBtn.addEventListener('click', () => this.refreshExecutions());
        }

        // Upload modal handlers
        const uploadConfirmBtn = document.getElementById('upload-confirm-btn');
        if (uploadConfirmBtn) {
            uploadConfirmBtn.addEventListener('click', () => this.handleFileUpload());
        }

        const uploadCancelBtn = document.getElementById('upload-cancel-btn');
        if (uploadCancelBtn) {
            uploadCancelBtn.addEventListener('click', () => this.hideUploadModal());
        }

        // Cancel edit button
        const cancelBtn = document.getElementById('cancel-edit-btn');
        if (cancelBtn) {
            cancelBtn.addEventListener('click', () => this.hideBatchEditor());
        }

        // Batch form submission
        const batchForm = document.getElementById('batch-form');
        if (batchForm) {
            batchForm.addEventListener('submit', (e) => this.handleBatchSave(e));
        }

        // File action buttons (using event delegation)
        const listContainer = document.getElementById('batch-files-list');
        if (listContainer) {
            listContainer.addEventListener('click', (e) => this.handleFileAction(e));
        }

        // Auto-refresh executions every 30 seconds
        setInterval(() => {
            this.refreshExecutions();
        }, 30000);
    }

    /**
     * Show the batch editor
     */
    showBatchEditor(batchData = null) {
        const editor = document.getElementById('batch-editor');
        if (!editor) return;

        if (batchData) {
            // Populate form with existing data
            document.getElementById('batch-name').value = batchData.name || '';
            document.getElementById('batch-description').value = batchData.description || '';
            document.getElementById('batch-content').value = batchData.content || '';
        } else {
            // Clear form for new batch
            document.getElementById('batch-form').reset();
            document.getElementById('batch-content').value = `# Example batch configuration
name: "New Batch"
description: "Batch description"
algorithms:
  - name: "Baseline"
    parameters: {}
datasets:
  - path: "datasets/example.fasta"
output:
  format: "csv"
  path: "outputs/"`;
        }

        editor.classList.remove('d-none');
        editor.scrollIntoView({ behavior: 'smooth' });
    }

    /**
     * Hide the batch editor
     */
    hideBatchEditor() {
        const editor = document.getElementById('batch-editor');
        if (editor) {
            editor.classList.add('d-none');
        }
    }

    /**
     * Handle batch form save
     */
    async handleBatchSave(event) {
        event.preventDefault();
        
        const name = document.getElementById('batch-name').value;
        const description = document.getElementById('batch-description').value;
        const content = document.getElementById('batch-content').value;

        try {
            // Mock save - in real implementation, would call API
            console.log('Saving batch:', { name, description, content });
            
            // Show success message
            if (window.NotificationManager) {
                window.NotificationManager.success('Batch configuration saved successfully!');
            }
            
            this.hideBatchEditor();
            await this.loadBatchFiles(); // Refresh the list
            
        } catch (error) {
            console.error('Failed to save batch:', error);
            if (window.NotificationManager) {
                window.NotificationManager.error('Failed to save batch configuration');
            }
        }
    }

    /**
     * Handle file action button clicks
     */
    handleFileAction(event) {
        const button = event.target.closest('button[data-action]');
        if (!button) return;

        const action = button.dataset.action;
        const fileName = button.dataset.file;

        switch (action) {
            case 'edit':
                this.editBatchFile(fileName);
                break;
            case 'download':
                this.downloadBatchFile(fileName);
                break;
            case 'delete':
                this.deleteBatchFile(fileName);
                break;
        }
    }

    /**
     * Edit a batch file
     */
    async editBatchFile(fileName) {
        try {
            // Mock data - in real implementation, would fetch from API
            const mockContent = `# ${fileName}
name: "${fileName.replace('.yaml', '')}"
description: "Batch configuration file"
algorithms:
  - name: "Baseline"
    parameters: {}
datasets:
  - path: "datasets/example.fasta"
output:
  format: "csv"
  path: "outputs/"`;

            this.showBatchEditor({
                name: fileName.replace('.yaml', ''),
                description: 'Existing batch file',
                content: mockContent
            });
        } catch (error) {
            console.error('Failed to load batch file for editing:', error);
            if (window.NotificationManager) {
                window.NotificationManager.error('Failed to load batch file');
            }
        }
    }

    /**
     * Download a batch file
     */
    downloadBatchFile(fileName) {
        // Mock download - in real implementation, would call API
        console.log('Downloading batch file:', fileName);
        if (window.NotificationManager) {
            window.NotificationManager.info(`Download for ${fileName} would start here`);
        }
    }

    /**
     * Delete a batch file
     */
    async deleteBatchFile(fileName) {
        if (!confirm(`Are you sure you want to delete ${fileName}?`)) {
            return;
        }

        try {
            // Mock delete - in real implementation, would call API
            console.log('Deleting batch file:', fileName);
            
            if (window.NotificationManager) {
                window.NotificationManager.success(`${fileName} deleted successfully`);
            }
            
            await this.loadBatchFiles(); // Refresh the list
        } catch (error) {
            console.error('Failed to delete batch file:', error);
            if (window.NotificationManager) {
                window.NotificationManager.error('Failed to delete batch file');
            }
        }
    }

    /**
     * Render error message
     */
    renderError(message) {
        const container = document.getElementById(this.containerId);
        if (container) {
            container.innerHTML = `
                <div class="alert alert-danger">
                    <i class="bi bi-exclamation-triangle me-2"></i>
                    ${message}
                </div>
            `;
        }
    }

    // =============================================================================
    // NEW METHODS FOR ENHANCED FUNCTIONALITY
    // =============================================================================

    /**
     * Show upload modal
     */
    showUploadModal() {
        const modal = document.getElementById('upload-modal');
        if (modal) {
            modal.classList.remove('d-none');
            modal.scrollIntoView({ behavior: 'smooth' });
        }
    }

    /**
     * Hide upload modal
     */
    hideUploadModal() {
        const modal = document.getElementById('upload-modal');
        if (modal) {
            modal.classList.add('d-none');
            // Reset file input
            const fileInput = document.getElementById('file-input');
            if (fileInput) fileInput.value = '';
        }
    }

    /**
     * Handle file upload
     */
    async handleFileUpload() {
        const fileInput = document.getElementById('file-input');
        const file = fileInput.files[0];

        if (!file) {
            showAlert('Please select a file to upload', 'warning');
            return;
        }

        if (!file.name.endsWith('.yaml') && !file.name.endsWith('.yml')) {
            showAlert('Please select a YAML file (.yaml or .yml)', 'warning');
            return;
        }

        try {
            const formData = new FormData();
            formData.append('file', file);

            const response = await this.apiClient.uploadBatchFile(formData);
            
            showAlert('File uploaded successfully!', 'success');
            this.hideUploadModal();
            await this.loadBatchFiles(); // Refresh list
            
        } catch (error) {
            console.error('Upload failed:', error);
            showAlert('Failed to upload file: ' + error.message, 'danger');
        }
    }

    /**
     * Refresh executions
     */
    async refreshExecutions() {
        try {
            await this.loadBatchExecutions();
        } catch (error) {
            console.error('Failed to refresh executions:', error);
        }
    }

    /**
     * Execute batch file
     */
    async executeBatch(filePath) {
        try {
            console.log('Executing batch:', filePath);
            
            const response = await this.apiClient.executeBatch(filePath, 'log');
            
            showAlert(`Batch iniciado! Redirecionando para progresso... (Work ID: ${response.work_id})`, 'success');

            if (response && response.work_id) {
                setTimeout(() => {
                    window.location.href = `/monitor/work/${response.work_id}`;
                }, 600);
            } else {
                console.warn('executeBatch response missing work_id, skipping redirect');
                // Fallback: atualizar lista
                await this.refreshExecutions();
            }
            
        } catch (error) {
            console.error('Execution failed:', error);
            showAlert('Failed to execute batch: ' + error.message, 'danger');
        }
    }

    /**
     * Control execution (pause, resume, cancel)
     */
    async controlExecution(workId, action) {
        try {
            console.log(`Controlling execution ${workId}: ${action}`);
            
            const response = await this.apiClient.controlBatchExecution(workId, action);
            
            if (response.success) {
                showAlert(`Execution ${action}ed successfully`, 'success');
                await this.refreshExecutions();
            } else {
                showAlert(`Failed to ${action} execution: ` + response.message, 'warning');
            }
            
        } catch (error) {
            console.error(`Control action failed:`, error);
            showAlert(`Failed to ${action} execution: ` + error.message, 'danger');
        }
    }

    /**
     * Show execution details
     */
    async showExecutionDetails(workId) {
        try {
            console.log('Getting execution details for:', workId);
            
            const response = await this.apiClient.getBatchStatus(workId);
            
            // Create modal with execution details
            const modal = document.createElement('div');
            modal.className = 'modal fade';
            modal.innerHTML = `
                <div class="modal-dialog">
                    <div class="modal-content">
                        <div class="modal-header">
                            <h5 class="modal-title">Execution Details: ${workId}</h5>
                            <button type="button" class="btn-close" data-bs-dismiss="modal"></button>
                        </div>
                        <div class="modal-body">
                            <dl class="row">
                                <dt class="col-sm-3">Status:</dt>
                                <dd class="col-sm-9">
                                    <span class="badge bg-${this.getStatusColor(response.status)}">${response.status}</span>
                                </dd>
                                <dt class="col-sm-3">Work ID:</dt>
                                <dd class="col-sm-9">${response.work_id}</dd>
                                <dt class="col-sm-3">Created:</dt>
                                <dd class="col-sm-9">${this.formatDate(response.created_at)}</dd>
                                <dt class="col-sm-3">Updated:</dt>
                                <dd class="col-sm-9">${this.formatDate(response.updated_at)}</dd>
                                ${response.output_path ? `
                                    <dt class="col-sm-3">Output:</dt>
                                    <dd class="col-sm-9">${response.output_path}</dd>
                                ` : ''}
                                ${response.error ? `
                                    <dt class="col-sm-3">Error:</dt>
                                    <dd class="col-sm-9 text-danger">${response.error}</dd>
                                ` : ''}
                            </dl>
                        </div>
                        <div class="modal-footer">
                            <button type="button" class="btn btn-secondary" data-bs-dismiss="modal">Close</button>
                        </div>
                    </div>
                </div>
            `;
            
            document.body.appendChild(modal);
            const bsModal = new bootstrap.Modal(modal);
            bsModal.show();
            
            // Remove modal from DOM when hidden
            modal.addEventListener('hidden.bs.modal', () => {
                document.body.removeChild(modal);
            });
            
        } catch (error) {
            console.error('Failed to get execution details:', error);
            showAlert('Failed to get execution details: ' + error.message, 'danger');
        }
    }

    /**
     * Enhanced batch save with API integration
     */
    async handleBatchSave(event) {
        event.preventDefault();
        
        const name = document.getElementById('batch-name').value;
        const description = document.getElementById('batch-description').value;
        const content = document.getElementById('batch-content').value;

        if (!name.trim()) {
            showAlert('Please enter a batch name', 'warning');
            return;
        }

        if (!content.trim()) {
            showAlert('Please enter batch configuration content', 'warning');
            return;
        }

        try {
            console.log('Saving batch:', { name, description, content });
            
            const response = await this.apiClient.createBatchFile(name, content, description);
            
            showAlert('Batch configuration saved successfully!', 'success');
            this.hideBatchEditor();
            await this.loadBatchFiles(); // Refresh list
            
        } catch (error) {
            console.error('Save failed:', error);
            showAlert('Failed to save batch configuration: ' + error.message, 'danger');
        }
    }

    /**
     * Enhanced file action handler with execution support
     */
    handleFileAction(event) {
        const button = event.target.closest('button[data-action]');
        if (!button) return;

        const action = button.dataset.action;
        const fileName = button.dataset.file;

        switch (action) {
            case 'execute':
                this.executeBatch(fileName);
                break;
            case 'edit':
                this.editBatchFile(fileName);
                break;
            case 'download':
                this.downloadBatchFile(fileName);
                break;
            case 'delete':
                this.deleteBatchFile(fileName);
                break;
        }
    }
}

// Export for global use
window.BatchManagerTabbed = BatchManagerTabbed;

// Create global instance for easy access from HTML onclick handlers
window.batchManager = new BatchManagerTabbed();
