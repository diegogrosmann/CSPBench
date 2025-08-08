/**
 * Batch Manager Component with Tabs
 * 
 * Handles the batch execution interface with separate tabs for:
 * - Batch Files management
 * - Execution monitoring
 */

class BatchManagerTabbed {
    constructor() {
        this.currentPage = 1;
        this.pageSize = 20;
        this.selectedBatches = new Set();
        this.selectedExecutions = new Set();
        this.batches = [];
        this.executions = [];
        this.filteredBatches = [];
        this.filteredExecutions = [];
        this.batchFilters = {
            search: '',
            dateFrom: '',
            dateTo: ''
        };
        this.executionFilters = {
            search: '',
            status: '',
            dateFrom: '',
            dateTo: ''
        };
        
        this.apiClient = window.apiClient || new APIClient();
        this.bindEvents();
    }

    initialize() {
        this.loadBatches();
        this.loadExecutions();
    }

    bindEvents() {
        // Tab change events
        document.getElementById('batch-files-tab')?.addEventListener('shown.bs.tab', () => {
            this.loadBatches();
        });
        document.getElementById('executions-tab')?.addEventListener('shown.bs.tab', () => {
            this.loadExecutions();
        });

        // Batch Files Events
        document.getElementById('batch-search-query')?.addEventListener('input', 
            this.debounce((e) => this.updateBatchFilter('search', e.target.value), 300));
        document.getElementById('batch-filter-date-from')?.addEventListener('change', 
            (e) => this.updateBatchFilter('dateFrom', e.target.value));
        document.getElementById('batch-filter-date-to')?.addEventListener('change', 
            (e) => this.updateBatchFilter('dateTo', e.target.value));
        document.getElementById('apply-batch-filters')?.addEventListener('click', 
            () => this.applyBatchFilters());

        document.getElementById('run-selected')?.addEventListener('click', 
            () => this.runSelectedBatches());
        document.getElementById('delete-selected-batches')?.addEventListener('click', 
            () => this.deleteSelectedBatches());
        document.getElementById('select-all-batches')?.addEventListener('change', 
            (e) => e.target.checked ? this.selectAllBatches() : this.clearBatchSelection());
        document.getElementById('refresh-batches')?.addEventListener('click', 
            () => this.loadBatches());
        document.getElementById('upload-batch')?.addEventListener('click', 
            () => this.showUploadModal());
        document.getElementById('upload-submit')?.addEventListener('click', 
            () => this.uploadBatchFile());

        // Executions Events  
        document.getElementById('execution-search-query')?.addEventListener('input', 
            this.debounce((e) => this.updateExecutionFilter('search', e.target.value), 300));
        document.getElementById('execution-filter-status')?.addEventListener('change', 
            (e) => this.updateExecutionFilter('status', e.target.value));
        document.getElementById('execution-filter-date-from')?.addEventListener('change', 
            (e) => this.updateExecutionFilter('dateFrom', e.target.value));
        document.getElementById('execution-filter-date-to')?.addEventListener('change', 
            (e) => this.updateExecutionFilter('dateTo', e.target.value));
        document.getElementById('apply-execution-filters')?.addEventListener('click', 
            () => this.applyExecutionFilters());

        document.getElementById('stop-selected')?.addEventListener('click', 
            () => this.stopSelectedExecutions());
        document.getElementById('select-all-executions')?.addEventListener('change', 
            (e) => e.target.checked ? this.selectAllExecutions() : this.clearExecutionSelection());
        document.getElementById('refresh-executions')?.addEventListener('click', 
            () => this.loadExecutions());

        // Filter toggle events
        document.getElementById('toggle-batch-filters')?.addEventListener('click', 
            () => this.toggleFiltersIcon('batch-filters'));
        document.getElementById('toggle-execution-filters')?.addEventListener('click', 
            () => this.toggleFiltersIcon('execution-filters'));

        // Retry events
        document.getElementById('retry-batch-load')?.addEventListener('click', 
            () => this.loadBatches());
        document.getElementById('retry-execution-load')?.addEventListener('click', 
            () => this.loadExecutions());
    }

    async loadBatches() {
        this.showBatchLoading();
        
        try {
            const response = await this.apiClient.get('/execution/api/batches', 'load_batches');
            this.batches = response.batch_files || [];
            this.updateStatistics();
            this.applyBatchFilters();
            this.hideBatchLoading();
        } catch (error) {
            console.error('Error loading batches:', error);
            this.showBatchError('Failed to load batch files: ' + (error.message || 'Unknown error'));
        }
    }

    async loadExecutions() {
        this.showExecutionLoading();
        
        try {
            const response = await this.apiClient.get('/execution/api/execution/sessions', 'load_executions');
            
            // Transform sessions data to match our interface
            this.executions = (response.sessions || []).map(session => ({
                id: session.session_id,
                batch_name: session.config?.batch_file || 'Unknown',
                started: session.started_at || new Date().toISOString(),
                status: session.status || 'unknown',
                progress: session.progress || 0,
                duration: session.duration || 0
            }));
            
            this.updateStatistics();
            this.applyExecutionFilters();
            this.hideExecutionLoading();
        } catch (error) {
            console.error('Error loading executions:', error);
            this.showExecutionError('Failed to load executions: ' + (error.message || 'Unknown error'));
        }
    }

    updateStatistics() {
        // Batch files statistics
        const totalBatches = this.batches.length;
        
        // Execution statistics
        const runningExecutions = this.executions.filter(e => e.status === 'running').length;
        const completedExecutions = this.executions.filter(e => e.status === 'completed').length;
        const failedExecutions = this.executions.filter(e => e.status === 'failed').length;

        // Update DOM elements
        const totalBatchesEl = document.getElementById('total-batches');
        const runningExecutionsEl = document.getElementById('running-executions');
        const completedExecutionsEl = document.getElementById('completed-executions');
        const failedExecutionsEl = document.getElementById('failed-executions');

        if (totalBatchesEl) totalBatchesEl.textContent = totalBatches;
        if (runningExecutionsEl) runningExecutionsEl.textContent = runningExecutions;
        if (completedExecutionsEl) completedExecutionsEl.textContent = completedExecutions;
        if (failedExecutionsEl) failedExecutionsEl.textContent = failedExecutions;
    }

    // Batch Files Methods
    updateBatchFilter(filterName, value) {
        this.batchFilters[filterName] = value;
        this.applyBatchFilters();
    }

    applyBatchFilters() {
        this.filteredBatches = this.batches.filter(batch => {
            // Search filter
            if (this.batchFilters.search) {
                const searchLower = this.batchFilters.search.toLowerCase();
                const searchableText = [
                    batch.name,
                    batch.description,
                    batch.filename
                ].filter(Boolean).join(' ').toLowerCase();
                
                if (!searchableText.includes(searchLower)) {
                    return false;
                }
            }

            // Date filters
            if (this.batchFilters.dateFrom) {
                const fromDate = new Date(this.batchFilters.dateFrom);
                const batchDate = new Date(batch.modified * 1000);
                if (batchDate < fromDate) {
                    return false;
                }
            }

            if (this.batchFilters.dateTo) {
                const toDate = new Date(this.batchFilters.dateTo);
                toDate.setHours(23, 59, 59, 999);
                const batchDate = new Date(batch.modified * 1000);
                if (batchDate > toDate) {
                    return false;
                }
            }

            return true;
        });

        this.renderBatches();
    }

    renderBatches() {
        const tbody = document.getElementById('batch-table-body');
        const tableContainer = document.getElementById('batch-table-container');
        const emptyState = document.getElementById('batch-empty-state');

        if (!tbody) return;

        if (this.filteredBatches.length === 0) {
            if (tableContainer) tableContainer.style.display = 'none';
            if (emptyState) emptyState.style.display = 'block';
            return;
        }

        if (tableContainer) tableContainer.style.display = 'block';
        if (emptyState) emptyState.style.display = 'none';

        tbody.innerHTML = this.filteredBatches.map(batch => this.renderBatchRow(batch)).join('');
        this.bindBatchRowEvents();
    }

    renderBatchRow(batch) {
        const createdDate = new Date(batch.modified * 1000).toLocaleString();
        const size = this.formatFileSize(batch.size);
        const isSelected = this.selectedBatches.has(batch.filename);

        return `
            <tr class="batch-row" data-batch-id="${batch.filename}">
                <td>
                    <input type="checkbox" class="form-check-input batch-checkbox" 
                           ${isSelected ? 'checked' : ''} 
                           data-batch-id="${batch.filename}">
                </td>
                <td>
                    <strong>${batch.name || batch.filename}</strong><br>
                    <small class="text-muted">${batch.filename}</small>
                </td>
                <td>
                    <small class="text-muted">${batch.description || 'No description'}</small>
                </td>
                <td>
                    <small>${createdDate}</small>
                </td>
                <td>
                    <small>${size}</small>
                </td>
                <td>
                    <div class="btn-group btn-group-sm" role="group">
                        <button class="btn btn-outline-success action-btn run-batch" 
                                data-batch-id="${batch.filename}" 
                                title="Run Batch">
                            <i class="bi bi-play"></i>
                        </button>
                        <button class="btn btn-outline-primary action-btn view-batch" 
                                data-batch-id="${batch.filename}" 
                                title="View/Edit">
                            <i class="bi bi-eye"></i>
                        </button>
                        <button class="btn btn-outline-secondary action-btn download-batch" 
                                data-batch-id="${batch.filename}" 
                                title="Download">
                            <i class="bi bi-download"></i>
                        </button>
                        <button class="btn btn-outline-danger action-btn delete-batch" 
                                data-batch-id="${batch.filename}" 
                                title="Delete">
                            <i class="bi bi-trash"></i>
                        </button>
                    </div>
                </td>
            </tr>
        `;
    }

    // Execution Methods
    updateExecutionFilter(filterName, value) {
        this.executionFilters[filterName] = value;
        this.applyExecutionFilters();
    }

    applyExecutionFilters() {
        this.filteredExecutions = this.executions.filter(execution => {
            // Search filter
            if (this.executionFilters.search) {
                const searchLower = this.executionFilters.search.toLowerCase();
                const searchableText = [
                    execution.id,
                    execution.batch_name
                ].filter(Boolean).join(' ').toLowerCase();
                
                if (!searchableText.includes(searchLower)) {
                    return false;
                }
            }

            // Status filter
            if (this.executionFilters.status && execution.status !== this.executionFilters.status) {
                return false;
            }

            // Date filters
            if (this.executionFilters.dateFrom) {
                const fromDate = new Date(this.executionFilters.dateFrom);
                const executionDate = new Date(execution.started);
                if (executionDate < fromDate) {
                    return false;
                }
            }

            if (this.executionFilters.dateTo) {
                const toDate = new Date(this.executionFilters.dateTo);
                toDate.setHours(23, 59, 59, 999);
                const executionDate = new Date(execution.started);
                if (executionDate > toDate) {
                    return false;
                }
            }

            return true;
        });

        this.renderExecutions();
    }

    renderExecutions() {
        const tbody = document.getElementById('execution-table-body');
        const tableContainer = document.getElementById('execution-table-container');
        const emptyState = document.getElementById('execution-empty-state');

        if (!tbody) return;

        if (this.filteredExecutions.length === 0) {
            if (tableContainer) tableContainer.style.display = 'none';
            if (emptyState) emptyState.style.display = 'block';
            return;
        }

        if (tableContainer) tableContainer.style.display = 'block';
        if (emptyState) emptyState.style.display = 'none';

        tbody.innerHTML = this.filteredExecutions.map(execution => this.renderExecutionRow(execution)).join('');
        this.bindExecutionRowEvents();
    }

    renderExecutionRow(execution) {
        const startedDate = new Date(execution.started).toLocaleString();
        const duration = this.formatDuration(execution.duration || 0);
        const progress = execution.progress || 0;
        const isSelected = this.selectedExecutions.has(execution.id);

        // Determine if execution is completed or still running
        const isCompleted = ['completed', 'failed', 'stopped'].includes(execution.status);
        const executionUrl = isCompleted ? `/results/${execution.id}` : `/execution/progress/${execution.id}`;

        return `
            <tr class="execution-row" data-execution-id="${execution.id}">
                <td>
                    <input type="checkbox" class="form-check-input execution-checkbox" 
                           ${isSelected ? 'checked' : ''} 
                           data-execution-id="${execution.id}">
                </td>
                <td>
                    <a href="${executionUrl}" class="text-decoration-none execution-link" 
                       data-execution-id="${execution.id}">
                        <code class="text-primary">${execution.id}</code>
                    </a>
                </td>
                <td>
                    <small>${execution.batch_name}</small>
                </td>
                <td>
                    <small>${startedDate}</small>
                </td>
                <td>
                    <span class="badge status-badge status-${execution.status}">
                        ${execution.status}
                    </span>
                </td>
                <td>
                    <div class="progress progress-sm">
                        <div class="progress-bar" role="progressbar" 
                             style="width: ${progress}%" 
                             aria-valuenow="${progress}" 
                             aria-valuemin="0" 
                             aria-valuemax="100">
                        </div>
                    </div>
                    <small class="text-muted">${progress}%</small>
                </td>
                <td>
                    <small>${duration}</small>
                </td>
                <td>
                    ${execution.status === 'running' ? `
                        <button class="btn btn-outline-warning btn-sm action-btn stop-execution" 
                                data-execution-id="${execution.id}" 
                                title="Stop">
                            <i class="bi bi-stop"></i>
                        </button>
                    ` : ''}
                </td>
            </tr>
        `;
    }

    // Upload functionality
    async uploadBatchFile() {
        const fileInput = document.getElementById('batch-file');
        const nameInput = document.getElementById('batch-name');
        
        if (!fileInput || !fileInput.files[0]) {
            alert('Please select a file to upload');
            return;
        }

        const formData = new FormData();
        formData.append('file', fileInput.files[0]);
        if (nameInput && nameInput.value) {
            formData.append('name', nameInput.value);
        }

        try {
            const response = await fetch('/execution/api/batches/upload', {
                method: 'POST',
                body: formData
            });

            if (!response.ok) {
                throw new Error(`Upload failed: ${response.statusText}`);
            }

            const result = await response.json();
            
            // Close modal and refresh
            const modalEl = document.getElementById('uploadModal');
            if (modalEl) {
                const modal = bootstrap.Modal.getInstance(modalEl);
                if (modal) modal.hide();
            }
            
            // Clear form
            if (fileInput) fileInput.value = '';
            if (nameInput) nameInput.value = '';
            
            // Refresh batch list
            this.loadBatches();
            
            // Show success message
            this.showSuccessMessage(`Batch file uploaded successfully: ${result.filename}`);
            
        } catch (error) {
            console.error('Upload error:', error);
            this.showErrorMessage('Failed to upload batch file: ' + error.message);
        }
    }

    // Utility methods
    formatFileSize(bytes) {
        if (!bytes || bytes === 0) return '0 Bytes';
        const k = 1024;
        const sizes = ['Bytes', 'KB', 'MB', 'GB'];
        const i = Math.floor(Math.log(bytes) / Math.log(k));
        return parseFloat((bytes / Math.pow(k, i)).toFixed(2)) + ' ' + sizes[i];
    }

    formatDuration(seconds) {
        if (!seconds || seconds < 60) return `${(seconds || 0).toFixed(1)}s`;
        const minutes = Math.floor(seconds / 60);
        const remainingSeconds = seconds % 60;
        return `${minutes}m ${remainingSeconds.toFixed(1)}s`;
    }

    showSuccessMessage(message) {
        console.log('Success:', message);
        
        // Create and show a temporary alert
        const alertDiv = document.createElement('div');
        alertDiv.className = 'alert alert-success alert-dismissible fade show position-fixed';
        alertDiv.style.cssText = 'top: 20px; right: 20px; z-index: 9999; min-width: 300px;';
        alertDiv.innerHTML = `
            <i class="bi bi-check-circle me-2"></i>${message}
            <button type="button" class="btn-close" data-bs-dismiss="alert"></button>
        `;
        
        document.body.appendChild(alertDiv);
        
        // Auto remove after 5 seconds
        setTimeout(() => {
            if (alertDiv.parentNode) {
                alertDiv.remove();
            }
        }, 5000);
    }

    showErrorMessage(message) {
        console.error('Error:', message);
        
        // Create and show a temporary alert
        const alertDiv = document.createElement('div');
        alertDiv.className = 'alert alert-danger alert-dismissible fade show position-fixed';
        alertDiv.style.cssText = 'top: 20px; right: 20px; z-index: 9999; min-width: 300px;';
        alertDiv.innerHTML = `
            <i class="bi bi-exclamation-triangle me-2"></i>${message}
            <button type="button" class="btn-close" data-bs-dismiss="alert"></button>
        `;
        
        document.body.appendChild(alertDiv);
        
        // Auto remove after 8 seconds (longer for errors)
        setTimeout(() => {
            if (alertDiv.parentNode) {
                alertDiv.remove();
            }
        }, 8000);
    }

    // Loading states
    showBatchLoading() {
        const loading = document.getElementById('batch-loading-indicator');
        const table = document.getElementById('batch-table-container');
        const empty = document.getElementById('batch-empty-state');
        const error = document.getElementById('batch-error-state');

        if (loading) loading.style.display = 'block';
        if (table) table.style.display = 'none';
        if (empty) empty.style.display = 'none';
        if (error) error.style.display = 'none';
    }

    hideBatchLoading() {
        const loading = document.getElementById('batch-loading-indicator');
        if (loading) loading.style.display = 'none';
    }

    showBatchError(message) {
        const loading = document.getElementById('batch-loading-indicator');
        const table = document.getElementById('batch-table-container');
        const empty = document.getElementById('batch-empty-state');
        const error = document.getElementById('batch-error-state');
        const errorMsg = document.getElementById('batch-error-message');

        if (loading) loading.style.display = 'none';
        if (table) table.style.display = 'none';
        if (empty) empty.style.display = 'none';
        if (error) error.style.display = 'block';
        if (errorMsg) errorMsg.textContent = message;
    }

    showExecutionLoading() {
        const loading = document.getElementById('execution-loading-indicator');
        const table = document.getElementById('execution-table-container');
        const empty = document.getElementById('execution-empty-state');
        const error = document.getElementById('execution-error-state');

        if (loading) loading.style.display = 'block';
        if (table) table.style.display = 'none';
        if (empty) empty.style.display = 'none';
        if (error) error.style.display = 'none';
    }

    hideExecutionLoading() {
        const loading = document.getElementById('execution-loading-indicator');
        if (loading) loading.style.display = 'none';
    }

    showExecutionError(message) {
        const loading = document.getElementById('execution-loading-indicator');
        const table = document.getElementById('execution-table-container');
        const empty = document.getElementById('execution-empty-state');
        const error = document.getElementById('execution-error-state');
        const errorMsg = document.getElementById('execution-error-message');

        if (loading) loading.style.display = 'none';
        if (table) table.style.display = 'none';
        if (empty) empty.style.display = 'none';
        if (error) error.style.display = 'block';
        if (errorMsg) errorMsg.textContent = message;
    }

    debounce(func, wait) {
        let timeout;
        return function executedFunction(...args) {
            const later = () => {
                clearTimeout(timeout);
                func(...args);
            };
            clearTimeout(timeout);
            timeout = setTimeout(later, wait);
        };
    }

    // Event binding methods
    bindBatchRowEvents() {
        // Checkbox events
        document.querySelectorAll('.batch-checkbox').forEach(checkbox => {
            checkbox.addEventListener('change', (e) => {
                const batchId = e.target.dataset.batchId;
                if (e.target.checked) {
                    this.selectedBatches.add(batchId);
                } else {
                    this.selectedBatches.delete(batchId);
                }
                this.updateBatchSelectionUI();
            });
        });

        // Action button events for batches
        document.querySelectorAll('.run-batch').forEach(btn => {
            btn.addEventListener('click', (e) => {
                const batchId = e.target.closest('[data-batch-id]').dataset.batchId;
                this.runSingleBatch(batchId);
            });
        });

        document.querySelectorAll('.delete-batch').forEach(btn => {
            btn.addEventListener('click', (e) => {
                const batchId = e.target.closest('[data-batch-id]').dataset.batchId;
                this.deleteSingleBatch(batchId);
            });
        });

        document.querySelectorAll('.view-batch').forEach(btn => {
            btn.addEventListener('click', (e) => {
                const batchId = e.target.closest('[data-batch-id]').dataset.batchId;
                this.viewBatchDetails(batchId);
            });
        });

        document.querySelectorAll('.download-batch').forEach(btn => {
            btn.addEventListener('click', (e) => {
                const batchId = e.target.closest('[data-batch-id]').dataset.batchId;
                this.downloadBatch(batchId);
            });
        });
    }

    bindExecutionRowEvents() {
        // Checkbox events
        document.querySelectorAll('.execution-checkbox').forEach(checkbox => {
            checkbox.addEventListener('change', (e) => {
                const executionId = e.target.dataset.executionId;
                if (e.target.checked) {
                    this.selectedExecutions.add(executionId);
                } else {
                    this.selectedExecutions.delete(executionId);
                }
                this.updateExecutionSelectionUI();
            });
        });

        // Stop execution button events (only for running executions)
        document.querySelectorAll('.stop-execution').forEach(btn => {
            btn.addEventListener('click', (e) => {
                e.preventDefault(); // Prevent any default action
                const executionId = e.target.closest('[data-execution-id]').dataset.executionId;
                this.stopSingleExecution(executionId);
            });
        });
    }

    updateBatchSelectionUI() {
        const count = this.selectedBatches.size;
        const countEl = document.getElementById('selected-batch-count');
        if (countEl) countEl.textContent = `${count} selected`;
        
        const batchActions = document.getElementById('batch-actions');
        if (batchActions) {
            batchActions.style.display = count > 0 ? 'block' : 'none';
        }
    }

    updateExecutionSelectionUI() {
        const count = this.selectedExecutions.size;
        const countEl = document.getElementById('selected-execution-count');
        if (countEl) countEl.textContent = `${count} selected`;
        
        const executionActions = document.getElementById('execution-actions');
        if (executionActions) {
            executionActions.style.display = count > 0 ? 'block' : 'none';
        }
    }

    // Action methods
    async runSingleBatch(batchId) {
        try {
            console.log('Running batch:', batchId);
            const response = await fetch(`/execution/api/batches/${batchId}/execute`, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                }
            });

            if (!response.ok) {
                throw new Error(`Failed to start batch execution: ${response.statusText}`);
            }

            const result = await response.json();
            this.showSuccessMessage(`Batch execution started: ${result.session_id}`);
            
            // Navigate to execution progress page
            setTimeout(() => {
                window.location.href = `/execution/progress/${result.session_id}`;
            }, 1000);
            
        } catch (error) {
            console.error('Error running batch:', error);
            this.showErrorMessage('Failed to start batch execution: ' + error.message);
        }
    }

    async deleteSingleBatch(batchId) {
        if (confirm(`Are you sure you want to delete batch "${batchId}"?`)) {
            try {
                console.log('Deleting batch:', batchId);
                const response = await fetch(`/execution/api/batches/${batchId}`, {
                    method: 'DELETE'
                });

                if (!response.ok) {
                    throw new Error(`Failed to delete batch: ${response.statusText}`);
                }

                const result = await response.json();
                this.showSuccessMessage(result.message);
                
                // Refresh batch list
                this.loadBatches();
                
            } catch (error) {
                console.error('Error deleting batch:', error);
                this.showErrorMessage('Failed to delete batch: ' + error.message);
            }
        }
    }

    async stopSingleExecution(executionId) {
        if (confirm(`Are you sure you want to stop execution "${executionId}"?`)) {
            try {
                console.log('Stopping execution:', executionId);
                const response = await fetch(`/execution/api/execution/sessions/${executionId}/cancel`, {
                    method: 'POST'
                });

                if (!response.ok) {
                    throw new Error(`Failed to stop execution: ${response.statusText}`);
                }

                const result = await response.json();
                this.showSuccessMessage(result.message);
                
                // Refresh executions
                setTimeout(() => this.loadExecutions(), 1000);
                
            } catch (error) {
                console.error('Error stopping execution:', error);
                this.showErrorMessage('Failed to stop execution: ' + error.message);
            }
        }
    }

    async deleteSingleExecution(executionId) {
        // Note: Backend doesn't have delete execution endpoint yet
        this.showErrorMessage('Delete execution functionality not yet implemented in backend');
    }

    async runSelectedBatches() {
        if (this.selectedBatches.size === 0) {
            alert('Please select at least one batch to run');
            return;
        }
        
        try {
            const batchIds = Array.from(this.selectedBatches);
            console.log('Running selected batches:', batchIds);
            
            // Run first batch and navigate to its progress
            if (batchIds.length === 1) {
                await this.runSingleBatch(batchIds[0]);
                return;
            }
            
            // For multiple batches, run the first one and navigate to it
            const firstBatchId = batchIds[0];
            const response = await fetch(`/execution/api/batches/${firstBatchId}/execute`, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json'
                }
            });

            if (!response.ok) {
                throw new Error(`Failed to start batch execution: ${response.statusText}`);
            }

            const result = await response.json();
            
            // Start remaining batches in background
            for (let i = 1; i < batchIds.length; i++) {
                try {
                    await fetch(`/execution/api/batches/${batchIds[i]}/execute`, {
                        method: 'POST',
                        headers: { 'Content-Type': 'application/json' }
                    });
                } catch (error) {
                    console.error(`Failed to start batch ${batchIds[i]}:`, error);
                }
            }
            
            this.showSuccessMessage(`Started ${batchIds.length} batch executions. Navigating to first execution...`);
            
            // Clear selection and navigate to first execution
            this.clearBatchSelection();
            
            setTimeout(() => {
                window.location.href = `/execution/progress/${result.session_id}`;
            }, 1500);
            
        } catch (error) {
            console.error('Error running selected batches:', error);
            this.showErrorMessage('Failed to run selected batches: ' + error.message);
        }
    }

    async deleteSelectedBatches() {
        if (this.selectedBatches.size === 0) {
            alert('Please select at least one batch to delete');
            return;
        }
        
        if (confirm(`Are you sure you want to delete ${this.selectedBatches.size} selected batches?`)) {
            try {
                const batchIds = Array.from(this.selectedBatches);
                console.log('Deleting selected batches:', batchIds);
                
                // Delete batches sequentially
                for (const batchId of batchIds) {
                    const response = await fetch(`/execution/api/batches/${batchId}`, {
                        method: 'DELETE'
                    });
                    
                    if (!response.ok) {
                        console.error(`Failed to delete ${batchId}: ${response.statusText}`);
                    }
                }
                
                this.showSuccessMessage(`${batchIds.length} batch files deleted successfully`);
                
                // Clear selection and refresh
                this.clearBatchSelection();
                this.loadBatches();
                
            } catch (error) {
                console.error('Error deleting selected batches:', error);
                this.showErrorMessage('Failed to delete selected batches: ' + error.message);
            }
        }
    }

    async stopSelectedExecutions() {
        if (this.selectedExecutions.size === 0) {
            alert('Please select at least one execution to stop');
            return;
        }
        
        if (confirm(`Are you sure you want to stop ${this.selectedExecutions.size} selected executions?`)) {
            try {
                const executionIds = Array.from(this.selectedExecutions);
                console.log('Stopping selected executions:', executionIds);
                
                // Stop executions sequentially
                for (const executionId of executionIds) {
                    const response = await fetch(`/execution/api/execution/sessions/${executionId}/cancel`, {
                        method: 'POST'
                    });
                    
                    if (!response.ok) {
                        console.error(`Failed to stop ${executionId}: ${response.statusText}`);
                    }
                }
                
                this.showSuccessMessage(`${executionIds.length} executions stopped successfully`);
                
                // Clear selection and refresh
                this.clearExecutionSelection();
                setTimeout(() => this.loadExecutions(), 1000);
                
            } catch (error) {
                console.error('Error stopping selected executions:', error);
                this.showErrorMessage('Failed to stop selected executions: ' + error.message);
            }
        }
    }

    selectAllBatches() {
        this.selectedBatches = new Set(this.filteredBatches.map(b => b.filename));
        this.renderBatches();
        this.updateBatchSelectionUI();
    }

    clearBatchSelection() {
        this.selectedBatches.clear();
        this.renderBatches();
        this.updateBatchSelectionUI();
    }

    selectAllExecutions() {
        this.selectedExecutions = new Set(this.filteredExecutions.map(e => e.id));
        this.renderExecutions();
        this.updateExecutionSelectionUI();
    }

    clearExecutionSelection() {
        this.selectedExecutions.clear();
        this.renderExecutions();
        this.updateExecutionSelectionUI();
    }

    toggleFiltersIcon(type) {
        const toggleBtn = document.getElementById(`toggle-${type}`);
        if (toggleBtn) {
            const icon = toggleBtn.querySelector('i');
            if (icon) {
                const section = document.getElementById(`${type}-section`);
                const isExpanded = section?.classList.contains('show');
                icon.className = isExpanded ? 'bi bi-funnel-fill me-1' : 'bi bi-funnel me-1';
            }
        }
    }

    showUploadModal() {
        // Reset form when opening modal
        const fileInput = document.getElementById('batch-file');
        const nameInput = document.getElementById('batch-name');
        if (fileInput) fileInput.value = '';
        if (nameInput) nameInput.value = '';
    }

    async viewBatchDetails(batchId) {
        try {
            const response = await fetch(`/execution/api/batches/${batchId}`);
            if (!response.ok) {
                throw new Error(`Failed to load batch details: ${response.statusText}`);
            }
            
            const result = await response.json();
            
            // Show details in a modal or navigate to details page
            console.log('Batch details:', result);
            this.showSuccessMessage(`Loaded details for batch: ${batchId}`);
            
        } catch (error) {
            console.error('Error loading batch details:', error);
            this.showErrorMessage('Failed to load batch details: ' + error.message);
        }
    }

    async downloadBatch(batchId) {
        try {
            const response = await fetch(`/execution/api/batches/${batchId}/download`);
            if (!response.ok) {
                throw new Error(`Failed to download batch: ${response.statusText}`);
            }
            
            // Create download link
            const blob = await response.blob();
            const url = window.URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = batchId;
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
            window.URL.revokeObjectURL(url);
            
            this.showSuccessMessage(`Downloaded batch file: ${batchId}`);
            
        } catch (error) {
            console.error('Error downloading batch:', error);
            this.showErrorMessage('Failed to download batch: ' + error.message);
        }
    }

    async downloadExecutionResults(executionId) {
        try {
            // This endpoint might not exist yet in the backend
            this.showErrorMessage('Download results functionality not yet implemented in backend');
            
        } catch (error) {
            console.error('Error downloading execution results:', error);
            this.showErrorMessage('Failed to download execution results: ' + error.message);
        }
    }
}

// Make available globally
window.BatchManagerTabbed = BatchManagerTabbed;
