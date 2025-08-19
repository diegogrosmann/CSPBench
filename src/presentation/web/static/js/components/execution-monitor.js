/**
 * Real-time Execution Monitor Component
 * Displays active executions with real-time updates
 */

class ExecutionMonitor {
    constructor(containerId) {
        this.containerId = containerId;
        this.executions = new Map(); // work_id -> execution data
        this.updateInterval = null;
        this.selectedWorkId = null;
        this.currentFilter = 'all'; // Track current filter
        
        // Initialize API client and WebSocket client
        this.apiClient = null;
        this.wsClient = null;
        
        // Bind methods
        this.initialize = this.initialize.bind(this);
        this.loadActiveExecutions = this.loadActiveExecutions.bind(this);
        this.renderExecutions = this.renderExecutions.bind(this);
        this.renderExecutionTable = this.renderExecutionTable.bind(this);
        this.handleWebSocketMessage = this.handleWebSocketMessage.bind(this);
        this.showExecutionDetails = this.showExecutionDetails.bind(this);
        this.showExecutionResults = this.showExecutionResults.bind(this);
        this.showLogs = this.showLogs.bind(this);
        this.showConnectionStatus = this.showConnectionStatus.bind(this);
        this.showError = this.showError.bind(this);
        this.showModal = this.showModal.bind(this);
        this.renderExecutionCard = this.renderExecutionCard.bind(this);
        this.updateExecutionCard = this.updateExecutionCard.bind(this);
        this.renderExecutionDetailsModal = this.renderExecutionDetailsModal.bind(this);
        this.renderExecutionResultsModal = this.renderExecutionResultsModal.bind(this);
        this.getStatusColor = this.getStatusColor.bind(this);
        this.getStatusIcon = this.getStatusIcon.bind(this);
        this.getTimeAgo = this.getTimeAgo.bind(this);
        this.handleExecutionClick = this.handleExecutionClick.bind(this);
        this.handleFilterChange = this.handleFilterChange.bind(this);
        this.filterExecutions = this.filterExecutions.bind(this);
        this.isActiveStatus = this.isActiveStatus.bind(this);
        this.downloadFile = this.downloadFile.bind(this);
    }
    
    async initialize() {
        console.log('Initializing Execution Monitor...');
        
        try {
            // Initialize API client and WebSocket client
            this.apiClient = new APIClient();
            this.wsClient = new WebSocketClient();
            
            // Setup WebSocket connection
            this.wsClient.on('connected', () => {
                console.log('Execution Monitor: WebSocket connected');
                this.showConnectionStatus('connected');
            });
            
            this.wsClient.on('disconnected', () => {
                console.log('Execution Monitor: WebSocket disconnected');
                this.showConnectionStatus('disconnected');
            });
            
            this.wsClient.on('work_update', this.handleWebSocketMessage);
            
            await this.wsClient.connect();
            
            // Initial load
            await this.loadActiveExecutions();
            
            // Setup filter change listener
            const filterSelect = document.getElementById('status-filter');
            if (filterSelect) {
                filterSelect.addEventListener('change', this.handleFilterChange);
            }
            
            // Setup periodic refresh (fallback for WebSocket)
            this.updateInterval = setInterval(() => {
                this.loadActiveExecutions();
            }, 30000); // Every 30 seconds
            
            console.log('Execution Monitor initialized successfully');
            
        } catch (error) {
            console.error('Error initializing Execution Monitor:', error);
            this.showError('Failed to initialize execution monitor: ' + error.message);
        }
    }
    
    async loadActiveExecutions() {
        try {
            console.log('Loading active executions...');
            
            const response = await this.apiClient.request('/api/monitoring/active-executions');
            const data = response.executions || [];
            
            // Update executions map
            this.executions.clear();
            data.forEach(execution => {
                this.executions.set(execution.work_id, execution);
                
                // Subscribe to updates for running executions
                if (execution.status === 'running' && this.wsClient.isConnected()) {
                    this.wsClient.subscribeToWork(execution.work_id);
                }
            });
            
            this.renderExecutions();
            
        } catch (error) {
            console.error('Error loading active executions:', error);
            this.showError('Failed to load executions: ' + error.message);
        }
    }
    
    renderExecutions() {
        const container = document.getElementById(this.containerId);
        if (!container) {
            console.error('Execution monitor container not found:', this.containerId);
            return;
        }
        
        const allExecutions = Array.from(this.executions.values());
        const filteredExecutions = this.filterExecutions(allExecutions);
        
        if (allExecutions.length === 0) {
            container.innerHTML = `
                <div class="no-executions">
                    <i class="bi bi-inbox"></i>
                    <h5>No executions found</h5>
                    <p class="text-muted">Execute a batch to see monitoring data here</p>
                </div>
            `;
            return;
        }
        
        container.innerHTML = this.renderExecutionTable(filteredExecutions, allExecutions.length);
    }
    
    renderExecutionTable(executions, totalCount) {
        const filterSelect = document.getElementById('status-filter');
        const currentFilter = filterSelect ? filterSelect.value : 'all';
        
        if (executions.length === 0) {
            return `
                <div class="no-executions">
                    <i class="bi bi-filter"></i>
                    <h5>No executions match the current filter</h5>
                    <p class="text-muted">Try changing the filter or check back later</p>
                    <p class="small text-muted">Total executions: ${totalCount}</p>
                </div>
            `;
        }
        
        const tableRows = executions.map(execution => {
            const statusColor = this.getStatusColor(execution.status);
            const statusIcon = this.getStatusIcon(execution.status);
            const createdAt = new Date(execution.created_at);
            const timeAgo = this.getTimeAgo(createdAt);
            const isActive = this.isActiveStatus(execution.status);
            
            return `
                <tr data-work-id="${execution.work_id}" 
                    onclick="executionMonitor.handleExecutionClick('${execution.work_id}', '${execution.status}')"
                    title="Click to ${isActive ? 'monitor execution' : 'view results'}">
                    <td class="work-id-cell">${execution.work_id}</td>
                    <td>
                        <span class="status-badge ${execution.status}">
                            <i class="bi bi-${statusIcon} status-icon"></i>
                            ${execution.status.charAt(0).toUpperCase() + execution.status.slice(1)}
                        </span>
                    </td>
                    <td class="time-cell" title="${createdAt.toLocaleString()}">
                        ${createdAt.toLocaleString()}
                    </td>
                    <td class="time-cell" title="${new Date(execution.last_modified).toLocaleString()}">
                        ${new Date(execution.last_modified).toLocaleString()}
                    </td>
                    <td>
                        <div class="d-flex gap-1">
                            <button class="btn btn-outline-primary btn-sm" 
                                    onclick="event.stopPropagation(); window.open('/execution/${execution.work_id}', '_blank')"
                                    title="Detailed monitoring">
                                <i class="bi bi-graph-up"></i>
                            </button>
                            
                            ${isActive ? `
                                <button class="btn btn-outline-info btn-sm" 
                                        onclick="event.stopPropagation(); executionMonitor.showExecutionDetails('${execution.work_id}')"
                                        title="Quick view">
                                    <i class="bi bi-eye"></i>
                                </button>
                            ` : `
                                <button class="btn btn-outline-success btn-sm" 
                                        onclick="event.stopPropagation(); executionMonitor.showExecutionResults('${execution.work_id}')"
                                        title="View results">
                                    <i class="bi bi-file-earmark-text"></i>
                                </button>
                            `}
                            
                            ${execution.log_files && execution.log_files.length > 0 ? `
                                <button class="btn btn-outline-secondary btn-sm" 
                                        onclick="event.stopPropagation(); executionMonitor.showLogs('${execution.work_id}')"
                                        title="View logs">
                                    <i class="bi bi-file-text"></i>
                                </button>
                            ` : ''}
                        </div>
                    </td>
                </tr>
            `;
        }).join('');
        
        return `
            <table class="table executions-table">
                <thead>
                    <tr>
                        <th>Work ID</th>
                        <th>Status</th>
                        <th>Started</th>
                        <th>Last Modified</th>
                        <th>Actions</th>
                    </tr>
                </thead>
                <tbody>
                    ${tableRows}
                </tbody>
            </table>
            <div class="d-flex justify-content-between align-items-center px-3 py-2 bg-light">
                <small class="text-muted">
                    Showing ${executions.length} of ${totalCount} executions
                    ${currentFilter !== 'all' ? `(filtered by: ${currentFilter})` : ''}
                </small>
                <small class="text-muted">
                    <i class="bi bi-info-circle me-1"></i>
                    Click row to ${this.isActiveStatus('running') ? 'monitor' : 'view results'}
                </small>
            </div>
        `;
    }
    
    filterExecutions(executions) {
        const filterSelect = document.getElementById('status-filter');
        const filter = filterSelect ? filterSelect.value : this.currentFilter;
        
        switch (filter) {
            case 'active':
                return executions.filter(ex => this.isActiveStatus(ex.status));
            case 'running':
                return executions.filter(ex => ex.status === 'running');
            case 'finished':
                return executions.filter(ex => ex.status === 'finished');
            case 'failed':
                return executions.filter(ex => ex.status === 'failed');
            default:
                return executions;
        }
    }
    
    isActiveStatus(status) {
        return ['running', 'paused', 'starting'].includes(status);
    }
    
    handleFilterChange(event) {
        this.currentFilter = event.target.value;
        this.renderExecutions();
    }
    
    handleExecutionClick(workId, status) {
        if (this.isActiveStatus(status)) {
            // Navigate to monitoring view
            this.showExecutionDetails(workId);
        } else {
            // Navigate to results view
            this.showExecutionResults(workId);
        }
    }
    
    renderExecutionCard(execution) {
        const statusColor = this.getStatusColor(execution.status);
        const statusIcon = this.getStatusIcon(execution.status);
        
        const createdAt = new Date(execution.created_at);
        const timeAgo = this.getTimeAgo(createdAt);
        
        return `
            <div class="col-md-6 col-lg-4 mb-3">
                <div class="card h-100 execution-card" data-work-id="${execution.work_id}">
                    <div class="card-header d-flex justify-content-between align-items-center">
                        <small class="text-muted">${execution.work_id}</small>
                        <span class="badge bg-${statusColor}">
                            <i class="bi bi-${statusIcon} me-1"></i>${execution.status}
                        </span>
                    </div>
                    <div class="card-body">
                        <div class="mb-2">
                            <small class="text-muted">Started: ${timeAgo}</small>
                        </div>
                        
                        ${execution.status === 'running' ? `
                            <div class="mb-2">
                                <div class="progress" style="height: 6px;">
                                    <div class="progress-bar progress-bar-striped progress-bar-animated" 
                                         style="width: 100%"></div>
                                </div>
                            </div>
                        ` : ''}
                        
                        <div class="d-flex justify-content-between align-items-center">
                            <button class="btn btn-outline-primary btn-sm" 
                                    onclick="executionMonitor.showExecutionDetails('${execution.work_id}')">
                                <i class="bi bi-eye me-1"></i>Details
                            </button>
                            
                            ${execution.log_files && execution.log_files.length > 0 ? `
                                <button class="btn btn-outline-info btn-sm" 
                                        onclick="executionMonitor.showLogs('${execution.work_id}')">
                                    <i class="bi bi-file-text me-1"></i>Logs
                                </button>
                            ` : ''}
                        </div>
                    </div>
                </div>
            </div>
        `;
    }
    
    getStatusColor(status) {
        const colors = {
            'running': 'primary',
            'finished': 'success',
            'failed': 'danger',
            'cancelled': 'secondary',
            'unknown': 'warning'
        };
        return colors[status] || 'secondary';
    }
    
    getStatusIcon(status) {
        const icons = {
            'running': 'play-circle',
            'finished': 'check-circle',
            'failed': 'x-circle',
            'cancelled': 'stop-circle',
            'unknown': 'question-circle'
        };
        return icons[status] || 'circle';
    }
    
    getTimeAgo(date) {
        const now = new Date();
        const diffMs = now - date;
        const diffMins = Math.floor(diffMs / 60000);
        const diffHours = Math.floor(diffMins / 60);
        const diffDays = Math.floor(diffHours / 24);
        
        if (diffMins < 1) return 'Just now';
        if (diffMins < 60) return `${diffMins}m ago`;
        if (diffHours < 24) return `${diffHours}h ago`;
        return `${diffDays}d ago`;
    }
    
    handleWebSocketMessage(data) {
        console.log('Received WebSocket message:', data);
        
        if (data.work_id) {
            // Update execution in map
            this.executions.set(data.work_id, data);
            
            // Update display
            this.updateExecutionCard(data.work_id, data);
        }
    }
    
    updateExecutionCard(workId, execution) {
        const card = document.querySelector(`[data-work-id="${workId}"]`);
        if (card) {
            // Re-render executions to update the display
            this.renderExecutions();
        }
    }
    
    async showExecutionDetails(workId) {
        console.log('Showing execution details for:', workId);
        this.selectedWorkId = workId;
        
        try {
            // Load detailed execution data
            const execution = await this.apiClient.request(`/api/monitoring/execution/${workId}`);
            
            // Show modal with details
            this.showModal('Execution Monitoring', this.renderExecutionDetailsModal(execution));
            
        } catch (error) {
            console.error('Error loading execution details:', error);
            this.showError('Failed to load execution details: ' + error.message);
        }
    }
    
    async showExecutionResults(workId) {
        console.log('Showing execution results for:', workId);
        this.selectedWorkId = workId;
        
        try {
            // Load detailed execution data
            const execution = await this.apiClient.request(`/api/monitoring/execution/${workId}`);
            
            // Show modal with results and download options
            this.showModal('Execution Results', this.renderExecutionResultsModal(execution));
            
        } catch (error) {
            console.error('Error loading execution results:', error);
            this.showError('Failed to load execution results: ' + error.message);
        }
    }
    
    renderExecutionDetailsModal(execution) {
        const metrics = execution.metrics || {};
        
        return `
            <div class="execution-details">
                <div class="row">
                    <div class="col-md-6">
                        <h6>Basic Information</h6>
                        <table class="table table-sm">
                            <tr>
                                <td><strong>Work ID:</strong></td>
                                <td><code>${execution.work_id}</code></td>
                            </tr>
                            <tr>
                                <td><strong>Status:</strong></td>
                                <td>
                                    <span class="badge bg-${this.getStatusColor(execution.status)}">
                                        ${execution.status}
                                    </span>
                                </td>
                            </tr>
                            <tr>
                                <td><strong>Created:</strong></td>
                                <td>${new Date(execution.created_at).toLocaleString()}</td>
                            </tr>
                            <tr>
                                <td><strong>Duration:</strong></td>
                                <td>${metrics.duration ? Math.round(metrics.duration) + 's' : 'N/A'}</td>
                            </tr>
                        </table>
                    </div>
                    <div class="col-md-6">
                        <h6>Files & Output</h6>
                        <table class="table table-sm">
                            <tr>
                                <td><strong>File Count:</strong></td>
                                <td>${metrics.file_count || 0}</td>
                            </tr>
                            <tr>
                                <td><strong>Total Size:</strong></td>
                                <td>${metrics.total_size_mb || 0} MB</td>
                            </tr>
                            <tr>
                                <td><strong>Output Path:</strong></td>
                                <td><code>${execution.output_path}</code></td>
                            </tr>
                        </table>
                    </div>
                </div>
                
                <div class="mt-3">
                    <h6>Recent Log Output</h6>
                    <pre class="bg-dark text-light p-3 rounded" style="max-height: 300px; overflow-y: auto;"><code>${execution.latest_log || 'No log content available'}</code></pre>
                </div>
                
                <div class="mt-3 d-flex gap-2">
                    <button class="btn btn-outline-info" onclick="executionMonitor.showLogs('${execution.work_id}')">
                        <i class="bi bi-file-text me-1"></i>View Full Logs
                    </button>
                    <button class="btn btn-outline-primary" onclick="executionMonitor.loadActiveExecutions()">
                        <i class="bi bi-arrow-clockwise me-1"></i>Refresh
                    </button>
                </div>
            </div>
        `;
    }
    
    renderExecutionResultsModal(execution) {
        const metrics = execution.metrics || {};
        const files = execution.files || {};
        
        return `
            <div class="execution-results">
                <div class="row">
                    <div class="col-md-6">
                        <h6>Execution Summary</h6>
                        <table class="table table-sm">
                            <tr>
                                <td><strong>Work ID:</strong></td>
                                <td><code>${execution.work_id}</code></td>
                            </tr>
                            <tr>
                                <td><strong>Final Status:</strong></td>
                                <td>
                                    <span class="badge bg-${this.getStatusColor(execution.status)}">
                                        <i class="bi bi-${this.getStatusIcon(execution.status)} me-1"></i>
                                        ${execution.status.charAt(0).toUpperCase() + execution.status.slice(1)}
                                    </span>
                                </td>
                            </tr>
                            <tr>
                                <td><strong>Completed:</strong></td>
                                <td>${new Date(execution.last_modified).toLocaleString()}</td>
                            </tr>
                            <tr>
                                <td><strong>Duration:</strong></td>
                                <td>${metrics.duration ? Math.round(metrics.duration) + 's' : 'N/A'}</td>
                            </tr>
                        </table>
                    </div>
                    <div class="col-md-6">
                        <h6>Output Statistics</h6>
                        <table class="table table-sm">
                            <tr>
                                <td><strong>Files Generated:</strong></td>
                                <td>${metrics.file_count || 0}</td>
                            </tr>
                            <tr>
                                <td><strong>Total Size:</strong></td>
                                <td>${metrics.total_size_mb || 0} MB</td>
                            </tr>
                            <tr>
                                <td><strong>Output Directory:</strong></td>
                                <td><code>${execution.output_path}</code></td>
                            </tr>
                        </table>
                    </div>
                </div>
                
                <div class="mt-3">
                    <h6>Output Files</h6>
                    <div class="row">
                        ${files.outputs && files.outputs.length > 0 ? `
                            <div class="col-md-6">
                                <h6 class="small">Result Files</h6>
                                <div class="list-group list-group-flush">
                                    ${files.outputs.map(file => `
                                        <div class="list-group-item d-flex justify-content-between align-items-center">
                                            <span class="font-monospace small">${file.split('/').pop()}</span>
                                            <button class="btn btn-outline-primary btn-sm" 
                                                    onclick="executionMonitor.downloadFile('${file}')"
                                                    title="Download file">
                                                <i class="bi bi-download"></i>
                                            </button>
                                        </div>
                                    `).join('')}
                                </div>
                            </div>
                        ` : ''}
                        
                        ${files.logs && files.logs.length > 0 ? `
                            <div class="col-md-6">
                                <h6 class="small">Log Files</h6>
                                <div class="list-group list-group-flush">
                                    ${files.logs.map(file => `
                                        <div class="list-group-item d-flex justify-content-between align-items-center">
                                            <span class="font-monospace small">${file.split('/').pop()}</span>
                                            <div class="btn-group" role="group">
                                                <button class="btn btn-outline-info btn-sm" 
                                                        onclick="executionMonitor.showLogs('${execution.work_id}')"
                                                        title="View logs">
                                                    <i class="bi bi-eye"></i>
                                                </button>
                                                <button class="btn btn-outline-primary btn-sm" 
                                                        onclick="executionMonitor.downloadFile('${file}')"
                                                        title="Download file">
                                                    <i class="bi bi-download"></i>
                                                </button>
                                            </div>
                                        </div>
                                    `).join('')}
                                </div>
                            </div>
                        ` : ''}
                        
                        ${files.all && files.all.length > 0 && (!files.outputs || files.outputs.length === 0) && (!files.logs || files.logs.length === 0) ? `
                            <div class="col-12">
                                <h6 class="small">All Files</h6>
                                <div class="list-group list-group-flush">
                                    ${files.all.map(file => `
                                        <div class="list-group-item d-flex justify-content-between align-items-center">
                                            <span class="font-monospace small">${file.split('/').pop()}</span>
                                            <button class="btn btn-outline-primary btn-sm" 
                                                    onclick="executionMonitor.downloadFile('${file}')"
                                                    title="Download file">
                                                <i class="bi bi-download"></i>
                                            </button>
                                        </div>
                                    `).join('')}
                                </div>
                            </div>
                        ` : ''}
                    </div>
                    
                    ${(!files.outputs || files.outputs.length === 0) && (!files.logs || files.logs.length === 0) && (!files.all || files.all.length === 0) ? `
                        <div class="text-center text-muted py-3">
                            <i class="bi bi-file-earmark-x display-4"></i>
                            <p class="mt-2">No output files available</p>
                        </div>
                    ` : ''}
                </div>
                
                <div class="mt-3 d-flex gap-2">
                    <button class="btn btn-outline-success" onclick="executionMonitor.downloadAllFiles('${execution.work_id}')">
                        <i class="bi bi-download me-1"></i>Download All Files
                    </button>
                    <button class="btn btn-outline-info" onclick="executionMonitor.showLogs('${execution.work_id}')">
                        <i class="bi bi-file-text me-1"></i>View Logs
                    </button>
                </div>
            </div>
        `;
    }
    
    async downloadFile(filePath) {
        try {
            // Create a temporary link and trigger download
            const response = await fetch(`/api/files/download?path=${encodeURIComponent(filePath)}`);
            
            if (!response.ok) {
                throw new Error('Failed to download file');
            }
            
            const blob = await response.blob();
            const fileName = filePath.split('/').pop();
            
            const url = window.URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = fileName;
            document.body.appendChild(a);
            a.click();
            window.URL.revokeObjectURL(url);
            document.body.removeChild(a);
            
        } catch (error) {
            console.error('Error downloading file:', error);
            this.showError('Failed to download file: ' + error.message);
        }
    }
    
    async downloadAllFiles(workId) {
        try {
            // Request a zip of all files for this execution
            const response = await fetch(`/api/files/download-zip/${workId}`);
            
            if (!response.ok) {
                throw new Error('Failed to create download archive');
            }
            
            const blob = await response.blob();
            const fileName = `execution_${workId}_results.zip`;
            
            const url = window.URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = fileName;
            document.body.appendChild(a);
            a.click();
            window.URL.revokeObjectURL(url);
            document.body.removeChild(a);
            
        } catch (error) {
            console.error('Error downloading all files:', error);
            this.showError('Failed to download files: ' + error.message);
        }
    }
    
    async showLogs(workId) {
        console.log('Showing logs for execution:', workId);
        
        try {
            const response = await this.apiClient.request(`/api/monitoring/execution/${workId}/logs?lines=200`);
            
            this.showModal('Execution Logs', `
                <div class="logs-viewer">
                    <div class="d-flex justify-content-between align-items-center mb-3">
                        <h6 class="mb-0">Work ID: <code>${workId}</code></h6>
                        <span class="text-muted">Showing last ${response.showing_lines} lines</span>
                    </div>
                    <pre class="bg-dark text-light p-3 rounded" style="height: 400px; overflow-y: auto;"><code>${response.logs || 'No log content available'}</code></pre>
                    <div class="mt-2">
                        <button class="btn btn-outline-primary btn-sm" onclick="executionMonitor.showLogs('${workId}')">
                            <i class="bi bi-arrow-clockwise me-1"></i>Refresh
                        </button>
                    </div>
                </div>
            `);
            
        } catch (error) {
            console.error('Error loading logs:', error);
            this.showError('Failed to load logs: ' + error.message);
        }
    }
    
    showConnectionStatus(status) {
        const statusElement = document.getElementById('ws-status');
        if (statusElement) {
            if (status === 'connected') {
                statusElement.className = 'badge bg-success me-2';
                statusElement.innerHTML = '<i class="bi bi-wifi me-1"></i>Live';
            } else {
                statusElement.className = 'badge bg-warning me-2';
                statusElement.innerHTML = '<i class="bi bi-wifi-off me-1"></i>Offline';
            }
        }
    }
    
    showModal(title, content) {
        // Create or update modal
        let modal = document.getElementById('execution-details-modal');
        if (!modal) {
            modal = document.createElement('div');
            modal.id = 'execution-details-modal';
            modal.className = 'modal fade';
            modal.innerHTML = `
                <div class="modal-dialog modal-lg">
                    <div class="modal-content">
                        <div class="modal-header">
                            <h5 class="modal-title"></h5>
                            <button type="button" class="btn-close" data-bs-dismiss="modal" aria-label="Close"></button>
                        </div>
                        <div class="modal-body"></div>
                    </div>
                </div>
            `;
            document.body.appendChild(modal);
        }
        
        modal.querySelector('.modal-title').textContent = title;
        modal.querySelector('.modal-body').innerHTML = content;
        
        // Ensure close button works
        const closeButton = modal.querySelector('.btn-close');
        if (closeButton) {
            closeButton.addEventListener('click', () => {
                const bsModal = bootstrap.Modal.getInstance(modal);
                if (bsModal) {
                    bsModal.hide();
                }
            });
        }
        
        const bsModal = new bootstrap.Modal(modal);
        bsModal.show();
    }
    
    showError(message) {
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
    
    destroy() {
        console.log('Destroying Execution Monitor...');
        
        if (this.updateInterval) {
            clearInterval(this.updateInterval);
        }
        
        if (this.wsClient) {
            this.wsClient.disconnect();
        }
    }
}

// Export for global use
if (typeof window !== 'undefined') {
    window.ExecutionMonitor = ExecutionMonitor;
}
