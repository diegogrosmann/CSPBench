/**
 * Execution Monitor Component
 * 
 * Real-time monitoring and progress tracking for algorithm execution
 */

class ExecutionMonitor {
    constructor(container, options = {}) {
        this.container = typeof container === 'string' ? document.querySelector(container) : container;
        this.options = {
            updateInterval: 1000,
            showLogs: true,
            showProgress: true,
            showMetrics: true,
            autoScroll: true,
            onComplete: () => {},
            onError: () => {},
            ...options
        };
        
        this.executionId = null;
        this.status = 'idle';
        this.startTime = null;
        this.progress = 0;
        this.logs = [];
        this.metrics = {};
        this.updateTimer = null;
        
        this.init();
    }

    init() {
        this.container.classList.add('execution-monitor');
        this.render();
    }

    render() {
        this.container.innerHTML = `
            <div class="execution-header">
                <div class="execution-status">
                    <span class="status-indicator ${this.status}" id="status-indicator"></span>
                    <span class="status-text" id="status-text">${this.getStatusText()}</span>
                </div>
                
                <div class="execution-time" id="execution-time">
                    ${this.getExecutionTime()}
                </div>
                
                <div class="execution-actions">
                    <button class="btn btn-danger" id="cancel-btn" ${this.status !== 'running' ? 'disabled' : ''}>
                        Cancel
                    </button>
                    <button class="btn btn-outline" id="download-logs" ${this.logs.length === 0 ? 'disabled' : ''}>
                        Download Logs
                    </button>
                </div>
            </div>
            
            ${this.options.showProgress ? this.renderProgress() : ''}
            ${this.options.showMetrics ? this.renderMetrics() : ''}
            ${this.options.showLogs ? this.renderLogs() : ''}
        `;
        
        this.setupEventListeners();
    }

    renderProgress() {
        return `
            <div class="progress-section">
                <div class="progress-header">
                    <h4>Progress</h4>
                    <span class="progress-percentage" id="progress-percentage">${this.progress}%</span>
                </div>
                
                <div class="progress-bar">
                    <div class="progress-fill" id="progress-fill" style="width: ${this.progress}%"></div>
                </div>
                
                <div class="progress-details">
                    <div class="progress-info" id="progress-info">
                        ${this.getProgressInfo()}
                    </div>
                </div>
            </div>
        `;
    }

    renderMetrics() {
        return `
            <div class="metrics-section">
                <h4>Real-time Metrics</h4>
                <div class="metrics-grid" id="metrics-grid">
                    ${this.renderMetricsGrid()}
                </div>
            </div>
        `;
    }

    renderMetricsGrid() {
        if (Object.keys(this.metrics).length === 0) {
            return '<div class="no-metrics">No metrics available yet...</div>';
        }

        return Object.entries(this.metrics).map(([key, value]) => `
            <div class="metric-card">
                <div class="metric-label">${this.formatMetricName(key)}</div>
                <div class="metric-value">${this.formatMetricValue(value)}</div>
            </div>
        `).join('');
    }

    renderLogs() {
        return `
            <div class="logs-section">
                <div class="logs-header">
                    <h4>Execution Logs</h4>
                    <div class="logs-controls">
                        <label class="checkbox-label">
                            <input type="checkbox" id="auto-scroll" ${this.options.autoScroll ? 'checked' : ''} />
                            Auto-scroll
                        </label>
                        <button class="btn btn-outline btn-small" id="clear-logs" ${this.logs.length === 0 ? 'disabled' : ''}>
                            Clear
                        </button>
                    </div>
                </div>
                
                <div class="logs-container" id="logs-container">
                    ${this.renderLogEntries()}
                </div>
            </div>
        `;
    }

    renderLogEntries() {
        if (this.logs.length === 0) {
            return '<div class="no-logs">No logs yet...</div>';
        }

        return this.logs.map(log => `
            <div class="log-entry ${log.level}" data-timestamp="${log.timestamp}">
                <span class="log-time">${this.formatLogTime(log.timestamp)}</span>
                <span class="log-level">${log.level.toUpperCase()}</span>
                <span class="log-message">${this.escapeHtml(log.message)}</span>
            </div>
        `).join('');
    }

    setupEventListeners() {
        // Cancel button
        document.getElementById('cancel-btn')?.addEventListener('click', () => {
            this.cancelExecution();
        });

        // Download logs
        document.getElementById('download-logs')?.addEventListener('click', () => {
            this.downloadLogs();
        });

        // Auto-scroll toggle
        document.getElementById('auto-scroll')?.addEventListener('change', (e) => {
            this.options.autoScroll = e.target.checked;
        });

        // Clear logs
        document.getElementById('clear-logs')?.addEventListener('click', () => {
            this.clearLogs();
        });
    }

    startExecution(executionData) {
        this.executionId = executionData.id || Date.now().toString();
        this.status = 'running';
        this.startTime = new Date();
        this.progress = 0;
        this.logs = [];
        this.metrics = {};
        
        this.addLog('info', `Starting execution of ${executionData.algorithm} on ${executionData.dataset}`);
        this.addLog('info', `Parameters: ${JSON.stringify(executionData.params, null, 2)}`);
        
        this.render();
        this.startPolling();
    }

    startPolling() {
        if (this.updateTimer) {
            clearInterval(this.updateTimer);
        }
        
        this.updateTimer = setInterval(() => {
            this.updateStatus();
        }, this.options.updateInterval);
    }

    stopPolling() {
        if (this.updateTimer) {
            clearInterval(this.updateTimer);
            this.updateTimer = null;
        }
    }

    async updateStatus() {
        if (!this.executionId || this.status !== 'running') {
            this.stopPolling();
            return;
        }

        try {
            const response = await window.apiClient.getExecutionStatus(this.executionId);
            this.handleStatusUpdate(response);
        } catch (error) {
            console.error('Failed to update execution status:', error);
            this.addLog('error', `Failed to get status update: ${error.message}`);
        }
    }

    handleStatusUpdate(data) {
        // Update progress
        if (data.progress !== undefined) {
            this.updateProgress(data.progress, data.progressMessage);
        }

        // Update metrics
        if (data.metrics) {
            this.updateMetrics(data.metrics);
        }

        // Add new logs
        if (data.logs && Array.isArray(data.logs)) {
            data.logs.forEach(log => {
                this.addLog(log.level, log.message, log.timestamp);
            });
        }

        // Check if completed
        if (data.status === 'completed') {
            this.handleCompletion(data.result);
        } else if (data.status === 'failed') {
            this.handleError(data.error);
        }
    }

    updateProgress(progress, message = '') {
        this.progress = Math.max(0, Math.min(100, progress));
        
        const progressFill = document.getElementById('progress-fill');
        const progressPercentage = document.getElementById('progress-percentage');
        const progressInfo = document.getElementById('progress-info');
        
        if (progressFill) progressFill.style.width = `${this.progress}%`;
        if (progressPercentage) progressPercentage.textContent = `${this.progress}%`;
        if (progressInfo && message) progressInfo.textContent = message;
    }

    updateMetrics(newMetrics) {
        this.metrics = { ...this.metrics, ...newMetrics };
        
        const metricsGrid = document.getElementById('metrics-grid');
        if (metricsGrid) {
            metricsGrid.innerHTML = this.renderMetricsGrid();
        }
    }

    addLog(level, message, timestamp = null) {
        const logEntry = {
            level: level.toLowerCase(),
            message,
            timestamp: timestamp || new Date().toISOString()
        };
        
        this.logs.push(logEntry);
        
        // Update logs container
        const logsContainer = document.getElementById('logs-container');
        if (logsContainer) {
            const logElement = document.createElement('div');
            logElement.className = `log-entry ${logEntry.level}`;
            logElement.dataset.timestamp = logEntry.timestamp;
            logElement.innerHTML = `
                <span class="log-time">${this.formatLogTime(logEntry.timestamp)}</span>
                <span class="log-level">${logEntry.level.toUpperCase()}</span>
                <span class="log-message">${this.escapeHtml(logEntry.message)}</span>
            `;
            
            logsContainer.appendChild(logElement);
            
            // Auto-scroll if enabled
            if (this.options.autoScroll) {
                logElement.scrollIntoView({ behavior: 'smooth' });
            }
        }
        
        // Update clear button state
        const clearBtn = document.getElementById('clear-logs');
        if (clearBtn) clearBtn.disabled = false;
        
        // Update download button state
        const downloadBtn = document.getElementById('download-logs');
        if (downloadBtn) downloadBtn.disabled = false;
    }

    handleCompletion(result) {
        this.status = 'completed';
        this.progress = 100;
        this.stopPolling();
        
        this.addLog('success', 'Execution completed successfully!');
        if (result) {
            this.addLog('info', `Best string: ${result.best_string}`);
            this.addLog('info', `Max distance: ${result.max_distance}`);
            if (result.metadata) {
                this.addLog('info', `Metadata: ${JSON.stringify(result.metadata, null, 2)}`);
            }
        }
        
        this.updateStatusIndicator();
        this.options.onComplete(result);
    }

    handleError(error) {
        this.status = 'failed';
        this.stopPolling();
        
        this.addLog('error', `Execution failed: ${error}`);
        this.updateStatusIndicator();
        this.options.onError(error);
    }

    async cancelExecution() {
        if (!this.executionId || this.status !== 'running') return;
        
        try {
            await window.apiClient.cancelExecution(this.executionId);
            this.status = 'cancelled';
            this.stopPolling();
            this.addLog('warning', 'Execution cancelled by user');
            this.updateStatusIndicator();
        } catch (error) {
            this.addLog('error', `Failed to cancel execution: ${error.message}`);
        }
    }

    updateStatusIndicator() {
        const indicator = document.getElementById('status-indicator');
        const statusText = document.getElementById('status-text');
        const cancelBtn = document.getElementById('cancel-btn');
        
        if (indicator) indicator.className = `status-indicator ${this.status}`;
        if (statusText) statusText.textContent = this.getStatusText();
        if (cancelBtn) cancelBtn.disabled = this.status !== 'running';
        
        // Update execution time
        const timeEl = document.getElementById('execution-time');
        if (timeEl) timeEl.textContent = this.getExecutionTime();
    }

    getStatusText() {
        const statusTexts = {
            idle: 'Ready',
            running: 'Running',
            completed: 'Completed',
            failed: 'Failed',
            cancelled: 'Cancelled'
        };
        return statusTexts[this.status] || 'Unknown';
    }

    getExecutionTime() {
        if (!this.startTime) return 'Not started';
        
        const elapsed = Date.now() - this.startTime.getTime();
        const seconds = Math.floor(elapsed / 1000);
        const minutes = Math.floor(seconds / 60);
        const hours = Math.floor(minutes / 60);
        
        if (hours > 0) {
            return `${hours}h ${minutes % 60}m ${seconds % 60}s`;
        } else if (minutes > 0) {
            return `${minutes}m ${seconds % 60}s`;
        } else {
            return `${seconds}s`;
        }
    }

    getProgressInfo() {
        if (this.status === 'idle') return 'Waiting to start...';
        if (this.status === 'running') return `${this.progress}% complete`;
        if (this.status === 'completed') return 'Execution completed successfully';
        if (this.status === 'failed') return 'Execution failed';
        if (this.status === 'cancelled') return 'Execution cancelled';
        return '';
    }

    formatMetricName(name) {
        return name.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase());
    }

    formatMetricValue(value) {
        if (typeof value === 'number') {
            return value % 1 === 0 ? value.toString() : value.toFixed(3);
        }
        return String(value);
    }

    formatLogTime(timestamp) {
        return new Date(timestamp).toLocaleTimeString();
    }

    escapeHtml(text) {
        const div = document.createElement('div');
        div.textContent = text;
        return div.innerHTML;
    }

    downloadLogs() {
        if (this.logs.length === 0) return;
        
        const logText = this.logs.map(log => 
            `[${log.timestamp}] ${log.level.toUpperCase()}: ${log.message}`
        ).join('\n');
        
        const blob = new Blob([logText], { type: 'text/plain' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = `execution-logs-${this.executionId || 'unknown'}.txt`;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
    }

    clearLogs() {
        this.logs = [];
        const logsContainer = document.getElementById('logs-container');
        if (logsContainer) {
            logsContainer.innerHTML = '<div class="no-logs">No logs yet...</div>';
        }
        
        const clearBtn = document.getElementById('clear-logs');
        if (clearBtn) clearBtn.disabled = true;
        
        const downloadBtn = document.getElementById('download-logs');
        if (downloadBtn) downloadBtn.disabled = true;
    }

    reset() {
        this.stopPolling();
        this.executionId = null;
        this.status = 'idle';
        this.startTime = null;
        this.progress = 0;
        this.logs = [];
        this.metrics = {};
        this.render();
    }

    getExecutionData() {
        return {
            id: this.executionId,
            status: this.status,
            startTime: this.startTime,
            progress: this.progress,
            logs: this.logs,
            metrics: this.metrics
        };
    }
}

// Export para uso global
window.ExecutionMonitor = ExecutionMonitor;
