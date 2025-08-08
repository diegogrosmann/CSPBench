/**
 * Modern API Client Component for CSPBench Web Interface
 * Centralized HTTP client with modern fetch API, error handling, and loading states
 */

class APIClient {
    constructor(baseURL = '') {
        this.baseURL = baseURL;
        this.defaultHeaders = {
            'Content-Type': 'application/json',
            'Accept': 'application/json'
        };
        this.loadingStates = new Map();
    }

    /**
     * Get current loading state for a request ID
     */
    isLoading(requestId) {
        return this.loadingStates.get(requestId) || false;
    }

    /**
     * Set loading state for a request
     */
    setLoading(requestId, isLoading) {
        this.loadingStates.set(requestId, isLoading);
        this.dispatchLoadingEvent(requestId, isLoading);
    }

    /**
     * Dispatch loading state change event
     */
    dispatchLoadingEvent(requestId, isLoading) {
        const event = new CustomEvent('api:loading', {
            detail: { requestId, isLoading }
        });
        document.dispatchEvent(event);
    }

    /**
     * Make HTTP request with modern error handling and loading states
     */
    async request(endpoint, options = {}, requestId = null) {
        const url = `${this.baseURL}${endpoint}`;
        const finalRequestId = requestId || `${options.method || 'GET'}_${endpoint}`;
        
        try {
            this.setLoading(finalRequestId, true);

            const requestOptions = {
                headers: { ...this.defaultHeaders, ...options.headers },
                ...options
            };

            const response = await fetch(url, requestOptions);

            if (!response.ok) {
                throw new APIError(
                    `HTTP ${response.status}: ${response.statusText}`,
                    response.status,
                    endpoint
                );
            }

            const contentType = response.headers.get('content-type');
            let data;

            if (contentType && contentType.includes('application/json')) {
                data = await response.json();
            } else if (contentType && contentType.includes('text/')) {
                data = await response.text();
            } else {
                data = await response.blob();
            }

            this.dispatchSuccessEvent(finalRequestId, data);
            return data;

        } catch (error) {
            this.dispatchErrorEvent(finalRequestId, error);
            throw error;
        } finally {
            this.setLoading(finalRequestId, false);
        }
    }

    /**
     * Dispatch success event
     */
    dispatchSuccessEvent(requestId, data) {
        const event = new CustomEvent('api:success', {
            detail: { requestId, data }
        });
        document.dispatchEvent(event);
    }

    /**
     * Dispatch error event
     */
    dispatchErrorEvent(requestId, error) {
        const event = new CustomEvent('api:error', {
            detail: { requestId, error }
        });
        document.dispatchEvent(event);
    }

    // Convenience methods for common HTTP verbs

    // Convenience methods for common HTTP verbs

    /**
     * GET request
     */
    async get(endpoint, requestId = null) {
        return this.request(endpoint, { method: 'GET' }, requestId);
    }

    /**
     * POST request
     */
    async post(endpoint, data = null, requestId = null) {
        const options = {
            method: 'POST',
            body: data ? JSON.stringify(data) : null
        };
        return this.request(endpoint, options, requestId);
    }

    /**
     * DELETE request
     */
    async delete(endpoint, requestId = null) {
        return this.request(endpoint, { method: 'DELETE' }, requestId);
    }

    // Algorithm endpoints
    async getAlgorithms() {
        return this.get('/api/algorithms', 'algorithms');
    }

    async executeAlgorithm(request) {
        return this.post('/api/execute', request, 'execute');
    }

    // Dataset endpoints
    async getDatasets() {
        return this.get('/api/datasets', 'datasets');
    }

    async generateDataset(config) {
        return this.post('/api/generate-dataset', config, 'generate_dataset');
    }

    // Dataset Generator endpoints
    async generateSyntheticDataset(params) {
        return this.post('/api/datasets-simple/generate/synthetic', params, 'generate_synthetic');
    }

    async downloadNCBIDataset(params) {
        return this.post('/api/datasets/generate/ncbi', params, 'download_ncbi');
    }

    async getGenerationStatus(sessionId) {
        return this.get(`/api/generation/status/${sessionId}`, `status_${sessionId}`);
    }

    async getSavedDatasets() {
        const response = await this.get('/api/datasets/saved', 'saved_datasets');
        return response.saved_datasets || [];
    }

    async getDatasetDownloadUrl(datasetId) {
        // For saved datasets, return the download URL directly
        return `${this.baseURL}/api/datasets/download/${datasetId}`;
    }

    async deleteDataset(filename) {
        return this.request(`/api/datasets/saved/${filename}`, {
            method: 'DELETE'
        }, `delete_dataset_${filename}`);
    }

    async cancelGeneration(sessionId) {
        return this.post(`/api/generation/cancel/${sessionId}`, null, `cancel_${sessionId}`);
    }

    // Execution management
    async getExecutionStatus(sessionId) {
        return this.get(`/api/status/${sessionId}`, `status_${sessionId}`);
    }

    async getExecutionResults(sessionId) {
        return this.get(`/api/results/${sessionId}`, `results_${sessionId}`);
    }

    async cancelExecution(sessionId) {
        return this.post(`/api/cancel/${sessionId}`, null, `cancel_${sessionId}`);
    }

    /**
     * Download results as blob
     */
    async downloadResults(sessionId) {
        const response = await fetch(`${this.baseURL}/api/download/${sessionId}`);
        if (!response.ok) {
            throw new APIError(`Failed to download results: ${response.statusText}`, response.status);
        }
        return response.blob();
    }

    /**
     * Upload file with progress tracking
     */
    async uploadFile(endpoint, file, onProgress = null, requestId = null) {
        const formData = new FormData();
        formData.append('file', file);

        const finalRequestId = requestId || `UPLOAD_${endpoint}`;
        
        try {
            this.setLoading(finalRequestId, true);

            const xhr = new XMLHttpRequest();
            
            return new Promise((resolve, reject) => {
                xhr.upload.addEventListener('progress', (event) => {
                    if (event.lengthComputable && onProgress) {
                        const progress = (event.loaded / event.total) * 100;
                        onProgress(progress);
                    }
                });

                xhr.addEventListener('load', () => {
                    if (xhr.status >= 200 && xhr.status < 300) {
                        try {
                            const data = JSON.parse(xhr.responseText);
                            this.dispatchSuccessEvent(finalRequestId, data);
                            resolve(data);
                        } catch (e) {
                            resolve(xhr.responseText);
                        }
                    } else {
                        const error = new APIError(
                            `HTTP ${xhr.status}: ${xhr.statusText}`,
                            xhr.status,
                            endpoint
                        );
                        this.dispatchErrorEvent(finalRequestId, error);
                        reject(error);
                    }
                });

                xhr.addEventListener('error', () => {
                    const error = new APIError('Network error', 0, endpoint);
                    this.dispatchErrorEvent(finalRequestId, error);
                    reject(error);
                });

                xhr.open('POST', `${this.baseURL}${endpoint}`);
                xhr.send(formData);
            });

        } finally {
            this.setLoading(finalRequestId, false);
        }
    }
}

/**
 * Custom API Error class
 */
class APIError extends Error {
    constructor(message, status = 0, endpoint = '') {
        super(message);
        this.name = 'APIError';
        this.status = status;
        this.endpoint = endpoint;
        this.timestamp = new Date().toISOString();
    }

    toString() {
        return `APIError: ${this.message} (${this.status}) at ${this.endpoint}`;
    }
}

// Global API client instance
const apiClient = new APIClient();

// Global loading indicator management
document.addEventListener('api:loading', (event) => {
    const { requestId, isLoading } = event.detail;
    
    // Update global loading state in UI
    const loadingIndicators = document.querySelectorAll(`[data-loading="${requestId}"]`);
    loadingIndicators.forEach(indicator => {
        if (isLoading) {
            indicator.classList.add('loading');
        } else {
            indicator.classList.remove('loading');
        }
    });
});

// Global error handling
document.addEventListener('api:error', (event) => {
    const { requestId, error } = event.detail;
    
    console.error(`API Error for ${requestId}:`, error);
    
    // Show user-friendly error message
    if (error instanceof APIError) {
        showAlert(`API Error: ${error.message}`, 'danger');
    } else {
        showAlert('An unexpected error occurred. Please try again.', 'danger');
    }
});

/**
 * Utility function to show Bootstrap alerts
 */
function showAlert(message, type = 'info', duration = 5000) {
    const alertContainer = document.getElementById('alert-container');
    if (!alertContainer) return;

    const alertId = `alert-${Date.now()}`;
    const alertHTML = `
        <div class="alert alert-${type} alert-dismissible fade show" role="alert" id="${alertId}">
            <i class="bi bi-${getAlertIcon(type)} me-2"></i>
            ${message}
            <button type="button" class="btn-close" data-bs-dismiss="alert"></button>
        </div>
    `;

    alertContainer.insertAdjacentHTML('beforeend', alertHTML);

    // Auto-dismiss after duration
    if (duration > 0) {
        setTimeout(() => {
            const alertElement = document.getElementById(alertId);
            if (alertElement) {
                const alert = new bootstrap.Alert(alertElement);
                alert.close();
            }
        }, duration);
    }
}

/**
 * Get appropriate icon for alert type
 */
function getAlertIcon(type) {
    const icons = {
        success: 'check-circle-fill',
        danger: 'exclamation-triangle-fill',
        warning: 'exclamation-triangle-fill',
        info: 'info-circle-fill'
    };
    return icons[type] || 'info-circle-fill';
}

// Global API client instance
window.apiClient = new APIClient();
