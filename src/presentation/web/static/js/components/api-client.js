/**
 * Modern API Client Component for CSPBench Web Interface
 * Centralized HTTP client with modern fetch API, error handling, and loading states
 *
 * This file is idempotent: it won't redeclare classes or rebind listeners
 * if loaded multiple times.
 */

(function(global){
    // Define APIError only if not already defined
    if (!global.APIError) {
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
        global.APIError = APIError;
    }
    const APIError = global.APIError;

    // Define APIClient only if not already defined
    if (!global.APIClient) {
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

    // =============================================================================
    // BATCH EXECUTION API METHODS
    // =============================================================================

    /**
     * Execute a batch configuration file
     */
    async executeBatch(batchFile, monitorType = 'log') {
        return this.request('/api/batch/execute', {
            method: 'POST',
            body: JSON.stringify({
                batch_file: batchFile,
                monitor_type: monitorType
            })
        }, 'execute_batch');
    }

    /**
     * Get batch execution status
     */
    async getBatchStatus(workId) {
        return this.request(`/api/batch/${workId}/status`, {
            method: 'GET'
        }, `status_${workId}`);
    }

    /**
     * Get batch execution results
     */
    async getBatchResults(workId) {
        return this.request(`/api/batch/${workId}/results`, {
            method: 'GET'
        }, `results_${workId}`);
    }

    /**
     * List all batch executions with optional status filter
     */
    async listBatchExecutions(status = null) {
        const params = status ? `?status=${status}` : '';
        return this.request(`/api/batch/list${params}`, {
            method: 'GET'
        }, 'list_batches');
    }

    /**
     * Control batch execution (cancel, pause, resume)
     */
    async controlBatchExecution(workId, action) {
        return this.request(`/api/batch/${workId}/control`, {
            method: 'POST',
            body: JSON.stringify({ action })
        }, `control_${workId}`);
    }

    /**
     * Delete batch execution
     */
    async deleteBatchExecution(workId) {
        return this.request(`/api/batch/${workId}`, {
            method: 'DELETE'
        }, `delete_${workId}`);
    }

    // =============================================================================
    // BATCH FILE MANAGEMENT API METHODS  
    // =============================================================================

    /**
     * Get list of batch configuration files
     */
    async getBatchFiles() {
        return this.request('/api/batches', {
            method: 'GET'
        }, 'get_batch_files');
    }

    /**
     * Upload a new batch configuration file
     */
    async uploadBatchFile(formData) {
        // Note: FormData automatically sets Content-Type with boundary
        return this.request('/api/batches/upload', {
            method: 'POST',
            body: formData,
            headers: {} // Let browser set Content-Type for FormData
        }, 'upload_batch');
    }

    /**
     * Create a new batch configuration file
     */
    async createBatchFile(name, content, description = '') {
        return this.request('/api/batches', {
            method: 'POST',
            body: JSON.stringify({
                name,
                content,
                description
            })
        }, 'create_batch');
    }

    /**
     * Get batch file content
     */
    async getBatchFile(fileName) {
        return this.request(`/api/batches/${fileName}`, {
            method: 'GET'
        }, `get_batch_${fileName}`);
    }

    /**
     * Update batch file content
     */
    async updateBatchFile(fileName, content, description = '') {
        return this.request(`/api/batches/${fileName}`, {
            method: 'PUT',
            body: JSON.stringify({
                content,
                description
            })
        }, `update_batch_${fileName}`);
    }

    /**
     * Delete batch file
     */
    async deleteBatchFile(fileName) {
        return this.request(`/api/batches/${fileName}`, {
            method: 'DELETE'
        }, `delete_batch_${fileName}`);
    }   
        }
        global.APIClient = APIClient;
    }
    const APIClient = global.APIClient;

    // Create a single global apiClient instance if not existing
    if (!global.apiClient) {
        global.apiClient = new APIClient();
    }

    // Register global listeners once
    if (!global.__apiClientListenersRegistered) {
        // Global loading indicator management
        document.addEventListener('api:loading', (event) => {
            const { requestId, isLoading } = event.detail;
            const loadingIndicators = document.querySelectorAll(`[data-loading="${requestId}"]`);
            loadingIndicators.forEach(indicator => {
                indicator.classList.toggle('loading', !!isLoading);
            });
        });

        // Global error handling
        document.addEventListener('api:error', (event) => {
            const { requestId, error } = event.detail;
            console.error(`API Error for ${requestId}:`, error);
            if (error instanceof APIError) {
                if (typeof global.showAlert === 'function') {
                    global.showAlert(`API Error: ${error.message}`, 'danger');
                } else if (global.notifications && typeof global.notifications.error === 'function') {
                    global.notifications.error(`API Error: ${error.message}`);
                }
            } else {
                if (typeof global.showAlert === 'function') {
                    global.showAlert('An unexpected error occurred. Please try again.', 'danger');
                } else if (global.notifications && typeof global.notifications.error === 'function') {
                    global.notifications.error('An unexpected error occurred. Please try again.');
                }
            }
        });

        global.__apiClientListenersRegistered = true;
    }
})(window);
