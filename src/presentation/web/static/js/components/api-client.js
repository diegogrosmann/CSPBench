/**
 * API Client Component
 * 
 * Centralized API communication for all components
 */

class APIClient {
    constructor(baseURL = '') {
        this.baseURL = baseURL;
    }

    async request(endpoint, options = {}) {
        const url = `${this.baseURL}${endpoint}`;
        const config = {
            headers: {
                'Content-Type': 'application/json',
                ...options.headers
            },
            ...options
        };

        try {
            const response = await fetch(url, config);
            
            if (!response.ok) {
                throw new Error(`HTTP error! status: ${response.status}`);
            }
            
            return await response.json();
        } catch (error) {
            console.error(`API request failed: ${endpoint}`, error);
            throw error;
        }
    }

    // Algorithm endpoints
    async getAlgorithms() {
        return this.request('/api/algorithms');
    }

    async executeAlgorithm(request) {
        return this.request('/api/execute', {
            method: 'POST',
            body: JSON.stringify(request)
        });
    }

    // Dataset endpoints
    async getDatasets() {
        return this.request('/api/datasets');
    }

    async generateDataset(config) {
        return this.request('/api/generate-dataset', {
            method: 'POST',
            body: JSON.stringify(config)
        });
    }

    // Dataset Generator endpoints
    async generateSyntheticDataset(params) {
        return this.request('/api/dataset/generate/synthetic', {
            method: 'POST',
            body: JSON.stringify(params)
        });
    }

    async downloadNCBIDataset(params) {
        return this.request('/api/dataset/generate/ncbi', {
            method: 'POST',
            body: JSON.stringify(params)
        });
    }

    async getGenerationStatus(sessionId) {
        return this.request(`/api/generation/status/${sessionId}`);
    }

    async getSavedDatasets() {
        return this.request('/api/datasets/saved');
    }

    async getDatasetDownloadUrl(datasetId) {
        // For saved datasets, return the download URL
        return `${this.baseURL}/api/dataset/download/${datasetId}`;
    }

    async deleteDataset(datasetId) {
        return this.request(`/api/datasets/${datasetId}`, {
            method: 'DELETE'
        });
    }

    async cancelGeneration(sessionId) {
        return this.request(`/api/generation/cancel/${sessionId}`, {
            method: 'POST'
        });
    }

    // Session endpoints
    async getSessions() {
        return this.request('/api/sessions');
    }

    async getSession(sessionId) {
        return this.request(`/api/sessions/${sessionId}`);
    }

    // Download endpoint
    getDownloadURL(sessionId) {
        return `${this.baseURL}/api/download/${sessionId}`;
    }
}

// Global API client instance
window.apiClient = new APIClient();
