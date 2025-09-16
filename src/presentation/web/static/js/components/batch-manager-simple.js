/**
 * Batch Manager Component for CSPBench Web Interface - Simple Test Version
 */

class BatchManagerTabbed {
    constructor() {
        this.containerId = 'batch-manager-container';
        console.log('BatchManagerTabbed: Constructor called');
    }

    /**
     * Initialize the batch manager component
     */
    async initialize() {
        console.log('BatchManagerTabbed: Initialize called');
        
        try {
            console.log('BatchManagerTabbed: Looking for container...');
            const container = document.getElementById(this.containerId);
            
            if (!container) {
                console.error('BatchManagerTabbed: Container not found!');
                return;
            }
            
            console.log('BatchManagerTabbed: Container found, updating content...');
            
            container.innerHTML = `
                <div class="alert alert-success">
                    <h5><i class="bi bi-check-circle me-2"></i>Batch Manager Loaded Successfully!</h5>
                    <p>This is a test version to verify the component is working.</p>
                    <button class="btn btn-primary" onclick="alert('Test button clicked!')">
                        Test Button
                    </button>
                </div>
                
                <div class="card mt-3">
                    <div class="card-header">
                        <h6><i class="bi bi-files me-2"></i>Batch Files (Test Data)</h6>
                    </div>
                    <div class="card-body">
                        <div class="list-group">
                            <div class="list-group-item">
                                <div class="d-flex justify-content-between align-items-center">
                                    <div>
                                        <h6 class="mb-1">example_batch.yaml</h6>
                                        <p class="mb-1 small text-muted">Example batch configuration file</p>
                                        <small>Size: 1.2 KB | Created: ${new Date().toLocaleDateString()}</small>
                                    </div>
                                    <div class="btn-group">
                                        <button class="btn btn-outline-primary btn-sm" onclick="alert('Edit clicked')">
                                            <i class="bi bi-pencil"></i>
                                        </button>
                                        <button class="btn btn-outline-info btn-sm" onclick="alert('Download clicked')">
                                            <i class="bi bi-download"></i>
                                        </button>
                                        <button class="btn btn-outline-danger btn-sm" onclick="alert('Delete clicked')">
                                            <i class="bi bi-trash"></i>
                                        </button>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
            `;
            
            console.log('BatchManagerTabbed: Content updated successfully!');
            
        } catch (error) {
            console.error('BatchManagerTabbed: Error during initialization:', error);
            
            const container = document.getElementById(this.containerId);
            if (container) {
                container.innerHTML = `
                    <div class="alert alert-danger">
                        <h5><i class="bi bi-exclamation-triangle me-2"></i>Initialization Error</h5>
                        <p>Error: ${error.message}</p>
                    </div>
                `;
            }
        }
    }
}

// Export for global use
window.BatchManagerTabbed = BatchManagerTabbed;

console.log('BatchManagerTabbed: Script loaded, class available as window.BatchManagerTabbed');
