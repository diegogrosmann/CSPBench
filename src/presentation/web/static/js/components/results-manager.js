/**
 * Results Manager Component
 * 
 * Handles the results management interface including:
 * - Loading and displaying results
 * - Filtering and searching
 * - Batch operations
 * - Downloads
 */

class ResultsManager {
    constructor() {
        this.currentPage = 1;
        this.pageSize = 20;
        this.selectedSessions = new Set();
        this.results = [];
        this.filteredResults = [];
        this.filters = {
            search: '',
            dateFrom: '',
            dateTo: ''
        };
        
        this.apiClient = window.apiClient || new APIClient();
        this.bindEvents();
    }

    initialize() {
        this.loadResults();
    }

    bindEvents() {
        // Filter events
        document.getElementById('search-query')?.addEventListener('input', 
            this.debounce((e) => this.updateFilter('search', e.target.value), 300));
        document.getElementById('filter-date-from')?.addEventListener('change', 
            (e) => this.updateFilter('dateFrom', e.target.value));
        document.getElementById('filter-date-to')?.addEventListener('change', 
            (e) => this.updateFilter('dateTo', e.target.value));
        document.getElementById('apply-filters')?.addEventListener('click', 
            () => this.applyFilters());

        // Batch action events
        document.getElementById('select-all')?.addEventListener('click', 
            () => this.selectAll());
        document.getElementById('clear-selection')?.addEventListener('click', 
            () => this.clearSelection());
        document.getElementById('download-selected')?.addEventListener('click', 
            () => this.downloadSelected());
        document.getElementById('delete-selected')?.addEventListener('click', 
            () => this.deleteSelected());

        // Select all checkbox
        document.getElementById('select-all-checkbox')?.addEventListener('change', 
            (e) => e.target.checked ? this.selectAll() : this.clearSelection());

        // Refresh
        document.getElementById('refresh-results')?.addEventListener('click', 
            () => this.loadResults());

        // Toggle filters button
        document.getElementById('toggle-filters')?.addEventListener('click', 
            () => this.toggleFiltersIcon());

        // Filter section collapse events
        const filtersSection = document.getElementById('filters-section');
        if (filtersSection) {
            filtersSection.addEventListener('shown.bs.collapse', () => this.updateFiltersIcon(true));
            filtersSection.addEventListener('hidden.bs.collapse', () => this.updateFiltersIcon(false));
        }

        // Retry
        document.getElementById('retry-load')?.addEventListener('click', 
            () => this.loadResults());

        // Delete confirmation
        document.getElementById('confirm-delete')?.addEventListener('click', 
            () => this.confirmDelete());
    }

    async loadResults() {
        this.showLoading();
        
        try {
            const response = await this.apiClient.get('/api/results/', {
                limit: 1000, // Load all for client-side filtering
                offset: 0
            });

            this.results = response.sessions || [];
            this.updateSummary(response);
            this.applyFilters();
            this.hideLoading();

        } catch (error) {
            console.error('Failed to load results:', error);
            this.showError('Failed to load results: ' + error.message);
        }
    }

    async loadAlgorithms() {
        try {
            const algorithms = await this.apiClient.get('/api/algorithms');
            const select = document.getElementById('filter-algorithm');
            
            if (select) {
                // Clear existing options except first
                select.innerHTML = '<option value="">All algorithms</option>';
                
                algorithms.forEach(algorithm => {
                    const option = document.createElement('option');
                    option.value = algorithm.name;
                    option.textContent = algorithm.name;
                    select.appendChild(option);
                });
            }
        } catch (error) {
            console.warn('Failed to load algorithms for filter:', error);
        }
    }

    updateFilter(key, value) {
        this.filters[key] = value;
    }

    applyFilters() {
        this.filteredResults = this.results.filter(session => {
            // Search filter
            if (this.filters.search) {
                const searchLower = this.filters.search.toLowerCase();
                const searchableText = [
                    session.session_id,
                    session.algorithm,
                    session.dataset_name
                ].filter(Boolean).join(' ').toLowerCase();
                
                if (!searchableText.includes(searchLower)) {
                    return false;
                }
            }

            // Date filters
            if (this.filters.dateFrom) {
                const fromDate = new Date(this.filters.dateFrom);
                const sessionDate = new Date(session.timestamp);
                if (sessionDate < fromDate) {
                    return false;
                }
            }

            if (this.filters.dateTo) {
                const toDate = new Date(this.filters.dateTo);
                toDate.setHours(23, 59, 59, 999); // End of day
                const sessionDate = new Date(session.timestamp);
                if (sessionDate > toDate) {
                    return false;
                }
            }

            return true;
        });

        this.currentPage = 1;
        this.renderResults();
        this.renderPagination();
    }

    renderResults() {
        const tbody = document.getElementById('results-table-body');
        const tableContainer = document.getElementById('results-table-container');
        const emptyState = document.getElementById('empty-state');

        if (!tbody) return;

        if (this.filteredResults.length === 0) {
            tableContainer.style.display = 'none';
            emptyState.style.display = 'block';
            this.hidePagination();
            return;
        }

        tableContainer.style.display = 'block';
        emptyState.style.display = 'none';

        // Pagination
        const startIndex = (this.currentPage - 1) * this.pageSize;
        const endIndex = startIndex + this.pageSize;
        const pageResults = this.filteredResults.slice(startIndex, endIndex);

        tbody.innerHTML = pageResults.map(session => this.renderResultRow(session)).join('');

        // Bind row events
        this.bindRowEvents();
        this.showPagination();
    }

    renderResultRow(session) {
        const timestamp = new Date(session.timestamp).toLocaleString();
        const size = session.size_mb ? `${session.size_mb} MB` : '-';
        const isSelected = this.selectedSessions.has(session.session_id);

        return `
            <tr class="session-row" data-session-id="${session.session_id}">
                <td>
                    <input type="checkbox" class="form-check-input session-checkbox" 
                           ${isSelected ? 'checked' : ''} 
                           data-session-id="${session.session_id}">
                </td>
                <td>
                    <code class="text-primary">
                        <a href="/results/${session.session_id}" class="text-decoration-none">${session.session_id}</a>
                    </code>
                </td>
                <td>
                    <small>${timestamp}</small>
                </td>
                <td>
                    <span class="badge status-badge status-${session.status}">
                        ${this.getStatusIcon(session.status)} ${session.status}
                    </span>
                </td>
                <td>${size}</td>
                <td>
                    <div class="btn-group btn-group-sm" role="group">
                        <button class="btn btn-outline-primary action-btn view-details" 
                                data-session-id="${session.session_id}" 
                                title="View Details">
                            <i class="bi bi-eye"></i>
                        </button>
                        <button class="btn btn-outline-success action-btn download-session" 
                                data-session-id="${session.session_id}" 
                                title="Download ZIP">
                            <i class="bi bi-download"></i>
                        </button>
                        <div class="btn-group btn-group-sm" role="group">
                            <button class="btn btn-outline-secondary action-btn dropdown-toggle" 
                                    data-bs-toggle="dropdown" 
                                    title="Download specific files">
                                <i class="bi bi-file-earmark"></i>
                            </button>
                            <ul class="dropdown-menu">
                                <li><a class="dropdown-item download-file" 
                                       data-session-id="${session.session_id}" 
                                       data-file-type="json">JSON Results</a></li>
                                <li><a class="dropdown-item download-file" 
                                       data-session-id="${session.session_id}" 
                                       data-file-type="csv">CSV Results</a></li>
                                <li><a class="dropdown-item download-file" 
                                       data-session-id="${session.session_id}" 
                                       data-file-type="log">Log File</a></li>
                                <li><a class="dropdown-item download-file" 
                                       data-session-id="${session.session_id}" 
                                       data-file-type="plots">Plots</a></li>
                            </ul>
                        </div>
                        <button class="btn btn-outline-danger action-btn delete-session" 
                                data-session-id="${session.session_id}" 
                                title="Delete">
                            <i class="bi bi-trash"></i>
                        </button>
                    </div>
                </td>
            </tr>
        `;
    }

    bindRowEvents() {
        // Checkbox events
        document.querySelectorAll('.session-checkbox').forEach(checkbox => {
            checkbox.addEventListener('change', (e) => {
                const sessionId = e.target.dataset.sessionId;
                if (e.target.checked) {
                    this.selectedSessions.add(sessionId);
                } else {
                    this.selectedSessions.delete(sessionId);
                }
                this.updateSelectionUI();
            });
        });

        // Action button events
        document.querySelectorAll('.view-details').forEach(btn => {
            btn.addEventListener('click', (e) => {
                const sessionId = e.target.closest('[data-session-id]').dataset.sessionId;
                this.viewDetails(sessionId);
            });
        });

        document.querySelectorAll('.download-session').forEach(btn => {
            btn.addEventListener('click', (e) => {
                const sessionId = e.target.closest('[data-session-id]').dataset.sessionId;
                this.downloadSession(sessionId);
            });
        });

        document.querySelectorAll('.download-file').forEach(btn => {
            btn.addEventListener('click', (e) => {
                e.preventDefault();
                const sessionId = e.target.dataset.sessionId;
                const fileType = e.target.dataset.fileType;
                this.downloadFile(sessionId, fileType);
            });
        });

        document.querySelectorAll('.delete-session').forEach(btn => {
            btn.addEventListener('click', (e) => {
                const sessionId = e.target.closest('[data-session-id]').dataset.sessionId;
                this.deleteSingleSession(sessionId);
            });
        });
    }

    updateSelectionUI() {
        const count = this.selectedSessions.size;
        document.getElementById('selected-count').textContent = `${count} selected`;
        
        const downloadBtn = document.getElementById('download-selected');
        const deleteBtn = document.getElementById('delete-selected');
        const batchActions = document.getElementById('batch-actions');
        
        if (downloadBtn) downloadBtn.disabled = count === 0;
        if (deleteBtn) deleteBtn.disabled = count === 0;
        
        // Show/hide batch actions based on selection
        if (batchActions) {
            batchActions.style.display = count > 0 ? 'block' : 'none';
        }

        // Update select all checkbox
        const selectAllCheckbox = document.getElementById('select-all-checkbox');
        if (selectAllCheckbox) {
            const visibleSessions = this.getVisibleSessionIds();
            const allSelected = visibleSessions.length > 0 && 
                                visibleSessions.every(id => this.selectedSessions.has(id));
            selectAllCheckbox.checked = allSelected;
            selectAllCheckbox.indeterminate = !allSelected && 
                                            visibleSessions.some(id => this.selectedSessions.has(id));
        }
    }

    getVisibleSessionIds() {
        const startIndex = (this.currentPage - 1) * this.pageSize;
        const endIndex = startIndex + this.pageSize;
        return this.filteredResults.slice(startIndex, endIndex).map(s => s.session_id);
    }

    selectAll() {
        this.getVisibleSessionIds().forEach(id => this.selectedSessions.add(id));
        this.updateSelectionUI();
        this.renderResults(); // Re-render to update checkboxes
    }

    clearSelection() {
        this.selectedSessions.clear();
        this.updateSelectionUI();
        this.renderResults(); // Re-render to update checkboxes
    }

    async downloadSelected() {
        if (this.selectedSessions.size === 0) return;

        try {
            const sessionIds = Array.from(this.selectedSessions);
            
            // Show loading indication
            const btn = document.getElementById('download-selected');
            const originalText = btn.innerHTML;
            btn.innerHTML = '<i class="bi bi-hourglass-split me-1"></i>Preparing...';
            btn.disabled = true;

            const response = await fetch('/api/results/download/batch', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify(sessionIds)
            });

            if (!response.ok) {
                throw new Error(`HTTP ${response.status}: ${response.statusText}`);
            }

            // Download the file
            const blob = await response.blob();
            const url = window.URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = `cspbench_batch_results_${new Date().toISOString().slice(0, 19).replace(/:/g, '-')}.zip`;
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
            window.URL.revokeObjectURL(url);

            // Reset button
            btn.innerHTML = originalText;
            btn.disabled = false;

        } catch (error) {
            console.error('Failed to download selected results:', error);
            this.showToast('Failed to download selected results: ' + error.message, 'error');
            
            // Reset button
            const btn = document.getElementById('download-selected');
            btn.innerHTML = '<i class="bi bi-download me-1"></i>Download Selected';
            btn.disabled = false;
        }
    }

    deleteSelected() {
        if (this.selectedSessions.size === 0) return;
        
        const modal = new bootstrap.Modal(document.getElementById('deleteConfirmModal'));
        modal.show();
    }

    async confirmDelete() {
        const sessionIds = Array.from(this.selectedSessions);
        
        try {
            const deletePromises = sessionIds.map(sessionId => 
                this.apiClient.delete(`/api/results/${sessionId}`)
            );
            
            await Promise.all(deletePromises);
            
            // Remove from local data
            this.results = this.results.filter(s => !sessionIds.includes(s.session_id));
            this.clearSelection();
            this.applyFilters();
            
            this.showToast(`Successfully deleted ${sessionIds.length} session(s)`, 'success');
            
            // Hide modal
            const modal = bootstrap.Modal.getInstance(document.getElementById('deleteConfirmModal'));
            modal.hide();

        } catch (error) {
            console.error('Failed to delete sessions:', error);
            this.showToast('Failed to delete some sessions: ' + error.message, 'error');
        }
    }

    async deleteSingleSession(sessionId) {
        if (!confirm('Are you sure you want to delete this session? This action cannot be undone.')) {
            return;
        }

        try {
            await this.apiClient.delete(`/api/results/${sessionId}`);
            
            // Remove from local data
            this.results = this.results.filter(s => s.session_id !== sessionId);
            this.selectedSessions.delete(sessionId);
            this.applyFilters();
            
            this.showToast('Session deleted successfully', 'success');

        } catch (error) {
            console.error('Failed to delete session:', error);
            this.showToast('Failed to delete session: ' + error.message, 'error');
        }
    }

    async downloadSession(sessionId) {
        try {
            window.open(`/api/results/${sessionId}/download`, '_blank');
        } catch (error) {
            console.error('Failed to download session:', error);
            this.showToast('Failed to download session: ' + error.message, 'error');
        }
    }

    async downloadFile(sessionId, fileType) {
        try {
            window.open(`/api/results/${sessionId}/files/${fileType}`, '_blank');
        } catch (error) {
            console.error('Failed to download file:', error);
            this.showToast('Failed to download file: ' + error.message, 'error');
        }
    }

    async viewDetails(sessionId) {
        try {
            const details = await this.apiClient.get(`/api/results/${sessionId}`);
            this.showDetailsModal(details);
        } catch (error) {
            console.error('Failed to load details:', error);
            this.showToast('Failed to load details: ' + error.message, 'error');
        }
    }

    showDetailsModal(details) {
        const modalBody = document.getElementById('result-details-content');
        if (!modalBody) return;

        modalBody.innerHTML = `
            <div class="row">
                <div class="col-md-6">
                    <h6>Session Information</h6>
                    <table class="table table-sm">
                        <tr><td><strong>Session ID:</strong></td><td><code>${details.session_id}</code></td></tr>
                        <tr><td><strong>Algorithm:</strong></td><td>${details.algorithm || '-'}</td></tr>
                        <tr><td><strong>Dataset:</strong></td><td>${details.dataset_name || '-'}</td></tr>
                        <tr><td><strong>Status:</strong></td><td><span class="badge status-${details.status}">${details.status}</span></td></tr>
                        <tr><td><strong>Timestamp:</strong></td><td>${new Date(details.timestamp).toLocaleString()}</td></tr>
                        <tr><td><strong>Execution Time:</strong></td><td>${details.execution_time ? details.execution_time.toFixed(2) + 's' : '-'}</td></tr>
                        <tr><td><strong>Size:</strong></td><td>${details.size_mb ? details.size_mb + ' MB' : '-'}</td></tr>
                    </table>
                </div>
                <div class="col-md-6">
                    <h6>Results</h6>
                    <table class="table table-sm">
                        <tr><td><strong>Best String:</strong></td><td><code>${details.best_string || '-'}</code></td></tr>
                        <tr><td><strong>Max Distance:</strong></td><td>${details.max_distance || '-'}</td></tr>
                    </table>
                    
                    <h6>Files</h6>
                    <ul class="list-group list-group-flush">
                        ${details.files.map(file => `<li class="list-group-item py-1"><small>${file}</small></li>`).join('')}
                    </ul>
                </div>
            </div>
            
            ${details.detailed_results ? `
                <div class="mt-3">
                    <h6>Detailed Results</h6>
                    <pre class="bg-light p-3 rounded"><code>${JSON.stringify(details.detailed_results, null, 2)}</code></pre>
                </div>
            ` : ''}
        `;

        const modal = new bootstrap.Modal(document.getElementById('resultDetailsModal'));
        modal.show();
    }

    renderPagination() {
        const pagination = document.getElementById('pagination');
        if (!pagination) return;

        const totalPages = Math.ceil(this.filteredResults.length / this.pageSize);
        
        if (totalPages <= 1) {
            this.hidePagination();
            return;
        }

        let paginationHTML = '';

        // Previous button
        paginationHTML += `
            <li class="page-item ${this.currentPage === 1 ? 'disabled' : ''}">
                <a class="page-link" href="#" data-page="${this.currentPage - 1}">Previous</a>
            </li>
        `;

        // Page numbers
        const startPage = Math.max(1, this.currentPage - 2);
        const endPage = Math.min(totalPages, this.currentPage + 2);

        if (startPage > 1) {
            paginationHTML += `<li class="page-item"><a class="page-link" href="#" data-page="1">1</a></li>`;
            if (startPage > 2) {
                paginationHTML += `<li class="page-item disabled"><span class="page-link">...</span></li>`;
            }
        }

        for (let i = startPage; i <= endPage; i++) {
            paginationHTML += `
                <li class="page-item ${i === this.currentPage ? 'active' : ''}">
                    <a class="page-link" href="#" data-page="${i}">${i}</a>
                </li>
            `;
        }

        if (endPage < totalPages) {
            if (endPage < totalPages - 1) {
                paginationHTML += `<li class="page-item disabled"><span class="page-link">...</span></li>`;
            }
            paginationHTML += `<li class="page-item"><a class="page-link" href="#" data-page="${totalPages}">${totalPages}</a></li>`;
        }

        // Next button
        paginationHTML += `
            <li class="page-item ${this.currentPage === totalPages ? 'disabled' : ''}">
                <a class="page-link" href="#" data-page="${this.currentPage + 1}">Next</a>
            </li>
        `;

        pagination.innerHTML = paginationHTML;

        // Bind pagination events
        pagination.querySelectorAll('a.page-link').forEach(link => {
            link.addEventListener('click', (e) => {
                e.preventDefault();
                const page = parseInt(e.target.dataset.page);
                if (page && page !== this.currentPage) {
                    this.currentPage = page;
                    this.renderResults();
                    this.renderPagination();
                }
            });
        });

        this.showPagination();
    }

    updateSummary(response) {
        document.getElementById('total-sessions').textContent = response.total_count || 0;
        document.getElementById('total-size').textContent = `${response.total_size_mb || 0} MB`;
        
        const completed = this.results.filter(s => s.status === 'completed').length;
        const failed = this.results.filter(s => s.status === 'failed').length;
        
        document.getElementById('completed-sessions').textContent = completed;
        document.getElementById('failed-sessions').textContent = failed;
    }

    getStatusIcon(status) {
        const icons = {
            completed: '<i class="bi bi-check-circle"></i>',
            failed: '<i class="bi bi-x-circle"></i>',
            timeout: '<i class="bi bi-clock"></i>',
            running: '<i class="bi bi-hourglass-split"></i>'
        };
        return icons[status] || '<i class="bi bi-question-circle"></i>';
    }

    showLoading() {
        document.getElementById('loading-indicator').style.display = 'block';
        document.getElementById('results-table-container').style.display = 'none';
        document.getElementById('empty-state').style.display = 'none';
        document.getElementById('error-state').style.display = 'none';
    }

    hideLoading() {
        document.getElementById('loading-indicator').style.display = 'none';
    }

    showError(message) {
        document.getElementById('loading-indicator').style.display = 'none';
        document.getElementById('results-table-container').style.display = 'none';
        document.getElementById('empty-state').style.display = 'none';
        document.getElementById('error-state').style.display = 'block';
        document.getElementById('error-message').textContent = message;
    }

    showPagination() {
        document.getElementById('pagination-container').style.display = 'block';
    }

    hidePagination() {
        document.getElementById('pagination-container').style.display = 'none';
    }

    showToast(message, type = 'info') {
        // Create toast element
        const toastHTML = `
            <div class="toast align-items-center text-white bg-${type === 'error' ? 'danger' : type === 'success' ? 'success' : 'primary'} border-0" role="alert">
                <div class="d-flex">
                    <div class="toast-body">${message}</div>
                    <button type="button" class="btn-close btn-close-white me-2 m-auto" data-bs-dismiss="toast"></button>
                </div>
            </div>
        `;

        // Add to container
        let container = document.getElementById('toast-container');
        if (!container) {
            container = document.createElement('div');
            container.id = 'toast-container';
            container.className = 'toast-container position-fixed top-0 end-0 p-3';
            document.body.appendChild(container);
        }

        container.insertAdjacentHTML('beforeend', toastHTML);
        const toastElement = container.lastElementChild;
        const toast = new bootstrap.Toast(toastElement);
        toast.show();

        // Auto remove after hiding
        toastElement.addEventListener('hidden.bs.toast', () => {
            toastElement.remove();
        });
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

    updateFiltersIcon(isExpanded) {
        const toggleBtn = document.getElementById('toggle-filters');
        if (toggleBtn) {
            const icon = toggleBtn.querySelector('i');
            if (icon) {
                icon.className = isExpanded ? 'bi bi-funnel-fill me-1' : 'bi bi-funnel me-1';
            }
        }
    }

    toggleFiltersIcon() {
        const filtersSection = document.getElementById('filters-section');
        const isExpanded = filtersSection?.classList.contains('show');
        this.updateFiltersIcon(!isExpanded);
    }
}

// Make available globally
window.ResultsManager = ResultsManager;
