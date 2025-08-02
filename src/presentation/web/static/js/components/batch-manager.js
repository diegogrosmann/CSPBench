/**
 * Batch Manager Component
 * 
 * Handles the batch execution interface including:
 * - Loading and displaying batch files
 * - Filtering and searching
 * - Batch operations
 * - File uploads
 */

class BatchManager {
    constructor() {
        this.currentPage = 1;
        this.pageSize = 20;
        this.selectedBatches = new Set();
        this.batches = [];
        this.filteredBatches = [];
        this.filters = {
            search: '',
            status: '',
            dateFrom: ''
        };
        
        this.apiClient = window.apiClient || new APIClient();
        this.bindEvents();
    }

    initialize() {
        this.loadBatches();
    }

    bindEvents() {
        // Filter events
        document.getElementById('search-query')?.addEventListener('input', 
            this.debounce((e) => this.updateFilter('search', e.target.value), 300));
        document.getElementById('filter-status')?.addEventListener('change', 
            (e) => this.updateFilter('status', e.target.value));
        document.getElementById('filter-date-from')?.addEventListener('change', 
            (e) => this.updateFilter('dateFrom', e.target.value));
        document.getElementById('apply-filters')?.addEventListener('click', 
            () => this.applyFilters());

        // Batch action events
        document.getElementById('run-selected')?.addEventListener('click', 
            () => this.runSelected());
        document.getElementById('delete-selected')?.addEventListener('click', 
            () => this.deleteSelected());

        // Select all checkbox
        document.getElementById('select-all-checkbox')?.addEventListener('change', 
            (e) => e.target.checked ? this.selectAll() : this.clearSelection());

        // Refresh
        document.getElementById('refresh-batches')?.addEventListener('click', 
            () => this.loadBatches());

        // Toggle filters button
        document.getElementById('toggle-filters')?.addEventListener('click', 
            () => this.toggleFiltersIcon());

        // Filter section collapse events
        const filtersSection = document.getElementById('filters-section');
        if (filtersSection) {
            filtersSection.addEventListener('shown.bs.collapse', () => this.updateFiltersIcon(true));
            filtersSection.addEventListener('hidden.bs.collapse', () => this.updateFiltersIcon(false));
        }

        // Upload events
        document.getElementById('upload-submit')?.addEventListener('click', 
            () => this.uploadBatch());

        // Modal events
        document.getElementById('confirm-run')?.addEventListener('click', 
            () => this.confirmRun());
        document.getElementById('confirm-delete')?.addEventListener('click', 
            () => this.confirmDelete());

        // Retry
        document.getElementById('retry-load')?.addEventListener('click', 
            () => this.loadBatches());
    }

    async loadBatches() {
        this.showLoading();
        
        try {
            const response = await this.apiClient.get('/execution/api/batches', 'load_batches');
            this.batches = response.batch_files || [];
            this.updateStatistics();
            this.applyFilters();
            this.hideLoading();
        } catch (error) {
            console.error('Error loading batches:', error);
            this.showError('Failed to load batch files: ' + (error.message || 'Unknown error'));
        }
    }

    updateStatistics() {
        const total = this.batches.length;
        const running = this.batches.filter(b => b.status === 'running').length;
        const completed = this.batches.filter(b => b.status === 'completed').length;
        const failed = this.batches.filter(b => b.status === 'failed').length;

        document.getElementById('total-batches').textContent = total;
        document.getElementById('running-batches').textContent = running;
        document.getElementById('completed-batches').textContent = completed;
        document.getElementById('failed-batches').textContent = failed;
    }

    updateFilter(filterName, value) {
        this.filters[filterName] = value;
        this.applyFilters();
    }

    applyFilters() {
        this.filteredBatches = this.batches.filter(batch => {
            // Search filter
            if (this.filters.search) {
                const searchLower = this.filters.search.toLowerCase();
                const searchableText = [
                    batch.name,
                    batch.description,
                    batch.filename
                ].filter(Boolean).join(' ').toLowerCase();
                
                if (!searchableText.includes(searchLower)) {
                    return false;
                }
            }

            // Status filter
            if (this.filters.status && batch.status !== this.filters.status) {
                return false;
            }

            // Date filters
            if (this.filters.dateFrom) {
                const fromDate = new Date(this.filters.dateFrom);
                const batchDate = new Date(batch.modified * 1000);
                if (batchDate < fromDate) {
                    return false;
                }
            }

            return true;
        });

        this.currentPage = 1;
        this.renderBatches();
        this.renderPagination();
    }

    renderBatches() {
        const tbody = document.getElementById('batch-table-body');
        const tableContainer = document.getElementById('batch-table-container');
        const emptyState = document.getElementById('empty-state');

        if (!tbody) return;

        if (this.filteredBatches.length === 0) {
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
        const pageBatches = this.filteredBatches.slice(startIndex, endIndex);

        tbody.innerHTML = pageBatches.map(batch => this.renderBatchRow(batch)).join('');

        // Bind row events
        this.bindRowEvents();
        this.showPagination();
    }

    renderBatchRow(batch) {
        const created = new Date(batch.modified * 1000).toLocaleString();
        const size = this.formatFileSize(batch.size);
        const isSelected = this.selectedBatches.has(batch.filename);
        const status = batch.status || 'pending';
        const progress = batch.progress || 0;

        return `
            <tr class="batch-row" data-batch-id="${batch.filename}">
                <td>
                    <input type="checkbox" class="form-check-input batch-checkbox" 
                           ${isSelected ? 'checked' : ''} 
                           data-batch-id="${batch.filename}">
                </td>
                <td>
                    <div class="d-flex flex-column">
                        <code class="text-primary mb-1">${batch.name}</code>
                        <small class="text-muted">${batch.description || 'No description'}</small>
                    </div>
                </td>
                <td>
                    <small>${created}</small>
                </td>
                <td>
                    <span class="badge status-badge status-${status}">
                        ${this.getStatusIcon(status)} ${status}
                    </span>
                </td>
                <td>
                    <div class="progress progress-batch">
                        <div class="progress-bar ${this.getProgressBarClass(status)}" 
                             role="progressbar" 
                             style="width: ${progress}%" 
                             aria-valuenow="${progress}" 
                             aria-valuemin="0" 
                             aria-valuemax="100">
                            ${progress > 0 ? progress + '%' : ''}
                        </div>
                    </div>
                </td>
                <td>
                    <div class="btn-group btn-group-sm" role="group">
                        <button class="btn btn-outline-primary action-btn view-details" 
                                data-batch-id="${batch.filename}" 
                                title="View Details">
                            <i class="bi bi-eye"></i>
                        </button>
                        <button class="btn btn-outline-success action-btn run-batch" 
                                data-batch-id="${batch.filename}" 
                                title="Run Batch"
                                ${status === 'running' ? 'disabled' : ''}>
                            <i class="bi bi-play"></i>
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

    bindRowEvents() {
        // Checkbox events
        document.querySelectorAll('.batch-checkbox').forEach(checkbox => {
            checkbox.addEventListener('change', (e) => {
                const batchId = e.target.dataset.batchId;
                if (e.target.checked) {
                    this.selectedBatches.add(batchId);
                } else {
                    this.selectedBatches.delete(batchId);
                }
                this.updateSelectionUI();
            });
        });

        // Action button events
        document.querySelectorAll('.view-details').forEach(btn => {
            btn.addEventListener('click', (e) => {
                const batchId = e.target.closest('button').dataset.batchId;
                this.viewBatchDetails(batchId);
            });
        });

        document.querySelectorAll('.run-batch').forEach(btn => {
            btn.addEventListener('click', (e) => {
                const batchId = e.target.closest('button').dataset.batchId;
                this.runSingleBatch(batchId);
            });
        });

        document.querySelectorAll('.download-batch').forEach(btn => {
            btn.addEventListener('click', (e) => {
                const batchId = e.target.closest('button').dataset.batchId;
                this.downloadBatch(batchId);
            });
        });

        document.querySelectorAll('.delete-batch').forEach(btn => {
            btn.addEventListener('click', (e) => {
                const batchId = e.target.closest('button').dataset.batchId;
                this.deleteSingleBatch(batchId);
            });
        });
    }

    updateSelectionUI() {
        const count = this.selectedBatches.size;
        document.getElementById('selected-count').textContent = `${count} selected`;
        
        const runBtn = document.getElementById('run-selected');
        const deleteBtn = document.getElementById('delete-selected');
        const batchActions = document.getElementById('batch-actions');
        
        if (runBtn) runBtn.disabled = count === 0;
        if (deleteBtn) deleteBtn.disabled = count === 0;
        
        // Show/hide batch actions based on selection
        if (batchActions) {
            batchActions.style.display = count > 0 ? 'block' : 'none';
        }

        // Update select all checkbox
        const selectAllCheckbox = document.getElementById('select-all-checkbox');
        if (selectAllCheckbox) {
            const visibleBatches = this.filteredBatches.length;
            if (visibleBatches === 0) {
                selectAllCheckbox.checked = false;
                selectAllCheckbox.indeterminate = false;
            } else if (count === visibleBatches) {
                selectAllCheckbox.checked = true;
                selectAllCheckbox.indeterminate = false;
            } else if (count > 0) {
                selectAllCheckbox.checked = false;
                selectAllCheckbox.indeterminate = true;
            } else {
                selectAllCheckbox.checked = false;
                selectAllCheckbox.indeterminate = false;
            }
        }
    }

    selectAll() {
        this.filteredBatches.forEach(batch => {
            this.selectedBatches.add(batch.filename);
        });
        this.updateSelectionUI();
        this.renderBatches();
    }

    clearSelection() {
        this.selectedBatches.clear();
        this.updateSelectionUI();
        this.renderBatches();
    }

    // Status and progress helpers
    getStatusIcon(status) {
        const icons = {
            pending: '<i class="bi bi-clock"></i>',
            running: '<i class="bi bi-play-circle"></i>',
            completed: '<i class="bi bi-check-circle"></i>',
            failed: '<i class="bi bi-x-circle"></i>'
        };
        return icons[status] || '<i class="bi bi-question-circle"></i>';
    }

    getProgressBarClass(status) {
        const classes = {
            pending: 'bg-secondary',
            running: 'bg-primary',
            completed: 'bg-success',
            failed: 'bg-danger'
        };
        return classes[status] || 'bg-secondary';
    }

    formatFileSize(bytes) {
        if (bytes === 0) return '0 B';
        const k = 1024;
        const sizes = ['B', 'KB', 'MB', 'GB'];
        const i = Math.floor(Math.log(bytes) / Math.log(k));
        return parseFloat((bytes / Math.pow(k, i)).toFixed(2)) + ' ' + sizes[i];
    }

    // UI State management
    showLoading() {
        document.getElementById('loading-indicator').style.display = 'block';
        document.getElementById('batch-table-container').style.display = 'none';
        document.getElementById('empty-state').style.display = 'none';
        document.getElementById('error-state').style.display = 'none';
    }

    hideLoading() {
        document.getElementById('loading-indicator').style.display = 'none';
    }

    showError(message) {
        document.getElementById('loading-indicator').style.display = 'none';
        document.getElementById('batch-table-container').style.display = 'none';
        document.getElementById('empty-state').style.display = 'none';
        document.getElementById('error-state').style.display = 'block';
        document.getElementById('error-message').textContent = message;
    }

    showPagination() {
        if (this.filteredBatches.length > this.pageSize) {
            document.getElementById('pagination-container').style.display = 'block';
        }
    }

    hidePagination() {
        document.getElementById('pagination-container').style.display = 'none';
    }

    renderPagination() {
        const totalPages = Math.ceil(this.filteredBatches.length / this.pageSize);
        if (totalPages <= 1) {
            this.hidePagination();
            return;
        }

        const pagination = document.getElementById('pagination');
        if (!pagination) return;

        let paginationHTML = '';
        
        // Previous button
        paginationHTML += `
            <li class="page-item ${this.currentPage === 1 ? 'disabled' : ''}">
                <a class="page-link" href="#" data-page="${this.currentPage - 1}">&laquo;</a>
            </li>
        `;

        // Page numbers
        const startPage = Math.max(1, this.currentPage - 2);
        const endPage = Math.min(totalPages, this.currentPage + 2);

        for (let i = startPage; i <= endPage; i++) {
            paginationHTML += `
                <li class="page-item ${i === this.currentPage ? 'active' : ''}">
                    <a class="page-link" href="#" data-page="${i}">${i}</a>
                </li>
            `;
        }

        // Next button
        paginationHTML += `
            <li class="page-item ${this.currentPage === totalPages ? 'disabled' : ''}">
                <a class="page-link" href="#" data-page="${this.currentPage + 1}">&raquo;</a>
            </li>
        `;

        pagination.innerHTML = paginationHTML;

        // Bind pagination events
        pagination.querySelectorAll('.page-link').forEach(link => {
            link.addEventListener('click', (e) => {
                e.preventDefault();
                const page = parseInt(e.target.dataset.page);
                if (page && page !== this.currentPage && page >= 1 && page <= totalPages) {
                    this.currentPage = page;
                    this.renderBatches();
                    this.renderPagination();
                }
            });
        });

        this.showPagination();
    }

    // Filter icon management
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

    // Placeholder methods for functionality to be implemented
    async uploadBatch() {
        console.log('Upload batch - to be implemented');
    }

    async runSelected() {
        console.log('Run selected batches:', Array.from(this.selectedBatches));
    }

    async deleteSelected() {
        console.log('Delete selected batches:', Array.from(this.selectedBatches));
    }

    async runSingleBatch(batchId) {
        console.log('Run single batch:', batchId);
    }

    async downloadBatch(batchId) {
        window.open(`/execution/api/batches/${batchId}/download`, '_blank');
    }

    async deleteSingleBatch(batchId) {
        console.log('Delete single batch:', batchId);
    }

    async viewBatchDetails(batchId) {
        console.log('View batch details:', batchId);
    }

    async confirmRun() {
        console.log('Confirm run - to be implemented');
    }

    async confirmDelete() {
        console.log('Confirm delete - to be implemented');
    }

    // Utility function for debouncing
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
}

// Make available globally
window.BatchManager = BatchManager;
