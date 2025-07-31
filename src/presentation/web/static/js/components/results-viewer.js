/**
 * Results Viewer Component
 * 
 * Component for displaying and analyzing execution results
 */

class ResultsViewer {
    constructor(container, options = {}) {
        this.container = typeof container === 'string' ? document.querySelector(container) : container;
        this.options = {
            showDownload: true,
            showVisualization: true,
            showComparison: false,
            showExport: true,
            ...options
        };
        
        this.results = null;
        this.comparisonResults = [];
        
        this.init();
    }

    init() {
        this.container.classList.add('results-viewer');
    }

    loadResults(results) {
        this.results = results;
        this.render();
    }

    loadComparisonResults(resultsArray) {
        this.comparisonResults = resultsArray;
        this.renderComparison();
    }

    render() {
        if (!this.results) {
            this.container.innerHTML = `
                <div class="no-results">
                    <div class="no-results-icon">📊</div>
                    <h3>No Results Yet</h3>
                    <p>Execute an algorithm to see results here.</p>
                </div>
            `;
            return;
        }

        this.container.innerHTML = `
            <div class="results-header">
                <h2>Execution Results</h2>
                <div class="results-actions">
                    ${this.options.showDownload ? this.renderDownloadButtons() : ''}
                    ${this.options.showExport ? this.renderExportButtons() : ''}
                </div>
            </div>
            
            <div class="results-content">
                ${this.renderSummary()}
                ${this.renderStringResult()}
                ${this.renderMetrics()}
                ${this.options.showVisualization ? this.renderVisualization() : ''}
                ${this.renderExecutionDetails()}
            </div>
        `;

        this.setupEventListeners();
        
        if (this.options.showVisualization) {
            this.initializeVisualization();
        }
    }

    renderDownloadButtons() {
        return `
            <div class="download-buttons">
                <button class="btn btn-primary" id="download-zip">
                    <i class="icon-download"></i> Download ZIP
                </button>
                <button class="btn btn-outline" id="download-fasta">
                    <i class="icon-file"></i> Download FASTA
                </button>
                <button class="btn btn-outline" id="download-report">
                    <i class="icon-report"></i> Download Report
                </button>
            </div>
        `;
    }

    renderExportButtons() {
        return `
            <div class="export-buttons">
                <button class="btn btn-outline" id="export-json">JSON</button>
                <button class="btn btn-outline" id="export-csv">CSV</button>
                <button class="btn btn-outline" id="export-pdf">PDF</button>
            </div>
        `;
    }

    renderSummary() {
        const { algorithm, dataset, execution_time, status } = this.results;
        
        return `
            <div class="results-summary">
                <div class="summary-cards">
                    <div class="summary-card">
                        <div class="card-label">Algorithm</div>
                        <div class="card-value">${algorithm}</div>
                    </div>
                    
                    <div class="summary-card">
                        <div class="card-label">Dataset</div>
                        <div class="card-value">${dataset}</div>
                    </div>
                    
                    <div class="summary-card">
                        <div class="card-label">Execution Time</div>
                        <div class="card-value">${this.formatTime(execution_time)}</div>
                    </div>
                    
                    <div class="summary-card ${status.toLowerCase()}">
                        <div class="card-label">Status</div>
                        <div class="card-value">${status}</div>
                    </div>
                </div>
            </div>
        `;
    }

    renderStringResult() {
        const { best_string, max_distance } = this.results;
        
        return `
            <div class="string-result">
                <h3>Best String Found</h3>
                
                <div class="string-display">
                    <div class="string-header">
                        <span class="string-label">Result String</span>
                        <span class="string-distance">Distance: ${max_distance}</span>
                        <button class="btn btn-small copy-btn" id="copy-string">Copy</button>
                    </div>
                    
                    <div class="string-content">
                        <div class="string-text" id="result-string">${this.formatString(best_string)}</div>
                    </div>
                    
                    <div class="string-analysis">
                        <div class="analysis-item">
                            <span class="analysis-label">Length:</span>
                            <span class="analysis-value">${best_string.length}</span>
                        </div>
                        
                        <div class="analysis-item">
                            <span class="analysis-label">Unique Characters:</span>
                            <span class="analysis-value">${new Set(best_string).size}</span>
                        </div>
                        
                        <div class="analysis-item">
                            <span class="analysis-label">Character Distribution:</span>
                            <div class="char-distribution" id="char-distribution">
                                ${this.renderCharacterDistribution(best_string)}
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        `;
    }

    renderMetrics() {
        const { metadata } = this.results;
        if (!metadata || Object.keys(metadata).length === 0) {
            return '';
        }

        return `
            <div class="metrics-section">
                <h3>Algorithm Metrics</h3>
                
                <div class="metrics-grid">
                    ${Object.entries(metadata).map(([key, value]) => `
                        <div class="metric-item">
                            <div class="metric-label">${this.formatMetricName(key)}</div>
                            <div class="metric-value">${this.formatMetricValue(value)}</div>
                        </div>
                    `).join('')}
                </div>
                
                ${this.renderMetricsChart()}
            </div>
        `;
    }

    renderMetricsChart() {
        // Only render chart for numeric metrics
        const numericMetrics = Object.entries(this.results.metadata || {})
            .filter(([key, value]) => typeof value === 'number')
            .slice(0, 8); // Limit to 8 metrics for readability

        if (numericMetrics.length === 0) return '';

        return `
            <div class="metrics-chart">
                <h4>Metrics Overview</h4>
                <canvas id="metrics-chart" width="400" height="200"></canvas>
            </div>
        `;
    }

    renderVisualization() {
        return `
            <div class="visualization-section">
                <h3>Visualization</h3>
                
                <div class="visualization-tabs">
                    <button class="tab-btn active" data-tab="distance-matrix">Distance Matrix</button>
                    <button class="tab-btn" data-tab="convergence">Convergence</button>
                    <button class="tab-btn" data-tab="diversity">Diversity</button>
                </div>
                
                <div class="visualization-content">
                    <div class="tab-content active" id="distance-matrix">
                        <canvas id="distance-matrix-canvas"></canvas>
                    </div>
                    
                    <div class="tab-content" id="convergence">
                        <canvas id="convergence-canvas"></canvas>
                    </div>
                    
                    <div class="tab-content" id="diversity">
                        <canvas id="diversity-canvas"></canvas>
                    </div>
                </div>
            </div>
        `;
    }

    renderExecutionDetails() {
        const { parameters, start_time, end_time } = this.results;
        
        return `
            <div class="execution-details">
                <h3>Execution Details</h3>
                
                <div class="details-sections">
                    <div class="details-section">
                        <h4>Parameters Used</h4>
                        <pre class="parameters-display">${JSON.stringify(parameters, null, 2)}</pre>
                    </div>
                    
                    <div class="details-section">
                        <h4>Timing Information</h4>
                        <div class="timing-info">
                            <div class="timing-item">
                                <span class="timing-label">Started:</span>
                                <span class="timing-value">${this.formatDateTime(start_time)}</span>
                            </div>
                            <div class="timing-item">
                                <span class="timing-label">Finished:</span>
                                <span class="timing-value">${this.formatDateTime(end_time)}</span>
                            </div>
                            <div class="timing-item">
                                <span class="timing-label">Duration:</span>
                                <span class="timing-value">${this.formatTime(this.results.execution_time)}</span>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        `;
    }

    renderComparison() {
        if (this.comparisonResults.length === 0) return;

        const comparisonHTML = `
            <div class="comparison-section">
                <h3>Algorithm Comparison</h3>
                
                <div class="comparison-table">
                    <table>
                        <thead>
                            <tr>
                                <th>Algorithm</th>
                                <th>Max Distance</th>
                                <th>Execution Time</th>
                                <th>String Length</th>
                                <th>Parameters</th>
                                <th>Actions</th>
                            </tr>
                        </thead>
                        <tbody>
                            ${this.comparisonResults.map(result => this.renderComparisonRow(result)).join('')}
                        </tbody>
                    </table>
                </div>
                
                <div class="comparison-charts">
                    <canvas id="comparison-chart"></canvas>
                </div>
            </div>
        `;

        this.container.insertAdjacentHTML('beforeend', comparisonHTML);
        this.initializeComparisonChart();
    }

    renderComparisonRow(result) {
        return `
            <tr class="comparison-row">
                <td class="algorithm-name">${result.algorithm}</td>
                <td class="distance-value">${result.max_distance}</td>
                <td class="time-value">${this.formatTime(result.execution_time)}</td>
                <td class="length-value">${result.best_string.length}</td>
                <td class="params-value">
                    <button class="btn btn-small show-params" data-params='${JSON.stringify(result.parameters)}'>
                        View
                    </button>
                </td>
                <td class="actions">
                    <button class="btn btn-small download-result" data-result='${JSON.stringify(result)}'>
                        Download
                    </button>
                </td>
            </tr>
        `;
    }

    setupEventListeners() {
        // Download buttons
        document.getElementById('download-zip')?.addEventListener('click', () => {
            this.downloadZip();
        });

        document.getElementById('download-fasta')?.addEventListener('click', () => {
            this.downloadFasta();
        });

        document.getElementById('download-report')?.addEventListener('click', () => {
            this.downloadReport();
        });

        // Export buttons
        document.getElementById('export-json')?.addEventListener('click', () => {
            this.exportJson();
        });

        document.getElementById('export-csv')?.addEventListener('click', () => {
            this.exportCsv();
        });

        document.getElementById('export-pdf')?.addEventListener('click', () => {
            this.exportPdf();
        });

        // Copy string button
        document.getElementById('copy-string')?.addEventListener('click', () => {
            this.copyString();
        });

        // Visualization tabs
        this.container.querySelectorAll('.tab-btn').forEach(btn => {
            btn.addEventListener('click', (e) => {
                this.switchTab(e.target.dataset.tab);
            });
        });

        // Comparison events
        this.container.querySelectorAll('.show-params').forEach(btn => {
            btn.addEventListener('click', (e) => {
                const params = JSON.parse(e.target.dataset.params);
                this.showParametersModal(params);
            });
        });

        this.container.querySelectorAll('.download-result').forEach(btn => {
            btn.addEventListener('click', (e) => {
                const result = JSON.parse(e.target.dataset.result);
                this.downloadSingleResult(result);
            });
        });
    }

    initializeVisualization() {
        // Initialize visualization charts
        setTimeout(() => {
            if (this.results.metadata) {
                this.drawMetricsChart();
            }
        }, 100);
    }

    drawMetricsChart() {
        const canvas = document.getElementById('metrics-chart');
        if (!canvas) return;

        const ctx = canvas.getContext('2d');
        const numericMetrics = Object.entries(this.results.metadata || {})
            .filter(([key, value]) => typeof value === 'number');

        if (numericMetrics.length === 0) return;

        // Simple bar chart implementation
        const maxValue = Math.max(...numericMetrics.map(([key, value]) => value));
        const barWidth = canvas.width / numericMetrics.length;
        const barHeight = canvas.height - 40;

        ctx.clearRect(0, 0, canvas.width, canvas.height);
        ctx.fillStyle = '#007acc';

        numericMetrics.forEach(([key, value], index) => {
            const height = (value / maxValue) * barHeight;
            const x = index * barWidth;
            const y = canvas.height - height - 20;

            ctx.fillRect(x + 5, y, barWidth - 10, height);

            // Label
            ctx.fillStyle = '#333';
            ctx.font = '12px Arial';
            ctx.textAlign = 'center';
            ctx.fillText(key, x + barWidth / 2, canvas.height - 5);
            ctx.fillText(value.toString(), x + barWidth / 2, y - 5);
            ctx.fillStyle = '#007acc';
        });
    }

    switchTab(tabName) {
        // Update tab buttons
        this.container.querySelectorAll('.tab-btn').forEach(btn => {
            btn.classList.toggle('active', btn.dataset.tab === tabName);
        });

        // Update tab content
        this.container.querySelectorAll('.tab-content').forEach(content => {
            content.classList.toggle('active', content.id === tabName);
        });
    }

    // Utility methods
    formatString(str) {
        if (str.length <= 100) return str;
        return str.substring(0, 100) + '...';
    }

    formatTime(seconds) {
        if (seconds < 60) return `${seconds.toFixed(2)}s`;
        const minutes = Math.floor(seconds / 60);
        const remainingSeconds = seconds % 60;
        return `${minutes}m ${remainingSeconds.toFixed(2)}s`;
    }

    formatDateTime(timestamp) {
        return new Date(timestamp).toLocaleString();
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

    renderCharacterDistribution(str) {
        const charCounts = {};
        for (const char of str) {
            charCounts[char] = (charCounts[char] || 0) + 1;
        }

        return Object.entries(charCounts)
            .sort(([,a], [,b]) => b - a)
            .slice(0, 10) // Show top 10
            .map(([char, count]) => `
                <div class="char-item">
                    <span class="char">${char}</span>
                    <span class="count">${count}</span>
                </div>
            `).join('');
    }

    // Action methods
    async downloadZip() {
        if (!this.results.session_id) return;
        
        try {
            const response = await window.apiClient.downloadSession(this.results.session_id);
            const blob = await response.blob();
            this.downloadBlob(blob, `${this.results.session_id}.zip`);
        } catch (error) {
            console.error('Failed to download ZIP:', error);
        }
    }

    downloadFasta() {
        const fastaContent = `>${this.results.algorithm}_result\n${this.results.best_string}`;
        const blob = new Blob([fastaContent], { type: 'text/plain' });
        this.downloadBlob(blob, `${this.results.algorithm}_result.fasta`);
    }

    downloadReport() {
        const reportContent = this.generateTextReport();
        const blob = new Blob([reportContent], { type: 'text/plain' });
        this.downloadBlob(blob, `${this.results.algorithm}_report.txt`);
    }

    exportJson() {
        const blob = new Blob([JSON.stringify(this.results, null, 2)], { type: 'application/json' });
        this.downloadBlob(blob, `${this.results.algorithm}_results.json`);
    }

    exportCsv() {
        const csvContent = this.generateCsvReport();
        const blob = new Blob([csvContent], { type: 'text/csv' });
        this.downloadBlob(blob, `${this.results.algorithm}_results.csv`);
    }

    exportPdf() {
        // This would require a PDF library like jsPDF
        alert('PDF export feature coming soon!');
    }

    copyString() {
        navigator.clipboard.writeText(this.results.best_string).then(() => {
            const btn = document.getElementById('copy-string');
            const originalText = btn.textContent;
            btn.textContent = 'Copied!';
            setTimeout(() => {
                btn.textContent = originalText;
            }, 2000);
        });
    }

    downloadBlob(blob, filename) {
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = filename;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
    }

    generateTextReport() {
        return `
CSP Algorithm Execution Report
=============================

Algorithm: ${this.results.algorithm}
Dataset: ${this.results.dataset}
Execution Time: ${this.formatTime(this.results.execution_time)}
Status: ${this.results.status}

Results:
--------
Best String: ${this.results.best_string}
Max Distance: ${this.results.max_distance}
String Length: ${this.results.best_string.length}

Parameters:
-----------
${JSON.stringify(this.results.parameters, null, 2)}

Metadata:
---------
${JSON.stringify(this.results.metadata, null, 2)}

Generated on: ${new Date().toLocaleString()}
        `.trim();
    }

    generateCsvReport() {
        const rows = [
            ['Metric', 'Value'],
            ['Algorithm', this.results.algorithm],
            ['Dataset', this.results.dataset],
            ['Execution Time', this.results.execution_time],
            ['Max Distance', this.results.max_distance],
            ['String Length', this.results.best_string.length],
            ...Object.entries(this.results.metadata || {})
        ];

        return rows.map(row => row.join(',')).join('\n');
    }

    showParametersModal(params) {
        alert(`Parameters:\n\n${JSON.stringify(params, null, 2)}`);
    }

    downloadSingleResult(result) {
        const blob = new Blob([JSON.stringify(result, null, 2)], { type: 'application/json' });
        this.downloadBlob(blob, `${result.algorithm}_result.json`);
    }

    clear() {
        this.results = null;
        this.comparisonResults = [];
        this.render();
    }
}
