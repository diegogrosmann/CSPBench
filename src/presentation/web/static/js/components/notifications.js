(function(global){
/**
 * Modern Notification System for CSPBench
 * Bootstrap-based toast notifications with animations
 *
 * This file is idempotent and safe to load multiple times.
 */

// Define NotificationManager only once
if (!global.NotificationManager) {
class NotificationManager {
    constructor() {
        this.container = null;
        this.notifications = new Map();
        this.defaultDuration = 5000;
        this.maxNotifications = 5;
        
        this.init();
    }

    init() {
        this.createContainer();
        this.setupStyles();
    }

    createContainer() {
        // Create toast container if it doesn't exist
        this.container = document.querySelector('.toast-container');
        
        if (!this.container) {
            this.container = document.createElement('div');
            this.container.className = 'toast-container position-fixed top-0 end-0 p-3';
            this.container.style.zIndex = '9999';
            document.body.appendChild(this.container);
        }
    }

    setupStyles() {
        // Add custom styles for notifications
        const style = document.createElement('style');
        style.textContent = `
            .toast-custom {
                backdrop-filter: blur(10px);
                border: none;
                box-shadow: 0 0.5rem 1rem rgba(0, 0, 0, 0.15);
                border-radius: 0.75rem;
                overflow: hidden;
            }
            
            .toast-custom .toast-header {
                border-bottom: none;
                padding: 0.75rem 1rem;
                background: rgba(255, 255, 255, 0.95);
            }
            
            .toast-custom .toast-body {
                padding: 0.75rem 1rem;
                background: rgba(255, 255, 255, 0.95);
            }
            
            .toast-custom.toast-success {
                border-left: 4px solid #28a745;
            }
            
            .toast-custom.toast-error {
                border-left: 4px solid #dc3545;
            }
            
            .toast-custom.toast-warning {
                border-left: 4px solid #ffc107;
            }
            
            .toast-custom.toast-info {
                border-left: 4px solid #17a2b8;
            }
            
            .toast-progress {
                height: 3px;
                background: rgba(0, 0, 0, 0.1);
                position: relative;
                overflow: hidden;
            }
            
            .toast-progress-bar {
                height: 100%;
                background: currentColor;
                transition: width 0.1s linear;
                opacity: 0.7;
            }
            
            .notification-enter {
                animation: slideInRight 0.3s ease-out;
            }
            
            .notification-exit {
                animation: slideOutRight 0.3s ease-in;
            }
            
            @keyframes slideInRight {
                from {
                    transform: translateX(100%);
                    opacity: 0;
                }
                to {
                    transform: translateX(0);
                    opacity: 1;
                }
            }
            
            @keyframes slideOutRight {
                from {
                    transform: translateX(0);
                    opacity: 1;
                }
                to {
                    transform: translateX(100%);
                    opacity: 0;
                }
            }
        `;
        
        if (!document.querySelector('#notification-styles')) {
            style.id = 'notification-styles';
            document.head.appendChild(style);
        }
    }

    show(message, type = 'info', options = {}) {
        const config = {
            duration: this.defaultDuration,
            persistent: false,
            showProgress: true,
            title: null,
            actions: [],
            ...options
        };

        const id = this.generateId();
        const toast = this.createToast(id, message, type, config);
        
        // Limit number of notifications
        this.limitNotifications();
        
        this.container.appendChild(toast);
        this.notifications.set(id, { element: toast, config });

        // Show the toast
        const bsToast = new bootstrap.Toast(toast, {
            autohide: !config.persistent,
            delay: config.duration
        });

        bsToast.show();

        // Add progress bar if enabled
        if (config.showProgress && !config.persistent) {
            this.addProgressBar(toast, config.duration);
        }

        // Handle actions
        this.bindActions(toast, config.actions);

        // Clean up after hide
        toast.addEventListener('hidden.bs.toast', () => {
            this.remove(id);
        });

        return id;
    }

    createToast(id, message, type, config) {
        const toast = document.createElement('div');
        toast.className = `toast toast-custom toast-${type} notification-enter`;
        toast.id = `notification-${id}`;
        toast.setAttribute('role', 'alert');
        toast.setAttribute('aria-live', 'assertive');
        toast.setAttribute('aria-atomic', 'true');

        const icon = this.getIcon(type);
        const title = config.title || this.getDefaultTitle(type);

        toast.innerHTML = `
            <div class="toast-header">
                <i class="bi bi-${icon} me-2 text-${type}"></i>
                <strong class="me-auto">${title}</strong>
                <small class="text-muted">${this.getTimeStamp()}</small>
                <button type="button" class="btn-close" data-bs-dismiss="toast" aria-label="Close"></button>
            </div>
            <div class="toast-body">
                ${message}
                ${this.renderActions(config.actions)}
            </div>
            ${config.showProgress && !config.persistent ? '<div class="toast-progress"><div class="toast-progress-bar"></div></div>' : ''}
        `;

        return toast;
    }

    renderActions(actions) {
        if (!actions || actions.length === 0) return '';

        const actionsHtml = actions.map(action => {
            const btnClass = action.variant ? `btn-${action.variant}` : 'btn-outline-primary';
            return `
                <button type="button" class="btn btn-sm ${btnClass} me-2" data-action="${action.id}">
                    ${action.icon ? `<i class="bi bi-${action.icon} me-1"></i>` : ''}
                    ${action.label}
                </button>
            `;
        }).join('');

        return `<div class="mt-2">${actionsHtml}</div>`;
    }

    bindActions(toast, actions) {
        if (!actions || actions.length === 0) return;

        actions.forEach(action => {
            const button = toast.querySelector(`[data-action="${action.id}"]`);
            if (button && action.handler) {
                button.addEventListener('click', action.handler);
            }
        });
    }

    addProgressBar(toast, duration) {
        const progressBar = toast.querySelector('.toast-progress-bar');
        if (!progressBar) return;

        let width = 100;
        const interval = 50; // Update every 50ms
        const decrement = (interval / duration) * 100;

        const timer = setInterval(() => {
            width -= decrement;
            progressBar.style.width = `${Math.max(0, width)}%`;

            if (width <= 0) {
                clearInterval(timer);
            }
        }, interval);

        // Store timer for cleanup
        toast._progressTimer = timer;
    }

    limitNotifications() {
        const notifications = Array.from(this.notifications.entries());
        
        if (notifications.length >= this.maxNotifications) {
            const oldestId = notifications[0][0];
            this.remove(oldestId);
        }
    }

    remove(id) {
        const notification = this.notifications.get(id);
        if (!notification) return;

        const { element } = notification;

        // Clear progress timer if exists
        if (element._progressTimer) {
            clearInterval(element._progressTimer);
        }

        // Add exit animation
        element.classList.add('notification-exit');

        // Remove after animation
        setTimeout(() => {
            if (element.parentNode) {
                element.parentNode.removeChild(element);
            }
            this.notifications.delete(id);
        }, 300);
    }

    clear() {
        for (const [id] of this.notifications) {
            this.remove(id);
        }
    }

    // Convenience methods
    success(message, options = {}) {
        return this.show(message, 'success', options);
    }

    error(message, options = {}) {
        return this.show(message, 'error', { ...options, duration: 8000 });
    }

    warning(message, options = {}) {
        return this.show(message, 'warning', options);
    }

    info(message, options = {}) {
        return this.show(message, 'info', options);
    }

    // Utility methods
    generateId() {
        return Date.now().toString(36) + Math.random().toString(36).substr(2);
    }

    getIcon(type) {
        const icons = {
            success: 'check-circle-fill',
            error: 'exclamation-triangle-fill',
            warning: 'exclamation-triangle-fill',
            info: 'info-circle-fill'
        };
        return icons[type] || 'info-circle-fill';
    }

    getDefaultTitle(type) {
        const titles = {
            success: 'Success',
            error: 'Error',
            warning: 'Warning',
            info: 'Information'
        };
        return titles[type] || 'Notification';
    }

    getTimeStamp() {
        return new Date().toLocaleTimeString([], { 
            hour: '2-digit', 
            minute: '2-digit' 
        });
    }
}
global.NotificationManager = NotificationManager;
}

// Progress notification for long-running operations
if (!global.ProgressNotification) {
class ProgressNotification {
    constructor(title = 'Processing...', message = '') {
        this.title = title;
        this.message = message;
        this.id = null;
        this.toast = null;
        this.progressBar = null;
    }

    show() {
        const container = document.querySelector('.toast-container') || 
                        notifications.container;

        this.id = notifications.generateId();
        this.toast = document.createElement('div');
        this.toast.className = 'toast toast-custom notification-enter';
        this.toast.id = `progress-${this.id}`;
        this.toast.setAttribute('role', 'alert');
        this.toast.setAttribute('aria-live', 'polite');

        this.toast.innerHTML = `
            <div class="toast-header">
                <div class="spinner-border spinner-border-sm text-primary me-2" role="status">
                    <span class="visually-hidden">Loading...</span>
                </div>
                <strong class="me-auto">${this.title}</strong>
                <button type="button" class="btn-close" data-bs-dismiss="toast"></button>
            </div>
            <div class="toast-body">
                <div class="progress-message">${this.message}</div>
                <div class="progress mt-2" style="height: 8px;">
                    <div class="progress-bar bg-primary" role="progressbar" style="width: 0%"></div>
                </div>
            </div>
        `;

        container.appendChild(this.toast);
        this.progressBar = this.toast.querySelector('.progress-bar');

        const bsToast = new bootstrap.Toast(this.toast, {
            autohide: false
        });
        bsToast.show();

        return this;
    }

    updateProgress(percentage, message = null) {
        if (this.progressBar) {
            this.progressBar.style.width = `${Math.min(100, Math.max(0, percentage))}%`;
        }

        if (message && this.toast) {
            const messageElement = this.toast.querySelector('.progress-message');
            if (messageElement) {
                messageElement.textContent = message;
            }
        }
    }

    complete(message = 'Complete!', type = 'success') {
        if (this.toast) {
            // Convert to success toast
            const header = this.toast.querySelector('.toast-header');
            const icon = type === 'success' ? 'check-circle-fill' : 'exclamation-triangle-fill';
            header.innerHTML = `
                <i class="bi bi-${icon} me-2 text-${type}"></i>
                <strong class="me-auto">${message}</strong>
                <button type="button" class="btn-close" data-bs-dismiss="toast"></button>
            `;

            // Hide progress bar
            const progressContainer = this.toast.querySelector('.progress');
            if (progressContainer) {
                progressContainer.style.display = 'none';
            }

            // Auto-hide after 3 seconds
            setTimeout(() => {
                const bsToast = bootstrap.Toast.getInstance(this.toast);
                if (bsToast) {
                    bsToast.hide();
                }
            }, 3000);
        }
    }

    hide() {
        if (this.toast) {
            const bsToast = bootstrap.Toast.getInstance(this.toast);
            if (bsToast) {
                bsToast.hide();
            }
        }
    }
}
global.ProgressNotification = ProgressNotification;
}

// Create a single global notifications instance if not existing
if (!global.notifications) {
    global.notifications = new global.NotificationManager();
}

// Replace the global showAlert function
global.showAlert = function(message, type = 'info', duration = 5000) {
    if (global.notifications) {
        global.notifications.show(message, type, { duration });
    } else {
        console[type === 'error' ? 'error' : 'log'](message);
    }
};

})(window);
