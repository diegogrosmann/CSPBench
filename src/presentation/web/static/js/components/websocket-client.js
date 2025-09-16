/**
 * WebSocket Client for Real-time Monitoring
 * Handles WebSocket connections and real-time updates
 */

class WebSocketClient {
    constructor(clientId = null) {
        this.clientId = clientId || this.generateClientId();
        this.ws = null;
        this.reconnectAttempts = 0;
        this.maxReconnectAttempts = 5;
        this.reconnectDelay = 1000; // Start with 1 second
        this.subscribers = new Map(); // event -> [callbacks]
        this.connected = false;
        this.connecting = false;

        // Fila assíncrona de mensagens
        this.messageQueue = [];
        this.processingQueue = false;
        this.maxBatch = 50; // processar no máx. 50 msgs por frame

        // Bind methods
        this.connect = this.connect.bind(this);
        this.disconnect = this.disconnect.bind(this);
        this.send = this.send.bind(this);
        this.onMessage = this.onMessage.bind(this);
        this.onOpen = this.onOpen.bind(this);
        this.onClose = this.onClose.bind(this);
        this.onError = this.onError.bind(this);
        this.processMessages = this.processMessages.bind(this);
    }

    generateClientId() {
        return 'client_' + Math.random().toString(36).substr(2, 9);
    }

    async connect() {
        if (this.connected || this.connecting) {
            console.log('WebSocket already connected or connecting');
            return;
        }

        this.connecting = true;

        try {
            const protocol = window.location.protocol === 'https:' ? 'wss:' : 'ws:';
            const wsUrl = `${protocol}//${window.location.host}/ws/${this.clientId}`;

            console.log('Connecting to WebSocket:', wsUrl);

            this.ws = new WebSocket(wsUrl);
            this.ws.onopen = this.onOpen;
            this.ws.onmessage = this.onMessage;
            this.ws.onclose = this.onClose;
            this.ws.onerror = this.onError;
        } catch (error) {
            console.error('WebSocket connection error:', error);
            this.connecting = false;
            this.scheduleReconnect();
        }
    }

    disconnect() {
        console.log('Disconnecting WebSocket');
        this.reconnectAttempts = this.maxReconnectAttempts; // Prevent reconnection

        if (this.ws) {
            this.ws.close();
            this.ws = null;
        }

        this.connected = false;
        this.connecting = false;
    }

    onOpen(event) {
        console.log('WebSocket connected');
        this.connected = true;
        this.connecting = false;
        this.reconnectAttempts = 0;
        this.reconnectDelay = 1000;

        this.emit('connected', { clientId: this.clientId });

        // Send ping to confirm connection
        this.ping();
    }

    onMessage(event) {
        try {
            // Em vez de parsear e emitir sincronicamente, enfileira
            this.messageQueue.push(event.data);
            this.processMessages();
        } catch (error) {
            console.error('Error enqueuing WebSocket message:', error);
        }
    }

    processMessages() {
        if (this.processingQueue) return;
        this.processingQueue = true;

        const step = () => {
            let count = 0;
            while (this.messageQueue.length && count < this.maxBatch) {
                const raw = this.messageQueue.shift();
                try {
                    const data = typeof raw === 'string' ? JSON.parse(raw) : raw;
                    // Emitir eventos
                    this.emit(data.type, data);
                    this.emit('message', data);
                } catch (e) {
                    console.error('Error parsing WebSocket message:', e);
                }
                count++;
            }
            if (this.messageQueue.length) {
                // Ainda há mensagens: processa no próximo frame
                requestAnimationFrame(step);
            } else {
                this.processingQueue = false;
            }
        };

        requestAnimationFrame(step);
    }

    onClose(event) {
        console.log('WebSocket closed:', event.code, event.reason);
        this.connected = false;
        this.connecting = false;

        this.emit('disconnected', { code: event.code, reason: event.reason });

        // Attempt to reconnect if not intentionally closed
        if (event.code !== 1000 && this.reconnectAttempts < this.maxReconnectAttempts) {
            this.scheduleReconnect();
        }
    }

    onError(error) {
        console.error('WebSocket error:', error);
        this.emit('error', error);
    }

    scheduleReconnect() {
        if (this.reconnectAttempts >= this.maxReconnectAttempts) {
            console.log('Max reconnection attempts reached');
            this.emit('max_reconnect_attempts');
            return;
        }

        this.reconnectAttempts++;
        const delay = this.reconnectDelay * Math.pow(2, this.reconnectAttempts - 1); // Exponential backoff

        console.log(`Scheduling reconnection attempt ${this.reconnectAttempts} in ${delay}ms`);

        setTimeout(() => {
            console.log(`Reconnection attempt ${this.reconnectAttempts}`);
            this.connect();
        }, delay);
    }

    send(message) {
        if (!this.connected || !this.ws) {
            console.warn('WebSocket not connected, cannot send message:', message);
            return false;
        }

        try {
            const messageStr = typeof message === 'string' ? message : JSON.stringify(message);
            this.ws.send(messageStr);
            return true;
        } catch (error) {
            console.error('Error sending WebSocket message:', error);
            return false;
        }
    }

    // High-level methods
    subscribeToWork(workId) {
        return this.send({
            type: 'subscribe_work',
            work_id: workId
        });
    }

    unsubscribeFromWork(workId) {
        return this.send({
            type: 'unsubscribe_work',
            work_id: workId
        });
    }

    ping() {
        return this.send({
            type: 'ping',
            timestamp: new Date().toISOString()
        });
    }

    // Event system
    on(event, callback) {
        if (!this.subscribers.has(event)) {
            this.subscribers.set(event, []);
        }
        this.subscribers.get(event).push(callback);
    }

    off(event, callback) {
        if (this.subscribers.has(event)) {
            const callbacks = this.subscribers.get(event);
            const index = callbacks.indexOf(callback);
            if (index > -1) {
                callbacks.splice(index, 1);
            }
        }
    }

    emit(event, data) {
        if (this.subscribers.has(event)) {
            this.subscribers.get(event).forEach(callback => {
                try {
                    callback(data);
                } catch (error) {
                    console.error(`Error in event callback for ${event}:`, error);
                }
            });
        }
    }

    // Status methods
    isConnected() {
        return this.connected && this.ws && this.ws.readyState === WebSocket.OPEN;
    }

    getStatus() {
        return {
            connected: this.connected,
            connecting: this.connecting,
            reconnectAttempts: this.reconnectAttempts,
            clientId: this.clientId,
            readyState: this.ws ? this.ws.readyState : null
        };
    }
}

// Export for use in other modules
if (typeof window !== 'undefined') {
    window.WebSocketClient = WebSocketClient;
}
