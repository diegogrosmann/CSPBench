"""
Web Presentation Layer - FastAPI-based Interface for CSPBench.

This package provides a comprehensive web-based user interface for the CSPBench
framework, built using FastAPI with modern web technologies and real-time
capabilities through WebSocket connections.

Key Components:
    - FastAPI Application: Main web server with RESTful API endpoints
    - WebSocket Infrastructure: Real-time progress monitoring and communication
    - Route Handlers: Organized endpoint handlers by functional domain
    - Security Layer: Input validation, sanitization, and access control
    - Static Assets: Frontend resources and templates

Architecture:
    The web interface follows a modular, domain-driven design:
    
    - Routes: Domain-specific endpoint handlers (datasets, algorithms, etc.)
    - Core: Shared models, security utilities, and configuration
    - WebSocket: Real-time communication infrastructure
    - Templates: Jinja2 templates for server-side rendering
    - Static: CSS, JavaScript, and asset files

Features:
    - RESTful API with comprehensive OpenAPI documentation
    - Real-time progress monitoring via WebSocket connections
    - Secure file upload and download capabilities
    - Batch execution management and monitoring
    - Interactive dataset generation and management
    - Results visualization and export functionality
    - Responsive design for multiple device types

Security Measures:
    - Input validation through Pydantic models
    - File path sanitization and access control
    - CORS configuration for cross-origin security
    - Rate limiting and request size restrictions
    - Secure WebSocket connection handling

Development Features:
    - Hot reload capability for rapid development
    - Comprehensive error handling and logging
    - Environment-based configuration management
    - Docker support for consistent deployment
    - Development server with debugging capabilities

Example Usage:
    Start the development server::
    
        from src.presentation.web.app import app
        import uvicorn
        
        uvicorn.run(app, host="0.0.0.0", port=8000, reload=True)
        
    Or use the provided launcher::
    
        python -m src.presentation.web.run_web

WebSocket Endpoints:
    - `/ws/{client_id}`: General WebSocket for testing and development
    - `/ws/work/{work_id}`: Real-time work progress monitoring

API Endpoints:
    - `/api/algorithms`: Algorithm discovery and metadata
    - `/api/datasets`: Dataset management and generation
    - `/api/batches`: Batch configuration management
    - `/api/batch`: Batch execution lifecycle
    - `/api/monitor`: Execution monitoring and status
    - `/api/files`: File download and result access
    - `/health`: System health and status checks

Deployment:
    The web interface supports multiple deployment scenarios:
    
    - Development: Built-in uvicorn server with hot reload
    - Production: Docker container with optimized settings
    - Cloud: Kubernetes deployment with horizontal scaling
    - Edge: Single-node deployment with minimal resources
"""

__version__ = "1.0.0"
__all__ = ["app"]

# Import main application for easy access
from .app import app
