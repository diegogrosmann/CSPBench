"""
Web routes and handlers organized by domain.

This package contains all FastAPI route handlers for the CSPBench web interface,
organized by functional domains to maintain clean separation of concerns and
facilitate maintainability.

Route Modules:
    pages: HTML page routes for serving Jinja2 templates
    health: System health check and monitoring endpoints  
    algorithms: Algorithm discovery and metadata endpoints
    datasets: Dataset management CRUD operations
    batches: Batch configuration file management
    batch_execution: Batch execution lifecycle management
    monitoring: Work execution monitoring and status tracking
    files: File download and result access endpoints

Architecture:
    - Each module contains routes for a specific domain
    - All routes use proper error handling and logging
    - Security validation is applied consistently
    - Pydantic models ensure type safety and validation
    - Comprehensive documentation follows Google/NumPy style

Usage:
    Routes are automatically registered with the FastAPI application
    through the main application factory. Each module exports a
    router instance that is included in the main app configuration.

Security:
    - Input validation through Pydantic models
    - Path sanitization for file operations  
    - Authorization checks where applicable
    - Rate limiting and abuse prevention
"""
