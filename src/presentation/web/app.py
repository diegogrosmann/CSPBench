"""FastAPI Web Application for CSPBench.

This module provides the main FastAPI application instance with comprehensive
configuration, middleware setup, and route registration. It implements a minimal
but feature-complete web interface with deterministic lifecycle management
through unified service contexts.
"""

import logging
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles

from src.application.services.work_service import work_service_lifespan

from .routes import algorithms, datasets, health, pages

logger = logging.getLogger(__name__)


@asynccontextmanager
async def combined_lifespan(app):  # pragma: no cover - framework integration
    """Combined application lifecycle manager.
    
    This async context manager handles the complete application lifecycle,
    including service initialization, dependency setup, and graceful shutdown.
    It integrates with the WorkService lifecycle to ensure proper resource
    management throughout the application's runtime.
    
    Args:
        app: FastAPI application instance (unused in current implementation)
        
    Yields:
        None: Application is ready for serving requests
        
    Note:
        This function is marked as no-cover since it's framework integration
        code that's difficult to test in isolation.
    """
    # Enter WorkService lifespan first to initialize core services
    async with work_service_lifespan(app):
        # Web application is ready for requests
        logger.info("Web application lifecycle initialized")
        yield
        # Future shutdown hooks can be added here
        logger.info("Web application lifecycle completed")


# Create FastAPI application instance with comprehensive metadata
app = FastAPI(
    title="CSPBench - Closest String Problem Benchmark",
    description="Comprehensive web interface for running and comparing CSP algorithms with real-time monitoring",
    version="1.0.0",
    docs_url=None,  # Disabled for security in production
    redoc_url=None,  # Disabled for security in production  
    lifespan=combined_lifespan,
    # Additional metadata for API documentation
    contact={
        "name": "CSPBench Development Team",
        "url": "https://github.com/your-org/cspbench",
    },
    license_info={
        "name": "Academic Research License",
    },
)

# Configure CORS middleware for cross-origin requests
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # TODO: Restrict to specific origins in production
    allow_credentials=True,
    allow_methods=["GET", "POST", "PUT", "DELETE", "OPTIONS"],
    allow_headers=["*"],
    # Performance and security headers
    expose_headers=["X-Total-Count", "X-Pagination-Page", "X-Pagination-Pages"],
)

# Configure static file serving
app.mount(
    "/static", 
    StaticFiles(directory="src/presentation/web/static"), 
    name="static"
)

# Configure dataset directory mounting with environment-aware fallbacks
import os

# Determine dataset directory with multiple fallback options
datasets_path = Path(__file__).parent.parent.parent.parent / "datasets"
dataset_dir_env = os.getenv("DATASET_DIRECTORY", "/app/data/datasets")
data_datasets_path = Path(dataset_dir_env)

# Priority-based dataset directory resolution
if data_datasets_path.exists():
    # Environment-configured directory (production)
    app.mount(
        "/datasets", 
        StaticFiles(directory=str(data_datasets_path)), 
        name="datasets"
    )
    logger.info(f"Mounted datasets from configured directory: {data_datasets_path}")
elif datasets_path.exists():
    # Project-relative directory (development)
    app.mount("/datasets", StaticFiles(directory=str(datasets_path)), name="datasets")
    logger.info(f"Mounted datasets from project directory: {datasets_path}")
else:
    # Create fallback directory structure
    try:
        data_datasets_path.mkdir(parents=True, exist_ok=True)
        app.mount(
            "/datasets", 
            StaticFiles(directory=str(data_datasets_path)), 
            name="datasets"
        )
        logger.info(f"Created and mounted datasets directory: {data_datasets_path}")
    except (PermissionError, OSError) as e:
        logger.warning(f"Could not create configured datasets directory {data_datasets_path}: {e}")
        # Final fallback to local directory
        fallback_path = Path("./data/datasets")
        fallback_path.mkdir(parents=True, exist_ok=True)
        app.mount(
            "/datasets", 
            StaticFiles(directory=str(fallback_path)), 
            name="datasets"
        )
        logger.info(f"Using fallback datasets directory: {fallback_path}")

# Register core API routers with organized endpoint structure
app.include_router(algorithms.router, tags=["algorithms"])
app.include_router(datasets.router, tags=["datasets"]) 
app.include_router(health.router, tags=["health"])
app.include_router(pages.router, tags=["pages"])

# Register batch management routers
from .routes import batches  # noqa: E402
app.include_router(batches.router)

from .routes import batch_execution  # noqa: E402
app.include_router(batch_execution.router)

# Register file management and download routers
from .routes import files  # noqa: E402
app.include_router(files.router)

# Register WebSocket routes for real-time monitoring
from .websocket import routes as websocket_routes  # noqa: E402
app.include_router(websocket_routes.router)

# Register monitoring and progress tracking routes
from .routes import monitoring  # noqa: E402
app.include_router(monitoring.router)

# Development server entry point
if __name__ == "__main__":  # pragma: no cover
    """Development server entry point.
    
    Provides a simple way to run the application during development.
    For production deployment, use the run_web.py script or proper
    ASGI server configuration.
    """
    import uvicorn

    uvicorn.run(
        app, 
        host="0.0.0.0", 
        port=8000,
        reload=True,
        log_level="info"
    )
