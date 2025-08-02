"""
FastAPI Web Application for CSPBench

Refactored web interface following hexagonal architecture principles.
Modular design with proper separation of concerns.
"""

import logging

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from slowapi import Limiter, _rate_limit_exceeded_handler
from slowapi.errors import RateLimitExceeded
from slowapi.util import get_remote_address

from .core.config import web_config
from .routes import algorithms, batch_execution, datasets, execution, health, pages, results

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Create FastAPI app
app = FastAPI(
    title="CSPBench - Closest String Problem Benchmark",
    description="Web interface for running and comparing CSP algorithms",
    version="1.0.0",
)

# Middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Configure for production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Rate limiting
limiter = Limiter(key_func=get_remote_address)
app.state.limiter = limiter
app.add_exception_handler(RateLimitExceeded, _rate_limit_exceeded_handler)

# Static files
from pathlib import Path

datasets_path = Path(__file__).parent.parent.parent.parent / "datasets"
app.mount(
    "/static", StaticFiles(directory="src/presentation/web/static"), name="static"
)
app.mount("/datasets", StaticFiles(directory=str(datasets_path)), name="datasets")

# Include route modules
app.include_router(health.router)
app.include_router(pages.router)
app.include_router(algorithms.router)
app.include_router(datasets.router)
app.include_router(execution.router)
app.include_router(batch_execution.router)
app.include_router(results.router)

# Include simple datasets router for testing
from .routes import datasets_simple

app.include_router(datasets_simple.router)

# Include debug router for testing
from .routes import test_debug

app.include_router(test_debug.router)


@app.on_event("startup")
async def startup_event():
    """Initialize the application on startup."""
    try:
        logger.info("Starting CSPBench Web Interface initialization...")

        # Initialize configuration and services
        success = web_config.initialize_services()

        if success:
            logger.info("CSPBench Web Interface started successfully")
        else:
            logger.warning(
                "Application started with errors - some features may be unavailable"
            )

    except Exception as e:
        logger.error(f"Startup failed: {e}")
        logger.warning(
            "Application started with errors - some features may be unavailable"
        )


@app.on_event("shutdown")
async def shutdown_event():
    """Clean up on application shutdown."""
    try:
        logger.info("Shutting down CSPBench Web Interface...")
        # Add any cleanup logic here
        logger.info("Shutdown complete")
    except Exception as e:
        logger.error(f"Error during shutdown: {e}")


# Legacy compatibility endpoints (can be removed after migration)
@app.get("/health")
async def legacy_health():
    """Legacy health endpoint for backwards compatibility."""
    return await health.simple_health()


@app.get("/metrics")
async def legacy_metrics():
    """Legacy metrics endpoint for backwards compatibility."""
    return await health.get_metrics()


# For uvicorn direct execution
if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8000)
