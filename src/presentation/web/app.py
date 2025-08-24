"""
FastAPI Web Application for CSPBench

Refactored web interface following hexagonal architecture principles.
Modular design with proper separation of concerns and global WorkManager.
"""

import logging
from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles

# Import WorkService lifespan
from src.application.services.work_service import work_service_lifespan

Limiter = None
_rate_limit_exceeded_handler = None
RateLimitExceeded = Exception


def get_remote_address(request):
    """Extract remote address from request."""
    forwarded_for = request.headers.get("X-Forwarded-For")
    if forwarded_for:
        return forwarded_for.split(",")[0].strip()
    return request.client.host


from .core.config import web_config
from .routes import algorithms, datasets, health, pages, monitoring
from . import websocket_routes

# Get logger (logging is configured globally in main.py via LoggerConfig)
logger = logging.getLogger(__name__)

# Create FastAPI app with global WorkManager lifespan
app = FastAPI(
    title="CSPBench - Closest String Problem Benchmark",
    description="Web interface for running and comparing CSP algorithms",
    version="1.0.0",
    docs_url=None,  # Disable API docs
    redoc_url=None,  # Disable ReDoc
    lifespan=work_service_lifespan,  # Add WorkService lifecycle
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
# Rate limiting disabled in this build

# Static files
datasets_path = Path(__file__).parent.parent.parent.parent / "datasets"
app.mount(
    "/static", StaticFiles(directory="src/presentation/web/static"), name="static"
)
app.mount("/datasets", StaticFiles(directory=str(datasets_path)), name="datasets")

# Include route modules
# Include routers
app.include_router(algorithms.router, tags=["algorithms"])
app.include_router(datasets.router, tags=["datasets"])
app.include_router(health.router, tags=["health"])
app.include_router(pages.router, tags=["pages"])
app.include_router(monitoring.router, tags=["monitoring"])
app.include_router(websocket_routes.router, tags=["websockets"])

# Include progress monitoring routes
from .routes import progress_monitor
app.include_router(progress_monitor.router, tags=["progress-monitoring"])

# Import and include batch routes
from .routes import batches

app.include_router(batches.router)

# Import and include batch execution routes
from .routes import batch_execution

app.include_router(batch_execution.router)

# Import and include file routes
from .routes import files

app.include_router(files.router)


@app.get("/test-progress")
async def test_progress_page():
    """Serve test progress page for debugging."""
    try:
        with open("src/presentation/web/static/test-progress.html", "r") as f:
            content = f.read()
        return HTMLResponse(content=content)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail="Test progress page not found")


@app.get("/monitor")
async def progress_monitor_page():
    """Serve progress monitor page."""
    try:
        with open("src/presentation/web/templates/progress_monitor.html", "r") as f:
            content = f.read()
        return HTMLResponse(content=content)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail="Progress monitor page not found")


# Remove old startup/shutdown events as they're replaced by lifespan
# Global WorkManager is now managed by the lifespan context manager


# Legacy endpoints removidos


# For uvicorn direct execution
if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8000)


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


# Legacy endpoints removidos


# For uvicorn direct execution
if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8000)
