"""
FastAPI Web Application for CSPBench

Refactored web interface following hexagonal architecture principles.
Modular design with proper separation of concerns.
"""

import logging
from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles

Limiter = None
_rate_limit_exceeded_handler = None
RateLimitExceeded = Exception


def get_remote_address(request):
    return "0.0.0.0"


from .core.config import web_config
from .routes import algorithms, batch_execution, datasets, health, pages, results

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Create FastAPI app
app = FastAPI(
    title="CSPBench - Closest String Problem Benchmark",
    description="Web interface for running and comparing CSP algorithms",
    version="1.0.0",
    docs_url=None,  # Disable API docs
    redoc_url=None,  # Disable ReDoc
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
app.include_router(health.router)
app.include_router(pages.router)
app.include_router(algorithms.router)
app.include_router(datasets.router)
app.include_router(batch_execution.router)
app.include_router(results.router)

# Include simple datasets router for testing
from .routes import datasets_simple

app.include_router(datasets_simple.router)


@app.get("/test-progress")
async def test_progress_page():
    """Serve test progress page for debugging."""
    try:
        with open("/workspaces/CSPBench/test_progress.html") as f:
            html_content = f.read()
        return HTMLResponse(content=html_content)
    except Exception as e:
        raise HTTPException(
            status_code=500, detail=f"Error loading test page: {e}"
        ) from e


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
