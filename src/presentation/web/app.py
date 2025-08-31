"""FastAPI Web Application for CSPBench.

Provides a minimal deterministic web interface. Lifecycle is handled via a
single lifespan context combining WorkService lifecycle and web_config service
initialization (no deprecated @on_event handlers).
"""

from contextlib import asynccontextmanager
import logging
from pathlib import Path

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles

from src.application.services.work_service import work_service_lifespan
from .core.config import get_web_config
from .routes import algorithms, datasets, health, pages

logger = logging.getLogger(__name__)


@asynccontextmanager
async def combined_lifespan(app):  # pragma: no cover - framework integration
    # Enter WorkService lifespan first
    async with work_service_lifespan(app):
        # Additional startup (web_config)
        try:
            logger.info("Initializing web_config services...")
            success = web_config.initialize_services()
            if success:
                logger.info("web_config services initialized")
            else:
                logger.warning("web_config initialization reported issues")
        except Exception as e:  # noqa: BLE001
            logger.error(f"web_config initialization failed: {e}")
        yield
        # Optional shutdown hooks here in future


app = FastAPI(
    title="CSPBench - Closest String Problem Benchmark",
    description="Web interface for running and comparing CSP algorithms",
    version="1.0.0",
    docs_url=None,
    redoc_url=None,
    lifespan=combined_lifespan,
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # TODO: restrict in production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Static files
app.mount("/static", StaticFiles(directory="src/presentation/web/static"), name="static")

# Mount datasets directory only if it exists
import os
datasets_path = Path(__file__).parent.parent.parent.parent / "datasets"
# Use environment variable for dataset directory, fallback to Docker path for production
dataset_dir_env = os.getenv("DATASET_DIRECTORY", "/app/data/datasets")
data_datasets_path = Path(dataset_dir_env)

# Prefer the configured directory, fallback to project datasets
if data_datasets_path.exists():
    app.mount("/datasets", StaticFiles(directory=str(data_datasets_path)), name="datasets")
elif datasets_path.exists():
    app.mount("/datasets", StaticFiles(directory=str(datasets_path)), name="datasets")
else:
    # Create the configured datasets directory if it doesn't exist and we have permissions
    try:
        data_datasets_path.mkdir(parents=True, exist_ok=True)
        app.mount("/datasets", StaticFiles(directory=str(data_datasets_path)), name="datasets")
    except (PermissionError, OSError) as e:
        logger.warning(f"Could not create datasets directory {data_datasets_path}: {e}")
        # Fallback to a local datasets directory
        fallback_path = Path("./data/datasets")
        fallback_path.mkdir(parents=True, exist_ok=True)
        app.mount("/datasets", StaticFiles(directory=str(fallback_path)), name="datasets")

# Routers
app.include_router(algorithms.router, tags=["algorithms"])
app.include_router(datasets.router, tags=["datasets"])
app.include_router(health.router, tags=["health"])
app.include_router(pages.router, tags=["pages"])

from .routes import batches  # noqa: E402

app.include_router(batches.router)

from .routes import batch_execution  # noqa: E402

app.include_router(batch_execution.router)

from .routes import files  # noqa: E402

app.include_router(files.router)

from .websocket import routes as websocket_routes  # noqa: E402

app.include_router(websocket_routes.router)

from .routes import monitoring  # noqa: E402

app.include_router(monitoring.router)


if __name__ == "__main__":  # pragma: no cover
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8000)
