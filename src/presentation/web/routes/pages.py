"""
Web page routes (HTML responses).

This module provides FastAPI routes for serving HTML pages of the CSPBench
web interface. All routes return rendered Jinja2 templates with appropriate
context data for the frontend application.
"""

from fastapi import APIRouter, Request
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates

templates = Jinja2Templates(directory="src/presentation/web/templates")

router = APIRouter(tags=["pages"])


@router.get("/", response_class=HTMLResponse)
async def get_main_page(request: Request):
    """Serve the main page with execution type selection.

    This is the landing page of the CSPBench web interface where users
    can choose between different execution modes and access main features.

    Args:
        request (Request): FastAPI request object for template context.

    Returns:
        HTMLResponse: Rendered index.html template.
    """
    return templates.TemplateResponse("index.html", {"request": request})


@router.get("/dataset/manage", response_class=HTMLResponse)
async def get_dataset_manager_page(request: Request):
    """Serve the dataset management page.

    Provides the interface for managing datasets including viewing,
    editing, deleting, and organizing existing datasets.

    Args:
        request (Request): FastAPI request object for template context.

    Returns:
        HTMLResponse: Rendered dataset_manager.html template.
    """
    return templates.TemplateResponse("dataset_manager.html", {"request": request})


@router.get("/dataset-manager", response_class=HTMLResponse)
async def redirect_old_dataset_manager(request: Request):
    """Redirect old dataset-manager URL to new dataset/manage.

    Maintains backward compatibility by redirecting the old URL structure
    to the new organized route structure.

    Args:
        request (Request): FastAPI request object (unused in redirect).

    Returns:
        RedirectResponse: Permanent redirect to new URL structure.
    """
    from fastapi.responses import RedirectResponse

    return RedirectResponse(url="/dataset/manage", status_code=301)


@router.get("/dataset-generator", response_class=HTMLResponse)
async def get_dataset_generator_page(request: Request):
    """Serve the dataset generator page.

    Provides the interface for creating new datasets through various
    methods including synthetic generation, file upload, and NCBI download.

    Args:
        request (Request): FastAPI request object for template context.

    Returns:
        HTMLResponse: Rendered dataset_generator.html template.
    """
    return templates.TemplateResponse("dataset_generator.html", {"request": request})


@router.get("/batch/manage", response_class=HTMLResponse)
async def get_batch_manager_page(request: Request):
    """Serve the batch manager page.

    Provides the interface for managing batch configuration files
    including creation, editing, and organization of batch workflows.

    Args:
        request (Request): FastAPI request object for template context.

    Returns:
        HTMLResponse: Rendered batch_manager.html template.
    """
    return templates.TemplateResponse("batch_manager.html", {"request": request})


@router.get("/monitor/works", response_class=HTMLResponse)
async def get_works_monitor_page(request: Request):
    """Serve the works monitor list page.

    Provides the interface for viewing and monitoring all batch executions
    with status information, filtering, and management capabilities.

    Args:
        request (Request): FastAPI request object for template context.

    Returns:
        HTMLResponse: Rendered monitor_works.html template.
    """
    return templates.TemplateResponse("monitor_works.html", {"request": request})


@router.get("/monitor/work/{work_id}", response_class=HTMLResponse)
async def get_work_progress_page(request: Request, work_id: str):
    """Serve the work progress monitor page.

    Provides detailed monitoring interface for a specific batch execution
    including real-time progress, logs, and detailed status information.

    Args:
        request (Request): FastAPI request object for template context.
        work_id (str): Unique identifier of the work to monitor.

    Returns:
        HTMLResponse: Rendered monitor_progress.html template with work context.
    """
    return templates.TemplateResponse(
        "monitor_progress.html", {"request": request, "work_id": work_id}
    )


@router.get("/results", response_class=HTMLResponse)
async def get_results_page(request: Request):
    """Serve the results viewer page.

    Provides the interface for viewing, analyzing, and downloading
    results from completed batch executions.

    Args:
        request (Request): FastAPI request object for template context.

    Returns:
        HTMLResponse: Rendered results.html template.
    """
    return templates.TemplateResponse("results.html", {"request": request})
