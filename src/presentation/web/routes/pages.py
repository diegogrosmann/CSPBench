"""
Web page routes (HTML responses).
"""

from fastapi import APIRouter, Request
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates

templates = Jinja2Templates(directory="src/presentation/web/templates")

router = APIRouter(tags=["pages"])


@router.get("/", response_class=HTMLResponse)
async def get_main_page(request: Request):
    """Main page with execution type selection."""
    return templates.TemplateResponse("index.html", {"request": request})


@router.get("/dataset/manage", response_class=HTMLResponse)
async def get_dataset_manager_page(request: Request):
    """Dataset management page."""
    return templates.TemplateResponse("dataset_manager.html", {"request": request})


@router.get("/dataset-manager", response_class=HTMLResponse)
async def redirect_old_dataset_manager(request: Request):
    """Redirect old dataset-manager URL to new dataset/manage."""
    from fastapi.responses import RedirectResponse

    return RedirectResponse(url="/dataset/manage", status_code=301)


@router.get("/dataset-generator", response_class=HTMLResponse)
async def get_dataset_generator_page(request: Request):
    """Dataset generator page."""
    return templates.TemplateResponse("dataset_generator.html", {"request": request})


@router.get("/batch/manage", response_class=HTMLResponse)
async def get_batch_manager_page(request: Request):
    """Batch manager page."""
    return templates.TemplateResponse("batch_manager.html", {"request": request})


@router.get("/monitor/works", response_class=HTMLResponse)
async def get_works_monitor_page(request: Request):
    """Works monitor list page."""
    return templates.TemplateResponse("monitor_works.html", {"request": request})


@router.get("/monitor/work/{work_id}", response_class=HTMLResponse)
async def get_work_progress_page(request: Request, work_id: str):
    """Work progress monitor page."""
    return templates.TemplateResponse("monitor_progress.html", {
        "request": request, 
        "work_id": work_id
    })


@router.get("/results", response_class=HTMLResponse)
async def get_results_page(request: Request):
    """Results viewer page."""
    return templates.TemplateResponse("results.html", {"request": request})
