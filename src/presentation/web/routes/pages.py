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


@router.get("/monitoring", response_class=HTMLResponse)
async def get_monitoring_page(request: Request):
    """General monitoring dashboard page."""
    return templates.TemplateResponse("monitoring.html", {"request": request})


@router.get("/execution/{work_id}", response_class=HTMLResponse)
async def get_execution_detail_page(request: Request, work_id: str):
    """Detailed execution monitoring page."""
    return templates.TemplateResponse(
        "execution_detail.html", {"request": request, "work_id": work_id}
    )


@router.get("/batch/manage", response_class=HTMLResponse)
async def get_batch_manager_page(request: Request):
    """Batch management page."""
    return templates.TemplateResponse("batch_manager.html", {"request": request})


@router.get("/results", response_class=HTMLResponse)
async def get_results_page(request: Request):
    """Results viewer page."""
    return templates.TemplateResponse("results.html", {"request": request})
