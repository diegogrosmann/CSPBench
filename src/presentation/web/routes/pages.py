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


@router.get("/execution/{execution_type}", response_class=HTMLResponse)
async def execution_page(request: Request, execution_type: str):
    """Execution pages for different types."""
    # Map short names to full template names
    type_mapping = {
        "batch": "batch_execution",
        "comparison": "algorithm_comparison",
        "benchmark": "benchmark_suite",
        "custom": "custom_workflow",
        "generator": "dataset_generator",
    }

    template_name = type_mapping.get(execution_type, execution_type)
    template_file = f"{template_name}.html"

    try:
        return templates.TemplateResponse(
            template_file, {"request": request, "execution_type": execution_type}
        )
    except Exception:
        # Fallback to batch execution if template not found
        return templates.TemplateResponse(
            "batch_execution.html",
            {"request": request, "execution_type": execution_type},
        )


@router.get("/dataset-generator", response_class=HTMLResponse)
async def dataset_generator_page(request: Request):
    """Dataset generator page."""
    return templates.TemplateResponse("dataset_generator.html", {"request": request})


@router.get("/generator", response_class=HTMLResponse)
async def generator_redirect(request: Request):
    """Alternative route for dataset generator."""
    return templates.TemplateResponse("dataset_generator.html", {"request": request})


@router.get("/results", response_class=HTMLResponse)
async def results_page(request: Request):
    """Results management page."""
    return templates.TemplateResponse("results.html", {"request": request})


@router.get("/results/{session_id}", response_class=HTMLResponse)
async def result_details_page(request: Request, session_id: str):
    """Result details page for a specific session."""
    return templates.TemplateResponse(
        "result_details.html", {"request": request, "session_id": session_id}
    )


@router.get("/debug", response_class=HTMLResponse)
async def debug_execution_page(request: Request):
    """Debug execution page for testing batch execution with detailed logs."""
    return templates.TemplateResponse("debug_execution.html", {"request": request})
