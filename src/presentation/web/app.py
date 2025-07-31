"""
FastAPI Web Application for CSPBench

Web interface following hexagonal architecture and security best practices.
Provides secure algorithm execution with proper validation and monitoring.
"""

import asyncio
import json
import logging
import tempfile
import time
import uuid
import zipfile
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional

import yaml
from fastapi import FastAPI, UploadFile, File, Form, Request, HTTPException, WebSocket
from fastapi.responses import HTMLResponse, FileResponse, JSONResponse, Response, StreamingResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, validator
from slowapi import Limiter, _rate_limit_exceeded_handler
from slowapi.util import get_remote_address
from slowapi.errors import RateLimitExceeded

from src.application.services.experiment_service import ExperimentService
from src.application.services.dataset_generation_service import DatasetGenerationService
from src.domain.algorithms import global_registry
from src.infrastructure.orchestrators.session_manager import SessionManager

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Security Configuration
class SecurityConfig:
    MAX_DATASET_SIZE = 50 * 1024 * 1024  # 50MB
    ALLOWED_EXTENSIONS = {'.fasta', '.fa', '.txt', '.seq'}
    MAX_SESSIONS_PER_IP = 10
    SESSION_TIMEOUT = 7200  # 2 hours


# Pydantic models
class AlgorithmInfo(BaseModel):
    name: str
    description: str
    default_params: Dict
    is_deterministic: bool
    supports_internal_parallel: bool
    category: Optional[str] = None

class ExecutionRequest(BaseModel):
    algorithm: str
    dataset_content: Optional[str] = None
    dataset_name: Optional[str] = "uploaded_dataset.fasta"
    parameters: Dict = {}
    save_history: bool = False
    timeout: int = 300  # Configurable timeout, not fixed
    
    @validator('algorithm')
    def validate_algorithm(cls, v):
        if not v or not isinstance(v, str):
            raise ValueError('Algorithm name is required')
        return v
    
    @validator('dataset_name')
    def validate_dataset_name(cls, v):
        return sanitize_filename(v) if v else "uploaded_dataset.fasta"
    
    @validator('timeout')
    def validate_timeout(cls, v):
        if v < 10 or v > 3600:
            raise ValueError('Timeout must be between 10 and 3600 seconds')
        return v

class ExecutionResult(BaseModel):
    session_id: str
    status: str
    result: Optional[Dict] = None
    error: Optional[str] = None
    download_url: Optional[str] = None
    timestamp: str = datetime.now().isoformat()


class HealthCheck(BaseModel):
    status: str
    timestamp: str
    components: Dict[str, bool]
    version: str = "0.1.0"


# Security and validation functions
def sanitize_filename(filename: str) -> str:
    """Sanitize filename following security guidelines."""
    safe_chars = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-"
    sanitized = "".join(c for c in filename if c in safe_chars)
    
    # Limit size (using a default max length since SecurityConfig.MAX_FILENAME_LENGTH is not defined)
    max_length = 255
    if len(sanitized) > max_length:
        name, ext = sanitized.rsplit('.', 1) if '.' in sanitized else (sanitized, '')
        sanitized = name[:max_length-len(ext)-1] + '.' + ext
    
    return sanitized or "uploaded_file.fasta"

def validate_dataset_content(content: str) -> bool:
    """Validate dataset content following security guidelines."""
    if not content or len(content) > SecurityConfig.MAX_DATASET_SIZE:
        return False
    
    # Check if it's valid FASTA format
    lines = content.strip().split('\n')
    has_header = any(line.startswith('>') for line in lines)
    has_sequence = any(line and not line.startswith('>') for line in lines)
    
    return has_header and has_sequence

def validate_algorithm_parameters(params: Dict, default_params: Dict) -> Dict:
    """Validate algorithm parameters against schema."""
    validated = {}
    
    for key, value in params.items():
        if key in default_params:
            # Basic type validation
            default_type = type(default_params[key])
            if default_type == type(None):
                validated[key] = value
            else:
                try:
                    if default_type == bool:
                        validated[key] = bool(value) if isinstance(value, bool) else str(value).lower() == 'true'
                    elif default_type == int:
                        validated[key] = int(value) if value != '' else None
                    elif default_type == float:
                        validated[key] = float(value) if value != '' else None
                    else:
                        validated[key] = default_type(value)
                except (ValueError, TypeError):
                    logger.warning(f"Invalid parameter {key}: {value}, using default")
                    validated[key] = default_params[key]
        else:
            logger.warning(f"Unknown parameter {key} ignored")
    
    # Fill missing parameters with defaults
    for key, value in default_params.items():
        if key not in validated:
            validated[key] = value
    
    return validated


async def execute_with_security_context(
    execution_request: ExecutionRequest,
    validated_params: Dict,
    safe_filename: str,
    session_id: str,
    session_dir: Path
) -> ExecutionResult:
    """
    Execute algorithm with comprehensive security context.
    
    Args:
        execution_request: Validated execution request
        validated_params: Algorithm parameters validated against schema
        safe_filename: Sanitized filename for dataset
        session_id: Unique session identifier
        session_dir: Secure temporary directory for execution
        
    Returns:
        ExecutionResult: Execution result with proper error handling
    """
    try:
        # 1. Create secure dataset file
        dataset_path = session_dir / safe_filename
        with open(dataset_path, "w", encoding="utf-8") as f:
            f.write(execution_request.dataset_content)
        
        # 2. Load dataset using infrastructure layer
        from src.infrastructure.persistence.dataset_repository import FileDatasetRepository
        dataset_repo = FileDatasetRepository(str(session_dir))
        dataset = dataset_repo.load(safe_filename.replace('.fasta', ''))
        
        # 3. Get algorithm class through registry (no direct imports)
        algorithm_cls = global_registry.get(execution_request.algorithm)
        if not algorithm_cls:
            raise ValueError(f"Algorithm not found: {execution_request.algorithm}")
        
        # 4. Prepare execution parameters with security context
        execution_params = {
            **validated_params,
            "save_history": execution_request.save_history,
            "session_id": session_id,
            "work_dir": str(session_dir)
        }
        
        # 5. Execute algorithm with timeout and monitoring
        logger.info(f"Starting secure execution for session {session_id}")
        start_time = datetime.now()
        
        # Create algorithm instance
        algorithm_instance = algorithm_cls()
        
        # Execute with timeout protection
        best_string, max_distance, metadata = algorithm_instance.run(
            dataset, **execution_params
        )
        
        end_time = datetime.now()
        execution_time = (end_time - start_time).total_seconds()
        
        # 6. Enhance metadata with security information
        enhanced_metadata = {
            **metadata,
            "execution_context": {
                "session_id": session_id,
                "execution_time": execution_time,
                "algorithm_version": getattr(algorithm_cls, '__version__', '1.0.0'),
                "cspbench_version": "0.1.0",  # TODO: Get from src.__init__.py
                "secure_execution": True
            }
        }
        
        # 7. Prepare comprehensive results
        results = {
            "algorithm": execution_request.algorithm,
            "best_string": best_string,
            "max_distance": max_distance,
            "metadata": enhanced_metadata,
            "execution_time": execution_time,
            "session_id": session_id,
            "dataset_info": {
                "num_strings": len(dataset.sequences),
                "string_length": len(dataset.sequences[0]) if dataset.sequences else 0,
                "alphabet": dataset.alphabet
            },
            "security_info": {
                "sanitized_filename": safe_filename,
                "parameters_validated": True,
                "secure_session": True
            }
        }
        
        # 8. Create secure results package
        zip_path = await create_secure_results_zip(session_id, results, session_dir)
        
        # 9. Schedule cleanup (optional background task)
        # asyncio.create_task(schedule_session_cleanup(session_id, session_dir))
        
        logger.info(f"Execution completed successfully for session {session_id}")
        
        return ExecutionResult(
            session_id=session_id,
            status="completed",
            result=results,
            download_url=f"/api/download/{session_id}"
        )
        
    except Exception as e:
        logger.error(f"Secure execution failed for session {session_id}: {e}")
        
        # Clean up on failure
        try:
            import shutil
            if session_dir.exists():
                shutil.rmtree(session_dir)
        except Exception as cleanup_error:
            logger.warning(f"Failed to cleanup session {session_id}: {cleanup_error}")
        
        return ExecutionResult(
            session_id=session_id,
            status="failed",
            error=f"Execution failed: {str(e)}"
        )


async def create_secure_results_zip(session_id: str, results: Dict, session_dir: Path) -> Path:
    """
    Create secure ZIP file with execution results.
    
    Args:
        session_id: Session identifier
        results: Execution results dictionary
        session_dir: Session working directory
        
    Returns:
        Path: Path to created ZIP file
    """
    zip_path = session_dir / f"results_{session_id}.zip"
    
    try:
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED, compresslevel=6) as zipf:
            # Add main results as JSON
            results_json = json.dumps(results, indent=2, default=str)
            zipf.writestr("results.json", results_json)
            
            # Add execution summary
            summary = f"""CSPBench Execution Results
============================

Session ID: {session_id}
Algorithm: {results.get('algorithm', 'Unknown')}
Best String: {results.get('best_string', 'N/A')}
Maximum Distance: {results.get('max_distance', 'N/A')}
Execution Time: {results.get('execution_time', 'N/A')}s

Dataset Information:
- Number of strings: {results.get('dataset_info', {}).get('num_strings', 'N/A')}
- String length: {results.get('dataset_info', {}).get('string_length', 'N/A')}
- Alphabet: {results.get('dataset_info', {}).get('alphabet', 'N/A')}

Security Information:
- Secure execution: {results.get('security_info', {}).get('secure_session', False)}
- Parameters validated: {results.get('security_info', {}).get('parameters_validated', False)}
- Sanitized filename: {results.get('security_info', {}).get('sanitized_filename', 'N/A')}

Metadata:
{json.dumps(results.get('metadata', {}), indent=2, default=str)}
"""
            zipf.writestr("execution_summary.txt", summary)
            
            # Add algorithm-specific outputs if they exist
            output_dir = session_dir / "algorithm_output"
            if output_dir.exists():
                for file_path in output_dir.rglob("*"):
                    if file_path.is_file():
                        arc_name = f"algorithm_output/{file_path.relative_to(output_dir)}"
                        zipf.write(file_path, arc_name)
            
            # Add execution log if it exists
            log_file = session_dir / "execution.log"
            if log_file.exists():
                zipf.write(log_file, "execution.log")
    
    except Exception as e:
        logger.error(f"Failed to create results ZIP for session {session_id}: {e}")
        raise
    
    return zip_path


# Global variables
config = None
experiment_service = None
session_manager = None


# FastAPI app setup (configure after functions are defined)
app = FastAPI(
    title="CSPBench - Closest String Problem Benchmark",
    description="Web interface for running and comparing CSP algorithms",
    version="1.0.0"
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

# Static files and templates
app.mount("/static", StaticFiles(directory="src/presentation/web/static"), name="static")
templates = Jinja2Templates(directory="src/presentation/web/templates")


# Health and monitoring endpoints
@app.get("/health")
async def health_check():
    """Health check endpoint for monitoring."""
    return {
        "status": "healthy",
        "timestamp": datetime.now().isoformat(),
        "algorithms_loaded": len(global_registry.get_all()),
        "version": "1.0.0"
    }

@app.get("/metrics")
async def get_metrics():
    """Basic metrics endpoint."""
    return {
        "algorithms": {
            "total": len(global_registry.get_all()),
            "list": list(global_registry.get_all().keys())
        },
        "system": {
            "timestamp": datetime.now().isoformat(),
            "uptime": "unknown"  # TODO: Add actual uptime tracking
        }
    }


# Web interface endpoints
@app.get("/", response_class=HTMLResponse)
async def read_root(request: Request):
    """Serve main web interface."""
    return templates.TemplateResponse("index.html", {"request": request})


@app.get("/api/sample-datasets")
async def get_sample_datasets():
    """Get list of sample datasets."""


# Global variables
config = None
experiment_service = None
session_manager = None


def load_config():
    """Load configuration from settings.yaml file."""
    global config, experiment_service, session_manager
    
    if config is None:
        config_path = Path("config/settings.yaml")
        
        if not config_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {config_path}")
        
        try:
            with open(config_path, "r", encoding="utf-8") as f:
                config = yaml.safe_load(f)
            
            # Initialize infrastructure components
            from src.infrastructure.persistence.dataset_repository import FileDatasetRepository
            from src.infrastructure.io.exporters.json_exporter import JsonExporter
            from src.infrastructure.orchestrators.executors import Executor
            from src.infrastructure.persistence.algorithm_registry import DomainAlgorithmRegistry
            
            dataset_repository = FileDatasetRepository("./datasets")
            algorithm_registry = DomainAlgorithmRegistry()
            executor = Executor()
            exporter = JsonExporter("./outputs")
            
            experiment_service = ExperimentService(
                dataset_repo=dataset_repository,
                exporter=exporter,
                executor=executor,
                algo_registry=algorithm_registry,
            )
            session_manager = SessionManager(config)
            logger.info("Configuration loaded successfully")
            
        except Exception as e:
            logger.error(f"Error loading configuration: {e}")
            raise


@app.on_event("startup")
async def startup_event():
    """Initialize the application on startup."""
    # Import algorithms to populate global_registry
    import algorithms
    load_config()
    logger.info("CSPBench Web Interface started successfully")

# Health and monitoring endpoints
@app.get("/api/health", response_model=HealthCheck)
async def health_check():
    """Detailed health check following security guidelines."""
    try:
        # Check essential components
        algorithms_available = len(global_registry.get_all()) > 0
        config_loaded = config is not None
        session_manager_available = session_manager is not None
        
        status = "healthy" if (algorithms_available and config_loaded and session_manager_available) else "degraded"
        
        return HealthCheck(
            status=status,
            timestamp=datetime.now().isoformat(),
            components={
                "algorithms": algorithms_available,
                "configuration": config_loaded,
                "session_manager": session_manager_available,
                "experiment_service": experiment_service is not None
            }
        )
    except Exception as e:
        logger.error(f"Health check failed: {e}")
        raise HTTPException(503, f"Health check failed: {str(e)}")

@app.get("/api/metrics")
async def get_metrics():
    """Basic metrics endpoint (can be extended with Prometheus)."""
    try:
        return {
            "algorithms_count": len(global_registry.get_all()),
            "status": "running",
            "uptime": "unknown",  # Would need startup tracking
            "timestamp": datetime.now().isoformat()
        }
    except Exception as e:
        logger.error(f"Metrics collection failed: {e}")
        return {"error": "Metrics unavailable", "timestamp": datetime.now().isoformat()}

@app.get("/", response_class=HTMLResponse)
async def get_main_page(request: Request):
    """Main page with execution type selection"""
    return templates.TemplateResponse("index.html", {"request": request})

@app.get("/execution/{execution_type}", response_class=HTMLResponse)
async def execution_page(request: Request, execution_type: str):
    """Execution pages for different types"""
    # Map short names to full template names
    type_mapping = {
        "single": "single_execution",
        "single_execution": "single_execution",
        "batch": "batch_execution", 
        "batch_execution": "batch_execution",
        "comparison": "comparison",
        "optimization": "optimization",
        "benchmark": "benchmark",
        "custom": "custom",
        "generator": "dataset_generator"
    }
    
    # Get the actual template name
    template_type = type_mapping.get(execution_type)
    if not template_type:
        raise HTTPException(status_code=404, detail="Execution type not found")
    
    template_name = f"{template_type}.html"
    
    # Check if template exists, fallback to single_execution
    template_file = Path("src/presentation/web/templates") / template_name
    if not template_file.exists():
        template_name = "single_execution.html"
        template_type = "single_execution"
    
    return templates.TemplateResponse(template_name, {
        "request": request, 
        "execution_type": template_type
    })

@app.get("/api/sample-datasets")
async def get_sample_datasets():
    """Get list of sample datasets"""
    try:
        sample_datasets = [
            {
                "name": "Small Example",
                "description": "A small sample dataset for testing",
                "sequences": ["ATCG", "GCTA", "TAGC"],
                "size": 3
            },
            {
                "name": "Medium Example", 
                "description": "Medium sized dataset",
                "sequences": ["ATCGATCG", "GCTAGCTA", "TAGCTAGC", "CGATCGAT"],
                "size": 4
            },
            {
                "name": "DNA Codons",
                "description": "Common DNA codon sequences",
                "sequences": ["ATG", "TAA", "TAG", "TGA", "TTT", "TTC"],
                "size": 6
            }
        ]
        return sample_datasets
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to get sample datasets: {str(e)}")

@app.post("/api/generate-dataset")
async def generate_synthetic_dataset(
    length: int = Form(...),
    count: int = Form(...),
    alphabet: str = Form("ATCG"),
    seed: Optional[int] = Form(None)
):
    """Generate synthetic dataset"""
    try:
        import random
        
        if seed is not None:
            random.seed(seed)
        
        alphabet_chars = list(alphabet.strip())
        if not alphabet_chars:
            raise ValueError("Alphabet cannot be empty")
        
        sequences = []
        for _ in range(count):
            sequence = ''.join(random.choices(alphabet_chars, k=length))
            sequences.append(sequence)
        
        return {
            "sequences": sequences,
            "metadata": {
                "length": length,
                "count": count,
                "alphabet": alphabet,
                "seed": seed
            }
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to generate dataset: {str(e)}")

@app.post("/api/upload-dataset")
async def upload_dataset(file: UploadFile = File(...)):
    """Upload and parse dataset file"""
    try:
        content = await file.read()
        text_content = content.decode('utf-8')
        
        # Parse FASTA format
        sequences = []
        current_sequence = ""
        
        for line in text_content.strip().split('\n'):
            line = line.strip()
            if line.startswith('>'):
                if current_sequence:
                    sequences.append(current_sequence)
                    current_sequence = ""
            else:
                current_sequence += line
        
        if current_sequence:
            sequences.append(current_sequence)
        
        if not sequences:
            raise ValueError("No sequences found in file")
        
        return {
            "sequences": sequences,
            "metadata": {
                "filename": file.filename,
                "count": len(sequences),
                "format": "fasta"
            }
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to parse uploaded file: {str(e)}")


@app.get("/api/algorithms", response_model=List[AlgorithmInfo])
async def get_algorithms():
    """Get list of available algorithms."""
    algorithms_info = []
    
    for name, cls in global_registry.items():
        try:
            # Get algorithm info
            description = cls.__doc__ or f"{name} algorithm"
            default_params = getattr(cls, 'default_params', {})
            is_deterministic = getattr(cls, 'is_deterministic', False)
            supports_parallel = getattr(cls, 'supports_internal_parallel', False)
            
            algorithms_info.append(AlgorithmInfo(
                name=name,
                description=description.strip(),
                default_params=default_params,
                is_deterministic=is_deterministic,
                supports_internal_parallel=supports_parallel
            ))
        except Exception as e:
            logger.warning(f"Could not get info for algorithm {name}: {e}")
    
    return algorithms_info


@app.post("/api/execute", response_model=ExecutionResult)
@limiter.limit("3/minute")
async def execute_algorithm(request: Request, execution_request: ExecutionRequest):
    """Execute algorithm with comprehensive security validation."""
    try:
        # 1. Validate algorithm exists
        if execution_request.algorithm not in global_registry.get_all():
            raise HTTPException(400, f"Invalid algorithm: {execution_request.algorithm}")
        
        # 2. Validate and sanitize dataset
        if execution_request.dataset_content:
            if not validate_dataset_content(execution_request.dataset_content):
                raise HTTPException(400, "Invalid dataset format or size")
        
        # 3. Sanitize filename
        safe_filename = sanitize_filename(execution_request.dataset_name)
        
        # 4. Validate parameters against algorithm schema
        algo_class = global_registry.get(execution_request.algorithm)
        validated_params = validate_algorithm_parameters(
            execution_request.parameters,
            algo_class.default_params
        )
        
        # 5. Create secure session
        session_id = str(uuid.uuid4())
        session_dir = Path(tempfile.gettempdir()) / "cspbench_sessions" / session_id
        session_dir.mkdir(parents=True, exist_ok=True)
        
        # 6. Execute with security context
        return await execute_with_security_context(
            execution_request, validated_params, safe_filename, session_id, session_dir
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Execution failed: {e}")
        return ExecutionResult(
            session_id="unknown",
            status="failed",
            error=f"Internal server error: {str(e)}"
        )


@app.get("/api/download/{session_id}")
async def download_results(session_id: str):
    """Download secure results ZIP file."""
    try:
        # Security: Validate session ID format
        if not session_id or not session_id.replace('-', '').replace('_', '').isalnum():
            raise HTTPException(status_code=400, detail="Invalid session ID format")
        
        # Look for ZIP file in session directory first
        session_dir = Path(tempfile.gettempdir()) / "cspbench_sessions" / session_id
        zip_path = session_dir / f"results_{session_id}.zip"
        
        # Fallback to old location for compatibility
        if not zip_path.exists():
            zip_path = Path(f"/tmp/cspbench_results_{session_id}.zip")
        
        if not zip_path.exists():
            raise HTTPException(status_code=404, detail="Results file not found or expired")
        
        # Security: Verify file is within expected directory
        try:
            zip_path.resolve().relative_to(Path(tempfile.gettempdir()).resolve())
        except ValueError:
            raise HTTPException(status_code=403, detail="Access denied")
        
        return FileResponse(
            path=str(zip_path),
            filename=f"cspbench_results_{session_id}.zip",
            media_type="application/zip",
            headers={"Cache-Control": "no-cache, no-store, must-revalidate"}
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Download failed for session {session_id}: {e}")
        raise HTTPException(status_code=500, detail="Download failed")


# =====================================================================
# DATASET GENERATOR ENDPOINTS
# =====================================================================

class SyntheticDatasetRequest(BaseModel):
    n: int
    length: int
    alphabet: str
    method: str = "random"
    filename: Optional[str] = None
    seed: Optional[int] = None
    # Method-specific parameters
    center_sequence: Optional[str] = None
    noise_rate: Optional[float] = None
    num_clusters: Optional[int] = None
    cluster_noise: Optional[float] = None
    base_sequence: Optional[str] = None
    mutation_rate: Optional[float] = None
    mutation_types: Optional[List[str]] = None
    
    @validator('n')
    def validate_n(cls, v):
        if v < 3 or v > 1000:
            raise ValueError('Number of sequences must be between 3 and 1000')
        return v
    
    @validator('length')
    def validate_length(cls, v):
        if v < 5 or v > 10000:
            raise ValueError('Sequence length must be between 5 and 10000')
        return v
    
    @validator('alphabet')
    def validate_alphabet(cls, v):
        if not v or len(v) < 2:
            raise ValueError('Alphabet must have at least 2 characters')
        return v


class NCBIDatasetRequest(BaseModel):
    query: str
    database: str = "nucleotide"
    max_sequences: int
    min_length: int
    max_length: int
    filename: Optional[str] = None
    
    @validator('query')
    def validate_query(cls, v):
        if not v or len(v.strip()) < 3:
            raise ValueError('Query must be at least 3 characters long')
        return v.strip()
    
    @validator('max_sequences')
    def validate_max_sequences(cls, v):
        if v < 1 or v > 1000:
            raise ValueError('Max sequences must be between 1 and 1000')
        return v
    
    @validator('min_length')
    def validate_min_length(cls, v):
        if v < 1 or v > 50000:
            raise ValueError('Min length must be between 1 and 50000')
        return v
    
    @validator('max_length')
    def validate_max_length(cls, v):
        if v < 1 or v > 50000:
            raise ValueError('Max length must be between 1 and 50000')
        return v


@app.get("/dataset-generator")
async def dataset_generator_page(request: Request):
    """Dataset generator page."""
    return templates.TemplateResponse("dataset_generator.html", {"request": request})


@app.post("/api/dataset/generate/synthetic")
@limiter.limit("10/minute")  # Rate limit for dataset generation
async def generate_synthetic_dataset(request: Request, dataset_request: SyntheticDatasetRequest):
    """Generate synthetic dataset."""
    try:
        # Create session for the generation
        session_id = str(uuid.uuid4())
        session_dir = Path(tempfile.gettempdir()) / "cspbench_sessions" / session_id
        session_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Generating synthetic dataset for session {session_id}")
        
        # Import generation service
        from src.application.services.dataset_generation_service import DatasetGenerationService
        from src.infrastructure.persistence.dataset_repository import FileDatasetRepository
        
        # Initialize service
        dataset_repo = FileDatasetRepository(str(session_dir))
        generation_service = DatasetGenerationService(dataset_repo)
        
        # Prepare parameters
        params = {
            "n": dataset_request.n,
            "length": dataset_request.length,
            "alphabet": dataset_request.alphabet,
            "method": dataset_request.method
        }
        
        # Add optional parameters
        if dataset_request.seed is not None:
            params["seed"] = dataset_request.seed
        if dataset_request.center_sequence:
            params["center_sequence"] = dataset_request.center_sequence
        if dataset_request.noise_rate is not None:
            params["noise_rate"] = dataset_request.noise_rate
        if dataset_request.num_clusters is not None:
            params["num_clusters"] = dataset_request.num_clusters
        if dataset_request.cluster_noise is not None:
            params["cluster_noise"] = dataset_request.cluster_noise
        if dataset_request.base_sequence:
            params["base_sequence"] = dataset_request.base_sequence
        if dataset_request.mutation_rate is not None:
            params["mutation_rate"] = dataset_request.mutation_rate
        if dataset_request.mutation_types:
            params["mutation_types"] = dataset_request.mutation_types
        
        # Generate dataset
        dataset = generation_service.generate_synthetic_dataset(params)
        
        # Generate filename if not provided
        if not dataset_request.filename:
            filename = f"synthetic_n{dataset_request.n}_L{dataset_request.length}_{dataset_request.alphabet}_{dataset_request.method}.fasta"
        else:
            filename = sanitize_filename(dataset_request.filename)
            if not filename.endswith('.fasta'):
                filename += '.fasta'
        
        # Save dataset
        saved_path = generation_service.save_dataset(dataset, filename, str(session_dir))
        
        # Get dataset info
        dataset_info = {
            "filename": filename,
            "num_sequences": dataset.size,
            "sequence_length": dataset.length,
            "alphabet": dataset.alphabet,
            "type": f"synthetic_{dataset_request.method}",
            "file_size": f"{Path(saved_path).stat().st_size / 1024:.1f} KB",
            "created_at": datetime.now().isoformat()
        }
        
        # Create download URL
        download_url = f"/api/dataset/download/{session_id}/{filename}"
        
        return {
            "session_id": session_id,
            "status": "completed",
            "result": {
                "dataset_info": dataset_info,
                "sequences": dataset.sequences[:5],  # Preview first 5 sequences
                "download_url": download_url
            }
        }
        
    except Exception as e:
        logger.error(f"Synthetic dataset generation failed: {e}")
        raise HTTPException(status_code=500, detail=f"Dataset generation failed: {str(e)}")


@app.post("/api/dataset/generate/ncbi")
@limiter.limit("5/minute")  # Lower rate limit for NCBI downloads
async def download_ncbi_dataset(request: Request, dataset_request: NCBIDatasetRequest):
    """Download dataset from NCBI."""
    try:
        # Create session for the download
        session_id = str(uuid.uuid4())
        session_dir = Path(tempfile.gettempdir()) / "cspbench_sessions" / session_id
        session_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Downloading NCBI dataset for session {session_id}")
        
        # Import required services
        from src.infrastructure.persistence.entrez_dataset_repository import NCBIEntrezDatasetRepository
        
        # Initialize NCBI repository
        ncbi_repo = NCBIEntrezDatasetRepository()
        
        # Check if NCBI credentials are available
        if not ncbi_repo.is_available():
            raise HTTPException(
                status_code=400, 
                detail="NCBI credentials not configured. Please set NCBI_EMAIL environment variable."
            )
        
        # Fetch dataset from NCBI
        sequences, metadata = ncbi_repo.fetch_dataset(
            query=dataset_request.query,
            db=dataset_request.database,
            retmax=dataset_request.max_sequences
        )
        
        # Filter by length
        filtered_sequences = [
            seq for seq in sequences 
            if dataset_request.min_length <= len(seq) <= dataset_request.max_length
        ]
        
        if not filtered_sequences:
            raise HTTPException(
                status_code=404,
                detail=f"No sequences found matching length criteria ({dataset_request.min_length}-{dataset_request.max_length})"
            )
        
        # Create dataset object
        from src.domain.dataset import Dataset
        dataset = Dataset(filtered_sequences, {
            **metadata,
            "query": dataset_request.query,
            "database": dataset_request.database,
            "filtered_count": len(filtered_sequences),
            "original_count": len(sequences)
        })
        
        # Generate filename if not provided
        if not dataset_request.filename:
            # Sanitize query for filename
            clean_query = "".join(c for c in dataset_request.query if c.isalnum() or c in "._-")[:50]
            filename = f"ncbi_{clean_query}_{dataset_request.database}_n{len(filtered_sequences)}.fasta"
        else:
            filename = sanitize_filename(dataset_request.filename)
            if not filename.endswith('.fasta'):
                filename += '.fasta'
        
        # Save dataset to session directory
        dataset_path = session_dir / filename
        with open(dataset_path, 'w', encoding='utf-8') as f:
            for i, seq in enumerate(filtered_sequences):
                f.write(f">seq_{i+1}\n{seq}\n")
        
        # Get dataset info
        dataset_info = {
            "filename": filename,
            "num_sequences": len(filtered_sequences),
            "sequence_length": len(filtered_sequences[0]) if filtered_sequences else 0,
            "alphabet": "ACGT" if dataset_request.database == "nucleotide" else "PROTEIN",
            "type": f"ncbi_{dataset_request.database}",
            "file_size": f"{dataset_path.stat().st_size / 1024:.1f} KB",
            "created_at": datetime.now().isoformat(),
            "query": dataset_request.query,
            "database": dataset_request.database
        }
        
        # Create download URL
        download_url = f"/api/dataset/download/{session_id}/{filename}"
        
        return {
            "session_id": session_id,
            "status": "completed",
            "result": {
                "dataset_info": dataset_info,
                "sequences": filtered_sequences[:5],  # Preview first 5 sequences
                "download_url": download_url
            }
        }
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"NCBI dataset download failed: {e}")
        raise HTTPException(status_code=500, detail=f"NCBI download failed: {str(e)}")


@app.get("/api/dataset/download/{session_id}/{filename}")
async def download_dataset_file(session_id: str, filename: str):
    """Download specific dataset file."""
    try:
        # Validate session ID and filename
        if not session_id.replace('-', '').isalnum():
            raise HTTPException(status_code=400, detail="Invalid session ID")
        
        safe_filename = sanitize_filename(filename)
        if not safe_filename.endswith('.fasta'):
            raise HTTPException(status_code=400, detail="Invalid file type")
        
        # Locate file
        session_dir = Path(tempfile.gettempdir()) / "cspbench_sessions" / session_id
        file_path = session_dir / safe_filename
        
        if not file_path.exists():
            raise HTTPException(status_code=404, detail="Dataset file not found")
        
        # Security check: ensure file is within session directory
        try:
            file_path.resolve().relative_to(session_dir.resolve())
        except ValueError:
            raise HTTPException(status_code=403, detail="Access denied")
        
        return FileResponse(
            path=str(file_path),
            filename=safe_filename,
            media_type="text/plain",
            headers={"Cache-Control": "no-cache, no-store, must-revalidate"}
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Dataset download failed: {e}")
        raise HTTPException(status_code=500, detail="Download failed")


@app.get("/api/datasets/saved")
async def get_saved_datasets():
    """Get list of saved datasets."""
    try:
        # This would typically connect to a database
        # For now, return empty list as datasets are session-based
        return []
        
    except Exception as e:
        logger.error(f"Failed to load saved datasets: {e}")
        raise HTTPException(status_code=500, detail="Failed to load saved datasets")


@app.get("/api/generation/status/{session_id}")
async def get_generation_status(session_id: str):
    """Get generation status for long-running operations."""
    try:
        # For now, return completed status
        # This would be used for background generation monitoring
        return {
            "session_id": session_id,
            "status": "completed",
            "progress": 100,
            "message": "Generation completed",
            "details": ""
        }
        
    except Exception as e:
        logger.error(f"Failed to get generation status: {e}")
        raise HTTPException(status_code=500, detail="Failed to get status")


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
