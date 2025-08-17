from dataclasses import dataclass, field
from pathlib import Path
import re
from typing import Any, Dict, List, Optional, Literal, Union, Sequence

import yaml

from src.domain.errors import BatchConfigurationError
import os

# ====================================================
# NOTA IMPORTANTE
# ====================================================
# O campo root "tasks: TaskGroup" referencia TaskGroup do asyncio importado acima.
# Provavelmente a intenção era usar a classe dataclass TasksGroup definida mais abaixo
# (ou um Union dos grupos específicos). Não alterei a assinatura para evitar
# reescrever o modelo conforme solicitado; durante o carregamento faço um cast.
# Se quiser ajustar depois, basta trocar o tipo para Union[ExperimentTasks, OptimizationTasks, SensitivityTasks].


# ====================================================
# SECTION 1: METADATA
# ====================================================
@dataclass
class MetadataConfig:
    name: str
    description: str
    author: str
    version: str
    creation_date: str
    tags: List[str] = field(default_factory=list)


# ====================================================
# SECTION 2: DATASETS — herança (type fixo na subclasse)
# ====================================================
@dataclass
class DatasetConfig:
    """Base de configuração para datasets (evita colisão com domain.dataset.Dataset)."""

    id: str
    name: str


@dataclass
class SyntheticDatasetConfig(DatasetConfig):
    type: Literal["synthetic"] = "synthetic"
    # Modo de geração
    mode: Literal["random", "noise", "mutations", "clustered"] = "random"
    # Parâmetros comuns
    n: int = 10
    L: int = 20
    alphabet: str = "ACGT"
    seed: Optional[int] = None
    parameters_mode: Dict[str, Any] = field(default_factory=dict)


@dataclass
class FileDatasetConfig(DatasetConfig):
    filename: str
    type: Literal["file"] = "file"


@dataclass
class EntrezDatasetConfig(DatasetConfig):
    type: Literal["entrez"] = "entrez"
    query: str = ""
    db: str = "nucleotide"
    retmax: int = 10
    # Filtros opcionais: None = não filtra
    min_length: Optional[int] = None
    max_length: Optional[int] = None
    # Política de uniformização: None = não uniformiza (permite comprimentos variados)
    uniform_policy: Optional[Literal["strict", "majority", "pad", "trim"]] = None


DatasetAny = Union[SyntheticDatasetConfig, FileDatasetConfig, EntrezDatasetConfig]


# ====================================================
# SECTION 3: ALGORITHMS (SIMPLIFICADO)
# ====================================================
@dataclass
class AlgParamsConfig:
    name: str
    params: Dict[str, Any] = field(default_factory=dict)


@dataclass
class AlgorithmsPresetConfig:
    id: str
    name: str
    description: str = ""
    items: List[AlgParamsConfig] = field(default_factory=list)


# ====================================================
# SECTION 5: TASKS (objetos) + GRUPOS
# ====================================================
# ---- Base Task (usa objetos) ----
@dataclass
class TaskConfig:
    id: str
    name: str
    type: Literal["experiment", "optimization", "sensitivity"]
    datasets: List[DatasetAny] = field(default_factory=list)
    algorithms: List[AlgorithmsPresetConfig] = field(default_factory=list)


# ---- Subtipos de Task ----
@dataclass
class ExperimentTaskConfig(TaskConfig):
    type: Literal["experiment"] = "experiment"
    repetitions: int = 1


@dataclass
class OptimizationTaskConfig(TaskConfig):
    type: Literal["optimization"] = "optimization"
    study_name: str = ""
    parameters: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    config: Optional[Dict[str, Any]] = None


@dataclass
class SensitivityTaskConfig(TaskConfig):
    type: Literal["sensitivity"] = "sensitivity"
    method: Literal["morris", "sobol", "fast", "delta"] = "morris"
    parameters: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    config: Optional[Dict[str, Any]] = None


# ---- Campos de grupo (para todos os grupos) ----
@dataclass
class TasksGroupConfig:
    """Generic tasks group (single type per batch)."""

    type: Literal["experiment", "optimization", "sensitivity"]
    items: Sequence[TaskConfig] = field(default_factory=list)


# EXPERIMENT — inclui method para padronizar a interface
@dataclass
class ExperimentTasksConfig(TasksGroupConfig):
    pass


@dataclass
class OptimizationTasksConfig(TasksGroupConfig):
    framework: Literal["optuna"] = "optuna"


@dataclass
class SensitivityTasksConfig(TasksGroupConfig):
    framework: Literal["SALib"] = "SALib"


# ====================================================
# SECTION 7: OUTPUT SETTINGS
# ====================================================
@dataclass
class ResultsFormats:
    csv: bool
    json: bool
    parquet: bool
    pickle: bool


@dataclass
class ResultsConfig:
    formats: ResultsFormats
    partial_results: bool


@dataclass
class OutputConfig:
    logging: bool
    results: ResultsConfig


# ====================================================
# SECTION 8: RESOURCES
# ====================================================
@dataclass
class CPUConfig:
    exclusive_cores: bool
    max_workers: Optional[int]
    internal_jobs: int = 1


@dataclass
class MemoryConfig:
    max_memory_gb: Optional[float]


@dataclass
class TimeoutsConfig:
    timeout_per_item: int
    timeout_total_batch: int


@dataclass
class ResourcesConfig:
    cpu: CPUConfig
    memory: MemoryConfig
    timeouts: TimeoutsConfig


# ====================================================
# SECTION 9: SYSTEM
# ====================================================
@dataclass
class SystemConfig:
    global_seed: Optional[int]


# ====================================================
# ROOT CONFIG
# ====================================================
@dataclass
class CSPBenchConfig:
    metadata: MetadataConfig
    tasks: TasksGroupConfig
    output: Optional[OutputConfig] = None
    resources: Optional[ResourcesConfig] = None
    system: Optional[SystemConfig] = None


# ====================================================
# LOADER / PARSER UTILITIES (merge settings + batch)
# ====================================================
def _load_yaml(path: Union[str, Path]) -> Dict[str, Any]:
    p = Path(path)
    if not p.exists():
        raise BatchConfigurationError(f"YAML file not found: {p}")
    try:
        with p.open(encoding="utf-8") as fh:
            data = yaml.safe_load(fh) or {}
        if not isinstance(data, dict):
            raise BatchConfigurationError(f"Root of YAML must be a mapping: {p}")
        return data
    except BatchConfigurationError:
        raise
    except Exception as e:  # noqa: BLE001
        raise BatchConfigurationError(f"Failed loading YAML {p}: {e}") from e


def _deep_merge(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    """Shallow+recursive dict merge returning new dict (override wins)."""
    result: Dict[str, Any] = {}
    keys = set(base) | set(override)
    for k in keys:
        if (
            k in base
            and k in override
            and isinstance(base[k], dict)
            and isinstance(override[k], dict)
        ):
            result[k] = _deep_merge(base[k], override[k])
        elif k in override:
            result[k] = override[k]
        else:
            result[k] = base[k]
    return result


def _slug(value: str) -> str:
    value = value.strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return re.sub(r"_+", "_", value).strip("_") or "task"


# ----------------------------- DATASETS -----------------------------
def _parse_datasets(
    batch: Dict[str, Any], global_seed: Optional[int] = None
) -> Dict[str, DatasetAny]:
    datasets_cfg = batch.get("datasets", [])
    if not isinstance(datasets_cfg, list):
        raise BatchConfigurationError("'datasets' deve ser lista")
    result: Dict[str, DatasetAny] = {}
    for d in datasets_cfg:
        if not isinstance(d, dict):
            raise BatchConfigurationError(
                "Entrada de dataset inválida (esperado mapeamento)"
            )
        did = d.get("id")
        dtype = d.get("type")
        name = d.get("name", did or "dataset")
        params = d.get("parameters", {})
        if not did or not dtype:
            raise BatchConfigurationError("Dataset precisa de campos 'id' e 'type'")
        if did in result:
            raise BatchConfigurationError(f"Dataset duplicado: {did}")
        if dtype == "synthetic":
            # Aplicar global_seed se definido, senão usar seed específico
            seed = global_seed if global_seed is not None else params.get("seed")
            ds = SyntheticDatasetConfig(
                id=did,
                name=name,
                mode=params.get("mode", "random"),
                n=params.get("n", 10),
                L=params.get("L", 20),
                alphabet=params.get("alphabet", "ACGT"),
                seed=seed,
                parameters_mode=d.get("parameters_mode")
                or params.get("parameters_mode")
                or {},
            )
        elif dtype == "file":
            fn = params.get("filename")
            if not fn:
                raise BatchConfigurationError(
                    "Dataset tipo 'file' requer parameters.filename"
                )
            ds = FileDatasetConfig(id=did, name=name, filename=fn)
        elif dtype == "entrez":
            ds = EntrezDatasetConfig(
                id=did,
                name=name,
                query=params.get("query", ""),
                db=params.get("db", "nucleotide"),
                retmax=params.get("retmax", 10),
                min_length=params.get("min_length"),
                max_length=params.get("max_length"),
                uniform_policy=params.get("uniform_policy"),
            )
        else:
            raise BatchConfigurationError(f"Tipo de dataset não suportado: {dtype}")

        result[did] = ds
    return result


# ----------------------------- ALGORITHMS PRESETS -----------------------------
def _parse_algorithm_presets(
    batch: Dict[str, Any], global_seed: Optional[int] = None
) -> Dict[str, AlgorithmsPresetConfig]:
    presets_cfg = batch.get("algorithms", [])
    if not isinstance(presets_cfg, list):
        raise BatchConfigurationError("'algorithms' deve ser lista")
    result: Dict[str, AlgorithmsPresetConfig] = {}
    for preset in presets_cfg:
        if not isinstance(preset, dict):
            raise BatchConfigurationError("Entrada de algorithms inválida (mapping)")
        pid = preset.get("id")
        if not pid:
            raise BatchConfigurationError("Algorithms preset sem 'id'")
        if pid in result:
            raise BatchConfigurationError(f"Algorithms preset duplicado: {pid}")
        alg_params: List[AlgParamsConfig] = []
        params_block = preset.get("algorithm_params", {}) or {}
        if not isinstance(params_block, dict):
            raise BatchConfigurationError("'algorithm_params' deve ser mapping")
        for alg_name, params in params_block.items():
            if params is None:
                params = {}
            if not isinstance(params, dict):
                raise BatchConfigurationError(
                    f"Parametros do algoritmo '{alg_name}' devem ser mapping"
                )
            # Aplicar global_seed se definido, sobrescrevendo seed existente
            final_params = params.copy()
            if global_seed is not None:
                final_params["seed"] = global_seed
            alg_params.append(AlgParamsConfig(name=alg_name, params=final_params))
        result[pid] = AlgorithmsPresetConfig(
            id=pid,
            name=preset.get("name") or pid,
            description=preset.get("description", ""),
            items=alg_params,
        )
    return result


# ----------------------------- TASK GROUPS -----------------------------
def _build_experiment_tasks(
    batch: Dict[str, Any],
    datasets: Dict[str, DatasetAny],
    presets: Dict[str, AlgorithmsPresetConfig],
) -> ExperimentTasksConfig:
    import logging
    logger = logging.getLogger("CSPBench.Config.Experiment")
    
    exp_section = batch.get("experiment") or {}
    tasks_list = exp_section.get("tasks") or exp_section.get("executions") or []
    if not isinstance(tasks_list, list):
        raise BatchConfigurationError("experiment.tasks deve ser lista")
        
    logger.debug(f"Processando {len(tasks_list)} tarefas de experimento")
    
    items: List[ExperimentTaskConfig] = []
    for _idx, t in enumerate(tasks_list):
        if not isinstance(t, dict):
            raise BatchConfigurationError("Item de task de experimento inválido")
        name = t.get("name")
        if not name:
            raise BatchConfigurationError("Task de experimento sem 'name'")
        tid = t.get("id", _slug(name))
        
        logger.debug(f"Processando task: id='{tid}', name='{name}'")
        
        ds_ids = t.get("datasets", []) or []
        alg_ids = t.get("algorithms", []) or []
        alg_config_ids = t.get("algorithm_config", []) or []
        
        # Support both 'algorithms' and 'algorithm_config' for backward compatibility
        if not alg_ids and alg_config_ids:
            alg_ids = alg_config_ids
            logger.debug(f"Task '{tid}': usando 'algorithm_config' em vez de 'algorithms' ({len(alg_ids)} itens)")
        
        if not ds_ids:
            raise BatchConfigurationError(f"Task '{tid}' sem datasets")
        if not alg_ids:
            logger.error(f"Task '{tid}': Erro - nem 'algorithms' nem 'algorithm_config' encontrados")
            logger.debug(f"Task '{tid}': Chaves disponíveis: {list(t.keys())}")
            raise BatchConfigurationError(
                f"Task '{tid}' sem algorithms (lista de ids de presets)"
            )
        ds_objs = []
        for d_id in ds_ids:
            if d_id not in datasets:
                raise BatchConfigurationError(
                    f"Dataset '{d_id}' não encontrado para task '{tid}'"
                )
            ds_objs.append(datasets[d_id])
        preset_objs: List[AlgorithmsPresetConfig] = []
        for p_id in alg_ids:
            if p_id not in presets:
                raise BatchConfigurationError(
                    f"Algorithms preset '{p_id}' não encontrado na task '{tid}'"
                )
            preset_objs.append(presets[p_id])
        repetitions = t.get("repetitions", 1)
        items.append(
            ExperimentTaskConfig(
                id=tid,
                name=name,
                datasets=ds_objs,
                algorithms=preset_objs,
                repetitions=repetitions,
            )
        )
    return ExperimentTasksConfig(type="experiment", items=items)


def _build_optimization_tasks(
    batch: Dict[str, Any],
    datasets: Dict[str, DatasetAny],
    presets: Dict[str, AlgorithmsPresetConfig],
) -> OptimizationTasksConfig:
    import logging
    logger = logging.getLogger("CSPBench.Config.Optimization")
    
    opt_section = batch.get("optimization") or {}
    tasks_list = opt_section.get("tasks") or opt_section.get("executions") or []
    
    if not isinstance(tasks_list, list):
        raise BatchConfigurationError("optimization.tasks deve ser lista")
    
    logger.debug(f"Processando {len(tasks_list)} tarefas de optimization")
    
    items: List[OptimizationTaskConfig] = []
    for t in tasks_list:
        if not isinstance(t, dict):
            raise BatchConfigurationError("Item de optimization inválido")
        name = t.get("name")
        if not name:
            raise BatchConfigurationError("Optimization task sem 'name'")
        tid = t.get("id", _slug(name))
        
        logger.debug(f"Processando task: id='{tid}', name='{name}'")
        
        ds_ids = t.get("datasets", []) or []
        alg_ids = t.get("algorithms", []) or []
        alg_config_ids = t.get("algorithm_config", []) or []
        
        # Support both 'algorithms' and 'algorithm_config' for backward compatibility
        if not alg_ids and alg_config_ids:
            alg_ids = alg_config_ids
            logger.debug(f"Task '{tid}': usando 'algorithm_config' em vez de 'algorithms' ({len(alg_ids)} itens)")
        
        if not ds_ids:
            raise BatchConfigurationError(f"Optimization task '{tid}' sem datasets")
        if not alg_ids:
            logger.error(f"Task '{tid}': Erro - nem 'algorithms' nem 'algorithm_config' encontrados")
            logger.debug(f"Task '{tid}': Chaves disponíveis: {list(t.keys())}")
            raise BatchConfigurationError(f"Optimization task '{tid}' sem algorithms")
        ds_objs = []
        for d_id in ds_ids:
            if d_id not in datasets:
                raise BatchConfigurationError(
                    f"Dataset '{d_id}' não encontrado para optimization '{tid}'"
                )
            ds_objs.append(datasets[d_id])
        preset_objs: List[AlgorithmsPresetConfig] = []
        for p_id in alg_ids:
            if p_id not in presets:
                raise BatchConfigurationError(
                    f"Algorithms preset '{p_id}' não encontrado em optimization '{tid}'"
                )
            preset_objs.append(presets[p_id])
        parameters = t.get("parameters", {}) or {}
        config = t.get("config") or {}
        study_name = t.get("study_name") or tid
        items.append(
            OptimizationTaskConfig(
                id=tid,
                name=name,
                datasets=ds_objs,
                algorithms=preset_objs,
                parameters=parameters,
                config=config,
                study_name=study_name,
            )
        )
    return OptimizationTasksConfig(type="optimization", items=items)  # type: ignore[arg-type]


def _build_sensitivity_tasks(
    batch: Dict[str, Any],
    datasets: Dict[str, DatasetAny],
    presets: Dict[str, AlgorithmsPresetConfig],
) -> SensitivityTasksConfig:
    import logging
    logger = logging.getLogger("CSPBench.Config.Sensitivity")
    
    sens_section = batch.get("sensitivity") or {}
    tasks_list = sens_section.get("tasks") or sens_section.get("executions") or []
    if not isinstance(tasks_list, list):
        raise BatchConfigurationError("sensitivity.tasks deve ser lista")
        
    logger.debug(f"Processando {len(tasks_list)} tarefas de análise de sensibilidade")
    
    items: List[SensitivityTaskConfig] = []
    for t in tasks_list:
        if not isinstance(t, dict):
            raise BatchConfigurationError("Item de sensitivity inválido")
        name = t.get("name")
        if not name:
            raise BatchConfigurationError("Sensitivity task sem 'name'")
        tid = t.get("id", _slug(name))
        
        logger.debug(f"Processando task: id='{tid}', name='{name}'")
        
        ds_ids = t.get("datasets", []) or []
        alg_ids = t.get("algorithms", []) or []
        alg_config_ids = t.get("algorithm_config", []) or []
        alg_singular_ids = t.get("algorithm", []) or []
        
        # Support 'algorithms', 'algorithm_config', and 'algorithm' (singular) for backward compatibility
        if not alg_ids and alg_config_ids:
            alg_ids = alg_config_ids
            logger.debug(f"Task '{tid}': usando 'algorithm_config' em vez de 'algorithms' ({len(alg_ids)} itens)")
        elif not alg_ids and alg_singular_ids:
            alg_ids = alg_singular_ids
            logger.debug(f"Task '{tid}': usando 'algorithm' (singular) em vez de 'algorithms' ({len(alg_ids)} itens)")
        
        if not ds_ids:
            raise BatchConfigurationError(f"Sensitivity task '{tid}' sem datasets")
        if not alg_ids:
            logger.error(f"Task '{tid}': Erro - nem 'algorithms' nem 'algorithm_config' nem 'algorithm' encontrados")
            logger.debug(f"Task '{tid}': Chaves disponíveis: {list(t.keys())}")
            raise BatchConfigurationError(f"Sensitivity task '{tid}' sem algorithms")
        ds_objs = []
        for d_id in ds_ids:
            if d_id not in datasets:
                raise BatchConfigurationError(
                    f"Dataset '{d_id}' não encontrado para sensitivity '{tid}'"
                )
            ds_objs.append(datasets[d_id])
        preset_objs: List[AlgorithmsPresetConfig] = []
        for p_id in alg_ids:
            if p_id not in presets:
                raise BatchConfigurationError(
                    f"Algorithms preset '{p_id}' não encontrado em sensitivity '{tid}'"
                )
            preset_objs.append(presets[p_id])
        parameters = t.get("parameters", {}) or {}
        config = t.get("config") or {}
        method = t.get("method", "morris")
        items.append(
            SensitivityTaskConfig(
                id=tid,
                name=name,
                datasets=ds_objs,
                algorithms=preset_objs,
                parameters=parameters,
                config=config,
                method=method,
            )
        )
    return SensitivityTasksConfig(type="sensitivity", items=items)  # type: ignore[arg-type]


def _build_tasks_group(
    batch: Dict[str, Any],
    datasets: Dict[str, DatasetAny],
    presets: Dict[str, AlgorithmsPresetConfig],
) -> TasksGroupConfig:
    task_section = batch.get("task", {})
    t_type = task_section.get("type")
    if t_type == "experiment":
        return _build_experiment_tasks(batch, datasets, presets)
    if t_type == "optimization":
        return _build_optimization_tasks(batch, datasets, presets)
    if t_type == "sensitivity":
        return _build_sensitivity_tasks(batch, datasets, presets)
    raise BatchConfigurationError(
        "Campo task.type ausente ou inválido (experiment|optimization|sensitivity)"
    )


# ----------------------------- OUTPUT / RESOURCES / SYSTEM -----------------------------
def _build_output(
    settings: Dict[str, Any], batch: Dict[str, Any]
) -> Optional[OutputConfig]:
    # Both structures identical per template: output.logging + output.results.{formats..., partial_results}
    merged = _deep_merge(settings.get("output", {}), batch.get("output", {}))
    if not merged:
        return None
    results = merged.get("results", {})
    formats_block = results.get("formats", {})
    rf = ResultsFormats(
        csv=bool(formats_block.get("csv", False)),
        json=bool(formats_block.get("json", False)),
        parquet=bool(formats_block.get("parquet", False)),
        pickle=bool(formats_block.get("pickle", False)),
    )
    rc = ResultsConfig(
        formats=rf, partial_results=bool(results.get("partial_results", False))
    )
    return OutputConfig(logging=bool(merged.get("logging", True)), results=rc)


def _build_resources(
    settings: Dict[str, Any], batch: Dict[str, Any]
) -> Optional[ResourcesConfig]:
    merged = _deep_merge(settings.get("resources", {}), batch.get("resources", {}))
    if not merged:
        return None
    cpu_block = merged.get("cpu") or {}
    mem_block = merged.get("memory") or {}
    timeouts_block = merged.get("timeouts") or {}
    cpu = CPUConfig(
        exclusive_cores=bool(cpu_block.get("exclusive_cores", True)),
        max_workers=cpu_block.get("max_workers") or cpu_block.get("max_cores"),  # Support both for compatibility
        internal_jobs=int(cpu_block.get("internal_jobs", 1)),
    )
    memory = MemoryConfig(max_memory_gb=mem_block.get("max_memory_gb"))
    if not timeouts_block:
        raise BatchConfigurationError(
            "timeouts obrigatorio em resources (timeout_per_item, timeout_total_batch)"
        )
    timeouts = TimeoutsConfig(
        timeout_per_item=int(timeouts_block.get("timeout_per_item", 3600)),
        timeout_total_batch=int(timeouts_block.get("timeout_total_batch", 86400)),
    )
    return ResourcesConfig(cpu=cpu, memory=memory, timeouts=timeouts)


def _build_system(
    settings: Dict[str, Any], batch: Dict[str, Any]
) -> Optional[SystemConfig]:
    merged = _deep_merge(settings.get("system", {}), batch.get("system", {}))
    if not merged:
        return None
    return SystemConfig(global_seed=merged.get("global_seed"))


def load_cspbench_config(
    batch_path: Union[str, Path], settings_path: Optional[Union[str, Path]] = None
) -> CSPBenchConfig:
    """Carrega configuração unificada (settings.yaml < batch.yaml).

    Precedência: settings (base) < batch (override).
    Faz deep-merge do root antes de parsear (listas são sobrescritas por inteiro).
    Aplica global_seed a todos os datasets synthetic quando definido.

    Parametros:
        settings_path: caminho para settings.yaml (defaults globais)
        batch_path: caminho para batch.yaml (override específico)

    Retorna:
        CSPBenchConfig populado com dataclasses de alto nível.
    """

    if settings_path is None:
        settings_path = os.environ.get("SETTINGS_PATH")
        if not settings_path:
            raise BatchConfigurationError(
                "settings_path não informado e SETTINGS_PATH não está definido no ambiente."
            )

    settings = _load_yaml(settings_path)
    batch = _load_yaml(batch_path)

    # Deep merge root: NOTE listas são substituídas (comportamento intencional aqui).
    merged_config = _deep_merge(settings, batch)

    # SYSTEM - carregar primeiro para obter global_seed
    system_cfg = _build_system({}, merged_config)
    global_seed = system_cfg.global_seed if system_cfg else None

    # METADATA (campos obrigatórios devem estar após merge)
    meta_block = merged_config.get("metadata") or {}
    required_meta = ["name", "description", "author", "version", "creation_date"]
    missing = [f for f in required_meta if f not in meta_block]
    if missing:
        raise BatchConfigurationError(
            f"Campos de metadata ausentes: {', '.join(missing)}"
        )
    metadata = MetadataConfig(
        name=meta_block["name"],
        description=meta_block["description"],
        author=meta_block["author"],
        version=meta_block["version"],
        creation_date=meta_block["creation_date"],
        tags=meta_block.get("tags", []),
    )

    # DATASETS & ALGORITHMS PRESETS (com global_seed aplicado aos datasets synthetic e algoritmos)
    datasets = _parse_datasets(merged_config, global_seed)
    presets = _parse_algorithm_presets(merged_config, global_seed)

    # TASKS GROUP (apenas experiment implementado)
    tasks_group = _build_tasks_group(merged_config, datasets, presets)

    # OUTPUT / RESOURCES
    output_cfg = _build_output({}, merged_config)
    resources_cfg = _build_resources({}, merged_config)

    return CSPBenchConfig(
        metadata=metadata,
        tasks=tasks_group,  # type: ignore[arg-type]
        output=output_cfg,
        resources=resources_cfg,
        system=system_cfg,
    )


# ---- Compat aliases (mantêm compatibilidade com nomes antigos) ----
# Classes já terminadas com Config permanecem inalteradas.
Metadata = MetadataConfig
SyntheticDataset = SyntheticDatasetConfig
FileDataset = FileDatasetConfig
EntrezDataset = EntrezDatasetConfig
AlgorithmsPreset = AlgorithmsPresetConfig
AlgParams = AlgParamsConfig
Task = TaskConfig
ExperimentTask = ExperimentTaskConfig
OptimizationTask = OptimizationTaskConfig
SensitivityTask = SensitivityTaskConfig
TasksGroup = TasksGroupConfig
ExperimentTasks = ExperimentTasksConfig
OptimizationTasks = OptimizationTasksConfig
SensitivityTasks = SensitivityTasksConfig

__all__ = [
    # Loader/root
    "CSPBenchConfig",
    "load_cspbench_config",
    # Aliases antigos (exportados)
    "Metadata",
    "DatasetAny",
    "SyntheticDataset",
    "FileDataset",
    "EntrezDataset",
    "AlgorithmsPreset",
    "AlgParams",
    "ExperimentTasks",
]
