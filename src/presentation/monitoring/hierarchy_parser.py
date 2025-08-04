"""Parser for extracting hierarchical structure from YAML configuration."""

from typing import Any, Dict, List, Tuple
import logging

from .interfaces import ExecutionHierarchy, HierarchyLevel, ExecutionLevel

logger = logging.getLogger(__name__)


class HierarchyParser:
    """Parser para extrair estrutura hierárquica de configurações YAML."""

    @staticmethod
    def parse_yaml_hierarchy(config: Dict[str, Any]) -> ExecutionHierarchy:
        """
        Extrai estrutura hierárquica completa do YAML de configuração.

        Args:
            config: Configuração YAML carregada como dicionário

        Returns:
            ExecutionHierarchy: Estrutura hierárquica inicializada com totais

        Raises:
            ValueError: Se configuração for inválida
        """
        hierarchy = ExecutionHierarchy()
        
        try:
            # Extrair estrutura hierárquica
            executions_data = HierarchyParser._extract_executions_structure(config)
            
            # Configurar níveis hierárquicos
            hierarchy.structure = executions_data
            
            # Calcular totais para cada nível
            total_executions = len(executions_data.get("executions", []))
            hierarchy.update_level(ExecutionLevel.EXECUTION, 0, total_executions)
            
            logger.info(f"HierarchyParser: Parsed {total_executions} executions")
            
            return hierarchy
            
        except Exception as e:
            logger.error(f"HierarchyParser: Error parsing hierarchy: {e}")
            raise ValueError(f"Erro ao extrair hierarquia do YAML: {e}")

    @staticmethod
    def _extract_executions_structure(config: Dict[str, Any]) -> Dict[str, Any]:
        """
        Extrai estrutura de execuções da configuração.

        Args:
            config: Configuração completa

        Returns:
            Dict com estrutura de execuções organizadas
        """
        task_type = config.get("task", {}).get("type", "execution")
        
        if task_type == "execution":
            return HierarchyParser._extract_execution_structure(config)
        elif task_type == "optimization":
            return HierarchyParser._extract_optimization_structure(config)
        elif task_type == "sensitivity":
            return HierarchyParser._extract_sensitivity_structure(config)
        else:
            raise ValueError(f"Tipo de tarefa não suportado: {task_type}")

    @staticmethod
    def _extract_execution_structure(config: Dict[str, Any]) -> Dict[str, Any]:
        """Extrai estrutura para tarefas de execução."""
        execution_config = config.get("execution", {})
        executions_list = execution_config.get("executions", [])
        
        datasets = config.get("datasets", [])
        algorithms_configs = config.get("algorithms", [])
        
        # Criar mapeamento de datasets
        dataset_map = {ds["id"]: ds for ds in datasets}
        
        # Criar mapeamento de configurações de algoritmos
        config_map = {alg["id"]: alg for alg in algorithms_configs}
        
        structured_executions = {}
        
        for exec_item in executions_list:
            exec_name = exec_item.get("name", "Unnamed")
            exec_datasets = exec_item.get("datasets", [])
            exec_algorithms = exec_item.get("algorithms", ["default_config"])
            exec_repetitions = exec_item.get("repetitions", 1)
            
            # Estruturar datasets para esta execução
            exec_structure = {
                "name": exec_name,
                "datasets": {},
                "total_datasets": len(exec_datasets),
                "total_configs": len(exec_algorithms),
                "repetitions": exec_repetitions
            }
            
            for dataset_id in exec_datasets:
                dataset_info = dataset_map.get(dataset_id, {"id": dataset_id, "name": dataset_id})
                
                # Estruturar configs para este dataset
                dataset_structure = {
                    "name": dataset_info.get("name", dataset_id),
                    "configs": {},
                    "total_configs": len(exec_algorithms)
                }
                
                for config_id in exec_algorithms:
                    config_info = config_map.get(config_id, {"id": config_id, "algorithms": []})
                    config_algorithms = config_info.get("algorithms", [])
                    
                    # Estruturar algoritmos para esta config
                    config_structure = {
                        "name": config_info.get("name", config_id),
                        "algorithms": {},
                        "total_algorithms": len(config_algorithms),
                        "repetitions": exec_repetitions
                    }
                    
                    # Estruturar cada algoritmo
                    for algo_name in config_algorithms:
                        config_structure["algorithms"][algo_name] = {
                            "name": algo_name,
                            "repetitions": exec_repetitions
                        }
                    
                    dataset_structure["configs"][config_id] = config_structure
                
                exec_structure["datasets"][dataset_id] = dataset_structure
            
            structured_executions[exec_name] = exec_structure
        
        return {
            "executions": structured_executions,
            "total_executions": len(structured_executions)
        }

    @staticmethod
    def _extract_optimization_structure(config: Dict[str, Any]) -> Dict[str, Any]:
        """Extrai estrutura para tarefas de otimização."""
        optimization_config = config.get("optimization", {})
        optimizations_list = optimization_config.get("optimizations", [])
        
        # Simplificado para otimização - pode ser expandido conforme necessário
        structured_optimizations = {}
        for i, opt in enumerate(optimizations_list):
            opt_name = opt.get("name", f"Optimization {i+1}")
            structured_optimizations[opt_name] = {"name": opt_name}
        
        return {
            "executions": structured_optimizations,
            "total_executions": len(optimizations_list)
        }

    @staticmethod
    def _extract_sensitivity_structure(config: Dict[str, Any]) -> Dict[str, Any]:
        """Extrai estrutura para tarefas de análise de sensibilidade."""
        sensitivity_config = config.get("sensitivity", {})
        analyses_list = sensitivity_config.get("analyses", [])
        
        # Simplificado para sensibilidade - pode ser expandido conforme necessário
        return {
            "executions": [{"name": analysis.get("name", f"Analysis {i+1}")} 
                          for i, analysis in enumerate(analyses_list)],
            "total_executions": len(analyses_list)
        }

    @staticmethod
    def extract_algorithm_counts(config: Dict[str, Any]) -> Dict[str, int]:
        """
        Extrai contagem de algoritmos por configuração.

        Args:
            config: Configuração completa

        Returns:
            Dict mapeando config_id para número de algoritmos
        """
        algorithms_configs = config.get("algorithms", [])
        
        counts = {}
        for alg_config in algorithms_configs:
            config_id = alg_config.get("id", "unknown")
            algorithms_list = alg_config.get("algorithms", [])
            counts[config_id] = len(algorithms_list)
        
        return counts

    @staticmethod
    def get_execution_totals(config: Dict[str, Any]) -> Tuple[int, int, int, int]:
        """
        Calcula totais para toda a configuração.

        Args:
            config: Configuração completa

        Returns:
            Tuple (total_executions, total_datasets, total_configs, total_algorithms)
        """
        try:
            structure = HierarchyParser._extract_executions_structure(config)
            executions = structure.get("executions", [])
            
            total_executions = len(executions)
            total_datasets = 0
            total_configs = 0
            total_algorithms = 0
            
            for execution in executions:
                if isinstance(execution, dict) and "datasets" in execution:
                    exec_datasets = execution["datasets"]
                    total_datasets += len(exec_datasets)
                    
                    for dataset_id, dataset_info in exec_datasets.items():
                        if isinstance(dataset_info, dict) and "configs" in dataset_info:
                            configs = dataset_info["configs"]
                            total_configs += len(configs)
                            
                            for config_id, config_info in configs.items():
                                if isinstance(config_info, dict) and "algorithms" in config_info:
                                    algorithms = config_info["algorithms"]
                                    total_algorithms += len(algorithms)
            
            return total_executions, total_datasets, total_configs, total_algorithms
            
        except Exception as e:
            logger.warning(f"Erro ao calcular totais: {e}")
            return 0, 0, 0, 0

    @staticmethod
    def validate_config(config: Dict[str, Any]) -> List[str]:
        """
        Valida configuração YAML para estrutura hierárquica.

        Args:
            config: Configuração para validar

        Returns:
            Lista de erros encontrados (vazia se válida)
        """
        errors = []
        
        # Verificar campos obrigatórios
        required_fields = ["task", "datasets", "algorithms"]
        for field in required_fields:
            if field not in config:
                errors.append(f"Campo obrigatório ausente: {field}")
        
        # Verificar tipo de tarefa
        task_config = config.get("task", {})
        task_type = task_config.get("type")
        if task_type not in ["execution", "optimization", "sensitivity"]:
            errors.append(f"Tipo de tarefa inválido: {task_type}")
        
        # Verificar estrutura específica por tipo
        if task_type == "execution":
            execution_config = config.get("execution", {})
            if "executions" not in execution_config:
                errors.append("Campo 'executions' ausente em configuração de execução")
        
        # Verificar datasets
        datasets = config.get("datasets", [])
        if not isinstance(datasets, list) or len(datasets) == 0:
            errors.append("Pelo menos um dataset deve ser definido")
        
        # Verificar algoritmos
        algorithms = config.get("algorithms", [])
        if not isinstance(algorithms, list) or len(algorithms) == 0:
            errors.append("Pelo menos uma configuração de algoritmos deve ser definida")
        
        return errors


def create_hierarchy_from_yaml(config: Dict[str, Any]) -> ExecutionHierarchy:
    """
    Função conveniente para criar hierarquia a partir de configuração YAML.

    Args:
        config: Configuração YAML carregada

    Returns:
        ExecutionHierarchy: Estrutura hierárquica inicializada

    Raises:
        ValueError: Se configuração for inválida
    """
    # Validar configuração
    errors = HierarchyParser.validate_config(config)
    if errors:
        raise ValueError(f"Configuração inválida: {'; '.join(errors)}")
    
    # Criar hierarquia
    return HierarchyParser.parse_yaml_hierarchy(config)
