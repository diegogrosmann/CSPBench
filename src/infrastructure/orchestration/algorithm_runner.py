"""AlgorithmRunner: executa uma unidade (repetition/trial/sample).
Assume retorno (center, distance, metadata) para CSPAlgorithm.
"""

from __future__ import annotations
from typing import Any, Dict, Tuple
import time
import traceback
from src.domain.algorithms import global_registry, CSPAlgorithm
from src.infrastructure.logging_config import get_logger

# Create module logger
logger = get_logger("CSPBench.AlgorithmRunner")


class AlgorithmExecutionError(Exception):
    pass


def run_algorithm(
    algorithm_name: str,
    strings: list[str],
    alphabet: str,
    params: Dict[str, Any],
) -> Dict[str, Any]:
    logger.info(f"Iniciando execução do algoritmo: {algorithm_name}")
    logger.debug(f"Parâmetros de entrada: strings={len(strings)}, alphabet='{alphabet}', params={params}")
    
    t0 = time.time()
    try:
        # Registry lookup with detailed logging
        logger.debug(f"Procurando algoritmo '{algorithm_name}' no registro global")
        cls = global_registry.get(algorithm_name)
        if cls is None:
            error_msg = f"Algorithm not registered: {algorithm_name}"
            logger.error(error_msg)
            logger.debug(f"Algoritmos disponíveis: {list(global_registry.keys())}")
            raise AlgorithmExecutionError(error_msg)
        
        logger.info(f"Algoritmo '{algorithm_name}' encontrado: {cls.__name__}")
        
        # Algorithm instantiation with detailed logging
        logger.debug(f"Criando instância do algoritmo com {len(strings)} strings")
        alg: CSPAlgorithm = cls(strings, alphabet, **params)
        logger.info(f"Instância do algoritmo '{algorithm_name}' criada com sucesso")
        
        # Algorithm execution with timing
        logger.info(f"Iniciando execução do algoritmo '{algorithm_name}'")
        start_time = time.time()
        
        center, dist, meta = alg.run()
        
        execution_time = time.time() - start_time
        total_duration = time.time() - t0
        
        logger.info(f"Algoritmo '{algorithm_name}' executado com sucesso")
        logger.info(f"Resultado: distância={dist}, tempo_execução={execution_time:.3f}s, tempo_total={total_duration:.3f}s")
        logger.debug(f"Centro encontrado: {center}")
        logger.debug(f"Metadata: {meta}")
        
        # Capture complete data with logging
        algorithm_metadata = alg.get_metadata() if hasattr(alg, 'get_metadata') else {}
        execution_history = alg.get_history() if hasattr(alg, 'get_history') else []
        actual_params = getattr(alg, 'params', params)
        
        if algorithm_metadata:
            logger.debug(f"Metadata do algoritmo capturada: {len(algorithm_metadata)} campos")
        if execution_history:
            logger.debug(f"Histórico de execução capturado: {len(execution_history)} entradas")
        
        result = {
            "status": "ok",
            "objective": float(dist),  # padrão: distância int -> float
            "center": center,
            "distance": dist,
            "duration_s": total_duration,
            "metadata": meta,
            # Extended data for complete storage
            "algorithm_metadata": algorithm_metadata,
            "execution_history": execution_history,
            "actual_params": actual_params,
        }
        
        logger.info(f"Execução do algoritmo '{algorithm_name}' finalizada com sucesso")
        return result
        
    except Exception as e:  # noqa: BLE001
        duration = time.time() - t0
        
        # Detect specific cases that should be treated as warnings, not errors
        error_message = str(e)
        is_resource_limitation = (
            "exceeds practical limit" in error_message or
            "requires all strings to have the same length" in error_message
        )
        
        if is_resource_limitation:
            # These are expected limitations, not programming errors
            if "exceeds practical limit" in error_message:
                log_msg = f"Algoritmo '{algorithm_name}' skipped: problema muito grande para recursos disponíveis"
                user_friendly_msg = f"Dataset muito grande para {algorithm_name} - considere usar um dataset menor"
            elif "requires all strings to have the same length" in error_message:
                log_msg = f"Algoritmo '{algorithm_name}' skipped: dataset contém strings de tamanhos diferentes"
                user_friendly_msg = f"Dataset incompatível com {algorithm_name} - strings devem ter o mesmo comprimento"
            
            logger.warning(log_msg)
            logger.info(f"Detalhes técnicos: {error_message}")
        else:
            # Genuine programming/runtime errors
            error_msg = f"Erro na execução do algoritmo '{algorithm_name}': {e}"
            user_friendly_msg = error_message
            logger.error(error_msg)
            logger.error(f"Tipo do erro: {e.__class__.__name__}")
            logger.debug(f"Traceback completo:", exc_info=True)
        
        error_result = {
            "status": "skipped" if is_resource_limitation else "error",
            "error_type": e.__class__.__name__,
            "error_message": user_friendly_msg,
            "technical_details": error_message,
            "traceback": traceback.format_exc(limit=5),
            "duration_s": duration,
            "objective": float("inf"),
            # Extended data even for errors
            "algorithm_metadata": {},
            "execution_history": [],
            "actual_params": params,
        }
        
        result_type = "skipped" if is_resource_limitation else "error"
        logger.info(f"Algoritmo '{algorithm_name}' finalizado com status '{result_type}' após {duration:.3f}s")
        return error_result
