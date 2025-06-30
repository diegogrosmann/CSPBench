"""
algorithm_executor.py
=====================

Executor de algoritmos em threads isoladas com timeout e monitoramento de recursos.
"""

import threading
import time
import queue
from typing import Any, Callable, Optional, Tuple, Union
import logging
from utils.resource_monitor import ResourceMonitor, ResourceLimits, force_garbage_collection

logger = logging.getLogger(__name__)

class TimeoutException(Exception):
    """Exceção lançada quando um algoritmo excede o tempo limite."""
    pass

class ResourceLimitException(Exception):
    """Exceção lançada quando recursos do sistema são excedidos."""
    pass

class AlgorithmExecutor:
    """Executa algoritmos em threads isoladas com timeout e monitoramento de recursos."""
    
    def __init__(self, timeout_seconds: int):
        self.timeout = timeout_seconds
        self.stop_event = threading.Event()
        
        # Configurar limites baseados no timeout e dataset
        limits = ResourceLimits(
            max_memory_mb=int(2048),  # 2GB máximo
            max_iterations=min(100000, timeout_seconds * 1000),  # Baseado no timeout
            check_interval=min(2.0, timeout_seconds / 10)  # Verificar mais frequentemente em timeouts curtos
        )
        self.resource_monitor = ResourceMonitor(limits)
        self.resource_violation = False
        
    def execute_with_timeout(
        self, 
        algorithm_instance: Any,
        progress_callback: Optional[Callable[[str], None]] = None,
        warning_callback: Optional[Callable[[str], None]] = None
    ) -> Tuple[Optional[Any], Union[int, float], dict]:
        """
        Executa um algoritmo em thread isolada com timeout e monitoramento de recursos.
        """
        result_queue = queue.Queue()
        self.stop_event.clear()
        self.resource_violation = False
        
        # Força garbage collection antes de começar
        force_garbage_collection()
        
        # Configurar monitor de recursos
        def resource_violation_handler(violation_msg: str):
            logger.warning(f"Violação de recursos detectada: {violation_msg}")
            self.resource_violation = True
            self.stop_event.set()
            if warning_callback:
                warning_callback(f"Recursos limitados: {violation_msg}")
        
        self.resource_monitor.set_violation_callback(resource_violation_handler)
        
        def run_algorithm():
            """Função executada na thread do algoritmo."""
            try:
                # Configurar callbacks se disponíveis
                if hasattr(algorithm_instance, 'set_progress_callback') and progress_callback is not None:
                    def monitored_progress(msg: str):
                        if self.stop_event.is_set():
                            if self.resource_violation:
                                raise ResourceLimitException("Algoritmo cancelado por limite de recursos")
                            else:
                                raise TimeoutException("Algoritmo cancelado por timeout")
                        if progress_callback is not None:
                            progress_callback(msg)
                    algorithm_instance.set_progress_callback(monitored_progress)
                    
                if hasattr(algorithm_instance, 'set_warning_callback') and warning_callback is not None:
                    algorithm_instance.set_warning_callback(warning_callback)
                
                # Executar algoritmo
                center, distance = algorithm_instance.run()
                
                # Coletar informações adicionais
                info = {'melhor_string': center}
                
                # Tentar obter número de iterações de diferentes atributos
                if hasattr(algorithm_instance, 'geracao'):
                    info['iteracoes'] = algorithm_instance.geracao
                elif hasattr(algorithm_instance, 'iterations'):
                    info['iteracoes'] = algorithm_instance.iterations
                elif hasattr(algorithm_instance, 'num_iteracoes'):
                    info['iteracoes'] = algorithm_instance.num_iteracoes
                else:
                    info['iteracoes'] = 1
                
                result_queue.put(('success', center, distance, info))
                
            except (TimeoutException, ResourceLimitException) as e:
                error_type = 'resource_limit' if isinstance(e, ResourceLimitException) else 'timeout'
                result_queue.put((error_type, None, float('inf'), {'erro': str(e)}))
            except Exception as e:
                error_msg = str(e)
                logger.error(f"Erro na execução do algoritmo: {error_msg}")
                result_queue.put(('error', None, float('inf'), {'erro': error_msg}))
        
        # Iniciar monitoramento de recursos
        self.resource_monitor.start_monitoring()
        
        # Iniciar thread do algoritmo
        algorithm_thread = threading.Thread(target=run_algorithm)
        algorithm_thread.daemon = True
        algorithm_thread.start()
        
        # Aguardar resultado ou timeout com verificações mais frequentes
        start_time = time.time()
        check_interval = 0.5  # Verificar a cada 0.5 segundos
        
        while algorithm_thread.is_alive():
            elapsed = time.time() - start_time
            
            if elapsed >= self.timeout:
                # Timeout atingido - sinalizar parada
                self.stop_event.set()
                
                # Aguardar um pouco para thread terminar graciosamente
                algorithm_thread.join(timeout=2.0)
                
                if algorithm_thread.is_alive():
                    logger.warning("Thread do algoritmo não terminou graciosamente após timeout")
                
                self.resource_monitor.stop_monitoring()
                raise TimeoutException(f"Algoritmo excedeu tempo limite de {self.timeout}s")
            
            # Verificar se há resultado disponível
            try:
                result = result_queue.get_nowait()
                status, center, distance, info = result
                
                self.resource_monitor.stop_monitoring()
                
                if status == 'timeout':
                    raise TimeoutException("Algoritmo cancelado por timeout durante execução")
                elif status == 'resource_limit':
                    raise ResourceLimitException("Algoritmo cancelado por limite de recursos")
                elif status == 'error':
                    return center, distance, info
                else:
                    return center, distance, info
                    
            except queue.Empty:
                time.sleep(check_interval)
        
        # Thread terminou - obter resultado final
        self.resource_monitor.stop_monitoring()
        
        try:
            result = result_queue.get_nowait()
            status, center, distance, info = result
            return center, distance, info
        except queue.Empty:
            return None, float('inf'), {'erro': 'Thread terminou sem resultado'}
    
    def cancel(self):
        """Cancela a execução atual."""
        self.stop_event.set()
        self.resource_monitor.stop_monitoring()
        self.resource_monitor.stop_monitoring()
