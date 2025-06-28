"""
algorithm_executor.py
=====================

Executor de algoritmos em threads isoladas com timeout e cancelamento.

Classes:
    AlgorithmExecutor: Executa algoritmos em threads separadas com controle de tempo.
    TimeoutException: Exceção lançada quando algoritmo excede tempo limite.
"""

import threading
import time
import queue
from typing import Any, Callable, Optional, Tuple, Union
import logging

logger = logging.getLogger(__name__)

class TimeoutException(Exception):
    """Exceção lançada quando um algoritmo excede o tempo limite."""
    pass

class AlgorithmExecutor:
    """Executa algoritmos em threads isoladas com timeout configurável."""
    
    def __init__(self, timeout_seconds: int):
        self.timeout = timeout_seconds
        self.stop_event = threading.Event()
        
    def execute_with_timeout(
        self, 
        algorithm_instance: Any,
        progress_callback: Optional[Callable[[str], None]] = None,
        warning_callback: Optional[Callable[[str], None]] = None
    ) -> Tuple[Optional[Any], Union[int, float], dict]:
        """
        Executa um algoritmo em thread isolada com timeout.
        
        Args:
            algorithm_instance: Instância do algoritmo a ser executado
            progress_callback: Callback para reportar progresso
            warning_callback: Callback para reportar avisos
            
        Returns:
            tuple: (center, distance, info) onde distance pode ser int ou float('inf')
            
        Raises:
            TimeoutException: Se o algoritmo exceder o tempo limite
        """
        result_queue = queue.Queue()
        self.stop_event.clear()
        
        def run_algorithm():
            """Função executada na thread do algoritmo."""
            try:
                # Configurar callbacks se disponíveis
                if hasattr(algorithm_instance, 'set_progress_callback') and progress_callback is not None:
                    # Wrapper para verificar cancelamento durante progresso
                    def monitored_progress(msg: str):
                        if self.stop_event.is_set():
                            raise TimeoutException("Algoritmo cancelado por timeout")
                        if progress_callback is not None:  # Verificação adicional para Pylance
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
                
            except TimeoutException:
                result_queue.put(('timeout', None, float('inf'), {'erro': 'timeout'}))
            except Exception as e:
                error_msg = str(e)
                logger.error(f"Erro na execução do algoritmo: {error_msg}")
                result_queue.put(('error', None, float('inf'), {'erro': error_msg}))
        
        # Iniciar thread do algoritmo
        algorithm_thread = threading.Thread(target=run_algorithm)
        algorithm_thread.daemon = True
        algorithm_thread.start()
        
        # Aguardar resultado ou timeout
        start_time = time.time()
        
        while algorithm_thread.is_alive():
            elapsed = time.time() - start_time
            
            if elapsed >= self.timeout:
                # Timeout atingido - sinalizar parada
                self.stop_event.set()
                
                # Aguardar um pouco para thread terminar graciosamente
                algorithm_thread.join(timeout=2.0)
                
                if algorithm_thread.is_alive():
                    logger.warning("Thread do algoritmo não terminou graciosamente após timeout")
                
                raise TimeoutException(f"Algoritmo excedeu tempo limite de {self.timeout}s")
            
            # Verificar se há resultado disponível
            try:
                result = result_queue.get_nowait()
                status, center, distance, info = result
                
                if status == 'timeout':
                    raise TimeoutException("Algoritmo cancelado por timeout durante execução")
                elif status == 'error':
                    return center, distance, info
                else:
                    return center, distance, info
                    
            except queue.Empty:
                time.sleep(0.1)  # Aguardar um pouco antes de verificar novamente
        
        # Thread terminou - obter resultado final
        try:
            result = result_queue.get_nowait()
            status, center, distance, info = result
            return center, distance, info
        except queue.Empty:
            # Caso raro onde thread terminou sem colocar resultado
            return None, float('inf'), {'erro': 'Thread terminou sem resultado'}
    
    def cancel(self):
        """Cancela a execução atual."""
        self.stop_event.set()
