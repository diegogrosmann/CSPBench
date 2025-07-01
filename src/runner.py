import threading
import time
import itertools
import signal
import sys
import logging

from utils.resource_monitor import check_algorithm_feasibility, get_safe_memory_limit, force_garbage_collection
from src.algorithm_executor import AlgorithmExecutor, TimeoutException, ResourceLimitException

from utils.config import ALGORITHM_TIMEOUT

"""
Utilitários para execução de algoritmos e exibição de progresso (spinner).

Classes:
    Spinner: Exibe animação de progresso durante execução de algoritmos.

Funções:
    execute_algorithm_runs(...): Executa múltiplas execuções de um algoritmo, coletando resultados.
"""

class Spinner:
    def __init__(self, prefix: str, console=None):
        self.prefix = prefix
        self.stop_event = threading.Event()
        self.thread = None
        self.spinner = itertools.cycle(['   ', '.  ', '.. ', '...'])
        self.console = console
        self.progress_override = False

    def start(self):
        if self.thread and self.thread.is_alive():
            return
        self.stop_event.clear()
        self.progress_override = False  # Reset override flag
        self.thread = threading.Thread(target=self._animate)
        self.thread.daemon = True
        self.thread.start()

    def _animate(self):
        while not self.stop_event.is_set():
            # Só mostrar spinner se não houver progresso específico
            if not self.progress_override:
                # Usar console.print_inline para ser thread-safe
                if self.console:
                    self.console.print_inline(f"{self.prefix:<25s}{next(self.spinner)}")
                else:
                    print(f"\r{self.prefix:<25s}{next(self.spinner)}", end="", flush=True)
            time.sleep(0.3)

    def set_progress_override(self, value: bool):
        self.progress_override = value

    def stop(self):
        self.stop_event.set()
        if self.thread and self.thread.is_alive():
            self.thread.join(timeout=0.5)  # Timeout mais curto para saída rápida
        # Limpar a linha após parar o spinner
        if self.console:
            self.console.print_inline(f"{' ' * 70}\r")
        else:
            print(f"\r{' ' * 70}\r", end="", flush=True)

def execute_algorithm_runs(alg_name, AlgClass, seqs, alphabet, num_execs, baseline_val, console=None, timeout=None):
    """
    Executa múltiplas execuções de um algoritmo em threads isoladas com timeout e monitoramento de recursos.
    """
    logger = logging.getLogger(__name__)
    logger.debug(f"[Runner] Iniciando {alg_name} (baseline_val={baseline_val})")
    
    if timeout is None:
        timeout = ALGORITHM_TIMEOUT
    
    # Força garbage collection antes de começar
    force_garbage_collection()
    
    # Verificar viabilidade do algoritmo antes da execução
    n, L = len(seqs), len(seqs[0])
    is_feasible, feasibility_msg = check_algorithm_feasibility(n, L, alg_name)
    
    if not is_feasible:
        if console:
            console.print(f"{alg_name:<25s}... INVIÁVEL: {feasibility_msg}")
        else:
            print(f"\r{alg_name:<25s}... INVIÁVEL: {feasibility_msg}")
        
        return [{
            'tempo': 0.0,
            'iteracoes': 0,
            'distancia': float('inf'),
            'melhor_string': '',
            'erro': f'Inviável: {feasibility_msg}'
        }]
    
    # Ajustar timeout baseado na complexidade estimada
    if alg_name == 'DP-CSP' and timeout < 60:
        adjusted_timeout = min(timeout * 3, 300)  # Máximo 5 minutos para DP-CSP
        if console:
            console.print(f"Timeout ajustado para {alg_name}: {adjusted_timeout}s")
        timeout = adjusted_timeout
        
    is_deterministic = getattr(AlgClass, 'is_deterministic', False)
    actual_execs = 1 if is_deterministic else num_execs
    executions = []
    
    for i in range(actual_execs):
        exec_prefix = f"{alg_name}"
        if actual_execs > 1:
            exec_prefix += f" ({i+1}/{actual_execs})"
        
        spinner = Spinner(exec_prefix, console)
        executor = AlgorithmExecutor(timeout)
        
        warning_holder = []
        
        def warning_callback(msg):
            warning_holder.append(msg)
        
        # Criar uma função de progresso que funciona com o spinner
        progress_shown = threading.Event()
        last_progress_time = [0.0]  # Lista para permitir modificação dentro da closure, usando float
        
        def progress_callback(msg: str):
            # Primeira chamada de progresso - parar o spinner
            if not progress_shown.is_set():
                spinner.set_progress_override(True)
                progress_shown.set()
            
            # Mostrar progresso na mesma linha usando console
            progress_text = f"{msg:<40s}"
            if console:
                console.print_inline(f"{exec_prefix:<25s}... {progress_text}")
            else:
                print(f"\r{exec_prefix:<25s}... {progress_text}", end="", flush=True)
                
            # Atualizar timestamp da última mensagem de progresso
            last_progress_time[0] = time.time()

        # Thread para monitorar inatividade do progresso e reativar o spinner
        def monitor_progress_activity():
            while not spinner.stop_event.is_set():
                if progress_shown.is_set() and time.time() - last_progress_time[0] > 0.5:
                    # Se passou mais de 0.5 segundos sem mensagem de progresso, voltar ao spinner
                    spinner.set_progress_override(False)
                    progress_shown.clear()
                time.sleep(0.5)
                
        # Iniciar thread de monitoramento
        monitor_thread = threading.Thread(target=monitor_progress_activity)
        monitor_thread.daemon = True
        monitor_thread.start()

        t0 = time.time()
        
        # Iniciar spinner antes da execução
        spinner.start()
        
        try:
            # Criar instância do algoritmo
            instance = AlgClass(seqs, alphabet)
            
            # Executar com timeout em thread isolada
            center, val, info = executor.execute_with_timeout(
                instance, 
                progress_callback=progress_callback,
                warning_callback=warning_callback
            )
            
            tempo_execucao = time.time() - t0
            
            if 'erro' in info:
                if info['erro'] == 'timeout' or 'timeout' in info['erro'].lower():
                    # Parar spinner e mostrar timeout
                    spinner.stop()
                    if console:
                        console.print(f"{exec_prefix:<25s}... TIMEOUT ({timeout}s)")
                    else:
                        print(f"\r{exec_prefix:<25s}... TIMEOUT ({timeout}s)")
                    executions.append({
                        'tempo': tempo_execucao,
                        'iteracoes': 0,
                        'distancia': float('inf'),
                        'melhor_string': '',
                        'erro': f'Timeout ({timeout}s)'
                    })
                elif 'recurso' in info['erro'].lower() or 'resource' in info['erro'].lower():
                    spinner.stop()
                    if console:
                        console.print(f"{exec_prefix:<25s}... RECURSOS LIMITADOS")
                    else:
                        print(f"\r{exec_prefix:<25s}... RECURSOS LIMITADOS")
                    executions.append({
                        'tempo': tempo_execucao,
                        'iteracoes': 0,
                        'distancia': float('inf'),
                        'melhor_string': '',
                        'erro': 'Limite de recursos atingido'
                    })
                else:
                    error_msg = info['erro'][:50]
                    spinner.stop()
                    if console:
                        console.print(f"{exec_prefix:<25s}... ERRO: {error_msg}")
                    else:
                        print(f"\r{exec_prefix:<25s}... ERRO: {error_msg}")
                    executions.append({
                        'tempo': tempo_execucao,
                        'iteracoes': 0,
                        'distancia': float('inf'),
                        'melhor_string': '',
                        'erro': error_msg
                    })
            else:
                iteracoes = info.get('iteracoes', 1)
                
                logger.debug(f"[Runner] {alg_name} exec {i+1}: dist={val}")
                
                # Parar spinner e mostrar resultado final
                spinner.stop()
                if console:
                    console.print(f"{exec_prefix:<25s}... dist={val}, tempo={tempo_execucao:.3f}s")
                else:
                    print(f"\r{exec_prefix:<25s}... dist={val}, tempo={tempo_execucao:.3f}s")
                
                for warning_msg in warning_holder:
                    if console:
                        console.print(f"  AVISO: {warning_msg}")
                        
                executions.append({
                    'tempo': tempo_execucao,
                    'iteracoes': iteracoes,
                    'distancia': val,
                    'melhor_string': center
                })
                
        except (TimeoutException, ResourceLimitException) as e:
            tempo_execucao = time.time() - t0
            spinner.stop()
            
            if isinstance(e, ResourceLimitException):
                if console:
                    console.print(f"{exec_prefix:<25s}... RECURSOS LIMITADOS")
                else:
                    print(f"\r{exec_prefix:<25s}... RECURSOS LIMITADOS")
                executions.append({
                    'tempo': tempo_execucao,
                    'iteracoes': 0,
                    'distancia': float('inf'),
                    'melhor_string': '',
                    'erro': 'Limite de recursos atingido'
                })
            else:
                if console:
                    console.print(f"{exec_prefix:<25s}... TIMEOUT ({timeout}s)")
                else:
                    print(f"\r{exec_prefix:<25s}... TIMEOUT ({timeout}s)")
                executions.append({
                    'tempo': tempo_execucao,
                    'iteracoes': 0,
                    'distancia': float('inf'),
                    'melhor_string': '',
                    'erro': f'Timeout ({timeout}s)'
                })
            
        except KeyboardInterrupt:
            # Cancelar executor e parar spinner
            executor.cancel()
            spinner.stop()
            if console:
                console.print("\nExecução interrompida pelo usuário. Encerrando.")
            else:
                print(f"\nExecução interrompida pelo usuário. Encerrando.")
            sys.exit(0)
            
        except Exception as e:
            tempo_execucao = time.time() - t0
            spinner.stop()
            error_msg = str(e)[:50]
            if console:
                console.print(f"{exec_prefix:<25s}... ERRO: {error_msg}")
            else:
                print(f"\r{exec_prefix:<25s}... ERRO: {error_msg}")
            executions.append({
                'tempo': tempo_execucao,
                'iteracoes': 0,
                'distancia': float('inf'),
                'melhor_string': '',
                'erro': error_msg
            })
        
        finally:
            # Garantir que spinner pare em qualquer caso
            if spinner.thread and spinner.thread.is_alive():
                spinner.stop()
            # Força garbage collection após cada execução
            force_garbage_collection()
    
    logger.debug(f"[Runner] {alg_name} finalizado: {len(executions)} execuções")
    return executions