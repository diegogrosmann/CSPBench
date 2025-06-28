import threading
import time
import itertools

def _gap(baseline: int, val):
    return 100 * (baseline - val) / baseline if (baseline and val is not None) else 0.0

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
        if self.console:
            self.console.set_spinner(self.prefix, self.spinner)
        self.thread = threading.Thread(target=self._animate)
        self.thread.daemon = True
        self.thread.start()

    def _animate(self):
        if self.console:
            self.console.start_spinner()
        while not self.stop_event.is_set():
            # Só mostrar spinner se não houver progresso específico
            if not self.progress_override:
                print(f"\r{self.prefix:<25s}{next(self.spinner)}", end="", flush=True)
            time.sleep(0.3)

    def set_progress_override(self, value: bool):
        self.progress_override = value

    def stop(self):
        self.stop_event.set()
        if self.thread and self.thread.is_alive():
            self.thread.join(timeout=1.0)
        if self.console:
            self.console.stop_spinner()

def execute_algorithm_runs(alg_name, AlgClass, seqs, alphabet, num_execs, baseline_val, console=None):
    is_deterministic = getattr(AlgClass, 'is_deterministic', False)
    actual_execs = 1 if is_deterministic else num_execs
    executions = []
    for i in range(actual_execs):
        exec_prefix = f"{alg_name}"
        if actual_execs > 1:
            exec_prefix += f" ({i+1}/{actual_execs})"
        spinner = Spinner(exec_prefix, console)
        result_holder = {}
        exc_holder = {}
        warning_holder = []

        def run_algorithm():
            try:
                def warning_callback(msg):
                    warning_holder.append(msg)
                instance = AlgClass(seqs, alphabet)
                if hasattr(instance, 'set_warning_callback'):
                    instance.set_warning_callback(warning_callback)
                if hasattr(instance, 'set_progress_callback'):
                    def progress_callback(msg: str):
                        # Parar spinner e mostrar progresso na mesma linha
                        spinner.set_progress_override(True)
                        progress_text = f"{msg:<40s}"
                        print(f"\r{exec_prefix:<25s}... {progress_text}", end="", flush=True)
                    instance.set_progress_callback(progress_callback)
                center, val = instance.run()
                result_holder['center'] = center
                result_holder['val'] = val
                if hasattr(instance, 'geracao'):
                    result_holder['iteracoes'] = instance.geracao
                elif hasattr(instance, 'iterations'):
                    result_holder['iteracoes'] = instance.iterations
                elif hasattr(instance, 'num_iteracoes'):
                    result_holder['iteracoes'] = instance.num_iteracoes
                else:
                    result_holder['iteracoes'] = 1
            except Exception as e:
                exc_holder['exc'] = e

        t0 = time.time()
        algo_thread = threading.Thread(target=run_algorithm)
        spinner.start()
        algo_thread.start()
        while algo_thread.is_alive():
            time.sleep(0.1)
        spinner.stop()
        tempo_execucao = time.time() - t0

        if 'exc' in exc_holder:
            error_msg = str(exc_holder['exc'])[:50]
            print(f"\r{exec_prefix:<25s}... ERRO: {error_msg}")
            executions.append({
                'tempo': tempo_execucao,
                'iteracoes': 0,
                'distancia': float('inf'),
                'melhor_string': '',
                'erro': error_msg
            })
        else:
            center = result_holder['center']
            val = result_holder['val']
            iteracoes = result_holder.get('iteracoes', 1)
            gap = _gap(baseline_val, val)
            # Substituir o conteúdo da linha atual com o resultado final
            print(f"\r{exec_prefix:<25s}... dist={val}, tempo={tempo_execucao:.3f}s, gap={gap:.1f}%")
            for warning_msg in warning_holder:
                if console:
                    console.print(f"  AVISO: {warning_msg}")
            executions.append({
                'tempo': tempo_execucao,
                'iteracoes': iteracoes,
                'distancia': val,
                'melhor_string': center,
                'gap': gap
            })
    return executions