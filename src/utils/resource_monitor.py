"""
Sistema de monitoramento de recursos para evitar travamento do sistema.

Classes:
    ResourceLimits: Limites de recursos do sistema.
    ResourceMonitor: Monitor de recursos com limites automáticos.

Funções:
    get_safe_memory_limit(): Retorna limite seguro de uso de memória.
    check_algorithm_feasibility(...): Verifica viabilidade de execução.
    force_garbage_collection(): Força coleta de lixo.
"""

import gc
import logging
import os
import threading
import time
from collections.abc import Callable
from dataclasses import dataclass

from .resource_limits_config import (
    get_merged_gc_config,
    get_merged_memory_safety,
    get_merged_resource_limits,
)

logger = logging.getLogger(__name__)


@dataclass
class ResourceLimits:
    """Limites de recursos do sistema."""

    max_memory_mb: int
    max_iterations: int
    check_interval: float
    gc_threshold: int

    # Configurações de segurança de memória
    memory_usage_ratio: float
    max_safe_memory_mb: int
    min_free_memory_mb: int
    memory_check_aggressive: bool

    # Configurações de garbage collection
    gc_auto_collect: bool
    gc_frequency: int
    gc_force_on_limit: bool
    gc_threshold_ratio: float

    @classmethod
    def from_config(cls) -> "ResourceLimits":
        """Cria ResourceLimits a partir das configurações."""
        limits = get_merged_resource_limits()
        memory_safety = get_merged_memory_safety()
        gc_config = get_merged_gc_config()

        return cls(
            # Limites básicos
            max_memory_mb=limits["max_memory_mb"],
            max_iterations=limits["max_iterations"],
            check_interval=limits["check_interval"],
            gc_threshold=limits["gc_threshold"],
            # Segurança de memória
            memory_usage_ratio=memory_safety["memory_usage_ratio"],
            max_safe_memory_mb=memory_safety["max_safe_memory_mb"],
            min_free_memory_mb=memory_safety["min_free_memory_mb"],
            memory_check_aggressive=memory_safety["memory_check_aggressive"],
            # Garbage collection
            gc_auto_collect=gc_config["gc_auto_collect"],
            gc_frequency=gc_config["gc_frequency"],
            gc_force_on_limit=gc_config["gc_force_on_limit"],
            gc_threshold_ratio=gc_config["gc_threshold_ratio"],
        )


class ResourceMonitor:
    """Monitor de recursos do sistema com limites automáticos."""

    def __init__(self, limits: ResourceLimits | None = None):
        self.limits = limits or ResourceLimits.from_config()
        self.monitoring = False
        self.monitor_thread: threading.Thread | None = None
        self.stop_event = threading.Event()
        self.violation_callback: Callable[[str], None] | None = None
        self.gc_counter = 0
        self.last_gc_time = 0.0
        self.memory_readings = []  # Histórico de leituras de memória

    def set_violation_callback(self, callback: Callable[[str], None]) -> None:
        """Define callback para quando limites são violados."""
        self.violation_callback = callback

    def start_monitoring(self) -> None:
        """Inicia monitoramento contínuo de recursos."""
        if self.monitoring:
            return

        self.monitoring = True
        self.stop_event.clear()
        self.monitor_thread = threading.Thread(target=self._monitor_loop, daemon=True)
        self.monitor_thread.start()
        logger.info("Monitoramento de recursos iniciado")

    def stop_monitoring(self) -> None:
        """Para monitoramento de recursos."""
        if not self.monitoring:
            return

        self.monitoring = False
        self.stop_event.set()

        if self.monitor_thread and self.monitor_thread.is_alive():
            self.monitor_thread.join(timeout=5.0)

        logger.info("Monitoramento de recursos parado")

    def get_process_memory_mb(self) -> float:
        """
        Retorna uso de memória do processo atual em MB.
        T10-2: Coleta GC em cada verificação de memória.
        """
        try:
            # T10-2: Força garbage collection para obter medição mais precisa
            if self.limits.gc_auto_collect:
                gc.collect()
                self.last_gc_time = time.time()

            # No Linux/Unix, tenta ler de /proc/self/status
            if os.path.exists("/proc/self/status"):
                with open("/proc/self/status") as f:
                    for line in f:
                        if line.startswith("VmRSS:"):
                            # VmRSS está em kB
                            kb = int(line.split()[1])
                            memory_mb = kb / 1024.0  # Converter para MB

                            # Manter histórico de leituras
                            self.memory_readings.append(memory_mb)
                            if len(self.memory_readings) > 10:
                                self.memory_readings.pop(0)

                            return memory_mb

            # Fallback: usar sys.getsizeof em objetos principais
            # Estimativa baseada no tamanho de objetos Python
            import sys

            total_size = 0
            for obj in gc.get_objects():
                try:
                    total_size += sys.getsizeof(obj)
                except (TypeError, AttributeError):
                    pass

            memory_mb = total_size / (1024 * 1024)  # Converter para MB

            # Manter histórico de leituras
            self.memory_readings.append(memory_mb)
            if len(self.memory_readings) > 10:
                self.memory_readings.pop(0)

            return memory_mb

        except Exception as e:
            logger.error(f"Erro ao obter uso de memória: {e}")
            return 0.0

    def check_limits(self) -> tuple[bool, list[str]]:
        """
        Verifica se os limites estão sendo respeitados.
        T10-2: Melhoria na verificação de limites com GC inteligente.
        """
        violations = []

        # Verificar memória do processo
        memory_mb = self.get_process_memory_mb()

        # Verificação básica de limite
        if memory_mb > self.limits.max_memory_mb:
            violations.append(f"Memória: {memory_mb:.1f}MB > {self.limits.max_memory_mb}MB")

        # Verificação de segurança de memória disponível
        if self.limits.memory_check_aggressive:
            available_mb = get_available_memory_mb()
            memory_usage_ratio = memory_mb / available_mb if available_mb > 0 else 1.0

            if memory_usage_ratio > self.limits.memory_usage_ratio:
                violations.append(f"Uso de memória: {memory_usage_ratio:.1%} > {self.limits.memory_usage_ratio:.1%}")

            # Verificar se há memória livre mínima
            free_memory = available_mb - memory_mb
            if free_memory < self.limits.min_free_memory_mb:
                violations.append(f"Memória livre baixa: {free_memory:.1f}MB < {self.limits.min_free_memory_mb}MB")

        # T10-2: GC preventivo se aproximando do limite
        if self.limits.gc_force_on_limit:
            threshold_memory = self.limits.max_memory_mb * self.limits.gc_threshold_ratio
            if memory_mb > threshold_memory:
                logger.info(f"Executando GC preventivo: {memory_mb:.1f}MB > {threshold_memory:.1f}MB")
                force_garbage_collection()

        return len(violations) == 0, violations

    def _monitor_loop(self) -> None:
        """
        Loop principal de monitoramento.
        T10-3: Lógica melhorada de GC com reset automático do contador.
        """
        while not self.stop_event.wait(self.limits.check_interval):
            try:
                # T10-3: Garbage collection periódico com lógica melhorada
                self.gc_counter += 1

                should_gc = False

                # Verificar se deve fazer GC por frequência
                if self.limits.gc_auto_collect and self.gc_counter >= self.limits.gc_frequency:
                    should_gc = True
                    logger.debug(f"GC por frequência: {self.gc_counter} >= {self.limits.gc_frequency}")

                # Verificar se deve fazer GC por tempo (a cada 30 segundos)
                current_time = time.time()
                if current_time - self.last_gc_time > 30.0:
                    should_gc = True
                    logger.debug("GC por tempo: > 30 segundos desde último GC")

                # Executar GC se necessário
                if should_gc:
                    force_garbage_collection()
                    self.gc_counter = 0  # T10-3: Reset do contador após GC
                    self.last_gc_time = current_time

                # Verificar limites
                is_safe, violations = self.check_limits()

                if not is_safe and self.violation_callback:
                    violation_msg = "; ".join(violations)
                    logger.warning(f"Limites de recursos violados: {violation_msg}")
                    self.violation_callback(violation_msg)

            except Exception as e:
                logger.error(f"Erro no loop de monitoramento: {e}")
                time.sleep(1.0)

    def get_memory_stats(self) -> dict:
        """
        Retorna estatísticas detalhadas de memória.
        """
        current_memory = self.get_process_memory_mb()
        available_memory = get_available_memory_mb()

        stats = {
            "current_memory_mb": current_memory,
            "available_memory_mb": available_memory,
            "memory_limit_mb": self.limits.max_memory_mb,
            "memory_usage_ratio": current_memory / available_memory if available_memory > 0 else 0,
            "gc_counter": self.gc_counter,
            "last_gc_time": self.last_gc_time,
            "memory_readings_count": len(self.memory_readings),
        }

        if self.memory_readings:
            stats.update(
                {
                    "memory_avg_mb": sum(self.memory_readings) / len(self.memory_readings),
                    "memory_max_mb": max(self.memory_readings),
                    "memory_min_mb": min(self.memory_readings),
                }
            )

        return stats


def get_available_memory_mb() -> float:
    """Estima memória disponível em MB."""
    try:
        # No Linux, tenta ler /proc/meminfo
        if os.path.exists("/proc/meminfo"):
            meminfo = {}
            with open("/proc/meminfo") as f:
                for line in f:
                    parts = line.split()
                    if len(parts) >= 2:
                        key = parts[0].rstrip(":")
                        value = int(parts[1])  # em kB
                        meminfo[key] = value

            # Calcular memória disponível
            if "MemAvailable" in meminfo:
                return meminfo["MemAvailable"] / 1024.0  # Converter para MB
            elif "MemFree" in meminfo and "Buffers" in meminfo and "Cached" in meminfo:
                available = meminfo["MemFree"] + meminfo["Buffers"] + meminfo["Cached"]
                return available / 1024.0

        # Fallback conservador: assume 2GB disponível
        return 2048.0

    except Exception as e:
        logger.error(f"Erro ao obter memória disponível: {e}")
        return 1024.0  # Fallback muito conservador


def get_safe_memory_limit() -> float:
    """Calcula limite seguro de memória baseado no sistema."""
    try:
        available_mb = get_available_memory_mb()
        # Usar no máximo 50% da memória disponível
        safe_limit = available_mb * 0.5
        return min(safe_limit, 2048.0)  # Máximo 2GB
    except (OSError, AttributeError):
        return 512.0  # Fallback conservador


def estimate_algorithm_memory(n: int, L: int, algorithm_name: str) -> float:
    """Estima uso de memória de um algoritmo em MB."""
    # Dados básicos: n strings de tamanho L
    base_data_mb = (n * L * 4) / (1024 * 1024)  # Assumindo 4 bytes por char

    # Estimativas por algoritmo (multiplicadores mais conservadores)
    multipliers = {
        "DP-CSP": min(1000, max(50, n**2)),  # Limitado mas ainda exponencial
        "BLF-GA": 20,  # População + estruturas auxiliares
        "H³-CSP": 10,  # Blocos + candidatos
        "CSC": 15,  # Clusters + matrizes
        "Baseline": 2,  # Apenas contadores
    }

    multiplier = multipliers.get(algorithm_name, 10)
    estimated = base_data_mb * multiplier

    logger.debug(f"Estimativa de memória para {algorithm_name}: {estimated:.2f}MB")
    return estimated


def check_algorithm_feasibility(n: int, L: int, algorithm_name: str) -> tuple[bool, str]:
    """
    Verifica se um algoritmo é viável para executar com os parâmetros dados.

    Returns:
        Tuple (is_viable, message)
    """
    if algorithm_name == "DP-CSP":
        # Estimativa mais realista para DP-CSP
        # Com d=baseline_distance, estados máximos = (d+1)^n

        # Estimar d baseado em n e L (heurística)
        estimated_d = min(15, max(3, int(L * 0.15)))  # Entre 3 e 15, ~15% do comprimento

        # Estados estimados
        estimated_states = (estimated_d + 1) ** n

        # Limites mais permissivos
        if n > 20:
            return False, f"Muitas strings (n={n} > 20)"

        if L > 200:
            return False, f"Strings muito longas (L={L} > 200)"

        if estimated_states > 10**9:  # 1 bilhão de estados
            return (
                False,
                f"Estados estimados muito altos: {estimated_states:.0e} (limite: 10^9)",
            )

        if estimated_states > 10**7:  # 10 milhões - dar aviso mas permitir
            return (
                True,
                f"Complexo mas viável: ~{estimated_states:.0e} estados estimados",
            )

        return (
            True,
            f"Viável: ~{estimated_states:.0e} estados estimados (d~{estimated_d})",
        )

    elif algorithm_name == "BLF-GA":
        # BLF-GA é geralmente escalável
        if n > 1000 or L > 5000:
            return False, f"Dataset muito grande: n={n}, L={L}"
        return True, "Escalável via metaheurística"

    elif algorithm_name == "H³-CSP":
        # H³-CSP com decomposição em blocos
        if n > 500 or L > 2000:
            return False, f"Dataset muito grande: n={n}, L={L}"
        return True, "Escalável via decomposição hierárquica"

    elif algorithm_name == "CSC":
        # CSC usando clustering
        if n > 200 or L > 1000:
            return False, f"Dataset muito grande: n={n}, L={L}"
        return True, "Escalável via clustering"

    elif algorithm_name == "Baseline":
        # Baseline é sempre viável (complexidade linear)
        return True, "Sempre viável (consenso ganancioso)"

    else:
        # Algoritmos desconhecidos - assumir viáveis
        return True, "Viável (algoritmo desconhecido)"


def force_garbage_collection():
    """
    Força garbage collection para liberar memória.
    T10-3: Versão melhorada com estatísticas e múltiplas passadas.
    """
    try:
        # Obter estatísticas antes da coleta
        objects_before = len(gc.get_objects())

        # Executar múltiplas passadas de GC para limpeza completa
        collected_total = 0
        for generation in range(3):  # GC tem 3 gerações (0, 1, 2)
            collected = gc.collect(generation)
            collected_total += collected

        # Estatísticas após coleta
        objects_after = len(gc.get_objects())
        objects_freed = objects_before - objects_after

        logger.debug(
            f"GC executado: {collected_total} objetos coletados, "
            f"{objects_freed} objetos liberados ({objects_before} → {objects_after})"
        )

        # Forçar limpeza adicional se muitos objetos foram coletados
        if collected_total > 1000:
            logger.debug("Executando GC adicional devido ao alto número de objetos coletados")
            gc.collect()

    except Exception as e:
        logger.error(f"Erro no garbage collection: {e}")
