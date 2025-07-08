"""
Módulo para execução em lote de algoritmos CSP.

Este módulo gerencia a execução de múltiplas configurações de dataset
usando a mesma estrutura da execução individual.
"""

import logging
import os
from typing import Any, Dict, List, Tuple

import yaml

from src.datasets.dataset_entrez import fetch_dataset
from src.datasets.dataset_file import load_dataset
from src.datasets.dataset_synthetic import generate_dataset as generate_synthetic
from src.datasets.dataset_synthetic import generate_dataset_from_params

logger = logging.getLogger(__name__)


def load_batch_config(config_path: str) -> Dict[str, Any]:
    """
    Carrega configuração de batch a partir de arquivo YAML.

    Args:
        config_path: Caminho para o arquivo de configuração

    Returns:
        Dicionário com a configuração carregada

    Raises:
        FileNotFoundError: Se o arquivo não for encontrado
        yaml.YAMLError: Se o arquivo não for um YAML válido
    """
    if not os.path.exists(config_path):
        raise FileNotFoundError(
            f"Arquivo de configuração não encontrado: {config_path}"
        )

    with open(config_path, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    # Validar estrutura básica
    if "execucoes" not in config:
        raise ValueError("Configuração deve conter seção 'execucoes'")

    return config


def generate_dataset_from_config(
    dataset_config: Dict[str, Any], silent: bool = True
) -> Tuple[List[str], Dict[str, Any]]:
    """
    Gera dataset baseado na configuração fornecida.

    Args:
        dataset_config: Configuração do dataset
        silent: Se True, não exibe mensagens

    Returns:
        Tupla contendo (sequências, parâmetros)

    Raises:
        ValueError: Se o tipo de dataset não for suportado
    """
    dataset_type = dataset_config.get("tipo", "synthetic")
    params = dataset_config.get("parametros", {})

    if dataset_type == "synthetic":
        # Usar nova função que não requer interação do usuário
        n = params.get("n", 20)
        L = params.get("L", 100)
        alphabet = params.get("alphabet", "ACGT")
        noise = params.get("noise", 0.1)
        fully_random = params.get("fully_random", False)
        seed = params.get("seed", None)

        seqs, generated_params = generate_dataset_from_params(
            n=n,
            L=L,
            alphabet=alphabet,
            noise=noise,
            fully_random=fully_random,
            seed=seed,
        )

        # Adicionar tipo de dataset aos parâmetros
        batch_params = {"dataset_source": "1"}  # 1 = sintético
        batch_params.update(generated_params)

        return seqs, batch_params

    elif dataset_type == "file":
        # Para arquivo, usar parâmetros específicos
        if "filename" in params:
            os.environ["DATASET_FILE"] = params["filename"]

        seqs, file_params = load_dataset(silent=silent)

        batch_params = {"dataset_source": "2"}  # 2 = arquivo
        batch_params.update(file_params)

        return seqs, batch_params

    elif dataset_type == "entrez":
        # Para Entrez, configurar parâmetros de busca
        if "query" in params:
            os.environ["ENTREZ_QUERY"] = params["query"]
        if "db" in params:
            os.environ["ENTREZ_DB"] = params["db"]
        if "retmax" in params:
            os.environ["ENTREZ_RETMAX"] = str(params["retmax"])

        seqs, entrez_params = fetch_dataset()

        batch_params = {"dataset_source": "3"}  # 3 = entrez
        batch_params.update(entrez_params)

        return seqs, batch_params

    else:
        raise ValueError(f"Tipo de dataset não suportado: {dataset_type}")


def execute_algorithms_core(
    algorithms: List[str],
    seqs: List[str],
    alphabet: str,
    num_execs: int,
    timeout: int,
    dataset_params: Dict[str, Any],
    silent: bool = False,
    use_curses: bool = False,
) -> Dict[str, Any]:
    """
    Função central para execução de algoritmos.

    Args:
        algorithms: Lista de algoritmos a executar
        seqs: Sequências do dataset
        alphabet: Alfabeto das sequências
        num_execs: Número de execuções por algoritmo
        timeout: Timeout em segundos
        dataset_params: Parâmetros do dataset
        silent: Se True, executa em modo silencioso
        use_curses: Se True, usa interface curses

    Returns:
        Dicionário com resultados da execução
    """
    if use_curses:
        from src.ui.curses_integration import CursesExecutionMonitor

        monitor = CursesExecutionMonitor(max_workers=4, timeout=timeout)

        results = monitor.execute_algorithms(
            algorithm_names=algorithms,
            seqs=seqs,
            alphabet=alphabet,
            num_execs=num_execs,
            dataset_params=dataset_params,
        )

        return results
    else:
        # Execução tradicional (fallback)
        import time

        from algorithms.base import global_registry
        from src.core.interfaces import TaskStatus, create_executor

        results = {}
        executor = create_executor(timeout_seconds=timeout, max_workers=4)

        try:
            for alg_name in algorithms:
                if alg_name not in global_registry:
                    if not silent:
                        print(f"❌ Algoritmo {alg_name} não encontrado!")
                    continue

                AlgClass = global_registry[alg_name]

                # Verificar se o algoritmo é determinístico
                is_deterministic = getattr(AlgClass, "is_deterministic", False)
                actual_num_execs = 1 if is_deterministic else num_execs

                if not silent:
                    if is_deterministic:
                        print(
                            f"  🔒 {alg_name} é determinístico - executando apenas 1 vez"
                        )
                    else:
                        print(
                            f"  🎲 {alg_name} é não-determinístico - executando {actual_num_execs} vezes"
                        )

                alg_results = []

                for i in range(actual_num_execs):
                    if not silent:
                        if actual_num_execs == 1:
                            print(f"  Executando {alg_name}")
                        else:
                            print(
                                f"  Executando {alg_name} - Run {i+1}/{actual_num_execs}"
                            )

                    instance = AlgClass(seqs, alphabet)
                    handle = executor.submit(instance)

                    # Aguardar conclusão
                    while executor.poll(handle) == TaskStatus.RUNNING:
                        time.sleep(0.1)

                    result = executor.result(handle)
                    alg_results.append(result)

                results[alg_name] = alg_results

        finally:
            if hasattr(executor, "shutdown"):
                executor.shutdown(wait=True)

        return results


def execute_batch_config(
    config_path: str, use_curses: bool = True, silent: bool = False
) -> Dict[str, Any]:
    """
    Executa configuração em lote.

    Args:
        config_path: Caminho para arquivo de configuração
        use_curses: Se True, usa interface curses
        silent: Se True, executa em modo silencioso

    Returns:
        Dicionário com resultados da execução
    """
    # Carregar configuração
    config = load_batch_config(config_path)
    batch_info = config.get("batch_info", {})
    execucoes = config["execucoes"]

    if not silent:
        print(f"📋 Executando batch: {batch_info.get('nome', 'Sem nome')}")
        print(f"📄 Descrição: {batch_info.get('descricao', 'N/A')}")
        print(f"🔧 Total de execuções: {len(execucoes)}")

    all_results = {}

    # Executar cada configuração
    for exec_idx, exec_config in enumerate(execucoes, 1):
        exec_name = exec_config.get("nome", f"Execução {exec_idx}")

        if not silent:
            print(f"\n{'='*60}")
            print(f"🚀 Execução {exec_idx}/{len(execucoes)}: {exec_name}")
            print(f"{'='*60}")

        try:
            # Configurar execução
            algorithms = exec_config.get("algoritmos", ["Baseline"])
            num_execs = exec_config.get("runs_per_algorithm_per_base", 1)
            num_bases = exec_config.get("num_bases", 1)
            timeout = exec_config.get("timeout", 300)

            if not silent:
                print(f"🧮 Algoritmos: {algorithms}")
                print(f"🔁 Execuções por algoritmo por base: {num_execs}")
                print(f"�️ Número de bases: {num_bases}")
                print(f"⏰ Timeout: {timeout}s")

            # Resultados desta configuração
            config_results = {}

            # Executar para cada base
            for base_idx in range(num_bases):
                if not silent:
                    print(f"\n📊 Gerando Base {base_idx + 1}/{num_bases}")

                # Gerar dataset para esta base (nova a cada iteração)
                dataset_config = exec_config["dataset"]
                seqs, dataset_params = generate_dataset_from_config(
                    dataset_config, silent=silent
                )

                if not silent:
                    print(f"📊 Dataset gerado: n={len(seqs)}, L={len(seqs[0])}")
                    print(f"� Tipo: {dataset_config.get('tipo', 'synthetic')}")

                # Extrair alfabeto
                alphabet = "".join(sorted(set("".join(seqs))))

                # Preparar parâmetros para esta base
                batch_dataset_params = dataset_params.copy()
                batch_dataset_params.update(
                    {
                        "batch_name": batch_info.get("nome", "Batch"),
                        "execution_name": exec_name,
                        "execution_index": exec_idx,
                        "total_executions": len(execucoes),
                        "base_index": base_idx + 1,
                        "total_bases": num_bases,
                    }
                )

                # Executar algoritmos para esta base
                base_results = execute_algorithms_core(
                    algorithms=algorithms,
                    seqs=seqs,
                    alphabet=alphabet,
                    num_execs=num_execs,
                    timeout=timeout,
                    dataset_params=batch_dataset_params,
                    silent=silent,
                    use_curses=use_curses,
                )

                # Armazenar resultados desta base
                config_results[f"base_{base_idx + 1}"] = {
                    "dataset_params": dataset_params,
                    "seqs": seqs,
                    "alphabet": alphabet,
                    "results": base_results,
                }

                if not silent:
                    print(f"✅ Base {base_idx + 1}/{num_bases} concluída!")

            # Armazenar resultados desta execução
            all_results[exec_name] = {
                "config": exec_config,
                "algorithms": algorithms,
                "num_bases": num_bases,
                "bases_results": config_results,
            }

            if not silent:
                print(f"✅ Execução {exec_idx} concluída com sucesso!")

        except Exception as e:
            logger.exception(f"Erro na execução {exec_idx}: {e}")
            if not silent:
                print(f"❌ Erro na execução {exec_idx}: {e}")

            all_results[exec_name] = {"error": str(e), "config": exec_config}

    if not silent:
        print(f"\n{'='*60}")
        print(f"🎉 Batch concluído!")
        print(
            f"✅ Execuções bem-sucedidas: {len([r for r in all_results.values() if 'error' not in r])}"
        )
        print(
            f"❌ Execuções com erro: {len([r for r in all_results.values() if 'error' in r])}"
        )
        print(f"{'='*60}")

    return all_results


def execute_synthetic_dataset(
    algorithms: List[str],
    num_execs: int,
    timeout: int,
    dataset_params: Dict[str, Any] | None = None,
    silent: bool = False,
    use_curses: bool = False,
) -> Dict[str, Any]:
    """
    Executa algoritmos em dataset sintético.

    Args:
        algorithms: Lista de algoritmos a executar
        num_execs: Número de execuções por algoritmo
        timeout: Timeout em segundos
        dataset_params: Parâmetros do dataset (se None, pergunta interativamente)
        silent: Se True, executa em modo silencioso
        use_curses: Se True, usa interface curses

    Returns:
        Dicionário com resultados da execução
    """
    if dataset_params is None:
        # Execução interativa - gera dataset perguntando ao usuário
        from src.datasets.dataset_synthetic import generate_dataset

        seqs, params = generate_dataset(silent=silent)
        params = {"dataset_source": "1", **params}
    else:
        # Execução batch - usa parâmetros fornecidos
        seqs, params = generate_dataset_from_params(
            n=dataset_params.get("n", 20),
            L=dataset_params.get("L", 100),
            alphabet=dataset_params.get("alphabet", "ACGT"),
            noise=dataset_params.get("noise", 0.1),
            fully_random=dataset_params.get("fully_random", False),
            seed=dataset_params.get("seed", None),
        )
        params = {"dataset_source": "1", **params}

    # Extrair alfabeto
    alphabet = "".join(sorted(set("".join(seqs))))

    # Executar algoritmos
    results = execute_algorithms_core(
        algorithms=algorithms,
        seqs=seqs,
        alphabet=alphabet,
        num_execs=num_execs,
        timeout=timeout,
        dataset_params=params,
        silent=silent,
        use_curses=use_curses,
    )

    return {
        "dataset_params": params,
        "results": results,
        "seqs": seqs,
        "alphabet": alphabet,
    }


if __name__ == "__main__":
    # Exemplo de uso
    import sys

    if len(sys.argv) > 1:
        config_path = sys.argv[1]
    else:
        config_path = "batch_configs/exemplo.yaml"

    try:
        results = execute_batch_config(config_path, use_curses=True, silent=False)
        print(f"\n📊 Resultados finais: {len(results)} execuções processadas")
    except Exception as e:
        print(f"❌ Erro: {e}")
        sys.exit(1)
