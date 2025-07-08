"""
Executor de otimização em lote baseado em arquivos YAML.

Este módulo permite executar otimização de hiperparâmetros usando
configurações definidas em arquivos YAML, similar aos batch_configs
mas focado em otimização.
"""

import json
import os
from datetime import datetime
from typing import Any, Dict, List

import yaml

from src.datasets.dataset_entrez import fetch_dataset_silent
from src.datasets.dataset_file import load_dataset_with_params
from src.datasets.dataset_synthetic import generate_dataset_from_params
from src.optimization.optuna_optimizer import optimize_algorithm
from src.ui.cli.console_manager import console


def load_optimization_config(config_file: str) -> Dict[str, Any]:
    """
    Carrega configuração de otimização de arquivo YAML.

    Args:
        config_file: Caminho para arquivo de configuração

    Returns:
        Dict com configuração carregada
    """
    if not os.path.exists(config_file):
        raise FileNotFoundError(
            f"Arquivo de configuração não encontrado: {config_file}"
        )

    with open(config_file, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    return config


def generate_dataset_from_config(dataset_config: Dict[str, Any]) -> tuple:
    """
    Gera dataset baseado na configuração.

    Args:
        dataset_config: Configuração do dataset

    Returns:
        tuple: (sequences, alphabet, dataset_info)
    """
    dataset_type = dataset_config.get("tipo", "synthetic")

    if dataset_type == "synthetic":
        params = dataset_config.get("parametros", {})
        from src.datasets.dataset_synthetic import generate_dataset_from_params

        seqs, dataset_params = generate_dataset_from_params(
            n=params.get("n", 20),
            L=params.get("L", 100),
            alphabet=params.get("alphabet", "ACGT"),
            noise=params.get("noise", 0.1),
            fully_random=params.get("fully_random", False),
            seed=params.get("seed", None),
        )
        alphabet = params.get("alphabet", "ACGT")
        dataset_info = {
            "type": "synthetic",
            "name": dataset_config.get("nome", "Sintético"),
            **params,
            **dataset_params,
        }

    elif dataset_type == "file":
        file_path = dataset_config.get("parametros", {}).get("file_path")
        if not file_path:
            raise ValueError("file_path é obrigatório para dataset tipo 'file'")
        from src.datasets.dataset_file import load_dataset_with_params

        seqs, dataset_params = load_dataset_with_params({"filepath": file_path})
        alphabet = "".join(sorted(set("".join(seqs))))
        dataset_info = {
            "type": "file",
            "name": dataset_config.get("nome", "Arquivo"),
            **dataset_params,
        }

    elif dataset_type == "entrez":
        entrez_params = dataset_config.get("parametros", {})
        from src.datasets.dataset_entrez import fetch_dataset_silent

        seqs, dataset_params = fetch_dataset_silent(entrez_params)
        alphabet = "".join(sorted(set("".join(seqs))))
        dataset_info = {
            "type": "entrez",
            "name": dataset_config.get("nome", "NCBI"),
            **dataset_params,
        }

    else:
        raise ValueError(f"Tipo de dataset não suportado: {dataset_type}")

    return seqs, alphabet, dataset_info


def execute_optimization_batch(config_file: str) -> Dict[str, Any]:
    """
    Executa otimização em lote baseada em arquivo de configuração.

    Args:
        config_file: Caminho para arquivo de configuração YAML

    Returns:
        Dict com resultados da otimização
    """
    console.print(f"📋 Carregando configuração: {config_file}")

    # Carregar configuração
    config = load_optimization_config(config_file)

    batch_info = config.get("batch_info", {})
    opt_config = config.get("optimization_config", {})
    algorithms = config.get("algorithms", ["BLF-GA"])
    datasets = config.get("datasets", [])
    algorithm_params = config.get("algorithm_params", {})
    output_config = config.get("output", {})

    console.print(f"🎯 Batch: {batch_info.get('nome', 'Sem nome')}")
    console.print(f"📝 Descrição: {batch_info.get('descricao', 'Sem descrição')}")
    console.print(f"🔬 Algoritmos: {', '.join(algorithms)}")
    console.print(f"📊 Datasets: {len(datasets)} configurações")

    # Inicializar resultados
    batch_results = {
        "batch_info": {
            "config_file": config_file,
            "timestamp": datetime.now().isoformat(),
            **batch_info,
        },
        "optimization_config": opt_config,
        "results": [],
    }

    # Executar para cada dataset
    for dataset_idx, dataset_config in enumerate(datasets):
        console.print(
            f"\n📊 Dataset {dataset_idx + 1}/{len(datasets)}: {dataset_config.get('nome', 'Sem nome')}"
        )

        try:
            # Gerar dataset
            seqs, alphabet, dataset_info = generate_dataset_from_config(dataset_config)
            console.print(
                f"✅ Dataset carregado: {len(seqs)} sequências de tamanho {len(seqs[0])}"
            )

            dataset_results = {
                "dataset_index": dataset_idx + 1,
                "dataset_info": dataset_info,
                "algorithms": {},
            }

            # Executar para cada algoritmo
            for alg_name in algorithms:
                console.print(f"\n🔬 Otimizando {alg_name}...")

                try:
                    # Configurar parâmetros específicos do algoritmo
                    alg_specific_config = algorithm_params.get(alg_name, {})

                    # Executar otimização
                    result = optimize_algorithm(
                        algorithm_name=alg_name,
                        sequences=seqs,
                        alphabet=alphabet,
                        n_trials=opt_config.get("n_trials", 50),
                        timeout_per_trial=opt_config.get("timeout_per_trial", 60),
                        direction=opt_config.get("direction", "minimize"),
                        sampler=opt_config.get("sampler", "TPE"),
                        pruner=opt_config.get("pruner", "Median"),
                        show_progress=True,
                    )

                    dataset_results["algorithms"][alg_name] = {
                        "success": True,
                        "best_value": result.best_value,
                        "best_params": result.best_params,
                        "n_trials": result.n_trials,
                        "optimization_time": result.optimization_time,
                        "study_name": result.study_name,
                    }

                    console.print(
                        f"✅ {alg_name}: Melhor valor = {result.best_value:.6f}"
                    )

                    # Salvar gráficos se solicitado
                    if opt_config.get("save_plots", False):
                        save_optimization_plots(
                            result, alg_name, dataset_idx, output_config
                        )

                except Exception as e:
                    console.print(f"❌ Erro em {alg_name}: {e}")
                    dataset_results["algorithms"][alg_name] = {
                        "success": False,
                        "error": str(e),
                        "best_value": float("inf"),
                    }

            batch_results["results"].append(dataset_results)

        except Exception as e:
            console.print(f"❌ Erro no dataset {dataset_idx + 1}: {e}")
            batch_results["results"].append(
                {"dataset_index": dataset_idx + 1, "error": str(e), "algorithms": {}}
            )

    # Salvar resultados
    if output_config.get("save_detailed_results", True):
        save_batch_results(batch_results, output_config)

    # Exibir resumo
    display_batch_summary(batch_results)

    return batch_results


def save_optimization_plots(
    result, algorithm_name: str, dataset_idx: int, output_config: Dict[str, Any]
):
    """Salva gráficos de otimização."""
    try:
        from src.optimization.visualization import OptimizationVisualizer

        visualizer = OptimizationVisualizer(result)
        plots_dir = output_config.get("results_dir", "outputs/batch_optimization")
        plots_dir = os.path.join(plots_dir, "plots")
        os.makedirs(plots_dir, exist_ok=True)

        # Gráfico de histórico
        history_path = os.path.join(
            plots_dir, f"{algorithm_name}_dataset{dataset_idx}_history.png"
        )
        visualizer.plot_optimization_history(save_path=history_path, interactive=False)

        # Gráfico de importância dos parâmetros
        importance_path = os.path.join(
            plots_dir, f"{algorithm_name}_dataset{dataset_idx}_importance.png"
        )
        visualizer.plot_parameter_importance(
            save_path=importance_path, interactive=False
        )

    except Exception as e:
        console.print(f"⚠️ Erro ao salvar gráficos: {e}")


def save_batch_results(batch_results: Dict[str, Any], output_config: Dict[str, Any]):
    """Salva resultados do batch."""
    results_dir = output_config.get("results_dir", "outputs/batch_optimization")
    os.makedirs(results_dir, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    results_file = os.path.join(results_dir, f"batch_optimization_{timestamp}.json")

    with open(results_file, "w", encoding="utf-8") as f:
        json.dump(batch_results, f, indent=2, default=str, ensure_ascii=False)

    console.print(f"💾 Resultados salvos em: {results_file}")


def display_batch_summary(batch_results: Dict[str, Any]):
    """Exibe resumo dos resultados."""
    console.print("\n📊 Resumo da Otimização em Lote:")
    console.print("=" * 50)

    for result in batch_results["results"]:
        if "error" in result:
            console.print(f"\n❌ Dataset {result['dataset_index']}: {result['error']}")
            continue

        dataset_idx = result["dataset_index"]
        dataset_name = result["dataset_info"].get("name", f"Dataset {dataset_idx}")

        console.print(f"\n📋 {dataset_name}:")
        console.print(f"   Tipo: {result['dataset_info'].get('type', 'N/A')}")
        console.print(f"   Sequências: {result['dataset_info'].get('n', 'N/A')}")
        console.print(f"   Tamanho: {result['dataset_info'].get('L', 'N/A')}")

        for alg_name, alg_result in result["algorithms"].items():
            if alg_result.get("success", False):
                console.print(
                    f"   ✅ {alg_name}: {alg_result['best_value']:.6f} ({alg_result['optimization_time']:.1f}s)"
                )
            else:
                console.print(
                    f"   ❌ {alg_name}: {alg_result.get('error', 'Erro desconhecido')}"
                )


# Função de conveniência para uso na interface
def run_yaml_optimization_batch():
    """Executa otimização em lote selecionando arquivo YAML."""
    console.print("\n=== Otimização em Lote (YAML) ===")

    # Listar arquivos de otimização disponíveis
    config_dir = "batch_configs"
    if os.path.exists(config_dir):
        optimization_files = [
            f
            for f in os.listdir(config_dir)
            if f.endswith(".yaml") and "otimizacao" in f.lower()
        ]

        if optimization_files:
            console.print("\nArquivos de otimização disponíveis:")
            for idx, file in enumerate(optimization_files, 1):
                console.print(f" {idx}) {file}")

            choice = input(f"Escolha o arquivo [1]: ").strip()
            if choice.isdigit() and 1 <= int(choice) <= len(optimization_files):
                config_file = os.path.join(
                    config_dir, optimization_files[int(choice) - 1]
                )
            else:
                config_file = os.path.join(config_dir, optimization_files[0])
        else:
            config_file = input("Caminho para arquivo de configuração: ").strip()
    else:
        config_file = input("Caminho para arquivo de configuração: ").strip()

    if not config_file:
        console.print("❌ Nenhum arquivo selecionado.")
        return

    try:
        execute_optimization_batch(config_file)
    except Exception as e:
        console.print(f"❌ Erro na execução: {e}")
        import traceback

        traceback.print_exc()


def run_batch_optimization(config_file: str) -> Dict[str, Any]:
    """
    Função de conveniência para executar otimização em lote.

    Args:
        config_file: Caminho para arquivo de configuração YAML

    Returns:
        Dict com resultados da otimização
    """
    return execute_optimization_batch(config_file)


if __name__ == "__main__":
    # Teste direto
    run_yaml_optimization_batch()
