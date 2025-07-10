"""
Executor de análise de sensibilidade em lote baseado em arquivos YAML.
"""

import os
from datetime import datetime

import yaml

from src.datasets.dataset_entrez import fetch_dataset
from src.datasets.dataset_file import load_dataset
from src.datasets.dataset_synthetic import generate_dataset_from_params
from src.optimization.sensitivity_analyzer import analyze_algorithm_sensitivity


def run_yaml_sensitivity_batch():
    """Executa análise de sensibilidade em lote selecionando arquivo YAML."""
    from src.ui.cli.menu import select_sensitivity_yaml_file

    print("\n=== Execução em lote de análise de sensibilidade (YAML) ===")

    # Selecionar arquivo usando a nova função
    config_file = select_sensitivity_yaml_file()
    if not config_file:
        print("❌ Nenhum arquivo selecionado.")
        return

    # Verificar se arquivo existe
    if not config_file or not os.path.exists(config_file):
        print(f"❌ Arquivo não encontrado: {config_file}")
        return

    print(f"\n📋 Carregando configuração: {os.path.basename(config_file)}")

    try:
        with open(config_file, "r", encoding="utf-8") as f:
            config = yaml.safe_load(f)
    except (FileNotFoundError, yaml.YAMLError, PermissionError) as e:
        print(f"❌ Erro ao carregar configuração: {e}")
        return

    # Validar se é um arquivo de sensibilidade válido
    if "sensitivity_config" not in config:
        print("❌ Arquivo inválido: não contém configuração 'sensitivity_config'")
        print(
            "   Este arquivo pode ser para otimização. Use um arquivo de sensibilidade."
        )
        return

    # Extrair configurações
    batch_info = config.get("batch_info", {})
    sensitivity_config = config.get("sensitivity_config", {})
    algorithms = config.get("algorithms", [])
    datasets = config.get("datasets", [])
    output_cfg = config.get("output", {})

    # Validações básicas
    if not algorithms:
        print("❌ Nenhum algoritmo especificado no arquivo YAML")
        return

    if not datasets:
        print("❌ Nenhum dataset especificado no arquivo YAML")
        return

    print(f"\n📊 Configuração carregada:")
    print(f"   Batch: {batch_info.get('nome', 'Sem nome')}")
    print(f"   Descrição: {batch_info.get('descricao', 'Sem descrição')}")
    print(f"   Algoritmos: {', '.join(algorithms)}")
    print(f"   Datasets: {len(datasets)} configurações")
    print(f"   Método: {sensitivity_config.get('method', 'morris')}")
    print(f"   Amostras: {sensitivity_config.get('n_samples', 100)}")

    # Confirmação antes de executar
    confirm = input("\n🚀 Executar análise de sensibilidade? [S/n]: ").strip().lower()
    if confirm and confirm not in ["s", "sim", "y", "yes"]:
        print("❌ Operação cancelada pelo usuário.")
        return

    print(f"\n{'='*60}")
    print("🔬 INICIANDO ANÁLISE DE SENSIBILIDADE EM LOTE")
    print(f"{'='*60}")

    results = []
    total_runs = len(datasets) * len(algorithms)
    current_run = 0

    for i, dataset in enumerate(datasets, 1):
        print(f"\n📊 Dataset {i}/{len(datasets)}: {dataset.get('nome', 'Sem nome')}")

        # Preparar dataset
        tipo = dataset.get("tipo", "synthetic")
        params = dataset.get("parametros", {})

        if tipo == "synthetic":
            seqs, _ = generate_dataset_from_params(**params)
            alphabet = params.get("alphabet", "ACGT")
        elif tipo == "file":
            # Usar load_dataset_with_params para carregar arquivo específico
            if "filename" in params:
                from src.datasets.dataset_file import load_dataset_with_params

                seqs, _ = load_dataset_with_params({"filepath": params["filename"]})
            else:
                seqs, _ = load_dataset(silent=True)
            alphabet = "".join(sorted(set("".join(seqs))))
        elif tipo == "entrez":
            seqs, _ = fetch_dataset()
            alphabet = "".join(sorted(set("".join(seqs))))
        else:
            print(f"❌ Tipo de dataset não suportado: {tipo}")
            continue

        print(
            f"   Sequências: {len(seqs)} | Comprimento: {len(seqs[0]) if seqs else 0} | Alfabeto: {alphabet}"
        )

        # Executar análise para cada algoritmo
        for alg in algorithms:
            current_run += 1
            print(f"\n🔬 [{current_run}/{total_runs}] Analisando {alg}...")

            try:
                result = analyze_algorithm_sensitivity(
                    algorithm_name=alg,
                    sequences=seqs,
                    alphabet=alphabet,
                    yaml_config=config,
                    **sensitivity_config,
                )

                print(
                    f"   ✅ Concluído: {result.method} | {result.n_samples} amostras | {result.analysis_time:.2f}s"
                )

                results.append(
                    {
                        "dataset": dataset.get("nome", ""),
                        "algorithm": alg,
                        "method": result.method,
                        "n_samples": result.n_samples,
                        "analysis_time": result.analysis_time,
                        "parameter_names": result.parameter_names,
                        "timestamp": datetime.now().isoformat(),
                    }
                )

            except (RuntimeError, ValueError, TypeError) as e:
                print(f"   ❌ Erro: {e}")
                results.append(
                    {
                        "dataset": dataset.get("nome", ""),
                        "algorithm": alg,
                        "error": str(e),
                        "timestamp": datetime.now().isoformat(),
                    }
                )

    print(f"\n{'='*60}")
    print("💾 SALVANDO RESULTADOS")
    print(f"{'='*60}")

    # Salvamento dos resultados
    results_dir = output_cfg.get("results_dir", "outputs/batch_sensitivity")
    os.makedirs(results_dir, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = os.path.join(results_dir, f"batch_sensitivity_{timestamp}.json")

    # Adicionar metadados ao arquivo de resultados
    output_data = {
        "metadata": {
            "config_file": os.path.basename(config_file),
            "timestamp": datetime.now().isoformat(),
            "batch_info": batch_info,
            "sensitivity_config": sensitivity_config,
            "total_runs": total_runs,
            "successful_runs": len([r for r in results if "error" not in r]),
            "failed_runs": len([r for r in results if "error" in r]),
        },
        "results": results,
    }

    import json

    with open(filename, "w", encoding="utf-8") as f:
        json.dump(output_data, f, indent=2, ensure_ascii=False)

    print(f"� Resultados salvos em: {filename}")

    # Resumo final
    successful = len([r for r in results if "error" not in r])
    failed = len([r for r in results if "error" in r])
    total_time = sum(r.get("analysis_time", 0) for r in results if "analysis_time" in r)

    print(f"\n{'='*60}")
    print("📊 RESUMO DA EXECUÇÃO")
    print(f"{'='*60}")
    print(f"✅ Execuções bem-sucedidas: {successful}")
    print(f"❌ Execuções com erro: {failed}")
    print(f"⏱️  Tempo total de análise: {total_time:.2f}s")
    print(f"📁 Arquivo de saída: {os.path.basename(filename)}")

    if failed > 0:
        print(
            f"\n⚠️  {failed} execuções falharam. Verifique os logs acima para detalhes."
        )

    print(f"\n🎉 Análise de sensibilidade em lote concluída!")
