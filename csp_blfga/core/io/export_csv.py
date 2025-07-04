import csv
from pathlib import Path


def export_results_to_csv(formatter, filename):
    """
    Exporta todos os dados de execuções detalhadas para um arquivo CSV.
    Cada linha corresponde a uma execução individual de um algoritmo.
    """
    file_path = Path(filename)
    file_path.parent.mkdir(parents=True, exist_ok=True)

    # Cabeçalhos padrão incluindo mais informações
    headers = [
        "algoritmo",
        "execucao",
        "melhor_string",
        "distancia",
        "tempo",
        "status",
        "erro",
        "distancia_string_base",
        "seed",
        "n",
        "L",
        "alphabet",
    ]

    # Obter informações extras se disponíveis
    extra_info = getattr(formatter, "extra_info", {})
    dataset_strings = extra_info.get("dataset_strings", [])
    params = extra_info.get("params", {})

    with open(file_path, "w", encoding="utf-8", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=headers)
        writer.writeheader()
        for alg_name, executions in formatter.results.items():
            for idx, exec_data in enumerate(executions, 1):
                row = {
                    "algoritmo": alg_name,
                    "execucao": idx,
                    "melhor_string": exec_data.get("melhor_string", ""),
                    "distancia": exec_data.get(
                        "distancia", exec_data.get("melhor_distancia", "")
                    ),
                    "tempo": exec_data.get("tempo", ""),
                    "status": exec_data.get(
                        "status",
                        (
                            "sucesso"
                            if "distancia" in exec_data
                            and exec_data["distancia"] != float("inf")
                            else "erro"
                        ),
                    ),
                    "erro": exec_data.get("erro", ""),
                    "distancia_string_base": exec_data.get(
                        "distancia_string_base", params.get("distancia_string_base", "")
                    ),
                    "seed": exec_data.get("seed", params.get("seed", "")),
                    "n": params.get(
                        "n", len(dataset_strings) if dataset_strings else ""
                    ),
                    "L": params.get(
                        "L",
                        (
                            len(dataset_strings[0])
                            if dataset_strings and len(dataset_strings) > 0
                            else ""
                        ),
                    ),
                    "alphabet": params.get("alphabet", ""),
                }
                writer.writerow(row)
