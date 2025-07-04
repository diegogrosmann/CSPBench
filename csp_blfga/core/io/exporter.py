"""
Sistema de exportação centralizado para resultados CSP.

Este módulo unifica todas as funcionalidades de exportação de resultados
para diferentes formatos (CSV, JSON, TXT).

Classes:
    CSPExporter: Exportador centralizado para resultados CSP.
"""

import csv
import json
import logging
from pathlib import Path
from typing import Any


class CSPExporter:
    """
    Exportador centralizado para resultados de algoritmos CSP.

    Suporta exportação para múltiplos formatos (CSV, JSON, TXT).
    """

    def __init__(self):
        self.logger = logging.getLogger(__name__)

    def export_to_csv(
        self,
        results: dict[str, list[dict[str, Any]]],
        filename: str,
        extra_info: dict[str, Any] | None = None,
    ) -> None:
        """
        Exporta resultados para arquivo CSV.

        Args:
            results: Resultados por algoritmo.
            filename: Caminho do arquivo CSV.
            extra_info: Informações extras do experimento.
        """
        file_path = Path(filename)
        file_path.parent.mkdir(parents=True, exist_ok=True)

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

        # Obter informações extras
        extra_info = extra_info or {}
        dataset_strings = extra_info.get("dataset_strings", [])
        params = extra_info.get("params", {})

        with open(file_path, "w", encoding="utf-8", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=headers)
            writer.writeheader()

            for alg_name, executions in results.items():
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
                            "distancia_string_base",
                            params.get("distancia_string_base", ""),
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

        self.logger.info(f"Resultados exportados para CSV: {file_path}")

    def export_batch_json_to_csv(self, json_path: str, csv_path: str) -> None:
        """
        Exporta resultados de batch (JSON) para CSV.

        Args:
            json_path: Caminho do arquivo JSON de entrada.
            csv_path: Caminho do arquivo CSV de saída.
        """
        with open(json_path, encoding="utf-8") as f:
            data = json.load(f)

        execucoes = data.get("execucoes", [])
        rows = []

        for execucao in execucoes:
            config_nome = execucao.get("config_nome", "")
            bases_info = execucao.get("bases_info", [])
            algoritmos_executados = execucao.get("algoritmos_executados", {})

            # Para cada algoritmo
            for alg, bases in algoritmos_executados.items():
                # Para cada base
                for base_key, result in bases.items():
                    base_idx = base_key.replace("base_", "")
                    base_params = {}

                    for b in bases_info:
                        if str(b.get("base_idx")) == str(base_idx):
                            base_params = b.get("params", {})
                            break

                    # Se houver execuções detalhadas, exportar todas
                    execucoes_detalhadas = result.get("execucoes_detalhadas")
                    if execucoes_detalhadas and isinstance(execucoes_detalhadas, list):
                        for exec_idx, exec_data in enumerate(execucoes_detalhadas, 1):
                            row = self._create_batch_row(
                                config_nome,
                                alg,
                                base_idx,
                                exec_idx,
                                exec_data,
                                result,
                                base_params,
                            )
                            rows.append(row)
                    else:
                        # Caso antigo: só resultado agregado
                        row = self._create_batch_row(
                            config_nome, alg, base_idx, 1, result, result, base_params
                        )
                        rows.append(row)

        self._write_batch_csv(csv_path, rows)
        self.logger.info(f"Batch exportado para CSV: {csv_path}")

    def _create_batch_row(
        self,
        config_nome: str,
        alg: str,
        base_idx: str,
        exec_idx: int,
        exec_data: dict[str, Any],
        result: dict[str, Any],
        base_params: dict[str, Any],
    ) -> dict[str, Any]:
        """Cria uma linha de dados para exportação batch."""
        return {
            "config_nome": config_nome,
            "algoritmo": alg,
            "base_idx": base_idx,
            "execucao_idx": exec_idx,
            "dist": exec_data.get(
                "distancia", exec_data.get("melhor_distancia", result.get("dist", ""))
            ),
            "dist_base": result.get("dist_base", ""),
            "tempo": exec_data.get("tempo", result.get("time", "")),
            "status": exec_data.get("status", result.get("status", "")),
            "erro": exec_data.get("erro", ""),
            "melhor_string": exec_data.get("melhor_string", ""),
            "n": base_params.get("n", ""),
            "L": base_params.get("L", ""),
            "alphabet": base_params.get("alphabet", ""),
            "noise": base_params.get("noise", ""),
            "seed": base_params.get("seed", ""),
            "distancia_string_base": base_params.get("distancia_string_base", ""),
        }

    def _write_batch_csv(self, csv_path: str, rows: list[dict[str, Any]]) -> None:
        """Escreve dados de batch para arquivo CSV."""
        headers = [
            "config_nome",
            "algoritmo",
            "base_idx",
            "execucao_idx",
            "dist",
            "dist_base",
            "tempo",
            "status",
            "erro",
            "melhor_string",
            "n",
            "L",
            "alphabet",
            "noise",
            "seed",
            "distancia_string_base",
        ]

        csv_path = Path(csv_path)
        csv_path.parent.mkdir(parents=True, exist_ok=True)

        with open(csv_path, "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=headers)
            writer.writeheader()
            for row in rows:
                writer.writerow(row)

    def export_to_json(
        self,
        data: dict[str, Any],
        filename: str,
        convert_keys: bool = True,
    ) -> None:
        """
        Exporta dados para arquivo JSON.

        Args:
            data: Dados a serem exportados.
            filename: Caminho do arquivo JSON.
            convert_keys: Se True, converte chaves de dict para string.
        """
        file_path = Path(filename)
        file_path.parent.mkdir(parents=True, exist_ok=True)

        if convert_keys:
            data = self._convert_dict_keys_to_str(data)

        with open(file_path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, ensure_ascii=False)

        self.logger.info(f"Dados exportados para JSON: {file_path}")

    def _convert_dict_keys_to_str(self, obj: Any) -> Any:
        """
        Recursivamente converte todas as chaves de dicionários para string.

        Args:
            obj: Objeto a ser convertido.
        Returns:
            Objeto com todas as chaves de dicionário como string.
        """
        if isinstance(obj, dict):
            return {str(k): self._convert_dict_keys_to_str(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [self._convert_dict_keys_to_str(i) for i in obj]
        else:
            return obj
