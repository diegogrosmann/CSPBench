"""
Formatação e salvamento de resultados de execuções de algoritmos CSP.

Classes:
    ResultsFormatter: Armazena, formata e salva resultados detalhados e comparativos.

Métodos:
    add_algorithm_results(...): Adiciona resultados de um algoritmo.
    format_detailed_results(): Formata resultados detalhados.
    format_quick_summary(): Formata resumo rápido dos resultados.
    print_quick_summary(...): Exibe resumo rápido no console.
    save_detailed_report(...): Salva relatório em arquivo.
    export_to_csv(...): Exporta resultados para CSV usando CSPExporter.
"""

import logging
import statistics
from datetime import datetime
from pathlib import Path
from typing import Any

from tabulate import tabulate

from src.core.io.exporter import CSPExporter


class ResultsFormatter:
    """
    Armazena, formata e salva resultados detalhados e comparativos de execuções de algoritmos CSP.

    Attributes:
        results (dict): Resultados por algoritmo.
        extra_info (dict): Informações extras do experimento.
    """

    def __init__(self):
        self.results = {}
        self.extra_info = {}
        self.exporter = CSPExporter()

    def add_algorithm_results(
        self, algorithm_name: str, executions: list[dict[str, Any]]
    ):
        """
        Adiciona resultados de um algoritmo.

        Args:
            algorithm_name (str): Nome do algoritmo.
            executions (list): Lista de execuções detalhadas.
        """
        self.results[algorithm_name] = executions

    def format_detailed_results(self) -> str:
        """
        Formata os resultados detalhados para apresentação.

        Returns:
            str: Resultados formatados.
        """
        logger = logging.getLogger(__name__)
        output = []
        output.append("=" * 80)
        output.append("RELATÓRIO DETALHADO DE RESULTADOS")
        output.append("=" * 80)

        # Adicionar informações do dataset
        output.append(self._format_dataset_info())

        # Tabelas individuais por algoritmo
        for algorithm_name in self.results:
            output.append(self._format_algorithm_executions_table(algorithm_name))
            output.append(self._format_algorithm_summary_statistics(algorithm_name))
            output.append(self._format_algorithm_best_strings(algorithm_name))
            output.append("\n" + "=" * 80 + "\n")

        # Tabela comparativa final
        output.append(self._format_algorithms_comparison())

        return "\n".join(output)

    def _format_algorithm_executions_table(self, algorithm_name: str) -> str:
        """Formata tabela individual de execuções de um algoritmo"""
        executions = self.results[algorithm_name]
        output = [f"\n📊 TABELA DE EXECUÇÕES - {algorithm_name.upper()}", "-" * 60]
        headers = ["Execução", "Tempo (s)", "Distância", "Status"]

        def format_row(i, exec_data):
            distancia = exec_data.get(
                "distancia", exec_data.get("melhor_distancia", "-")
            )
            if exec_data.get("erro"):
                return [i, f"{exec_data['tempo']:.4f}", "-", f"✗ {exec_data['erro']}"]
            if exec_data.get("timeout"):
                return [i, f"{exec_data['tempo']:.4f}", distancia, "⏰ Timeout"]
            if distancia == float("inf"):
                return [i, f"{exec_data['tempo']:.4f}", "∞", "∞ Sem solução"]
            return [i, f"{exec_data['tempo']:.4f}", distancia, "✓ OK"]

        table_data = [
            format_row(i, exec_data) for i, exec_data in enumerate(executions, 1)
        ]
        table = tabulate(
            table_data, headers=headers, tablefmt="grid", stralign="center"
        )
        output.append(table)
        return "\n".join(output)

    def _format_algorithm_summary_statistics(self, algorithm_name: str) -> str:
        """Formata estatísticas resumidas de um algoritmo"""
        executions = self.results[algorithm_name]
        valid = [
            e
            for e in executions
            if not e.get("erro")
            and e.get("distancia", e.get("melhor_distancia")) not in ["-", float("inf")]
        ]
        output = [f"\n📈 ESTATÍSTICAS DETALHADAS - {algorithm_name.upper()}", "-" * 60]
        if not valid:
            output.append("❌ Nenhuma execução válida para calcular estatísticas.")
            return "\n".join(output)
        tempos = [e["tempo"] for e in valid]
        distancias = [
            e.get("distancia", e.get("melhor_distancia", float("inf"))) for e in valid
        ]

        def stat_line(label, val, sufixo=""):
            return [label, f"{val}{sufixo}"]

        stats_data = [
            stat_line("Execuções Válidas", f"{len(valid)}/{len(executions)}"),
            stat_line("Média", f"{statistics.mean(tempos):.4f}", " s"),
            stat_line("Mediana", f"{statistics.median(tempos):.4f}", " s"),
            stat_line(
                "Desvio Padrão",
                f"{statistics.stdev(tempos) if len(tempos) > 1 else 0:.4f}",
                " s",
            ),
            stat_line("Mínimo", f"{min(tempos):.4f}", " s"),
            stat_line("Máximo", f"{max(tempos):.4f}", " s"),
            ["", ""],
            stat_line("Média Distância", f"{statistics.mean(distancias):.2f}"),
            stat_line("Mediana Distância", f"{statistics.median(distancias):.2f}"),
            stat_line(
                "Desvio Padrão Distância",
                f"{statistics.stdev(distancias) if len(distancias) > 1 else 0:.2f}",
            ),
            stat_line("Melhor (Mínima)", min(distancias)),
            stat_line("Pior (Máxima)", max(distancias)),
        ]
        table = tabulate(
            stats_data, headers=["Métrica", "Valor"], tablefmt="grid", stralign="left"
        )
        output.append(table)
        return "\n".join(output)

    def _format_algorithm_best_strings(self, algorithm_name: str) -> str:
        """Formata melhores strings encontradas para auditoria"""
        executions = self.results[algorithm_name]
        output = [f"\n🔍 STRINGS PARA AUDITORIA - {algorithm_name.upper()}", "-" * 60]
        for i, e in enumerate(executions, 1):
            distancia = e.get("distancia", e.get("melhor_distancia", "-"))
            string_result = e.get("melhor_string", e.get("string_resultado", ""))
            iteracoes = e.get("iteracoes", e.get("num_iteracoes", 0))
            output.append(f"Execução {i:2d}:")
            if e.get("params") is not None:
                output.append(f"  Parâmetros de geração: {e['params']}")
            if e.get("erro"):
                output.append(f"  ❌ Erro: {e['erro']}")
                output.append(f"  Tempo: {e['tempo']:.4f}s")
            else:
                output.append(f"  String Centro (Base): '{string_result}'")
                output.append(f"  Distância (algoritmo): {distancia}")
                output.append(f"  Iterações: {iteracoes}")
                output.append(f"  Tempo: {e['tempo']:.4f}s")
            output.append("")
        return "\n".join(output)

    def _format_algorithms_comparison(self) -> str:
        """Formata tabela comparativa entre algoritmos (renomeado para clareza)"""
        output = ["\n🏆 TABELA COMPARATIVA FINAL"]
        output.append("=" * 80)
        if not self.results:
            output.append("Nenhum resultado disponível para comparação.")
            return "\n".join(output)

        # Extrair dados comparativos
        comparative_data = []
        base_algorithms = {}

        for algorithm_full_name, executions in self.results.items():
            distancias = []
            tempos = []
            taxa_sucesso = 0.0

            # Simplificar parsing de algoritmo/base
            algo_name, base_display = self._parse_algorithm_name(algorithm_full_name)

            valid_executions = [
                exec_data
                for exec_data in executions
                if not exec_data.get("erro")
                and exec_data.get("distancia", exec_data.get("melhor_distancia"))
                not in ["-", float("inf")]
            ]

            if not valid_executions:
                tempos = [exec_data["tempo"] for exec_data in executions]
                row = [
                    base_display,
                    algo_name,
                    f"{statistics.mean(tempos) if tempos else 0:.4f}",
                    f"{statistics.stdev(tempos) if len(tempos) > 1 else 0:.4f}",
                    "ERRO",
                    "ERRO",
                    "ERRO",
                    "0.0",
                ]
            else:
                tempos = [exec_data["tempo"] for exec_data in valid_executions]
                distancias = [
                    exec_data.get(
                        "distancia", exec_data.get("melhor_distancia", float("inf"))
                    )
                    for exec_data in valid_executions
                ]
                taxa_sucesso = (
                    len(valid_executions) / len(executions) * 100
                    if len(executions) > 0
                    else 0
                )

                row = [
                    base_display,
                    algo_name,
                    f"{statistics.mean(tempos) if tempos else 0:.4f}",
                    f"{statistics.stdev(tempos) if len(tempos) > 1 else 0:.4f}",
                    f"{min(distancias) if distancias else 'ERRO'}",
                    f"{statistics.mean(distancias):.2f}" if distancias else "ERRO",
                    (
                        f"{statistics.stdev(distancias) if len(distancias) > 1 else 0:.2f}"
                        if distancias
                        else "ERRO"
                    ),
                    f"{taxa_sucesso:.1f}",
                ]

            comparative_data.append(row)

            if base_display not in base_algorithms:
                base_algorithms[base_display] = []

            if valid_executions:
                base_algorithms[base_display].append(
                    {
                        "name": algo_name,
                        "dist": min(distancias) if distancias else float("inf"),
                        "time": statistics.mean(tempos) if tempos else 0,
                        "success_rate": (
                            taxa_sucesso if isinstance(taxa_sucesso, int | float) else 0
                        ),
                    }
                )
            else:
                base_algorithms[base_display].append(
                    {
                        "name": algo_name,
                        "dist": float("inf"),
                        "time": statistics.mean(tempos) if tempos else 0,
                        "success_rate": 0,
                    }
                )

        headers = [
            "Base",
            "Algoritmo",
            "Tempo Médio (s)",
            "Desvio Tempo",
            "Melhor Distância",
            "Distância Média",
            "Desvio Distancia",
            "Taxa Sucesso (%)",
        ]

        def sort_key(row):
            base_name = row[0]
            if row[4] == "ERRO":
                dist = float("inf")
            else:
                try:
                    dist = float(row[4])
                except Exception:
                    dist = float("inf")
            return (base_name, dist, float(row[2]))

        comparative_data.sort(key=sort_key)
        table = tabulate(
            comparative_data, headers=headers, tablefmt="grid", stralign="center"
        )
        output.append(table)

        # Adicionar ranking POR BASE
        output.append("\n🥇 RANKING POR BASE:")
        output.append("-" * 60)

        for base_name, algorithms in sorted(base_algorithms.items()):
            output.append(f"\n💠 BASE: {base_name}")
            output.append("-" * 40)

            sorted_algorithms = sorted(algorithms, key=lambda x: (x["dist"], x["time"]))

            for i, alg in enumerate(sorted_algorithms, 1):
                if alg["dist"] == float("inf"):
                    medal = "❌"
                    info = f"Falha na execução | Tempo: {alg['time']:.4f}s"
                else:
                    medal = (
                        "🥇"
                        if i == 1
                        else "🥈" if i == 2 else "🥉" if i == 3 else f"{i}°"
                    )
                    info = f"Distância: {alg['dist']} | Tempo: {alg['time']:.4f}s | Sucesso: {alg['success_rate']:.1f}%"
                output.append(f"{medal} {alg['name']} - {info}")

        return "\n".join(output)

    def save_detailed_report(self, filename: str | None = None):
        """
        Salva o relatório detalhado em um arquivo na estrutura outputs/reports/.

        Args:
            filename: Nome do arquivo ou None para gerar automaticamente
        """
        if filename is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"detailed_report_{timestamp}.txt"

        # Usar estrutura padronizada outputs/reports/<timestamp>/
        if "/" not in filename:
            timestamp_dir = datetime.now().strftime("%Y%m%d_%H%M%S")
            file_path = Path("outputs") / "reports" / timestamp_dir / filename
        else:
            file_path = Path(filename)

        # Garantir que o diretório existe
        file_path.parent.mkdir(parents=True, exist_ok=True)

        with open(file_path, "w", encoding="utf-8") as f:
            f.write(self.format_detailed_results())

        logging.getLogger(__name__).info(f"Relatório salvo em: {file_path}")

    def get_detailed_report_data(self):
        """
        Obtém os dados do relatório detalhado.

        Returns:
            dict: Dados do relatório, incluindo informações extras e resultados.
        """
        return {"extra_info": self.extra_info, "results": self.results}

    def export_to_csv(self, filename: str) -> None:
        """
        Exporta resultados para arquivo CSV usando CSPExporter.

        Args:
            filename: Caminho do arquivo CSV de saída.
        """
        self.exporter.export_to_csv(self.results, filename, self.extra_info)

    def export_to_json(self, filename: str) -> None:
        """
        Exporta dados do relatório para arquivo JSON usando CSPExporter.

        Args:
            filename: Caminho do arquivo JSON de saída.
        """
        data = self.get_detailed_report_data()
        self.exporter.export_to_json(data, filename)

    def _format_dataset_info(self) -> str:
        """Formata informações de todos os datasets presentes em extra_info, iterando sobre cada um."""
        output = ["\n📊 INFORMAÇÕES DO DATASET", "-" * 60]

        if not hasattr(self, "extra_info") or not self.extra_info:
            output.append("❌ Nenhuma informação do dataset disponível.")
            return "\n".join(output)

        encontrou = False
        for nome, params in self.extra_info.items():
            if not isinstance(params, dict):
                continue
            encontrou = True
            output.append(f"\n📋 Parâmetros do Dataset: {nome}")
            if "n" in params:
                output.append(f"  • Número de strings: {params['n']}")
            if "L" in params:
                output.append(f"  • Comprimento das strings: {params['L']}")
            if "alphabet" in params:
                value = params["alphabet"]
                output.append(f"  • Alfabeto: '{value}' (|Σ| = {len(str(value))})")
            if "noise" in params:
                value = params["noise"]
                if isinstance(value, list):
                    output.append(
                        f"  • Taxa de ruído: variável (min: {min(value):.3f}, max: {max(value):.3f})"
                    )
                else:
                    output.append(f"  • Taxa de ruído: {value}")
            if "seed" in params:
                output.append(f"  • Semente: {params['seed']}")
            if "fully_random" in params:
                output.append(
                    f"  • Modo: {'Totalmente aleatório' if params['fully_random'] else 'Base + ruído'}"
                )
            if "base_string" in params or "distancia_string_base" in params:
                output.append("")
                output.append("🎯 STRING BASE PARA AUDITORIA:")
                if "base_string" in params:
                    output.append(f"  • String base: '{params['base_string']}'")
                if "distancia_string_base" in params:
                    output.append(
                        f"  • Distância da string base: {params['distancia_string_base']}"
                    )

        if not encontrou:
            output.append("❌ Nenhum parâmetro de dataset encontrado.")
        return "\n".join(output)

    def format_quick_summary(self) -> str:
        """
        Formata um resumo rápido dos resultados dos algoritmos.

        Returns:
            String formatada com resumo dos algoritmos
        """
        if not self.results:
            return "❌ Nenhum resultado disponível para resumo."

        output = []
        output.append("=" * 60)
        output.append("RESUMO RÁPIDO DOS RESULTADOS")
        output.append("=" * 60)

        # Cabeçalho
        header_format = f"{'Algoritmo':<15} {'Melhor Dist':<12} {'Dist. Base':<12} {'Tempo Médio':<12} {'Sucessos':<10}"
        output.append(header_format)
        output.append("-" * 70)

        for alg_name, executions in self.results.items():
            # Calcular estatísticas
            successful_execs = [e for e in executions if not e.get("erro")]

            if successful_execs:
                # Melhor distância
                distances = [e.get("distancia", float("inf")) for e in successful_execs]
                valid_distances = [d for d in distances if d != float("inf")]
                best_dist = min(valid_distances) if valid_distances else "∞"

                # Tempo médio
                times = [e.get("tempo", 0) for e in successful_execs]
                avg_time = statistics.mean(times) if times else 0

                # Distância da base (se disponível)
                dist_base = successful_execs[0].get("dist_base", "-")

                # Taxa de sucesso
                success_rate = f"{len(successful_execs)}/{len(executions)}"

                # Formatar valores
                best_dist_str = str(best_dist) if best_dist != "∞" else "∞"
                dist_base_str = str(dist_base) if dist_base != "-" else "-"
                time_str = f"{avg_time:.3f}s"

            else:
                best_dist_str = "FALHOU"
                dist_base_str = "-"
                time_str = "-"
                success_rate = f"0/{len(executions)}"

            # Linha do algoritmo
            row = f"{alg_name:<15} {best_dist_str:<12} {dist_base_str:<12} {time_str:<12} {success_rate:<10}"
            output.append(row)

        output.append("=" * 60)
        return "\n".join(output)

    def print_quick_summary(self, console=None):
        """
        Exibe um resumo rápido dos resultados no console.

        Args:
            console: Instância do gerenciador de console (opcional)
        """
        summary = self.format_quick_summary()

        if console:
            console.print("\n" + summary)
        else:
            print("\n" + summary)

    def _parse_algorithm_name(self, algorithm_full_name: str) -> tuple[str, str]:
        """
        Simplifica o parsing de nomes de algoritmos e bases.

        Args:
            algorithm_full_name: Nome completo do algoritmo (pode incluir base)

        Returns:
            Tupla (nome_algoritmo, nome_base_formatado)
        """
        # Parsing simplificado: buscar padrões comuns
        if "_Base" in algorithm_full_name:
            # Formato: "DatasetName_Base123_AlgorithmName"
            parts = algorithm_full_name.split("_")

            # Encontrar índice do "Base"
            base_index = -1
            for i, part in enumerate(parts):
                if "Base" in part:
                    base_index = i
                    break

            if base_index >= 0 and base_index < len(parts) - 1:
                # Dividir entre base e algoritmo
                base_parts = parts[: base_index + 1]
                algo_parts = parts[base_index + 1 :]

                base_name = "_".join(base_parts[:-1])
                base_num = base_parts[-1]
                algo_name = "_".join(algo_parts)

                base_display = f"{base_name} ({base_num})"
            else:
                # Fallback
                algo_name = algorithm_full_name
                base_display = "N/A"
        else:
            # Sem base identificável
            algo_name = algorithm_full_name
            base_display = "N/A"

        return algo_name, base_display

    def _simplify_algorithm_label(self, algorithm_name: str) -> str:
        """
        Simplifica rótulos de algoritmos para exibição.

        Args:
            algorithm_name: Nome do algoritmo

        Returns:
            Rótulo simplificado
        """
        # Remover prefixos comuns
        if algorithm_name.startswith("Algorithm"):
            algorithm_name = algorithm_name[9:]

        # Substituir underscores por espaços
        algorithm_name = algorithm_name.replace("_", " ")

        # Capitalizar primeira letra
        return algorithm_name.strip().title()
