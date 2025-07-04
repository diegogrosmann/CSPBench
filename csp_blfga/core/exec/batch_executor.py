"""
Sistema de execução em lote para múltiplas configurações de CSP.

Classes:
    BatchConfig: Representa uma configuração de execução
    BatchExecutor: Executa sequência de configurações
"""

import json
import logging
import time
import uuid
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path
from typing import Any

import yaml

from algorithms.base import global_registry
from csp_blfga.core.exec.runner import execute_algorithm_runs
from csp_blfga.core.io.results_formatter import ResultsFormatter
from csp_blfga.ui.cli.console_manager import console
from csp_blfga.utils.config import ALGORITHM_TIMEOUT, safe_input
from csp_blfga.utils.resource_monitor import check_algorithm_feasibility

logger = logging.getLogger(__name__)


def list_batch_configs() -> list[Path]:
    """Lista arquivos de configuração disponíveis."""
    config_dir = Path("batch_configs")
    if not config_dir.exists():
        config_dir.mkdir(exist_ok=True)

    # Procurar por arquivos JSON
    json_files = list(config_dir.glob("*.json"))
    return sorted(json_files)


def select_batch_config() -> str:
    """Menu para seleção de arquivo de configuração de lote."""
    batch_dir = Path("batch_configs")
    if not batch_dir.exists():
        batch_dir.mkdir(parents=True, exist_ok=True)
        console.print("📁 Diretório batch_configs criado.")

    # Buscar arquivos de configuração
    config_files = []
    for pattern in ["*.yaml", "*.yml", "*.json", "*.xml"]:
        config_files.extend(batch_dir.glob(pattern))

    if not config_files:
        console.print("❌ Nenhum arquivo de configuração encontrado em batch_configs/")
        console.print("Crie um arquivo .yaml, .json ou .xml neste diretório.")
        return ""

    console.print("\n📋 Arquivos de configuração disponíveis:")
    for i, config_file in enumerate(config_files, 1):
        console.print(f"  {i}) {config_file.name}")

    while True:
        choice = safe_input(f"\nEscolha um arquivo [1-{len(config_files)}]: ")
        if choice.isdigit() and 1 <= int(choice) <= len(config_files):
            selected_file = config_files[int(choice) - 1]
            console.print(f"✓ Selecionado: {selected_file.name}")
            return str(selected_file)
        else:
            console.print("❌ Opção inválida.")


def create_example_config():
    """Cria arquivo de configuração de exemplo."""
    config_dir = Path("batch_configs")
    config_dir.mkdir(exist_ok=True)

    example_config = {
        "batch_info": {
            "nome": "Experimento Exemplo CSP",
            "descricao": "Configuração de exemplo para testes comparativos",
            "timeout_global": 1800,
        },
        "execucoes": [
            {
                "nome": "Sintético Pequeno",
                "dataset": {
                    "tipo": "synthetic",
                    "parametros": {"n": 10, "L": 50, "alphabet": "ACGT", "noise": 0.1},
                },
                "algoritmos": ["Baseline", "BLF-GA", "H³-CSP"],
                "execucoes_por_algoritmo": 3,
                "timeout": 300,
            },
            {
                "nome": "Sintético Médio",
                "dataset": {
                    "tipo": "synthetic",
                    "parametros": {
                        "n": 20,
                        "L": 100,
                        "alphabet": "ACGT",
                        "noise": 0.15,
                    },
                },
                "algoritmos": ["Baseline", "BLF-GA", "CSC"],
                "execucoes_por_algoritmo": 3,
                "timeout": 600,
            },
            {
                "nome": "Dataset de Arquivo",
                "dataset": {
                    "tipo": "file",
                    "parametros": {"filepath": "saved_datasets/sequences.fasta"},
                },
                "algoritmos": ["Baseline", "BLF-GA"],
                "execucoes_por_algoritmo": 3,
                "timeout": 900,
            },
        ],
    }

    example_path = config_dir / "exemplo_batch.json"
    with open(example_path, "w", encoding="utf-8") as f:
        json.dump(example_config, f, indent=2, ensure_ascii=False)

    console.print(f"📄 Arquivo de exemplo criado: {example_path}")
    console.print("💡 Edite este arquivo para personalizar suas configurações")


class BatchConfig:
    """Representa uma configuração individual de execução."""

    def __init__(self, config_dict: dict[str, Any]):
        self.nome = config_dict.get("nome", "Execução Sem Nome")
        self.dataset_config = config_dict.get("dataset", {})
        self.algoritmos = config_dict.get("algoritmos", [])
        # Mudança de execucoes_por_algoritmo para execucoes_por_algoritmo_por_base
        self.execucoes_por_algoritmo_por_base = config_dict.get(
            "execucoes_por_algoritmo_por_base",
            config_dict.get("execucoes_por_algoritmo", 3),
        )  # Retro-compatibilidade
        self.num_bases = config_dict.get("num_bases", 1)
        self.timeout = config_dict.get("timeout", ALGORITHM_TIMEOUT)

    def __str__(self):
        return f"BatchConfig({self.nome}, {len(self.algoritmos)} algoritmos, {self.execucoes_por_algoritmo_por_base} exec por base, {self.num_bases} bases)"


class BatchExecutor:
    """Executa uma sequência de configurações em lote."""

    def __init__(self, config_file: str):
        self.config_file = Path(config_file)
        self.batch_info = {}
        self.execucoes = []
        self.results = {}
        self.consolidated_formatter = ResultsFormatter()

        # Identificador único para esta execução em lote
        self.batch_id = (
            f"batch_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:8]}"
        )

        # Criar diretório dedicado para os resultados deste lote
        self.results_dir = Path("results") / self.batch_id
        self.results_dir.mkdir(parents=True, exist_ok=True)
        console.print(f"📁 Diretório de resultados criado: {self.results_dir}")

        self._load_config()

    def _load_config(self):
        """Carrega configuração do arquivo YAML, JSON ou XML."""
        try:
            with open(self.config_file, encoding="utf-8") as f:
                if self.config_file.suffix.lower() in [".yaml", ".yml"]:
                    config = yaml.safe_load(f)
                elif self.config_file.suffix.lower() == ".json":
                    config = json.load(f)
                elif self.config_file.suffix.lower() == ".xml":
                    config = self._parse_xml(f)
                else:
                    raise ValueError(
                        f"Formato de arquivo não suportado: {self.config_file.suffix}"
                    )

            self.batch_info = config.get("batch_info", {})
            execucoes_raw = config.get("execucoes", [])

            # Converte para objetos BatchConfig
            self.execucoes = [BatchConfig(exec_config) for exec_config in execucoes_raw]

            console.print(
                f"✓ Configuração carregada: {len(self.execucoes)} execuções planejadas"
            )

        except Exception as e:
            console.print(f"❌ Erro ao carregar configuração: {e}")
            raise

    def _parse_xml(self, file_handle):
        """Converte XML para estrutura de dicionário compatível."""
        tree = ET.parse(file_handle)
        root = tree.getroot()

        config = {}

        # Parse batch_info
        batch_info_elem = root.find("batch_info")
        if batch_info_elem is not None:
            batch_info = {}
            for child in batch_info_elem:
                if child.tag in ["timeout_global"] and child.text is not None:
                    batch_info[child.tag] = int(child.text)
                else:
                    batch_info[child.tag] = child.text if child.text is not None else ""
            config["batch_info"] = batch_info

        # Parse execucoes
        execucoes_elem = root.find("execucoes")
        execucoes = []

        if execucoes_elem is not None:
            for exec_elem in execucoes_elem.findall("execucao"):
                exec_config = {}

                # Nome
                nome_elem = exec_elem.find("nome")
                if nome_elem is not None and nome_elem.text is not None:
                    exec_config["nome"] = nome_elem.text

                # Dataset
                dataset_elem = exec_elem.find("dataset")
                if dataset_elem is not None:
                    dataset_config = {}

                    tipo_elem = dataset_elem.find("tipo")
                    if tipo_elem is not None and tipo_elem.text is not None:
                        dataset_config["tipo"] = tipo_elem.text

                    params_elem = dataset_elem.find("parametros")
                    if params_elem is not None:
                        parametros = {}
                        for param in params_elem:
                            if param.text is None:
                                continue
                            if param.tag in ["n", "L"]:
                                parametros[param.tag] = int(param.text)
                            elif param.tag == "noise":
                                parametros[param.tag] = float(param.text)
                            else:
                                parametros[param.tag] = param.text
                        dataset_config["parametros"] = parametros

                    exec_config["dataset"] = dataset_config

                # Algoritmos
                algs_elem = exec_elem.find("algoritmos")
                if algs_elem is not None:
                    algoritmos = []
                    for alg_elem in algs_elem.findall("algoritmo"):
                        if alg_elem.text is not None:
                            algoritmos.append(alg_elem.text)
                    exec_config["algoritmos"] = algoritmos

                # Execuções por algoritmo por base
                exec_per_alg_elem = exec_elem.find("execucoes_por_algoritmo_por_base")
                if exec_per_alg_elem is None:
                    exec_per_alg_elem = exec_elem.find("execucoes_por_algoritmo")

                if exec_per_alg_elem is not None and exec_per_alg_elem.text is not None:
                    exec_config["execucoes_por_algoritmo_por_base"] = int(
                        exec_per_alg_elem.text
                    )

                # Timeout
                timeout_elem = exec_elem.find("timeout")
                if timeout_elem is not None and timeout_elem.text is not None:
                    exec_config["timeout"] = int(timeout_elem.text)

                execucoes.append(exec_config)

        config["execucoes"] = execucoes
        return config

    def _generate_dataset(
        self, dataset_config: dict[str, Any]
    ) -> tuple[list[str], dict[str, Any]]:
        """Gera ou carrega dataset baseado na configuração."""
        dataset_type = dataset_config.get("tipo", "synthetic")
        params = dataset_config.get("parametros", {})

        # Checar aleatorio_total global
        aleatorio_total_global = self.batch_info.get("aleatorio_total", None)
        if (
            aleatorio_total_global is not None
            and "fully_random" not in params
            and "aleatorio_total" not in params
        ):
            # Para compatibilidade, repassa como fully_random para o gerador
            params["fully_random"] = aleatorio_total_global

        if dataset_type == "synthetic":
            from datasets.dataset_synthetic import generate_dataset_with_params

            return generate_dataset_with_params(params)

        elif dataset_type == "file":
            from datasets.dataset_file import load_dataset_with_params

            return load_dataset_with_params(params)

        elif dataset_type == "entrez":
            from datasets.dataset_entrez import fetch_dataset_with_params

            return fetch_dataset_with_params(params)

        else:
            raise ValueError(f"Tipo de dataset não suportado: {dataset_type}")

    def _execute_single_config(
        self, config: BatchConfig, exec_index: int
    ) -> dict[str, Any]:
        """Executa uma configuração individual."""
        console.print(f"\n{'='*60}")
        console.print(f"EXECUÇÃO {exec_index + 1}/{len(self.execucoes)}: {config.nome}")
        console.print(f"{'='*60}")

        start_time = time.time()
        result = {
            "config_nome": config.nome,
            "inicio": datetime.now().isoformat(),
            "algoritmos_executados": {},
            "erro": None,
            "tempo_total": 0.0,
            "num_bases": config.num_bases,
            "execucoes_por_algoritmo_por_base": config.execucoes_por_algoritmo_por_base,
        }

        try:
            # Número total de bases a serem processadas
            num_bases = config.num_bases

            # Para cada base de dados
            for base_idx in range(num_bases):
                console.print(
                    f"\n📊 Gerando base de dados {base_idx + 1}/{num_bases}..."
                )

                # Gerar/carregar dataset para esta base
                seqs, dataset_params = self._generate_dataset(config.dataset_config)

                logger.debug(f"[BatchExecutor] Base {base_idx+1}: {len(seqs)} seqs")

                if not seqs:
                    raise ValueError(f"Dataset vazio gerado para base {base_idx + 1}")

                alphabet = "".join(sorted(set("".join(seqs))))
                n, L = len(seqs), len(seqs[0])
                seed = dataset_params.get("seed")
                distancia_string_base = dataset_params.get(
                    "distancia_string_base", None
                )

                # Armazenar informações da base de dados
                if "bases_info" not in result:
                    result["bases_info"] = []

                result["bases_info"].append(
                    {
                        "base_idx": base_idx + 1,
                        "n": n,
                        "L": L,
                        "alphabet_size": len(alphabet),
                        "params": dataset_params,
                    }
                )

                # Verificar se há algoritmos configurados
                if not config.algoritmos:
                    console.print(
                        f"⚠ Nenhum algoritmo configurado para a base {base_idx + 1} - pulando"
                    )
                    continue

                # Verificar viabilidade dos algoritmos
                viable_algs = []
                for alg_name in config.algoritmos:
                    is_viable, msg = check_algorithm_feasibility(n, L, alg_name)
                    if is_viable:
                        viable_algs.append(alg_name)
                        console.print(f"✓ {alg_name}: viável")
                    else:
                        console.print(f"⚠ {alg_name}: {msg} (pulado)")

                if not viable_algs:
                    console.print(
                        f"❌ Nenhum algoritmo viável para a base {base_idx + 1}"
                    )
                    continue

                # Criar formatter para esta execução
                exec_formatter = ResultsFormatter()

                # Definir chave única para o extra_info do formatter individual
                base_key_info = f"{config.nome}_Base{base_idx+1}"
                exec_formatter.extra_info = {base_key_info: dataset_params}

                logger.debug(
                    f"[BatchExecutor] extra_info definido para {base_key_info}"
                )

                # Atualizar informações extras no formatter consolidado
                if not hasattr(self.consolidated_formatter, "extra_info"):
                    self.consolidated_formatter.extra_info = {}
                # Armazenar os params diretamente na chave da base
                base_key = f"{config.nome}_Base{base_idx+1}"
                self.consolidated_formatter.extra_info[base_key] = dataset_params

                # Executar algoritmos
                for alg_name in viable_algs:
                    console.print(
                        f"\n🔄 Executando {alg_name} para base {base_idx + 1}..."
                    )
                    if alg_name not in global_registry:
                        console.print(f"❌ Algoritmo '{alg_name}' não encontrado!")
                        continue
                    AlgClass = global_registry[alg_name]
                    execucoes_por_algoritmo = config.execucoes_por_algoritmo_por_base
                    try:
                        executions = execute_algorithm_runs(
                            alg_name,
                            AlgClass,
                            seqs,
                            alphabet,
                            execucoes_por_algoritmo,
                            None,
                            console,
                            config.timeout,
                        )

                        exec_formatter.add_algorithm_results(alg_name, executions)
                        exec_key = f"{config.nome}_Base{base_idx+1}_{alg_name}"
                        self.consolidated_formatter.add_algorithm_results(
                            exec_key, executions
                        )
                        valid_results = [
                            e
                            for e in executions
                            if "distancia" in e and e["distancia"] != float("inf")
                        ]
                        if valid_results:
                            best_exec = min(valid_results, key=lambda e: e["distancia"])
                            # Usar distancia da string base dos params
                            dist_base = dataset_params.get("distancia_string_base", "-")
                            base_key_result = f"base_{base_idx+1}"
                            if alg_name not in result["algoritmos_executados"]:
                                result["algoritmos_executados"][alg_name] = {}
                            result["algoritmos_executados"][alg_name][
                                base_key_result
                            ] = {
                                "dist": best_exec["distancia"],
                                "dist_base": dist_base,
                                "time": best_exec["tempo"],
                                "status": "sucesso",
                                "execucoes_detalhadas": executions,  # Adicionar todas as execuções
                            }
                        else:
                            error_exec = next(
                                (e for e in executions if "erro" in e), executions[0]
                            )
                            base_key_result = f"base_{base_idx+1}"
                            if alg_name not in result["algoritmos_executados"]:
                                result["algoritmos_executados"][alg_name] = {}
                            result["algoritmos_executados"][alg_name][
                                base_key_result
                            ] = {
                                "dist": float("inf"),
                                "time": error_exec["tempo"],
                                "status": "erro",
                                "erro": error_exec.get("erro", "Erro desconhecido"),
                                "execucoes_detalhadas": executions,  # Adicionar todas as execuções mesmo em caso de erro
                            }
                    except Exception as e:
                        console.print(
                            f"❌ Erro executando {alg_name} na base {base_idx+1}: {e}"
                        )
                        base_key_result = f"base_{base_idx+1}"
                        if alg_name not in result["algoritmos_executados"]:
                            result["algoritmos_executados"][alg_name] = {}
                        # Criar execução de erro para manter consistência
                        error_execution = [
                            {
                                "distancia": float("inf"),
                                "tempo": 0.0,
                                "status": "erro",
                                "erro": str(e),
                            }
                        ]
                        result["algoritmos_executados"][alg_name][base_key_result] = {
                            "dist": float("inf"),
                            "time": 0.0,
                            "status": "erro",
                            "erro": str(e),
                            "execucoes_detalhadas": error_execution,
                        }

                # Salvar relatório individual para esta base
                report_filename = f"{exec_index+1}_base{base_idx+1}_{config.nome.replace(' ', '_')}.txt"
                report_path = self.results_dir / report_filename
                exec_formatter.save_detailed_report(str(report_path))
                console.print(
                    f"📄 Relatório para base {base_idx+1} salvo: {report_filename}"
                )

        except Exception as e:
            console.print(f"❌ Erro na execução: {e}")
            logger.exception(f"Erro na execução {config.nome}")
            result["erro"] = str(e)

        finally:
            result["tempo_total"] = time.time() - start_time
            result["fim"] = datetime.now().isoformat()

        return result

    def execute_batch(self) -> dict[str, Any]:
        """Executa todas as configurações em sequência."""
        console.print("\n🚀 INICIANDO EXECUÇÃO EM LOTE")
        console.print(f"Lote: {self.batch_info.get('nome', 'Sem nome')}")
        console.print(f"Descrição: {self.batch_info.get('descricao', 'Sem descrição')}")
        console.print(f"Total de execuções: {len(self.execucoes)}")

        batch_start = time.time()
        batch_result = {
            "batch_info": self.batch_info,
            "batch_id": self.batch_id,
            "inicio": datetime.now().isoformat(),
            "execucoes": [],
            "resumo": {},
            "tempo_total": 0.0,
        }

        # Executar cada configuração
        for i, config in enumerate(self.execucoes):
            try:
                exec_result = self._execute_single_config(config, i)
                batch_result["execucoes"].append(exec_result)

                # Mostrar progresso
                elapsed = time.time() - batch_start
                console.print(
                    f"\n⏱️ Execução {i+1} concluída em {exec_result['tempo_total']:.1f}s"
                )
                console.print(f"Tempo decorrido total: {elapsed:.1f}s")

            except KeyboardInterrupt:
                console.print("\n⚠️ Execução em lote interrompida pelo usuário")
                break
            except Exception as e:
                console.print(f"\n❌ Erro fatal na execução {i+1}: {e}")
                logger.exception(f"Erro fatal na execução {config.nome}")
                continue

        # Finalizar lote
        batch_result["tempo_total"] = time.time() - batch_start
        batch_result["fim"] = datetime.now().isoformat()

        # Gerar resumo
        self._generate_batch_summary(batch_result)

        # Salvar relatório consolidado
        self._save_batch_report(batch_result)

        return batch_result

    def _generate_batch_summary(
        self,
        batch_result: dict[str, Any],
        num_bases_sinteticas=None,
        execs_por_base_str=None,
    ):
        """Gera resumo consolidado do lote."""
        total_execucoes = len(batch_result["execucoes"])
        execucoes_com_sucesso = len(
            [e for e in batch_result["execucoes"] if not e.get("erro")]
        )

        console.print("\n📋 RESUMO DO LOTE")
        console.print(f"{'='*50}")
        console.print(f"Total de execuções: {total_execucoes}")
        console.print(f"Execuções com sucesso: {execucoes_com_sucesso}")
        console.print(
            f"Taxa de sucesso: {100 * execucoes_com_sucesso / total_execucoes:.1f}%"
        )
        console.print(f"Tempo total: {batch_result['tempo_total']:.1f}s")

        # Resumo por algoritmo
        alg_stats = {}

        # Para cada configuração de execução
        for exec_result in batch_result["execucoes"]:
            if exec_result.get("erro"):
                continue

            # Para cada algoritmo em cada base
            for alg_name, base_results in exec_result.get(
                "algoritmos_executados", {}
            ).items():
                if alg_name not in alg_stats:
                    alg_stats[alg_name] = {
                        "sucessos": 0,
                        "total": 0,
                        "tempo_total": 0.0,
                    }

                # Para cada base processada neste algoritmo
                for base_key, result in base_results.items():
                    alg_stats[alg_name]["total"] += 1
                    if result.get("status") == "sucesso":
                        alg_stats[alg_name]["sucessos"] += 1
                        alg_stats[alg_name]["tempo_total"] += result.get("time", 0.0)

        console.print("\n📈 Estatísticas por algoritmo:")
        for alg_name, stats in alg_stats.items():
            taxa = 100 * stats["sucessos"] / stats["total"] if stats["total"] > 0 else 0
            tempo_med = (
                stats["tempo_total"] / stats["sucessos"] if stats["sucessos"] > 0 else 0
            )
            console.print(
                f"  {alg_name}: {stats['sucessos']}/{stats['total']} sucessos ({taxa:.1f}%), tempo médio: {tempo_med:.3f}s"
            )

        batch_result["resumo"] = {
            "total_execucoes": total_execucoes,
            "execucoes_com_sucesso": execucoes_com_sucesso,
            "taxa_sucesso": (
                100 * execucoes_com_sucesso / total_execucoes
                if total_execucoes > 0
                else 0
            ),
            "algoritmo_stats": alg_stats,
        }

    def _save_batch_report(self, batch_result: dict[str, Any]):
        """Salva relatório consolidado do lote."""
        # Relatório JSON detalhado
        json_filename = "batch_results.json"
        json_path = self.results_dir / json_filename

        # Corrigir a serialização JSON do relatório do batch, convertendo chaves de dicionários para string antes de salvar
        def convert_keys_to_str(obj):
            if isinstance(obj, dict):
                return {str(k): convert_keys_to_str(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [convert_keys_to_str(i) for i in obj]
            else:
                return obj

        # Converter todas as chaves para string e salvar uma única vez
        batch_result_str_keys = convert_keys_to_str(batch_result)
        with open(json_path, "w", encoding="utf-8") as f:
            json.dump(batch_result_str_keys, f, indent=2, ensure_ascii=False)

        console.print(f"💾 Relatório JSON salvo: {self.results_dir / json_filename}")

        # Relatório consolidado de algoritmos
        consolidated_filename = "consolidated_algorithms.txt"
        consolidated_path = self.results_dir / consolidated_filename
        self.consolidated_formatter.save_detailed_report(str(consolidated_path))
        console.print(
            f"📄 Relatório consolidado salvo: {self.results_dir / consolidated_filename}"
        )

        # Salvamento de arquivo README com metadados do lote
        self._create_batch_readme(batch_result)

    def _create_batch_readme(self, batch_result: dict) -> None:
        """
        Cria arquivo README com metadados do lote.

        Args:
            batch_result: Resultado do lote para incluir no README.
        """
        readme_path = self.results_dir / "README.md"
        with open(readme_path, "w", encoding="utf-8") as f:
            f.write(
                f"# Resultados do Lote: {self.batch_info.get('nome', 'Sem nome')}\n\n"
            )
            f.write(
                f"**Descrição**: {self.batch_info.get('descricao', 'Sem descrição')}\n\n"
            )
            f.write(
                f"**Data de execução**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
            )
            f.write(f"**ID do lote**: {self.batch_id}\n")
            f.write(f"**Tempo total**: {batch_result['tempo_total']:.1f} segundos\n")
            f.write(
                f"**Configurações**: {len(batch_result['execucoes'])} configurações\n"
            )
            f.write(
                f"**Algoritmos**: {len(self.consolidated_formatter.results)} algoritmos\n"
            )

        console.print(f"📝 README criado: {self.results_dir / 'README.md'}")
