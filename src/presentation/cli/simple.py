"""
Interface CLI simples para CSPBench

Implementa comandos básicos sem formatação rich.
"""

from pathlib import Path
from typing import Optional

import typer

from src.application.services.experiment_service import ExperimentService
from src.infrastructure import (
    CsvExporter,
    DomainAlgorithmRegistry,
    FileDatasetRepository,
    JsonExporter,
    SequentialExecutor,
)

app = typer.Typer(
    help="CSPBench v0.1.0 - Framework para Closest String Problem",
    rich_markup_mode=None,
    add_completion=False,
)

# Configuração da infraestrutura
dataset_repository = FileDatasetRepository(
    "/home/diego_grosmann/CSPBench/saved_datasets"
)
algorithm_registry = DomainAlgorithmRegistry()
executor = SequentialExecutor()
exporter = JsonExporter("./outputs")

experiment_service = ExperimentService(
    dataset_repo=dataset_repository,
    exporter=exporter,
    executor=executor,
    algo_registry=algorithm_registry,
)


@app.command()
def run(
    algorithm: str = typer.Argument(..., help="Nome do algoritmo"),
    dataset: str = typer.Argument(..., help="Caminho ou nome do dataset"),
    output: Optional[str] = typer.Option(
        None, "--output", "-o", help="Arquivo de saída"
    ),
):
    """Executa um algoritmo em um dataset."""
    try:
        typer.echo(f"Executando {algorithm} em {dataset}...")

        # Carrega dataset
        ds = dataset_repository.load(dataset)

        # Executa algoritmo
        result = experiment_service.run_single_experiment(algorithm, dataset)

        # Mostra resultado
        typer.echo(f"Resultado: {result['result_string']}")
        typer.echo(f"Distância máxima: {result['max_distance']}")
        typer.echo(f"Tempo de execução: {result['execution_time']:.4f}s")

        # Salva resultado se especificado
        if output:
            if output.endswith(".json"):
                exporter = JsonExporter(str(Path(output).parent))
                exporter.export(result, Path(output).name)
            else:
                exporter = CsvExporter(str(Path(output).parent))
                exporter.export(result, Path(output).name)
            typer.echo(f"Resultado salvo em {output}")

    except Exception as e:
        typer.echo(f"Erro: {e}", err=True)
        raise typer.Exit(1)


@app.command()
def batch(config: str = typer.Argument(..., help="Arquivo de configuração do batch")):
    """Executa batch de experimentos."""
    try:
        typer.echo(f"Executando batch: {config}")

        # TODO: Implementar carregamento de configuração
        # Por agora, um exemplo simples
        ds = dataset_repository.load("synthetic_n10_L20_noise0.1_ACTG.fasta")

        results = experiment_service.run_batch(
            {"algorithms": ["Baseline"], "datasets": [ds], "params": {}}
        )

        typer.echo(f"Batch concluído com {len(results)} resultados")
        for result in results:
            if result.get("status") == "completed":
                typer.echo(f"  {result['algorithm']}: {result['max_distance']}")
            else:
                typer.echo(f"  {result['algorithm']}: ERRO")

    except Exception as e:
        typer.echo(f"Erro: {e}", err=True)
        raise typer.Exit(1)


@app.command()
def algorithms():
    """Lista algoritmos disponíveis."""
    try:
        algos = algorithm_registry.list_available()
        typer.echo("Algoritmos disponíveis:")
        for algo in algos:
            typer.echo(f"  - {algo}")
    except Exception as e:
        typer.echo(f"Erro: {e}", err=True)
        raise typer.Exit(1)


@app.command()
def datasets():
    """Lista datasets disponíveis."""
    try:
        datasets = dataset_repository.list_available()
        typer.echo("Datasets disponíveis:")
        for ds in datasets:
            typer.echo(f"  - {ds}")
    except Exception as e:
        typer.echo(f"Erro: {e}", err=True)
        raise typer.Exit(1)


if __name__ == "__main__":
    app()
