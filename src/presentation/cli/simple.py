"""
Simple CLI interface for CSPBench

Implements basic commands without rich formatting.
"""

from pathlib import Path
from typing import Optional

import typer

from src.application.services.experiment_service import ExperimentService
from src.infrastructure import (
    CsvExporter,
    DomainAlgorithmRegistry,
    Executor,
    FileDatasetRepository,
    JsonExporter,
)

app = typer.Typer(
    help="CSPBench v0.1.0 - Framework for Closest String Problem",
    rich_markup_mode=None,
    add_completion=False,
)

# Infrastructure setup
dataset_repository = FileDatasetRepository(
    "/home/diego_grosmann/CSPBench/saved_datasets"
)
algorithm_registry = DomainAlgorithmRegistry()
executor = Executor()
exporter = JsonExporter("./outputs")

experiment_service = ExperimentService(
    dataset_repo=dataset_repository,
    exporter=exporter,
    executor=executor,
    algo_registry=algorithm_registry,
)


@app.command()
def run(
    algorithm: str = typer.Argument(..., help="Algorithm name"),
    dataset: str = typer.Argument(..., help="Dataset path or name"),
    output: Optional[str] = typer.Option(None, "--output", "-o", help="Output file"),
):
    """Execute an algorithm on a dataset."""
    try:
        typer.echo(f"Executing {algorithm} on {dataset}...")

        # Load dataset
        ds = dataset_repository.load(dataset)

        # Execute algorithm
        result = experiment_service.run_single_experiment(algorithm, dataset)

        # Show result
        typer.echo(f"Result: {result['result_string']}")
        typer.echo(f"Maximum distance: {result['max_distance']}")
        typer.echo(f"Execution time: {result['execution_time']:.4f}s")

        # Save result if specified
        if output:
            if output.endswith(".json"):
                exporter = JsonExporter(str(Path(output).parent))
                exporter.export(result, Path(output).name)
            else:
                exporter = CsvExporter(str(Path(output).parent))
                exporter.export([result], Path(output).name)
            typer.echo(f"Result saved to {output}")

    except Exception as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(1)


@app.command()
def batch(config: str = typer.Argument(..., help="Batch configuration file")):
    """Execute batch of experiments."""
    try:
        typer.echo(f"Executing batch: {config}")

        # TODO: Implement configuration loading
        # For now, a simple example
        ds = dataset_repository.load("synthetic_n10_L20_noise0.1_ACTG.fasta")

        results = experiment_service.run_batch(
            {"algorithms": ["Baseline"], "datasets": [ds], "params": {}}
        )

        typer.echo(f"Batch completed with {len(results)} results")
        for result in results:
            if result.get("status") == "completed":
                typer.echo(f"  {result['algorithm']}: {result['max_distance']}")
            else:
                typer.echo(f"  {result['algorithm']}: ERROR")

    except Exception as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(1)


@app.command()
def algorithms():
    """List available algorithms."""
    try:
        algos = algorithm_registry.list_available()
        typer.echo("Available algorithms:")
        for algo in algos:
            typer.echo(f"  - {algo}")
    except Exception as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(1)


@app.command()
def datasets():
    """List available datasets."""
    try:
        datasets = dataset_repository.list_available()
        typer.echo("Available datasets:")
        for ds in datasets:
            typer.echo(f"  - {ds}")
    except Exception as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(1)


if __name__ == "__main__":
    app()
