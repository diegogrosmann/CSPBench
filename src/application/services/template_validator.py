"""
Template Validator for CSPBench

This module validates that batch configurations contain ONLY parameters
documented in the TEMPLATE.yaml file.
"""

from pathlib import Path
from typing import Any, Dict, List, Union

import yaml

from src.domain.errors import BatchConfigurationError


class TemplateValidator:
    """
    Validates batch configurations against the TEMPLATE.yaml specification.

    Ensures that ONLY documented parameters are used in configurations.
    """

    # Define valid parameters for each section based on TEMPLATE.yaml
    VALID_SECTIONS = {
        "metadata": {
            "required": {"name", "description", "author", "version", "creation_date"},
            "optional": {"tags"},
        },
        "infrastructure": {"required": set(), "optional": {"history", "result"}},
        "datasets": {
            "required": set(),  # List validation handled separately
            "optional": set(),
        },
        "algorithms": {
            "required": set(),  # List validation handled separately
            "optional": set(),
        },
        "task": {"required": {"type"}, "optional": set()},
        "execution": {
            "required": set(),  # Structure validation handled separately
            "optional": set(),
        },
        "optimization": {
            "required": set(),  # Structure validation handled separately
            "optional": set(),
        },
        "sensitivity": {
            "required": set(),  # Structure validation handled separately
            "optional": set(),
        },
        "export": {
            "required": set(),
            "optional": {
                "enabled",
                "destination",
                "formats",
                "format_options",
                "directory_structure",
                "include",
            },
        },
        "plots": {
            "required": set(),
            "optional": {
                "enabled",
                "convergence",
                "comparison",
                "boxplots",
                "scatter",
                "heatmap",
                "runtime",
                "success_rate",
                "optimization_history",
                "parameter_importance",
                "parallel_coordinate",
                "sensitivity_indices",
                "morris_trajectories",
                "interaction_effects",
                "formats",
            },
        },
        "monitoring": {
            "required": set(),
            "optional": {"enabled", "interface", "update_interval"},
        },
        "resources": {
            "required": set(),
            "optional": {"cpu", "memory", "parallel", "timeouts"},
        },
        "logging": {"required": set(), "optional": {"level", "output"}},
        "system": {"required": set(), "optional": {"reproducibility"}},
    }

    # Valid dataset types
    VALID_DATASET_TYPES = {"synthetic", "file", "entrez"}

    # Valid task types
    VALID_TASK_TYPES = {"execution", "optimization", "sensitivity"}

    # Valid algorithm names (from the algorithms registry)
    VALID_ALGORITHM_NAMES = {"Baseline", "BLF-GA", "CSC", "H³-CSP", "DP-CSP"}

    # Valid sensitivity methods
    VALID_SENSITIVITY_METHODS = {"morris", "sobol", "fast", "delta"}

    # Valid monitoring interfaces
    VALID_MONITORING_INTERFACES = {"simple", "tui"}

    # Valid logging levels
    VALID_LOGGING_LEVELS = {"DEBUG", "INFO", "WARNING", "ERROR"}

    @classmethod
    def validate_config(
        cls, config: Dict[str, Any], config_path: str = ""
    ) -> List[str]:
        """
        Validate a complete configuration against TEMPLATE.yaml.

        Args:
            config: Configuration dictionary to validate
            config_path: Path to config file (for error messages)

        Returns:
            List of validation errors (empty if valid)
        """
        errors = []

        # Validate top-level sections
        errors.extend(cls._validate_sections(config, config_path))

        # Validate specific sections
        if "metadata" in config:
            errors.extend(cls._validate_metadata(config["metadata"], config_path))
        else:
            errors.append(f"{config_path}: Missing required section 'metadata'")

        if "datasets" in config:
            errors.extend(cls._validate_datasets(config["datasets"], config_path))

        if "algorithms" in config:
            errors.extend(cls._validate_algorithms(config["algorithms"], config_path))

        if "task" in config:
            errors.extend(cls._validate_task(config["task"], config_path))
        else:
            errors.append(f"{config_path}: Missing required section 'task'")

        # Validate task-specific sections
        task_type = config.get("task", {}).get("type")
        if task_type == "execution" and "execution" not in config:
            errors.append(
                f"{config_path}: Task type 'execution' requires 'execution' section"
            )
        elif task_type == "optimization" and "optimization" not in config:
            errors.append(
                f"{config_path}: Task type 'optimization' requires 'optimization' section"
            )
        elif task_type == "sensitivity" and "sensitivity" not in config:
            errors.append(
                f"{config_path}: Task type 'sensitivity' requires 'sensitivity' section"
            )

        # Validate optional sections if present
        for section_name in [
            "infrastructure",
            "export",
            "plots",
            "monitoring",
            "resources",
            "logging",
            "system",
        ]:
            if section_name in config:
                errors.extend(
                    cls._validate_section_fields(
                        config[section_name], section_name, config_path
                    )
                )

        return errors

    @classmethod
    def _validate_sections(cls, config: Dict[str, Any], config_path: str) -> List[str]:
        """Validate that only known sections are present."""
        errors = []
        valid_sections = set(cls.VALID_SECTIONS.keys())
        present_sections = set(config.keys())

        invalid_sections = present_sections - valid_sections
        if invalid_sections:
            errors.append(
                f"{config_path}: Invalid sections found: {invalid_sections}. "
                f"Valid sections: {valid_sections}"
            )

        return errors

    @classmethod
    def _validate_section_fields(
        cls, section_data: Dict[str, Any], section_name: str, config_path: str
    ) -> List[str]:
        """Validate fields within a section."""
        errors = []

        if section_name not in cls.VALID_SECTIONS:
            return errors

        section_spec = cls.VALID_SECTIONS[section_name]
        required_fields = section_spec["required"]
        optional_fields = section_spec["optional"]
        valid_fields = required_fields | optional_fields

        present_fields = set(section_data.keys())

        # Check for missing required fields
        missing_required = required_fields - present_fields
        if missing_required:
            errors.append(
                f"{config_path}: Section '{section_name}' missing required fields: {missing_required}"
            )

        # Check for invalid fields
        invalid_fields = present_fields - valid_fields
        if invalid_fields:
            errors.append(
                f"{config_path}: Section '{section_name}' contains invalid fields: {invalid_fields}. "
                f"Valid fields: {valid_fields}"
            )

        return errors

    @classmethod
    def _validate_metadata(
        cls, metadata: Dict[str, Any], config_path: str
    ) -> List[str]:
        """Validate metadata section."""
        return cls._validate_section_fields(metadata, "metadata", config_path)

    @classmethod
    def _validate_datasets(
        cls, datasets: List[Dict[str, Any]], config_path: str
    ) -> List[str]:
        """Validate datasets section."""
        errors = []

        if not isinstance(datasets, list):
            errors.append(f"{config_path}: 'datasets' must be a list")
            return errors

        dataset_ids = set()
        for i, dataset in enumerate(datasets):
            # Check required fields
            required_fields = {"id", "name", "type", "parameters"}
            missing_fields = required_fields - set(dataset.keys())
            if missing_fields:
                errors.append(
                    f"{config_path}: Dataset {i} missing required fields: {missing_fields}"
                )
                continue

            # Check dataset ID uniqueness
            dataset_id = dataset.get("id")
            if dataset_id in dataset_ids:
                errors.append(f"{config_path}: Duplicate dataset ID: {dataset_id}")
            dataset_ids.add(dataset_id)

            # Validate dataset type
            dataset_type = dataset.get("type")
            if dataset_type not in cls.VALID_DATASET_TYPES:
                errors.append(
                    f"{config_path}: Invalid dataset type '{dataset_type}' in dataset {i}. "
                    f"Valid types: {cls.VALID_DATASET_TYPES}"
                )

        return errors

    @classmethod
    def _validate_algorithms(
        cls, algorithms: List[Dict[str, Any]], config_path: str
    ) -> List[str]:
        """Validate algorithms section."""
        errors = []

        if not isinstance(algorithms, list):
            errors.append(f"{config_path}: 'algorithms' must be a list")
            return errors

        algorithm_ids = set()
        for i, algorithm in enumerate(algorithms):
            # Check required fields
            required_fields = {"id", "name", "algorithms", "algorithm_params"}
            missing_fields = required_fields - set(algorithm.keys())
            if missing_fields:
                errors.append(
                    f"{config_path}: Algorithm config {i} missing required fields: {missing_fields}"
                )
                continue

            # Check algorithm ID uniqueness
            alg_id = algorithm.get("id")
            if alg_id in algorithm_ids:
                errors.append(f"{config_path}: Duplicate algorithm config ID: {alg_id}")
            algorithm_ids.add(alg_id)

            # Validate algorithm names
            algorithm_names = algorithm.get("algorithms", [])
            invalid_names = set(algorithm_names) - cls.VALID_ALGORITHM_NAMES
            if invalid_names:
                errors.append(
                    f"{config_path}: Invalid algorithm names in config {i}: {invalid_names}. "
                    f"Valid names: {cls.VALID_ALGORITHM_NAMES}"
                )

        return errors

    @classmethod
    def _validate_task(cls, task: Dict[str, Any], config_path: str) -> List[str]:
        """Validate task section."""
        errors = []

        # Validate task type
        task_type = task.get("type")
        if not task_type:
            errors.append(f"{config_path}: Task section missing required field 'type'")
        elif task_type not in cls.VALID_TASK_TYPES:
            errors.append(
                f"{config_path}: Invalid task type '{task_type}'. "
                f"Valid types: {cls.VALID_TASK_TYPES}"
            )

        # Check for unexpected fields
        valid_fields = {"type"}
        invalid_fields = set(task.keys()) - valid_fields
        if invalid_fields:
            errors.append(
                f"{config_path}: Task section contains invalid fields: {invalid_fields}. "
                f"Valid fields: {valid_fields}"
            )

        return errors

    @classmethod
    def validate_file(cls, config_path: Union[str, Path]) -> List[str]:
        """
        Validate a configuration file.

        Args:
            config_path: Path to the configuration file

        Returns:
            List of validation errors (empty if valid)
        """
        try:
            with open(config_path, encoding="utf-8") as file:
                config = yaml.safe_load(file)

            if not config:
                return [f"{config_path}: Empty configuration file"]

            return cls.validate_config(config, str(config_path))

        except FileNotFoundError:
            return [f"{config_path}: Configuration file not found"]
        except yaml.YAMLError as e:
            return [f"{config_path}: Invalid YAML syntax: {e}"]
        except Exception as e:
            return [f"{config_path}: Error loading configuration: {e}"]

    @classmethod
    def is_valid(cls, config: Union[Dict[str, Any], str, Path]) -> bool:
        """
        Check if a configuration is valid.

        Args:
            config: Configuration dict or path to config file

        Returns:
            True if valid, False otherwise
        """
        if isinstance(config, (str, Path)):
            errors = cls.validate_file(config)
        else:
            errors = cls.validate_config(config)

        return len(errors) == 0

    @classmethod
    def validate_and_raise(cls, config: Union[Dict[str, Any], str, Path]) -> None:
        """
        Validate configuration and raise exception if invalid.

        Args:
            config: Configuration dict or path to config file

        Raises:
            BatchConfigurationError: If configuration is invalid
        """
        if isinstance(config, (str, Path)):
            errors = cls.validate_file(config)
        else:
            errors = cls.validate_config(config)

        if errors:
            error_message = "Configuration validation failed:\n" + "\n".join(errors)
            raise BatchConfigurationError(error_message)


def validate_config_file(config_path: Union[str, Path]) -> List[str]:
    """Convenience function to validate a configuration file."""
    return TemplateValidator.validate_file(config_path)


def is_config_valid(config_path: Union[str, Path]) -> bool:
    """Convenience function to check if a configuration file is valid."""
    return TemplateValidator.is_valid(config_path)
