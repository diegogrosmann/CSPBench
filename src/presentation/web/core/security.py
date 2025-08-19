"""
Security utilities and validators for web interface.
"""

import logging
from pathlib import Path
from typing import Dict

logger = logging.getLogger(__name__)


class SecurityConfig:
    """Security configuration constants."""

    MAX_DATASET_SIZE = 50 * 1024 * 1024  # 50MB
    ALLOWED_EXTENSIONS = {".fasta", ".fa", ".txt", ".seq"}
    MAX_SESSIONS_PER_IP = 10
    SESSION_TIMEOUT = 7200  # 2 hours
    MAX_FILENAME_LENGTH = 255


class SecurityValidator:
    """Security validation utilities."""

    @staticmethod
    def sanitize_filename(filename: str) -> str:
        """Sanitize filename following security guidelines."""
        if not filename:
            return "uploaded_file.fasta"

        # Allow only safe characters
        safe_chars = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-"
        sanitized = "".join(c for c in filename if c in safe_chars)

        # Limit size
        max_length = SecurityConfig.MAX_FILENAME_LENGTH
        if len(sanitized) > max_length:
            name, ext = (
                sanitized.rsplit(".", 1) if "." in sanitized else (sanitized, "")
            )
            sanitized = name[: max_length - len(ext) - 1] + "." + ext

        return sanitized or "uploaded_file.fasta"

    @staticmethod
    def generate_unique_filename(base_filename: str, directory: Path) -> str:
        """Generate a unique filename by adding a counter if file already exists."""
        filename = base_filename
        counter = 1

        # Extract name and extension
        if "." in filename:
            name_part, ext_part = filename.rsplit(".", 1)
        else:
            name_part, ext_part = filename, ""

        # Keep trying until we find a unique name
        while (directory / filename).exists():
            if ext_part:
                filename = f"{name_part}_{counter}.{ext_part}"
            else:
                filename = f"{name_part}_{counter}"
            counter += 1

        return filename

    @staticmethod
    def validate_dataset_content(content: str) -> bool:
        """Validate dataset content following security guidelines."""
        if not content or len(content) > SecurityConfig.MAX_DATASET_SIZE:
            return False

        # Check if it's valid FASTA format
        lines = content.strip().split("\n")
        has_header = any(line.startswith(">") for line in lines)
        has_sequence = any(line and not line.startswith(">") for line in lines)

        return has_header and has_sequence

    @staticmethod
    def validate_algorithm_parameters(params: Dict, default_params: Dict) -> Dict:
        """Validate algorithm parameters against schema."""
        validated = {}

        for key, value in params.items():
            if key in default_params:
                # Basic type validation
                default_type = type(default_params[key])
                if default_type == type(None):
                    validated[key] = value
                else:
                    try:
                        if default_type == bool:
                            validated[key] = (
                                bool(value)
                                if isinstance(value, bool)
                                else str(value).lower() == "true"
                            )
                        elif default_type == int:
                            validated[key] = int(value) if value != "" else None
                        elif default_type == float:
                            validated[key] = float(value) if value != "" else None
                        else:
                            validated[key] = default_type(value)
                    except (ValueError, TypeError):
                        logger.warning(
                            f"Invalid parameter {key}: {value}, using default"
                        )
                        validated[key] = default_params[key]
            else:
                logger.warning(f"Unknown parameter {key} ignored")

        # Fill missing parameters with defaults
        for key, value in default_params.items():
            if key not in validated:
                validated[key] = value

        return validated

    @staticmethod
    def validate_file_extension(filename: str) -> bool:
        """Validate file extension against allowed list."""
        if not filename:
            return False

        extension = Path(filename).suffix.lower()
        return extension in SecurityConfig.ALLOWED_EXTENSIONS


# Convenience functions for easy import
def sanitize_filename(filename: str) -> str:
    """Sanitize filename using SecurityValidator."""
    return SecurityValidator.sanitize_filename(filename)
