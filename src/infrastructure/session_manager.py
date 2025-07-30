"""
CSPBench Session Manager

Manages creation of folders and paths for logs and results organized by session.
"""

import os
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional

from src.infrastructure.logging_config import get_logger


class SessionManager:
    """Manages execution sessions with folders organized by datetime."""

    def __init__(self, config: Dict[str, Any]):
        """
        Initialize the session manager.

        Args:
            config: System configuration
        """
        self._config = config
        self._logger = get_logger(__name__)
        self._session_folder: Optional[str] = None

        # Extrair configurações de logging e resultados
        logging_config = config.get("infrastructure", {}).get("logging", {})
        result_config = config.get("infrastructure", {}).get("result", {})

        self._base_log_dir = logging_config.get("base_log_dir", "./outputs/logs")
        self._base_result_dir = result_config.get(
            "base_result_dir", "./outputs/results"
        )
        self._session_folder_format = logging_config.get(
            "session_folder_format", "%Y%m%d_%H%M%S"
        )
        self._log_filename = logging_config.get("log_filename", "cspbench.log")
        self._result_filename = result_config.get("result_filename", "results.json")
        self._create_session_folders = logging_config.get(
            "create_session_folders", True
        )

    def create_session(self, session_name: Optional[str] = None) -> str:
        """
        Create a new session with specific folder.

        Args:
            session_name: Custom session name (optional)

        Returns:
            str: Created session folder name
        """
        if not self._create_session_folders:
            # If should not create session folders, use base folders
            self._session_folder = ""
            return ""

        if session_name:
            self._session_folder = session_name
        else:
            # Generate automatic name with timestamp
            timestamp = datetime.now().strftime(self._session_folder_format)
            self._session_folder = timestamp

        # Create directories if they don't exist
        log_dir = self.get_log_dir()
        result_dir = self.get_result_dir()

        Path(log_dir).mkdir(parents=True, exist_ok=True)
        Path(result_dir).mkdir(parents=True, exist_ok=True)

        self._logger.info(f"Session created: {self._session_folder}")
        self._logger.info(f"Log directory: {log_dir}")
        self._logger.info(f"Results directory: {result_dir}")

        return self._session_folder

    def get_session_folder(self) -> str:
        """Return the current session folder name."""
        return self._session_folder or ""

    def get_log_dir(self) -> str:
        """Return the complete directory for current session logs."""
        if not self._session_folder:
            return self._base_log_dir
        return os.path.join(self._base_log_dir, self._session_folder)

    def get_result_dir(self) -> str:
        """Return the complete directory for current session results."""
        if not self._session_folder:
            return self._base_result_dir
        return os.path.join(self._base_result_dir, self._session_folder)

    def get_log_dir(self) -> str:
        """Retorna o diretório completo para logs da sessão atual."""
        if not self._session_folder:
            return self._base_log_dir
        return os.path.join(self._base_log_dir, self._session_folder)

    def get_result_dir(self) -> str:
        """Retorna o diretório completo para resultados da sessão atual."""
        if not self._session_folder:
            return self._base_result_dir
        return os.path.join(self._base_result_dir, self._session_folder)

    def get_log_path(self, filename: Optional[str] = None) -> str:
        """
        Return the complete path for log file.

        Args:
            filename: File name (uses default if not specified)

        Returns:
            str: Complete path to log file
        """
        log_filename = filename or self._log_filename
        return os.path.join(self.get_log_dir(), log_filename)

    def get_result_path(self, filename: Optional[str] = None) -> str:
        """
        Return the complete path for result file.

        Args:
            filename: File name (uses default if not specified)

        Returns:
            str: Complete path to result file
        """
        result_filename = filename or self._result_filename
        return os.path.join(self.get_result_dir(), result_filename)

    def cleanup_old_sessions(self, keep_last: int = 10) -> None:
        """
        Remove sessões antigas, mantendo apenas as mais recentes.

        Args:
            keep_last: Número de sessões mais recentes para manter
        """
        try:
            # Listar pastas de sessão em logs
            log_base = Path(self._base_log_dir)
            if log_base.exists():
                log_sessions = [d for d in log_base.iterdir() if d.is_dir()]
                log_sessions.sort(key=lambda x: x.stat().st_mtime, reverse=True)

                # Remover sessões antigas de logs
                for session_dir in log_sessions[keep_last:]:
                    self._logger.info(
                        f"Removendo sessão antiga de logs: {session_dir.name}"
                    )
                    self._remove_directory(session_dir)

            # Listar pastas de sessão em resultados
            result_base = Path(self._base_result_dir)
            if result_base.exists():
                result_sessions = [d for d in result_base.iterdir() if d.is_dir()]
                result_sessions.sort(key=lambda x: x.stat().st_mtime, reverse=True)

                # Remover sessões antigas de resultados
                for session_dir in result_sessions[keep_last:]:
                    self._logger.info(
                        f"Removendo sessão antiga de resultados: {session_dir.name}"
                    )
                    self._remove_directory(session_dir)

        except Exception as e:
            self._logger.warning(f"Erro ao limpar sessões antigas: {e}")

    def _remove_directory(self, directory: Path) -> None:
        """Remove diretório e todo seu conteúdo."""
        import shutil

        try:
            shutil.rmtree(directory)
        except Exception as e:
            self._logger.warning(f"Erro ao remover diretório {directory}: {e}")

    def list_sessions(self) -> Dict[str, Dict[str, Any]]:
        """
        Lista todas as sessões disponíveis.

        Returns:
            Dict: Informações das sessões {nome: {logs: bool, results: bool, created: datetime}}
        """
        sessions = {}

        # Verificar sessões em logs
        log_base = Path(self._base_log_dir)
        if log_base.exists():
            for session_dir in log_base.iterdir():
                if session_dir.is_dir():
                    sessions[session_dir.name] = {
                        "logs": True,
                        "results": False,
                        "created": datetime.fromtimestamp(session_dir.stat().st_mtime),
                    }

        # Verificar sessões em resultados
        result_base = Path(self._base_result_dir)
        if result_base.exists():
            for session_dir in result_base.iterdir():
                if session_dir.is_dir():
                    if session_dir.name in sessions:
                        sessions[session_dir.name]["results"] = True
                    else:
                        sessions[session_dir.name] = {
                            "logs": False,
                            "results": True,
                            "created": datetime.fromtimestamp(
                                session_dir.stat().st_mtime
                            ),
                        }

        return sessions
