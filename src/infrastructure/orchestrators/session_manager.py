"""
Session Manager

Responsible for creating and managing session folders organized by datetime
for CSPBench logs and results.
"""

import shutil
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional


class SessionManager:
    """
    Manages sessions organized by datetime for logs and results.

    Creates folder structure in the format:
    - logs/YYYYMMDD_HHMMSS/cspbench.log
    - results/YYYYMMDD_HHMMSS/results.json
    """

    def __init__(self, config: Dict):
        """
        Initialize the session manager.

        Args:
            config: System configuration containing base paths
        """
        self.config = config
        self.session_format = config["infrastructure"]["logging"][
            "session_folder_format"
        ]
        self.create_session_folders = config["infrastructure"]["logging"][
            "create_session_folders"
        ]

        # Base directories for logs and results
        self.logs_base_dir = Path(config["infrastructure"]["logging"]["base_log_dir"])
        self.results_base_dir = Path(
            config["infrastructure"]["result"]["base_result_dir"]
        )

        # Current session
        self.current_session: Optional[str] = None

    def create_session(self) -> str:
        """
        Create a new session with current timestamp.

        Returns:
            str: Created session name (format: YYYYMMDD_HHMMSS)
        """
        if not self.create_session_folders:
            return "default"

        # Generate session name with current timestamp
        now = datetime.now()
        session_name = now.strftime(self.session_format)

        # Create session directories
        session_log_dir = self.logs_base_dir / session_name
        session_result_dir = self.results_base_dir / session_name

        session_log_dir.mkdir(parents=True, exist_ok=True)
        session_result_dir.mkdir(parents=True, exist_ok=True)

        # Set as current session
        self.current_session = session_name

        return session_name

    def get_log_path(self, session_name: Optional[str] = None) -> Path:
        """
        Get the log file path for a session.

        Args:
            session_name: Session name. If None, uses current session.

        Returns:
            Path: Complete path to log file
        """
        session_name = session_name or self.current_session
        if not session_name:
            raise ValueError("No active session")

        session_log_dir = self.logs_base_dir / session_name
        return session_log_dir / "cspbench.log"

    def get_result_path(self, session_name: Optional[str] = None) -> Path:
        """
        Get the result file path for a session.

        Args:
            session_name: Session name. If None, uses current session.

        Returns:
            Path: Complete path to result file
        """
        session_name = session_name or self.current_session
        if not session_name:
            raise ValueError("No active session")

        session_result_dir = self.results_base_dir / session_name
        return session_result_dir / "results.json"

    def get_result_dir(self) -> Path:
        """
        Get the results directory of current session.

        Returns:
            Path: Results directory of current session
        """
        if not self.current_session:
            raise ValueError("No active session")

        return self.results_base_dir / self.current_session

    def get_result_dir(self) -> Path:
        """
        Obtém o diretório de resultados da sessão atual.

        Returns:
            Path: Diretório de resultados da sessão atual
        """
        if not self.current_session:
            raise ValueError("Nenhuma sessão ativa")

        return self.results_base_dir / self.current_session

    def list_sessions(self) -> Dict[str, Dict]:
        """
        List all available sessions with their information.

        Returns:
            Dict containing information for each session:
            {
                "session_name": {
                    "created": datetime,
                    "logs": bool,
                    "results": bool
                }
            }
        """
        sessions = {}

        # Check log directories
        if self.logs_base_dir.exists():
            for log_dir in self.logs_base_dir.iterdir():
                if log_dir.is_dir() and self._is_valid_session_name(log_dir.name):
                    sessions[log_dir.name] = {
                        "created": self._parse_session_datetime(log_dir.name),
                        "logs": True,
                        "results": False,
                    }

        # Check result directories
        if self.results_base_dir.exists():
            for result_dir in self.results_base_dir.iterdir():
                if result_dir.is_dir() and self._is_valid_session_name(result_dir.name):
                    if result_dir.name in sessions:
                        sessions[result_dir.name]["results"] = True
                    else:
                        sessions[result_dir.name] = {
                            "created": self._parse_session_datetime(result_dir.name),
                            "logs": False,
                            "results": True,
                        }

        return sessions

    def cleanup_old_sessions(self, keep_last: int = 10) -> None:
        """
        Remove old sessions, keeping only the most recent ones.

        Args:
            keep_last: Number of most recent sessions to keep
        """
        sessions = self.list_sessions()
        if len(sessions) <= keep_last:
            return

        # Sort by creation date
        sorted_sessions = sorted(
            sessions.items(), key=lambda x: x[1]["created"], reverse=True
        )

        # Remove oldest sessions
        for session_name, _ in sorted_sessions[keep_last:]:
            self._remove_session(session_name)

    def cleanup_old_sessions(self, keep_last: int = 10) -> None:
        """
        Remove sessões antigas, mantendo apenas as mais recentes.

        Args:
            keep_last: Número de sessões mais recentes para manter
        """
        sessions = self.list_sessions()

        if len(sessions) <= keep_last:
            return  # Nada para remover

        # Ordenar por data de criação (mais antigas primeiro)
        sorted_sessions = sorted(sessions.items(), key=lambda x: x[1]["created"])

        # Determinar quais sessões remover
        sessions_to_remove = sorted_sessions[:-keep_last]  # Todas exceto as últimas N

        for session_name, _ in sessions_to_remove:
            self._remove_session(session_name)

    def _is_valid_session_name(self, name: str) -> bool:
        """Verifica se o nome segue o formato de sessão (YYYYMMDD_HHMMSS)."""
        try:
            self._parse_session_datetime(name)
            return True
        except ValueError:
            return False

    def _parse_session_datetime(self, session_name: str) -> datetime:
        """Converte nome da sessão em datetime."""
        return datetime.strptime(session_name, self.session_format)

    def _remove_session(self, session_name: str) -> None:
        """Remove uma sessão específica (logs e resultados)."""
        # Remover diretório de logs
        log_session_dir = self.logs_base_dir / session_name
        if log_session_dir.exists():
            shutil.rmtree(log_session_dir)

        # Remover diretório de resultados
        result_session_dir = self.results_base_dir / session_name
        if result_session_dir.exists():
            shutil.rmtree(result_session_dir)
