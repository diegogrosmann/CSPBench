"""
Gerenciador de Sessões

Responsável por criar e gerenciar pastas de sessão organizadas por datetime
para logs e resultados do CSPBench.
"""

import shutil
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional


class SessionManager:
    """
    Gerencia sessões organizadas por datetime para logs e resultados.

    Cria estrutura de pastas no formato:
    - logs/YYYYMMDD_HHMMSS/cspbench.log
    - results/YYYYMMDD_HHMMSS/results.json
    """

    def __init__(self, config: Dict):
        """
        Inicializa o gerenciador de sessões.

        Args:
            config: Configuração do sistema contendo paths base
        """
        self.config = config
        self.session_format = config["infrastructure"]["logging"][
            "session_folder_format"
        ]
        self.create_session_folders = config["infrastructure"]["logging"][
            "create_session_folders"
        ]

        # Diretórios base para logs e resultados
        self.logs_base_dir = Path(config["infrastructure"]["logging"]["base_log_dir"])
        self.results_base_dir = Path(
            config["infrastructure"]["result"]["base_result_dir"]
        )

        # Sessão atual
        self.current_session: Optional[str] = None

    def create_session(self) -> str:
        """
        Cria uma nova sessão com timestamp atual.

        Returns:
            str: Nome da sessão criada (formato: YYYYMMDD_HHMMSS)
        """
        if not self.create_session_folders:
            return ""

        # Gerar nome da sessão com timestamp atual
        now = datetime.now()
        session_name = now.strftime(self.session_format)

        # Criar diretórios da sessão
        session_log_dir = self.logs_base_dir / session_name
        session_result_dir = self.results_base_dir / session_name

        session_log_dir.mkdir(parents=True, exist_ok=True)
        session_result_dir.mkdir(parents=True, exist_ok=True)

        # Definir como sessão atual
        self.current_session = session_name

        return session_name

    def get_log_path(self, session_name: Optional[str] = None) -> Path:
        """
        Obtém o caminho do arquivo de log para uma sessão.

        Args:
            session_name: Nome da sessão. Se None, usa a sessão atual.

        Returns:
            Path: Caminho completo para o arquivo de log
        """
        session_name = session_name or self.current_session
        if not session_name:
            raise ValueError("Nenhuma sessão ativa ou especificada")

        session_log_dir = self.logs_base_dir / session_name
        return session_log_dir / "cspbench.log"

    def get_result_path(self, session_name: Optional[str] = None) -> Path:
        """
        Obtém o caminho do arquivo de resultado para uma sessão.

        Args:
            session_name: Nome da sessão. Se None, usa a sessão atual.

        Returns:
            Path: Caminho completo para o arquivo de resultado
        """
        session_name = session_name or self.current_session
        if not session_name:
            raise ValueError("Nenhuma sessão ativa ou especificada")

        session_result_dir = self.results_base_dir / session_name
        return session_result_dir / "results.json"

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
        Lista todas as sessões disponíveis com suas informações.

        Returns:
            Dict contendo informações de cada sessão:
            {
                "session_name": {
                    "created": datetime,
                    "logs": bool,
                    "results": bool
                }
            }
        """
        sessions = {}

        # Verificar diretórios de logs
        if self.logs_base_dir.exists():
            for session_dir in self.logs_base_dir.iterdir():
                if session_dir.is_dir() and self._is_valid_session_name(
                    session_dir.name
                ):
                    session_name = session_dir.name
                    if session_name not in sessions:
                        sessions[session_name] = {
                            "created": self._parse_session_datetime(session_name),
                            "logs": False,
                            "results": False,
                        }

                    # Verificar se tem arquivo de log
                    log_file = session_dir / "cspbench.log"
                    sessions[session_name]["logs"] = log_file.exists()

        # Verificar diretórios de resultados
        if self.results_base_dir.exists():
            for session_dir in self.results_base_dir.iterdir():
                if session_dir.is_dir() and self._is_valid_session_name(
                    session_dir.name
                ):
                    session_name = session_dir.name
                    if session_name not in sessions:
                        sessions[session_name] = {
                            "created": self._parse_session_datetime(session_name),
                            "logs": False,
                            "results": False,
                        }

                    # Verificar se tem arquivo de resultado
                    result_file = session_dir / "results.json"
                    sessions[session_name]["results"] = result_file.exists()

        return sessions

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
