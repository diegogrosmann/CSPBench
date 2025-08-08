import logging
from typing import Dict, Any, Optional
from datetime import datetime


class WebDisplay:
    """
    Versão de debug do WebDisplay - registra tudo que recebe nos logs
    para entender o que está sendo passado para os métodos.
    """

    def __init__(self, *args, **kwargs):
        """Inicializa o WebDisplay e registra todos os parâmetros recebidos."""
        self._logger = logging.getLogger(__name__)
        
        self._logger.info("=== WebDisplay.__init__ chamado ===")
        self._logger.info(f"Args recebidos: {args}")
        self._logger.info(f"Kwargs recebidos: {kwargs}")
        self._logger.info(f"")

    def handle_event(self, event):
        """Registra todos os detalhes do evento recebido."""
        self._logger.info("=== WebDisplay.handle_event chamado ===")
        
        # Log do evento principal
        self._logger.info(f"Event type: {type(event).__name__}")
        self._logger.info(f"Event value: {event}")