import logging
from pathlib import Path
from utils.config import DEBUG_DEFAULT, safe_input
from src.console_manager import console

def setup_logging(base_name: str, silent: bool = False) -> None:
    if silent:
        debug_mode = DEBUG_DEFAULT
    else:
        try:
            debug_input = safe_input("Habilitar modo debug? [s/N]: ", DEBUG_DEFAULT).lower()
        except KeyboardInterrupt:
            console.print("\nOperação cancelada pelo usuário.")
            import sys
            sys.exit(0)
        debug_mode = debug_input if debug_input else DEBUG_DEFAULT
    Path("logs").mkdir(exist_ok=True)
    log_file = Path("logs") / f"{base_name}.log"
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    if debug_mode == 's':
        level = logging.DEBUG
        if not silent:
            console.print(f"Debug ON -> {log_file}")
    else:
        level = logging.CRITICAL
    logging.basicConfig(
        level=level,
        filename=str(log_file),
        filemode='w',
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
