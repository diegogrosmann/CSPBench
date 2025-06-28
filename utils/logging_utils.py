# logging_utils.py
import logging

def setup_logging(debug_mode: bool = False, log_file: str = 'debug.log'):
    if debug_mode:
        logging.basicConfig(
            level=logging.DEBUG,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            filename=log_file,
            filemode='w'
        )
        print(f"Debug mode enabled. Log saved to {log_file}")
    else:
        logging.disable(logging.CRITICAL)
        print("Debug mode disabled.")
