import os
import sys
from loguru import logger

module_dir = os.path.dirname(os.path.abspath(__file__))


def get_updates_logger():
    updates_logger = logger
    updates_logger.remove()
    updates_log_path = os.path.join(module_dir, "db_updates.log")
    updates_logger.add(updates_log_path, rotation="10 MB", level="DEBUG", format="{level} {time:YYYY-MM-DD HH:mm:ss} {message}")
    updates_logger.add(sys.stdout, level="DEBUG")
    return updates_logger
    
def get_computation_logger():
    updates_logger = logger
    updates_logger.remove()
    updates_log_path = os.path.join(module_dir, "compute_updates.log")
    updates_logger.add(updates_log_path, rotation="10 MB", level="DEBUG", format="{level} {time:YYYY-MM-DD HH:mm:ss} {message}")
    updates_logger.add(sys.stdout, level="DEBUG")
    return updates_logger