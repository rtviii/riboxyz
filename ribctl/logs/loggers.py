import os
import sys
from loguru import logger

module_dir = os.path.dirname(os.path.abspath(__file__))



def get_etl_logger():
    etl_logger = logger
    etl_logger.remove()
    etl_log_path = os.path.join(module_dir, "etl.log")
    etl_logger.add(etl_log_path, rotation="10 MB", level="DEBUG", format="{level} {time:YYYY-MM-DD HH:mm:ss} {message}")
    etl_logger.add(sys.stdout, level="DEBUG")
    return etl_logger

def get_classification_logger():
    classification_logger = logger
    classification_logger.remove()
    etl_log_path = os.path.join(module_dir, "classification.log")
    classification_logger.add(etl_log_path, rotation="10 MB", level="DEBUG", format="{level} {time:YYYY-MM-DD HH:mm:ss} {message}")
    classification_logger.add(sys.stdout, level="DEBUG")
    return classification_logger