import sys
from loguru import logger

def configure_logging():
    # Remove default handler
    logger.remove()
    
    # Add our custom handler
    logger.add(
        sys.stdout,
        format="<blue>{time}</blue> | <level>{level: <8}</level> | <cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <white>{message}</white>",
        colorize=True,
        diagnose=False,
        backtrace=False,
        enqueue=True,
        catch=True  # This catches any error during logging
    )

    # Optionally, you can also set a specific error format
    logger.add(
        sys.stderr,
        format="<blue>{time}</blue> | <level>{level: <8}</level> | <cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <red>{message}</red>",
        colorize=True,
        diagnose=False,
        backtrace=False,
        enqueue=True,
        catch=True,
        level="ERROR"
    )