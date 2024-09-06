import logging
import os
import sys
from pathlib import Path
from typing import IO, Dict

LOGFILES: Dict[str, IO] = {}

def pytest_configure(config):
    """Configures logging for pytest workers when using pytest-xdist."""
    
    # Get the worker ID from the environment variable
    worker_id = os.environ.get("PYTEST_XDIST_WORKER")

    # Proceed if worker_id is defined and not already logged
    if worker_id and worker_id not in LOGFILES:
        logfile = Path(f"build/tests/{worker_id}.log")

        # Create parent directories if they don't exist
        logfile.parent.mkdir(exist_ok=True, parents=True)

        # Set up logging configuration
        logging.basicConfig(
            format=config.getini("log_file_format"),  # Logging format from pytest configuration
            filename=logfile,                         # Log file path for this worker
            level=config.getini("log_file_level"),    # Log level from pytest configuration
        )

        # Open the log file in write-binary mode
        with logfile.open("wb") as log_fh:
            LOGFILES[worker_id] = log_fh

            # Redirect stdout and stderr to the log file
            os.dup2(log_fh.fileno(), sys.stdout.fileno())
            os.dup2(log_fh.fileno(), sys.stderr.fileno())