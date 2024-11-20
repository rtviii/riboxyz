# update_manager.py
import subprocess
import logging
import time
from pathlib import Path
import sys
from datetime import datetime
import json
import os
from typing import Optional

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('update_manager.log'),
        logging.StreamHandler(sys.stdout)
    ]
)

logger = logging.getLogger('UpdateManager')

class UpdateManager:
    def __init__(self, base_dir: str):
        self.base_dir = Path(base_dir)
        self.status_file = self.base_dir / 'status.json'
        self.lock_file = self.base_dir / 'update.lock'
        
    def is_update_running(self) -> bool:
        return self.lock_file.exists()

    def create_lock(self):
        self.lock_file.touch()

    def remove_lock(self):
        if self.lock_file.exists():
            self.lock_file.unlink()

    def update_status(self, status: str, error: Optional[str] = None):
        current_status = {
            'last_update_attempt': datetime.now().isoformat(),
            'status': status,
            'error': error
        }
        with open(self.status_file, 'w') as f:
            json.dump(current_status, f)

    def run_update(self):
        if self.is_update_running():
            logger.warning("Update already in progress")
            return

        try:
            self.create_lock()
            self.update_status('starting')
            
            # Stop the services
            logger.info("Stopping services...")
            subprocess.run(['docker-compose', 'down'], check=True)

            # Run your ETL process
            logger.info("Running ETL process...")
            self.update_status('running_etl')
            # Replace with your actual ETL command
            subprocess.run(['python', 'your_etl_script.py'], check=True)

            # Start services back up
            logger.info("Starting services...")
            self.update_status('starting_services')
            subprocess.run(['docker-compose', 'up', '-d'], check=True)

            # Wait for Neo4j to be ready
            self._wait_for_neo4j()

            # Run database updates
            logger.info("Running database updates...")
            self.update_status('updating_database')
            subprocess.run([
                'docker', 'exec', 'django',
                'python', 'manage.py', 'update_database'
            ], check=True)

            self.update_status('completed')
            logger.info("Update completed successfully")

        except Exception as e:
            error_msg = str(e)
            logger.error(f"Update failed: {error_msg}")
            self.update_status('failed', error_msg)
            raise
        finally:
            self.remove_lock()

    def _wait_for_neo4j(self, timeout=300):
        """Wait for Neo4j to be ready, with timeout in seconds"""
        start_time = time.time()
        while time.time() - start_time < timeout:
            try:
                result = subprocess.run([
                    'docker', 'exec', 'neo',
                    'neo4j-admin', 'dbms', 'status'
                ], capture_output=True, text=True)
                if "Neo4j is running" in result.stdout:
                    return True
                time.sleep(5)
            except subprocess.CalledProcessError:
                time.sleep(5)
        raise TimeoutError("Neo4j failed to start within timeout period")

if __name__ == "__main__":
    manager = UpdateManager('/path/to/your/app')
    manager.run_update()