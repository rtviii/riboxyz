
import os


def ready(self):
    logs_dir = os.path.join(BASE_DIR, 'logs')
    if not os.path.exists(logs_dir):
        os.makedirs(logs_dir)
