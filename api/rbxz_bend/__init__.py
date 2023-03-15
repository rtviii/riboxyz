from django.apps import AppConfig
import os

from rbxz_bend.settings import BASE_DIR

class MyAppConfig(AppConfig):
    name = 'rbxz_bend'
    def ready(self):
        logs_dir = os.path.join(BASE_DIR, 'logs')
        if not os.path.exists(logs_dir):
            os.makedirs(logs_dir)
