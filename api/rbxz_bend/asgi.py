from asgiref.wsgi import WsgiToAsgi
from rbxz_bend.wsgi import application
from django.conf import settings
from django.contrib.staticfiles.urls import staticfiles_urlpatterns
from django.urls import re_path

app = WsgiToAsgi(application)
