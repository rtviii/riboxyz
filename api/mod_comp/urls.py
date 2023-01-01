from django.urls import path
from .views import * 


urlpatterns = [
    path('hello/', hello),
    path('edit_chain/', edit_chain)
]
app_name = 'mod_comp'