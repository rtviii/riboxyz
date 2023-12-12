from django.urls import path
from .views import * 


urlpatterns = [
    path('hello/', hello),
    path('edit_chain/', edit_chain),
    path('get_chain/', get_chain),
    path('get_profile/', get_profile)
]
app_name = 'mod_comp'