
from django.urls import path
from .views import * 


urlpatterns = [
    path('hello/'       , hello        ),
    path('get_struct/'  , get_struct   ),
    path('health/'      , check_health ),
    path('envs/'        , get_envs     ),
]
app_name = 'mod_db'