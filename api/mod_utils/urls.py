from django.urls import path
from .views import * 

#TODO:
# - map old endpoints to new ones
# - serve pfam mappings to protein families
# - "last updated" endpoint for fend serving : date, some new structure ids.

urlpatterns = [
    path('hello/' , hello),
    path('struct_commit_new_PDB/',struct_commit_new_PDB ),
    path('structs_get_ids', structs_get_ids),
]

app_name = 'mod_utils'