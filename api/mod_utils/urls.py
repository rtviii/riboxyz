from django.urls import path
from .views import * 

#TODO:
# - map old endpoints to new ones
# - serve pfam mappings to protein families
# - "last updated" endpoint for fend serving : date, some new structure ids.

urlpatterns = [
    path('hello/'                , hello                 ),
    path('struct_commit_new/'     , struct_commit_new     ),
    path('structs_all_ids/'       , structs_all_ids       ),
    path('structs_diff_pdb/'      , structs_diff_pdb      ),
    path('structs_sync_with_pdb/' , structs_sync_with_pdb ),
]

app_name = 'mod_utils'