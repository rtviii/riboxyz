from django.urls import path
from .views import * 

#TODO:
# - map old endpoints to new ones
# - serve pfam mappings to protein families
# - "last updated" endpoint for fend serving : date, some new structure ids.

urlpatterns = [
    path('struct_commit_new/'     , struct_commit_new     ),
    path('structs_all_ids/'       , structs_all_ids       ),
    path('structs_diff_pdb/'      , structs_diff_pdb      ),
    path('structs_sync_with_pdb/' , structs_sync_with_pdb ),


# Pre migration ednpoints 
## neo4j connector

#       path('get_struct/'              , get_struct              ),
#       path('get_full_struct/'         , get_full_structure      ),
#       path('get_RibosomeStructure/'   , get_RibosomeStructure      ),
#     # path('get_homologs/'            , get_homologs            ) ,
#       path('cypher/'                  , custom_cypher           ) ,
#       path('anything/'                , anything                ) ,
#       path('list_nom_classes/'        , list_nom_classes        ) ,
#       path('get_all_structs/'         , get_all_structs         ) ,
#       path('gmo_nom_class/'           , gmo_nom_class           ) ,
#       path('get_banclasses_metadata/' , get_banclasses_metadata ) ,
#       path('nomclass_visualize/'      , nomclass_visualize      ) ,
#       path('get_banclass_for_chain/', get_banclass_for_chain ) ,
#       path('banclass_annotation/'     , banclass_annotation     ) ,
#       path('get_rna_class/'           , get_rna_class           ) ,
#       path('get_individual_ligand/'   , get_individual_ligand   ) ,
#       path('get_rnas_by_struct/'      , get_rnas_by_struct      ) ,
#       path('get_ligands_by_struct/'   , get_ligands_by_struct   ) ,
#       path('get_all_ligands/'         , get_all_ligands         ) ,
#       path('get_all_ligandlike/'      , get_all_ligandlike      ) ,
#       path('match_structs/'           , match_structs           ) ,
#       path('proteins_number/'         , proteins_number         ),
#       path('tax_ids/'                 , tax_ids                 ),
#       path('nomenclature/'            , nomenclature            ),

# ## static_files
#     path('get_ligand_nbhd/'      , get_ligand_nbhd      ) ,
#     path('download_ligand_nbhd/' , download_ligand_nbhd ) ,
#     path('cif_chain/'            , cif_chain            ) ,
#     path('ligand_prediction/'    , ligand_prediction    ) ,
#     path('download_structure/'   , download_structure   ) ,
#     path('cif_chain_by_class/'   , cif_chain_by_class   ) ,
#     path('ranged_align/'         , ranged_align         ) ,

]

app_name = 'mod_utils'