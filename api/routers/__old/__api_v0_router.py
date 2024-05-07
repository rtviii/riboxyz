from io import BytesIO
from django.http import HttpResponse
from ninja import Router
from ribctl.lib.schema.types_ribosome import PolymerClass, CytosolicProteinClass, PolynucleotideClass, RibosomeStructure
from schema.v0 import BanClassMetadata, ExogenousRNAByStruct,LigandInstance, LigandlikeInstance, NeoStruct, NomenclatureClass, NomenclatureClassMember
from wsgiref.util import FileWrapper

apiv0 = Router()
@apiv0.get('/ranged_align',tags=['Static Files'])
def ranged_align(request,
                 range_start     : int,
                 range_end       : int,

                 src_rcsb_id     : str,
                 src_auth_asym_id: str,

                 tgt_rcsb_id     : str,
                 tgt_auth_asym_id: str): 

    print("Hit endpoint ranged align!!")
    # (src_auth_asym_id, src_path, src_range,
    #  tgt_auth_asym_id, tgt_path, tgt_range) = ranged_align_by_polyclass(src_rcsb_id, tgt_rcsb_id,res_range,poly_class)

    print("Got parameters")
    print(f"""      { range_start      }
                 { range_end        }
                 { src_rcsb_id      }
                 { src_auth_asym_id }
                 { tgt_rcsb_id      }
                 { tgt_auth_asym_id }""")
    (  _, _, src_range,
       _, _, tgt_range) = ranged_align_by_auth_asym_id(src_rcsb_id, src_auth_asym_id, tgt_rcsb_id, tgt_auth_asym_id, ( range_start,range_end ))

    cif_str:str = pymol_super(
        src_rcsb_id,
        src_range,
        src_auth_asym_id,

        tgt_rcsb_id,
        tgt_range,
        tgt_auth_asym_id,
    )

    response = HttpResponse(FileWrapper(BytesIO(bytes(cif_str,'utf-8'))), content_type='chemical/x-mmcif')
    response['Content-Disposition'] = 'attachment; filename="{}-{}_{}-{}.cif"'.format(src_rcsb_id,src_auth_asym_id,tgt_rcsb_id,tgt_auth_asym_id)
    return response

# ## ----- Old data endpoints ----- 

# @apiv0.get('/get_all_structures',tags=['Structure'], 
#         # response=list[NeoStruct]
#         # TODO: validate
#         )
# def get_all_structures(request,):
#     return db_connection.get_all_structures()

# @apiv0.get('/get_struct', 
#         # response=NeoStruct, # TODO: validate
#         tags=['Structure']
#         )

# def get_struct(request,rcsb_id:str):
#     return db_connection.get_struct(rcsb_id.upper())

# @apiv0.get('/get_full_structure', response=NeoStruct, tags=['Structure'])
# def get_full_structure(request,rcsb_id:str):
#     return db_connection.get_struct(rcsb_id.upper())

# @apiv0.get('/get_all_ligands', 
#         # response=list[NeoStruct],
#           tags=['Ligand'])
# def get_all_ligands(request,):
#     return db_connection.get_all_ligands()

# # @v0.get('/get_individual_ligand', response=list[LigandInstance], tags=['Ligand'])
# # def get_individual_ligand(request,chemicalId:str):
# #     return db_connection.get_individual_ligand(chemicalId)
    
# @apiv0.get('/get_all_ligandlike', response=list[LigandlikeInstance], tags=['Ligand'])
# def get_all_ligandlike(request,):
#     return db_connection.get_all_ligandlike()

# @apiv0.get('/get_RibosomeStructure', response=RibosomeStructure, tags=['Structure'])
# def get_RibosomeStructure(request,rcsb_id:str):
#     return db_connection.get_RibosomeStructure(rcsb_id.upper())

# @apiv0.get('/match_structs_w_proteins', response=RibosomeStructure, tags=['Structure'])
# def match_structs_w_proteins(request,has_proteins:list[CytosolicProteinClass]):
#     return db_connection.match_structs_w_proteins(has_proteins)

# @apiv0.get('/get_banclass_for_chain', response=list[CytosolicProteinClass], tags=['Classification'])
# def get_banclass_for_chain(request,rcsb_id:str, auth_asym_id:str):
#     return db_connection.get_banclass_for_chain(rcsb_id,auth_asym_id)

# @apiv0.get('/get_banclasses_metadata', response=list[BanClassMetadata], tags=['Classification'])
# def get_banclasses_metadata(request,family:typing.Literal['b','e','u'], subunit:typing.Literal['SSU', 'LSU']):
#     return db_connection.get_banclasses_metadata(family, subunit)
    
# @apiv0.get('/get_nom_classes', response=list[NomenclatureClass], tags=['Classification'])
# def get_nom_classes(request,):
#     return db_connection.list_nom_classes()

# @apiv0.get('/gmo_nom_class', response=list[ NomenclatureClassMember ], tags=['Classification'])
# def gmo_nom_class(request,class_id:CytosolicProteinClass):
#     return db_connection.gmo_nom_class(class_id)

# @apiv0.get('/proteins_number', response=int, tags=['Protein'])
# def proteins_number(request):
#     return db_connection.proteins_number()

# @apiv0.get('/number_of_structures', response=int, tags=['Structure'])
# def number_of_structures(request):
#     return db_connection.number_of_structures()

# @apiv0.get('/get_rnas_by_struct', response=list[ExogenousRNAByStruct], tags=['RNA'])
# def get_rnas_by_struct(request):
#     return db_connection.get_rnas_by_struct()

# @apiv0.get('/get_rna_class', response=list[NomenclatureClassMember], tags=['RNA'])
# def get_rna_class(request,rna_class:PolynucleotideClass):
#     return db_connection.get_rna_class(rna_class)


