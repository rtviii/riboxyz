from pprint import pprint
import typing
from django.forms import ValidationError
from ninja import Router, Schema
from ribctl.lib.types.types_polymer import RNAClass
from ribctl.lib.types.types_ribosome import ExogenousRNAByStruct, ProteinClass, RibosomeStructure
from schema.v0 import BanClassMetadata, LigandInstance, LigandlikeInstance, NeoStruct, NomenclatureClass, NomenclatureClassMember
from rbxz_bend.application import db_connection

v0 = Router()




@v0.get('/get_all_structures',tags=['Structure'], 
        # response=list[NeoStruct]
        # TODO: validate
        )
def get_all_structures(request,):
    return db_connection.get_all_structures()


@v0.get('/get_struct', 
        # response=NeoStruct, # TODO: validate
        tags=['Structure']
        )

def get_struct(request,rcsb_id:str):
    return db_connection.get_struct(rcsb_id.upper())

@v0.get('/get_full_structure', response=NeoStruct, tags=['Structure'])
def get_full_structure(request,rcsb_id:str):
    return db_connection.get_struct(rcsb_id.upper())

@v0.get('/get_all_ligands', 
        # response=list[NeoStruct],
          tags=['Ligand'])
def get_all_ligands(request,):
    return db_connection.get_all_ligands()

@v0.get('/get_individual_ligand', response=list[LigandInstance], tags=['Ligand'])
def get_individual_ligand(request,chemicalId:str):
    return db_connection.get_individual_ligand(chemicalId)
    
@v0.get('/get_all_ligandlike', response=list[LigandlikeInstance], tags=['Ligand'])
def get_all_ligandlike(request,):
    return db_connection.get_all_ligandlike()

@v0.get('/get_RibosomeStructure', response=RibosomeStructure, tags=['Structure'])
def get_RibosomeStructure(request,rcsb_id:str):
    return db_connection.get_RibosomeStructure(rcsb_id.upper())

@v0.get('/match_structs_w_proteins', response=RibosomeStructure, tags=['Structure'])
def match_structs_w_proteins(request,has_proteins:list[ProteinClass]):
    return db_connection.match_structs_w_proteins(has_proteins)

@v0.get('/get_banclass_for_chain', response=list[ProteinClass], tags=['Classification'])
def get_banclass_for_chain(request,rcsb_id:str, auth_asym_id:str):
    return db_connection.get_banclass_for_chain(rcsb_id,auth_asym_id)

@v0.get('/get_banclasses_metadata', response=list[BanClassMetadata], tags=['Classification'])
def get_banclasses_metadata(request,family:typing.Literal['b','e','u'], subunit:typing.Literal['SSU', 'LSU']):
    return db_connection.get_banclasses_metadata(family, subunit)
    
@v0.get('/get_nom_classes', response=list[NomenclatureClass], tags=['Classification'])
def get_nom_classes(request,):
    return db_connection.list_nom_classes()

@v0.get('/gmo_nom_class', response=list[ NomenclatureClassMember ], tags=['Classification'])
def gmo_nom_class(request,class_id:ProteinClass):
    return db_connection.gmo_nom_class(class_id)

@v0.get('/proteins_number', response=int, tags=['Protein'])
def proteins_number(request):
    return db_connection.proteins_number()

@v0.get('/number_of_structures', response=int, tags=['Structure'])
def number_of_structures(request):
    return db_connection.number_of_structures()

@v0.get('/get_rnas_by_struct', response=list[ExogenousRNAByStruct], tags=['RNA'])
def get_rnas_by_struct(request):
    return db_connection.get_rnas_by_struct()

@v0.get('/get_rna_class', response=list[NomenclatureClassMember], tags=['RNA'])
def get_rna_class(request,rna_class:RNAClass):
    return db_connection.get_rna_class(rna_class)


