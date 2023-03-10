import typing
from ninja import Router
from ribctl.lib.types.types_polymer import RNAClass
from ribctl.lib.types.types_ribosome import ExogenousRNAByStruct, ProteinClass, RibosomeStructure
from ribctl.db.data import QueryOps
from schema.v0 import BanClassMetadata, LigandInstance, LigandlikeInstance, NeoStruct, NomenclatureClass, NomenclatureClassMember

router = Router()
QO     = QueryOps()

@router.get('/v0/get_all_structures', response=list[NeoStruct])
def get_all_structures(request,):
    return QO.get_all_structures()

@router.get('/v0/get_struct', response=NeoStruct)
def get_struct(request,rcsb_id:str):
    return QO.get_struct(rcsb_id.upper())

@router.get('/v0/get_full_structure', response=NeoStruct)
def get_full_structure(request,rcsb_id:str):
    return QO.get_struct(rcsb_id.upper())

@router.get('/v0/get_all_ligands', response=list[NeoStruct])
def get_all_ligands(request,):
    return QO.get_all_ligands()

@router.get('/v0/get_individual_ligand', response=list[LigandInstance])
def get_individual_ligand(request,chemicalId:str):
    return QO.get_individual_ligand(chemicalId)
    
@router.get('/v0/get_all_ligandlike', response=list[LigandlikeInstance])
def get_all_ligandlike(request,):
    return QO.get_all_ligandlike()

@router.get('/v0/get_RibosomeStructure', response=RibosomeStructure)
def get_RibosomeStructure(request,rcsb_id:str):
    return QO.get_RibosomeStructure(rcsb_id.upper())

@router.get('/v0/match_structs_w_proteins', response=RibosomeStructure)
def match_structs_w_proteins(request,has_proteins:list[ProteinClass]):
    return QO.match_structs_w_proteins(has_proteins)
    

@router.get('/v0/get_banclass_for_chain', response=list[ProteinClass])
def get_banclass_for_chain(request,rcsb_id:str, auth_asym_id:str):
    return QO.get_banclass_for_chain(rcsb_id,auth_asym_id)
    

@router.get('/v0/get_banclasses_metadata', response=list[BanClassMetadata])
def get_banclasses_metadata(request,family:typing.Literal['b','e','u'], subunit:typing.Literal['SSU', 'LSU']):
    return QO.get_banclasses_metadata(family, subunit)
    
@router.get('/v0/get_nom_classes', response=list[NomenclatureClass])
def get_nom_classes(request,):
    return QO.list_nom_classes()

@router.get('/v0/gmo_nom_class', response=list[ NomenclatureClassMember ])
def gmo_nom_class(request,class_id:ProteinClass):
    return QO.gmo_nom_class(class_id)

@router.get('/v0/proteins_number', response=int)
def proteins_number(request):
    return QO.proteins_number()

@router.get('/v0/get_rnas_by_struct', response=list[ExogenousRNAByStruct])
def get_rnas_by_struct(request):
    return QO.get_rnas_by_struct()

@router.get('/v0/get_rna_class', response=list[NomenclatureClassMember])
def get_rna_class(request,class_id:RNAClass):
    return QO.get_rna_class(class_id)


