import asyncio
import random
import typing
from venv import logger
from ninja import Router
from rbxz_bend.settings import get_logger
from ribctl.lib.struct_rcsb_api import current_rcsb_structs
from ribctl.lib.types.types_ribosome_assets import RibosomeAssets
from ribctl.db.ribosomexyz import Neo4jDB
from ribctl.lib.types.types_polymer import RNAClass
from ribctl.lib.types.types_ribosome import ExogenousRNAByStruct, ProteinClass, RibosomeStructure
from ribctl.db.data import QueryOps
from schema.v0 import BanClassMetadata, LigandInstance, LigandlikeInstance, NeoStruct, NomenclatureClass, NomenclatureClassMember
from ninja import NinjaAPI
import concurrent.futures
import logging
import os

v0 = Router()
QO = QueryOps()



@v0.get('/log_test', tags=['0'], )
def log_test(request):

   get_logger('computations').info('test')
   get_logger('computations').critical('test')
   get_logger('custom').info('test')
   get_logger('custom2').debug('test')

# @router.get('/async_test')
# def async_test(request):

#     def testf():
#         # Create a logger with the same name as the file
#         script_name = os.path.splitext(os.path.basename(__file__))[0]
#         script_path = os.path.dirname(os.path.abspath(__file__))

#         logger = logging.getLogger(script_name)
#         logger.setLevel(logging.DEBUG)

#         console_handler = logging.StreamHandler()
#         console_handler.setLevel(logging.DEBUG)

#         file_handler = logging.FileHandler(os.path.join(script_path, f"{script_name}_log.txt"))
#         file_handler.setLevel(logging.INFO)

#         formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
#         console_handler.setFormatter(formatter)

#         logger.addHandler(console_handler)
#         logger.addHandler(file_handler)

#         import time

#         for i in range(10):
#             logger.debug("Passed:{}".format(i))
#             time.sleep(1)

#         logger.info("Done")

    # executor = concurrent.futures.ThreadPoolExecutor(max_workers=1)
    # executor.submit(testf)
    # return {"message": "Started work"}



@v0.get('/sync_with_rcsb', response=list[str], tags=['0-Operations'])
def sync_with_rcsb(request):
    D        = Neo4jDB()
    synced   = D.get_all_structs()
    unsynced = sorted(current_rcsb_structs())

    for rcsb_id in set(unsynced ) - set(synced):
        assets = RibosomeAssets(rcsb_id)
        try:
            assets._verify_json_profile(True)
            D.add_structure(assets)
        except Exception as e:
            print(e)
            logger.error("Exception occurred:", exc_info=True)


@v0.get('/get_all_structures', response=list[NeoStruct], tags=['Structure'])
def get_all_structures(request,):
    return QO.get_all_structures()

@v0.get('/get_struct', response=NeoStruct, tags=['Structure'])
def get_struct(request,rcsb_id:str):
    return QO.get_struct(rcsb_id.upper())

@v0.get('/get_full_structure', response=NeoStruct, tags=['Structure'])
def get_full_structure(request,rcsb_id:str):
    return QO.get_struct(rcsb_id.upper())

@v0.get('/get_all_ligands', response=list[NeoStruct], tags=['Ligand'])
def get_all_ligands(request,):
    return QO.get_all_ligands()

@v0.get('/get_individual_ligand', response=list[LigandInstance], tags=['Ligand'])
def get_individual_ligand(request,chemicalId:str):
    return QO.get_individual_ligand(chemicalId)
    
@v0.get('/get_all_ligandlike', response=list[LigandlikeInstance], tags=['Ligand'])
def get_all_ligandlike(request,):
    return QO.get_all_ligandlike()

@v0.get('/get_RibosomeStructure', response=RibosomeStructure, tags=['Structure'])
def get_RibosomeStructure(request,rcsb_id:str):
    return QO.get_RibosomeStructure(rcsb_id.upper())

@v0.get('/match_structs_w_proteins', response=RibosomeStructure, tags=['Structure'])
def match_structs_w_proteins(request,has_proteins:list[ProteinClass]):
    return QO.match_structs_w_proteins(has_proteins)
    

@v0.get('/get_banclass_for_chain', response=list[ProteinClass], tags=['Classification'])
def get_banclass_for_chain(request,rcsb_id:str, auth_asym_id:str):
    return QO.get_banclass_for_chain(rcsb_id,auth_asym_id)
    

@v0.get('/get_banclasses_metadata', response=list[BanClassMetadata], tags=['Classification'])
def get_banclasses_metadata(request,family:typing.Literal['b','e','u'], subunit:typing.Literal['SSU', 'LSU']):
    return QO.get_banclasses_metadata(family, subunit)
    
@v0.get('/get_nom_classes', response=list[NomenclatureClass], tags=['Classification'])
def get_nom_classes(request,):
    return QO.list_nom_classes()

@v0.get('/gmo_nom_class', response=list[ NomenclatureClassMember ], tags=['Classification'])
def gmo_nom_class(request,class_id:ProteinClass):
    return QO.gmo_nom_class(class_id)

@v0.get('/proteins_number', response=int, tags=['Protein'])
def proteins_number(request):
    return QO.proteins_number()

@v0.get('/get_rnas_by_struct', response=list[ExogenousRNAByStruct], tags=['RNA'])
def get_rnas_by_struct(request):
    return QO.get_rnas_by_struct()

@v0.get('/get_rna_class', response=list[NomenclatureClassMember], tags=['RNA'])
def get_rna_class(request,class_id:RNAClass):
    return QO.get_rna_class(class_id)

