from ninja import NinjaAPI
from .custom_router import v0

api              = NinjaAPI(
    docs_url='/',
    title='riboxyz-api',
    description='A programming interface for ribosome data.'
    
                            )
api.add_router('/v0', v0)


# ninja_api_router = Router()
# QO               = QueryOps()



# @ninja_api_router.get('/log_test')
# def log_test(request):

#    get_logger('computations').info('test')
#    get_logger('computations').critical('test')
#    get_logger('custom').info('test')
#    get_logger('custom2').debug('test')

# # @router.get('/async_test')
# # def async_test(request):

# #     def testf():
# #         # Create a logger with the same name as the file
# #         script_name = os.path.splitext(os.path.basename(__file__))[0]
# #         script_path = os.path.dirname(os.path.abspath(__file__))

# #         logger = logging.getLogger(script_name)
# #         logger.setLevel(logging.DEBUG)

# #         console_handler = logging.StreamHandler()
# #         console_handler.setLevel(logging.DEBUG)

# #         file_handler = logging.FileHandler(os.path.join(script_path, f"{script_name}_log.txt"))
# #         file_handler.setLevel(logging.INFO)

# #         formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
# #         console_handler.setFormatter(formatter)

# #         logger.addHandler(console_handler)
# #         logger.addHandler(file_handler)

# #         import time

# #         for i in range(10):
# #             logger.debug("Passed:{}".format(i))
# #             time.sleep(1)

# #         logger.info("Done")

#     # executor = concurrent.futures.ThreadPoolExecutor(max_workers=1)
#     # executor.submit(testf)
#     # return {"message": "Started work"}



# @ninja_api_router.get('/sync_with_rcsb', response=list[str])
# def sync_with_rcsb(request):
#     D        = Neo4jDB()
#     synced   = D.get_all_structs()
#     unsynced = sorted(current_rcsb_structs())

#     for rcsb_id in set(unsynced ) - set(synced):
#         assets = RibosomeAssets(rcsb_id)
#         try:
#             assets._verify_json_profile(True)
#             D.add_structure(assets)
#         except Exception as e:
#             print(e)
#             logger.error("Exception occurred:", exc_info=True)


# @ninja_api_router.get('/v0/get_all_structures', response=list[NeoStruct])
# def get_all_structures(request,):
#     return QO.get_all_structures()

# @ninja_api_router.get('/v0/get_struct', response=NeoStruct)
# def get_struct(request,rcsb_id:str):
#     return QO.get_struct(rcsb_id.upper())

# @ninja_api_router.get('/v0/get_full_structure', response=NeoStruct)
# def get_full_structure(request,rcsb_id:str):
#     return QO.get_struct(rcsb_id.upper())

# @ninja_api_router.get('/v0/get_all_ligands', response=list[NeoStruct])
# def get_all_ligands(request,):
#     return QO.get_all_ligands()

# @ninja_api_router.get('/v0/get_individual_ligand', response=list[LigandInstance])
# def get_individual_ligand(request,chemicalId:str):
#     return QO.get_individual_ligand(chemicalId)
    
# @ninja_api_router.get('/v0/get_all_ligandlike', response=list[LigandlikeInstance])
# def get_all_ligandlike(request,):
#     return QO.get_all_ligandlike()

# @ninja_api_router.get('/v0/get_RibosomeStructure', response=RibosomeStructure)
# def get_RibosomeStructure(request,rcsb_id:str):
#     return QO.get_RibosomeStructure(rcsb_id.upper())

# @ninja_api_router.get('/v0/match_structs_w_proteins', response=RibosomeStructure)
# def match_structs_w_proteins(request,has_proteins:list[ProteinClass]):
#     return QO.match_structs_w_proteins(has_proteins)
    

# @ninja_api_router.get('/v0/get_banclass_for_chain', response=list[ProteinClass])
# def get_banclass_for_chain(request,rcsb_id:str, auth_asym_id:str):
#     return QO.get_banclass_for_chain(rcsb_id,auth_asym_id)
    

# @ninja_api_router.get('/v0/get_banclasses_metadata', response=list[BanClassMetadata])
# def get_banclasses_metadata(request,family:typing.Literal['b','e','u'], subunit:typing.Literal['SSU', 'LSU']):
#     return QO.get_banclasses_metadata(family, subunit)
    
# @ninja_api_router.get('/v0/get_nom_classes', response=list[NomenclatureClass])
# def get_nom_classes(request,):
#     return QO.list_nom_classes()

# @ninja_api_router.get('/v0/gmo_nom_class', response=list[ NomenclatureClassMember ])
# def gmo_nom_class(request,class_id:ProteinClass):
#     return QO.gmo_nom_class(class_id)

# @ninja_api_router.get('/v0/proteins_number', response=int)
# def proteins_number(request):
#     return QO.proteins_number()

# @ninja_api_router.get('/v0/get_rnas_by_struct', response=list[ExogenousRNAByStruct])
# def get_rnas_by_struct(request):
#     return QO.get_rnas_by_struct()

# @ninja_api_router.get('/v0/get_rna_class', response=list[NomenclatureClassMember])
# def get_rna_class(request,class_id:RNAClass):
#     return QO.get_rna_class(class_id)

