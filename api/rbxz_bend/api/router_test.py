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

test = Router()
qo = QueryOps()



@test.get('/empty_db_query', tags=['0test'], )
def any_test(request):
    return qo.get_any()

@test.get('/log_test',tags=['0test']  )
def log_test(request):

   get_logger('computations').info('test')
   get_logger('computations').critical('test')
   get_logger('custom').info('test')
   get_logger('custom2').debug('test')



@test.get('/sync_with_rcsb', response=list[str], tags=['0-Operations'])
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