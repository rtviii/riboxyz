from concurrent.futures import ProcessPoolExecutor
import threading
import typing
from venv import logger
from ninja import Router
from rbxz_bend.settings import get_logger
from ribctl.lib.struct_rcsb_api import current_rcsb_structs
from ribctl.lib.types.types_ribosome_assets import RibosomeAssets
from rbxz_bend.db.ribosomexyz import ribosomexyzDB
from rbxz_bend.application import db_connection 
from rbxz_bend.settings import NEO4J_CURRENTDB,NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER
import asyncio

test = Router()

@test.get('/see_current_auth', tags=['0test'], )
def current_auth(request):
    return db_connection.see_current_auth()

@test.get('/empty_db_query', tags=['0test'], )
def any_test(request):
    return db_connection.get_any()

@test.get('/log_test',tags=['0test']  )
def log_test(request):

   get_logger('computations').info('test')
   get_logger('computations').critical('test')



@test.get('/sync_with_rcsb', tags=['0test'])
async def sync_with_rcsb(request):

    D        = ribosomexyzDB(NEO4J_URI, NEO4J_USER, NEO4J_PASSWORD)
    threading.Thread(target=D.sync_with_rcsb,args=(50,)).start()

    return {"message": "Started work"}






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