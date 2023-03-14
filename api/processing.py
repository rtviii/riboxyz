import concurrent.futures
import random
from venv import logger

# from api.ribctl.db.ribosomexyz import Neo4jDB
# from api.ribctl.lib.struct_rcsb_api import current_rcsb_structs
# from api.ribctl.lib.types.types_ribosome_assets import RibosomeAssets

import logging
import os

script_name = os.path.splitext(os.path.basename(__file__))[0]
script_path = os.path.dirname(os.path.abspath(__file__))

# Create a logger with the same name as the file
logger = logging.getLogger(script_name)
logger.setLevel(logging.DEBUG)

console_handler = logging.StreamHandler()
console_handler.setLevel(logging.DEBUG)

file_handler = logging.FileHandler(os.path.join(script_path, f"{script_name}_log.txt"))
file_handler.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)

# logger.addHandler(console_handler)
logger.addHandler(file_handler)


def testf():
    for i in range(1000000):
        number = random.randint(0,100100)
        if number > 10000:
            logger.debug("Exception occurred with number {}".format(number))



with concurrent.futures.ProcessPoolExecutor() as executor:
    future = executor.submit(testf)
    while future.running():
        ...
    print("Function started.")