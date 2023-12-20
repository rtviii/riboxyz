

import sys
from scripts.pymol_visualtion import by_chain, ray_picture, sload
rcsb_id = sys.argv[3].upper()
print("About to render pdbid: ",rcsb_id)
sload(rcsb_id)
by_chain(rcsb_id)
ray_picture(rcsb_id)