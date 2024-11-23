from ribctl.lib.libbsite import bsite_ligand, bsite_transpose

cid  ='PAR'

vsite= bsite_ligand(cid, '7K00', 10)
bsite_transpose("7K00",'5AFI',vsite, False, True)