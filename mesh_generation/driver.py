
from ribctl.lib.tunnel import encode_atoms, open_tunnel_csv, parse_struct_via_centerline


RCSB_ID = '6Z6K'

data = open_tunnel_csv(RCSB_ID)
atoms = parse_struct_via_centerline(RCSB_ID,data)

print(encode_atoms(RCSB_ID,atoms))