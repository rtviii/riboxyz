
from ribctl.lib.tunnel import open_tunnel_csv, parse_struct_via_centerline


RCSB_ID = '6Z6K'

data = open_tunnel_csv(RCSB_ID)
parse_struct_via_centerline(RCSB_ID,data)