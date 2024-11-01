from pprint import pprint
from ribctl.lib.landmarks.ptc_via_trna import PTC_location

# 7A5F,6HCQ,8RCS,7TM3
# loc, residues = PTC_location('4UG0')
# pprint(residues)
# loc, residues = PTC_location('8HKZ')
# pprint(residues)

loc, residues = PTC_location('8RCS')
pprint(loc)
pprint(residues)