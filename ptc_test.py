from pprint import pprint
from ribctl.lib.landmarks.ptc_via_trna import PTC_location


# loc, residues = PTC_location('4UG0')
# pprint(residues)
# loc, residues = PTC_location('8HKZ')
# pprint(residues)
loc, residues = PTC_location('5AFI')
pprint(residues)