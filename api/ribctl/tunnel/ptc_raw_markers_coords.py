#!/usr/bin/env python3
# kdd@math.ubc.ca
# rtkushner@gmail.com
# This script is courtesy of the ribosome.xyz and its authors.
# This relies on the following packages to run
# - gemmi   : https://gemmi.readthedocs.io/en/latest/install.html
# - bipython: https://biopython.org/
# And additionally "requests" to download missing structures: https://pypi.org/project/requests/

# Distribute freely.
from functools import reduce
import json
import os
import sys
from Bio import SeqRecord
from pprint import pprint
import subprocess
from typing import List, Tuple
import re
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from fuzzysearch import find_near_matches
import pandas
from Bio import pairwise2
from Bio import SeqIO
import gemmi
import argparse
from api.ribctl.etl.ribosome_assets import RibosomeAssets
from api.ribctl.lib.utils import open_structure
from api.scratch_tunnel_workflow import DORIS_ET_AL, pick_match
from api.taxonomy import filter_by_parent_tax

RIBETL_DATA = os.environ.get('RIBETL_DATA')
parser      = argparse.ArgumentParser(description='CLI for locating PTC residues of 23SrRNA in a given prokaryotic PDB file')

args = parser .parse_args()
argdict = vars(parser.parse_args())

bacteria_structs = filter_by_parent_tax(2)




if args.canon:
    domain         = 'bacteria'
    rcsb_id        = argdict["target"].upper()

    R = RibosomeAssets(rcsb_id)
    struct, profile = R.get_struct_and_profile()

    lsu_rna  = R.get_LSU_rRNA()
    chain_id = lsu_rna.auth_asym_id
    

    chain3d       : Chain = struct.child_dict[0].child_dict[chain_id]
    ress          : List[Residue] = chain3d.child_list
    _r            : Residue

    ress_sanitized: List[Residue] = [*filter(lambda r: r.get_resname() in ["A", "C", "G", "U", "PSU"], ress)]

    for _r in ress_sanitized:
        if _r.get_resname() == "PSU":
            _r.resname = "U"

    raw_seq         = reduce(lambda x,y: x + y.resname, ress_sanitized,'')
    matches         = find_near_matches(DORIS_ET_AL["subseq_9"], raw_seq, max_l_dist=1)
    m0              = pick_match(matches, len(raw_seq))
    sought_residues = [ ress_sanitized[i] for i in list(range(m0.start, m0.end))]
    PTC_MARKERS     = {}

    for res in sought_residues:
        if res.id[1] not in PTC_MARKERS:
            PTC_MARKERS[res.id[1]] = {}

        atom:Atom
        for atom in res.child_list:
            atom_name                              = atom.name
            atom_coords                            = atom.get_coord()
            PTC_MARKERS[res.id[1]][atom_name] = list(map(lambda x: float(x), list(atom_coords)))

    PTC_MARKERS = {
        chain_id: PTC_MARKERS
    }

    markers_dir = os.path.join(RIBETL_DATA, "PTC_MARKERS_RAW")
    outfile     = os.path.join(markers_dir, f"{rcsb_id.upper()}_PTC_MARKERS_RAW.json")
    with open(outfile, 'w') as outf:
        json.dump(PTC_MARKERS, outf, indent=4)
        print("Saved {} successfully.".format(outfile))

if args.generate:

    f     = open('rcsb_pdb_ids_20230106032038.txt', 'r')
    lines = f.readlines()

    f.close()
    bacteria_structs =[*map(lambda _: _.strip("\n"),lines)]

    struct_ids    = []
    parent_chain  = []
    residue_name  = []
    residue_seqid = []

    coord_x       = []
    coord_y       = []
    coord_z       = []

    i             = 0
    for struct in bacteria_structs:
        struct     = str.upper(struct).strip("\n")
        markerfile = os.path.join(RIBETL_DATA,"PTC_MARKERS_RAW", f"{struct}_PTC_MARKERS_RAW.json")

        if not os.path.isfile(markerfile):
            print(f"Could not find {markerfile} in RIBETL_DATA")
            continue
        try:
            with open(markerfile, 'r') as f:
                POSNS = json.load(f)

            chain = [*POSNS.keys()][0]
            if chain == None:
                exit("Could not identify chain")

            if "O4'" in [ *POSNS[chain].values() ][len(POSNS[chain]) - 2]:
                U_end_pos   = [ *POSNS[chain].values() ][len(POSNS[chain]) - 2]["O4'"] # pre  last residue of the comb
            else:
                U_end_pos   = [ *POSNS[chain].values() ][len(POSNS[chain]) - 2]["C4"] # pre  last residue of the comb

            if "O4'" in  [ *POSNS[chain].values() ][0]:
                U_start_pos = [ *POSNS[chain].values() ][0]["O4'"]                     # first residue of the comb
            else:
                U_start_pos = [ *POSNS[chain].values() ][0]["C4"]                     # first residue of the comb


            midpoint = [
                ( U_end_pos[0] + U_start_pos[0] ) / 2,
                ( U_end_pos[1] + U_start_pos[1] ) / 2,
                ( U_end_pos[2] + U_start_pos[2] ) / 2,
            ]


            struct_ids.append(struct)
            parent_chain.append(chain)
            coord_x.append(midpoint[0])
            coord_y.append(midpoint[1])
            coord_z.append(midpoint[2])

            i = i+1
            print(f"Processed structs  : {i}")
        except:
            i = i+1
            print("err")
            print(f"Processed structs  : {i}")
            continue
    # add this to the dataframe
    df = pandas.DataFrame.from_dict({
        'struct_ids'  : struct_ids,
        'parent_chain': parent_chain,
        'coord_x'     : coord_x,
        'coord_y'     : coord_y,
        'coord_z'     : coord_z
    })

    # save the dataframe to a csv
    df.to_csv(f"ptc_centroids.csv")
    exit(1)

def PTC_MARKERS_RAW(rcsb_id:str):


# if not args.display_all:
    # print("\nTo display more residues per target structure, use additional --display_all flag.")


# 5afi | A| resi 2610-2611 |~ C3
# 3j7z | sele c. A and resi 2609 and name C4 | O4 | n4
