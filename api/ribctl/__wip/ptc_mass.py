
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
from fuzzywuzzy import process
from fuzzysearch import find_near_matches
import re
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from ribctl.struct_extract_bsites import open_structure
import pandas
from Bio import pairwise2
from Bio import SeqIO
import gemmi
import argparse

RIBETL_DATA = os.environ.get('RIBETL_DATA')

# Change these two parameters to have a different "source" sequence to align *against*  ----|
#                                                                                           |
# List of conserved nucleotide sequences on 23s-28s can be found here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4574749/pdf/1719.pdf
# Some of the PTC residues in bacterial 23SrRNA                                             |


# -------------------------------------------------------------------------------------------|

parser = argparse.ArgumentParser(
    description='CLI for locating PTC residues of 23SrRNA in a given prokaryotic PDB file')

parser.add_argument("-t", "--target", type=str, required=False)
parser.add_argument("-d", "--domain", type=str,
                    required=False, choices=['e', 'b'])
parser.add_argument("-m", "--mode", type=str,  choices=['position', 'residue'], nargs='?', default='position',
                    help="Display the position of the PTC residue, or the residue itself. Default: position")
parser.add_argument("--display_all", action='store_true')
parser.add_argument("--generate", action='store_true')
parser.add_argument("--generate_from_raw", action='store_true')
parser.add_argument("--fuzzy", action='store_true')
parser.add_argument("--markers", action='store_true')
parser.add_argument("--canon", action='store_true')
parser.add_argument("--fasta_profile", action='store_true')
parser.add_argument("--ptc", action='store_true')

args = parser .parse_args()
argdict = vars(parser.parse_args())

bact_registry_file = open('rcsb_pdb_ids_20230106032038.txt', 'r')
line = bact_registry_file.readline()
bact_registry_file.close()
bacteria_structs = line.split(',')


def util__backwards_match(alntgt: str, aln_resid: int, verbose: bool = False) -> Tuple[int, str, int]:
    """
    returns (projected i.e. "ungapped" residue id, the residue itself residue)
    """
    if aln_resid > len(alntgt):
        raise IndexError(
            f"Passed residue with invalid index ({aln_resid}) to back-match to target. Seqlen:{len(alntgt)}")

    counter_proper = 0
    for i, char in enumerate(alntgt):
        if i == aln_resid:
            if verbose:
                print("[ {} ] <-----> id.[aligned: {} | orgiginal: {} ]".format(
                    alntgt[aln_resid], i, counter_proper))
            return (counter_proper,  alntgt[i], aln_resid)
        if char == '-':
            continue
        else:
            counter_proper += 1

    raise LookupError()


def util__forwards_match(string: str, resid: int):
    """Returns the index of a source-sequence residue in the (aligned) source sequence."""
    if resid >= len(string):
        raise IndexError(
            "Requested residue index({resid}) exceeds aligned(likely already gaps-extended) sequence. Something went wrong.")

    count_proper = 0
    for alignment_indx, char in enumerate(string):
        if count_proper == resid:
            return alignment_indx
        if char == '-':
            continue
        else:
            count_proper += 1




def get_sequence_by_nomclass(rcsb_id: str, nomenclature_class: str, canonical:bool=True,path:str=None) -> Tuple[str, str]:

    target       = gemmi.cif.read_file(path)
    block        = target.sole_block()
    model        = gemmi.read_structure(path)[0]

    STRAND = None
    SEQ    = None

    # Locate the chain of given nom. class
    for (strand, nomclass) in zip(
        block.find_loop('_ribosome_nomenclature.entity_poly.pdbx_strand_id'),
        block.find_loop('_ribosome_nomenclature.polymer_class')
    ):
        if nomclass == nomenclature_class:
            STRAND = strand
            break

    # Now find sequence of this class
    for (chain_id, one_letter_code) in zip(
        block.find_loop('_entity_poly.pdbx_strand_id'),
        block.find_loop('_entity_poly.pdbx_seq_one_letter_code_can') if canonical else block.find_loop('_entity_poly.pdbx_seq_one_letter_code')
    ):
        # X-RAY structures have 'dual' chains. Split on comma to check both.
        if STRAND in chain_id.split(','):
            SEQ = str(one_letter_code).strip(";").strip("\n")

    if SEQ == None:
        print("Could not locate {} sequence in {} CIF file".format(
            nomenclature_class, rcsb_id))
    return (STRAND, SEQ)

def retrieve_LSU_rRNA(rcsb_id, canonical:bool=True):
    annotated_cifpath = os.path.join(RIBETL_DATA, rcsb_id.upper(), f"{rcsb_id.upper()}_modified.cif")
    rna_type          = ""
    #--------------
    [chain_id, strand_target] = get_sequence_by_nomclass(
        rcsb_id,
        "23SrRNA",
        canonical,
        path=os.path.join(
            RIBETL_DATA, rcsb_id.upper(), f"{rcsb_id.upper()}_modified.cif")
    )
    rna_type = "23SrRNA"

    if chain_id == None or strand_target == None:
        [chain_id, strand_target] = get_sequence_by_nomclass(
            rcsb_id,
            "25SrRNA",
            canonical,
            path=annotated_cifpath
        )
        rna_type = "25SrRNA"

    if chain_id == None or strand_target == None:
        [chain_id, strand_target] = get_sequence_by_nomclass(
            rcsb_id,
            "28SrRNA",
            canonical,
            path=annotated_cifpath)
        rna_type = "28SrRNA"

    if chain_id == None or strand_target == None:
        print("Failed to locate either 23S or 25S or 28 rRNA in {}".format(rcsb_id))
        exit(1)

    return [chain_id, strand_target, rna_type]

# ※ ---------------------------- 23/25/28SrRNA PTC residue locations ---------------------------- ※
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4574749/pdf/1719.pdf
DORIS_ET_AL = {
    "subseq_6": "AAGACCC",
    "subseq_8": "GGAUAAC",
    "subseq_9": "GAGCUGGGUUUA",
    "e": {
        "site_6": [2400, 2401, 2402, 2403, 2404, 2405, 2406, 2407],
        "site_8": [2811, 2812, 2813, 2814, 2815, 2816, 2817, 2818],
        "site_9": [2941, 2942, 2943, 2944, 2945, 2946, 2947, 2948, 2949, 2950, 2951, 2952, 2953]
    },
    'b': {
        "site_6": [2059, 2060, 2061,
                   2062, 2063, 2064,
                   2065],
        "site_8": [2446, 2447, 2448,
                   2449, 2450, 2451, 2452],
        "site_9": [2576, 2577, 2578, 2579,
                   2580, 2581, 2582, 2583,
                   2584, 2585, 2586, 2587],
    }
}



class RibovisionAlignment:

    def __init__(self, mass_alignment_path: str) -> None:
        self.fasta_seqs: List[SeqIO.SeqRecord] = SeqIO.parse(
            open(mass_alignment_path), 'fasta')
        pass

    def find_aln_by_id(self, rcsb_id: str) -> SeqIO.SeqRecord:
        top_score = 0
        seq = None
        for fasta in self.fasta_seqs:
            rat = max([process.fuzz.ratio(rcsb_id, fasta.description),
                      process.fuzz.ratio(rcsb_id, fasta.id)])
            if rat > top_score:
                top_score = rat
                seq = fasta
        return seq


def add_target_to_domain_alignment(rcsb_id: str, domain: str):
    """
    @param rcsb_id: PDB ID of target structure
    @param domain:  Domain of target structure (euk or bac)
    """

    rcsb_id = argdict["target"]
    struct_profile = open_structure(rcsb_id, 'json')
    [chain_id, strand_target] = get_23SrRNA_strandseq(
        rcsb_id,
        custom_path=os.path.join(
            RIBETL_DATA, rcsb_id.upper(), f"{rcsb_id.upper()}_modified.cif")
    )

    fpath_23s = f'{rcsb_id.upper()}_{chain_id}_23SrRNA.fasta'
    domain_alignment: str = ''
    if domain == 'bacteria':
        domain_alignment = 'data/ribovision.bacteria.fasta'
    elif domain == 'eukarya':
        domain_alignment = 'data/ribovision.eukaryota.fasta'
    else:
        raise FileNotFoundError(
            "Domain misspecified. Must be either 'bacteria' or 'eukarya'.")

    seq_to_fasta(rcsb_id, strand_target, fpath_23s)
    muscle_combine_profile(domain_alignment, fpath_23s,
                           f'combined_{rcsb_id.upper()}_ribovision_{domain}.fasta')

def seq_to_fasta(rcsb_id: str, _seq: str, outfile: str):
    from Bio.Seq import Seq
    _seq          = _seq.replace("\n", "")
    seq_record    = SeqRecord.SeqRecord(Seq(_seq).upper())
    seq_record.id = seq_record.description = rcsb_id
    SeqIO.write(seq_record, outfile, 'fasta',)

def muscle_combine_profile(msa_path1: str, msa_path2: str, out_filepath: str):
    """Combine two MSA-profiles into a single one. Used here to "append" a target sequence two the ribovision alignment. """
    cmd = ['/home/rxz/dev/docker_ribxz/cli/scripts/muscle3.8', '-profile','-in1', msa_path1, '-in2', msa_path2, '-out', out_filepath]
    subprocess.Popen(cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE, env=os.environ.copy()).wait()
    sys.stdout.flush()

def retrieve_aligned_23s(rcsb_id: str, domain: str):
    if domain == 'b':
        domain_alignment_path = 'data/ribovision.bacteria.fasta'
    elif domain == 'e':
        domain_alignment_path = 'data/ribovision.eukaryota.fasta'
    else:
        raise FileNotFoundError(
            "Domain misspecified. Must be either 'bacteria' or 'eukarya'.")

    ribovision = RibovisionAlignment(domain_alignment_path)
    target_seq_record = ribovision.find_aln_by_id(rcsb_id)

    return target_seq_record

def pick_match(ms, rna_length: int):
    if len(ms) == 0:
        print("No matches found!!")
        exit(1)
    """Pick the match that is closest to the 3' end of the rRNA."""
    best = None
    farthest_dist = 0

    if len(ms) > 1:
        for m in ms:
            if rna_length - ((m.start + m.end) / 2) < farthest_dist:
                best = m
        return best
    else:
        return ms[0]

if args.fuzzy:

    domain = 'b'
    rcsb_id = args.target.upper()

    [chain_id, strand_target, rna_type] = retrieve_LSU_rRNA(rcsb_id)

    print("Target rRNA: {}".format(chain_id))
    ptc_projected = {"site_9": []}


    rnaS = SeqIO.SeqRecord(strand_target)
    matches_site9 = find_near_matches(
        DORIS_ET_AL["subseq_9"], strand_target, max_l_dist=1)
    print("Got {} matches for {}".format(len(matches_site9), rna_type))
    best_match9 = pick_match(matches_site9, len(rnaS))
    ptc_projected["site_9"] = [*range(best_match9.start, best_match9.end)]
    report = {
        chain_id: {
            "site_9": ptc_projected["site_9"]
        }
    }

    for res in [* map(lambda x: strand_target[x], report[chain_id]["site_9"])]:
        print("Res", res)

    ptcs_dir = os.path.join(RIBETL_DATA, "PTC_COORDINATES")
    ptc_fuzzy_path = os.path.join(
        ptcs_dir, f"{rcsb_id.upper()}_FUZZY_PTC.json")

    with open(ptc_fuzzy_path, 'w') as f:
        json.dump(report, f, indent=4)
        print("[Saved {} successfully.]".format(ptc_fuzzy_path))
        exit(1)

if args.canon:
    domain = 'bacteria'
    rcsb_id = argdict["target"].upper()
    struct_profile: Structure = open_structure(rcsb_id, 'cif')
    [chain_id, strand_target, rna_type] = retrieve_LSU_rRNA(rcsb_id, False)

    if chain_id in struct_profile.child_dict[0].child_dict:
        chain3d: Chain = struct_profile.child_dict[0].child_dict[chain_id]
    else:
        chain3d: Chain = struct_profile.child_dict[1].child_dict[chain_id]

    ress:List[Residue] = chain3d.child_list
    _r:Residue

    ress_sanitized: List[Residue] = [*filter(lambda r: r.get_resname() in ["A", "C", "G", "U", "PSU"], ress)]

    for _r in ress_sanitized:
        if _r.get_resname() == "PSU":
            _r.resname = "U"


    raw_seq         = reduce(lambda x,y: x + y.resname, ress_sanitized,'')
    matches         = find_near_matches(DORIS_ET_AL["subseq_9"], raw_seq, max_l_dist=1)
    m0              = pick_match(matches, len(raw_seq))
    sought_residues = [ ress_sanitized[i] for i in list(range(m0.start, m0.end))]

    PTC_MARKERKS_RAW = {}
    for res in sought_residues:
        if res.id[1] not in PTC_MARKERKS_RAW:
            PTC_MARKERKS_RAW[res.id[1]] = {}
        atom:Atom
        for atom in res.child_list:
            atom_name = atom.name
            atom_coords = atom.get_coord()
            PTC_MARKERKS_RAW[res.id[1]][atom_name] = list(map(lambda x: float(x), list(atom_coords)))

    PTC_MARKERKS_RAW = {
        chain_id: PTC_MARKERKS_RAW
    }

    markers_dir = os.path.join(RIBETL_DATA, "PTC_MARKERS_RAW")
    outfile     = os.path.join(markers_dir, f"{rcsb_id.upper()}_PTC_MARKERS_RAW.json")

    with open(outfile, 'w') as outf:
        json.dump(PTC_MARKERKS_RAW, outf, indent=4)
        print("Saved {} successfully.".format(outfile))

if args.markers:
    domain = 'bacteria'
    rcsb_id = argdict["target"].upper()
    struct_profile: Structure = open_structure(rcsb_id, 'cif')
    [chain_id, strand_target, rna_type] = retrieve_LSU_rRNA(rcsb_id)

    if chain_id in struct_profile.child_dict[0].child_dict:
        rnas: Chain = struct_profile.child_dict[0].child_dict[chain_id]
    else:
        rnas: Chain = struct_profile.child_dict[1].child_dict[chain_id]

    x: Residue
    SITE_DICT    = {}

    for res in rnas.child_list:
        res: Residue
        seq_id_raw = res.id[1]
        if seq_id_raw in LANDMARK_IDS:

            if seq_id_raw not in SITE_DICT:
                SITE_DICT[seq_id_raw] = {}
            for atom in res.child_dict.items():
                atom_name = atom[0]
                atom_coords = atom[1].get_coord()
                SITE_DICT[seq_id_raw][atom_name] = list(map(lambda x: float(x), list(atom_coords)))

    pprint(SITE_DICT)
    markers_dir = os.path.join(RIBETL_DATA, "PTC_MARKERS")
    outfile     = os.path.join(markers_dir, f"{rcsb_id.upper()}_PTC_MARKERS.json")

    with open(outfile, 'w') as outf:
        json.dump(SITE_DICT, outf, indent=4)
        print("Saved {} successfully.".format(outfile))

if args.ptc:

    domain                    = 'b'
    rcsb_id                   = args.target.upper()
    [chain_id, strand_target] = get_23SrRNA_strandseq(
        rcsb_id,
        custom_path=os.path.join(
            RIBETL_DATA, rcsb_id.upper(), f"{rcsb_id.upper()}_modified.cif")
    )

    tgt_seq = retrieve_aligned_23s(rcsb_id, domain)

    ptc_projected = {
        "site_6": [],
        "site_8": [],
        "site_9": []
    }

    report = {
        chain_id: {
            "site_6": [],
            "site_8": [],
            "site_9": []
        }
    }
    for (site, resids) in DORIS_ET_AL[domain].items():
        for resid in resids:
            report[chain_id][site].append(
                util__backwards_match(tgt_seq.seq, resid))

    pprint(report)

    pprint("Site 6:" +
           "".join([* map(lambda x: x[1], report[chain_id]["site_6"])]))
    pprint("Site 8:" +
           "".join([* map(lambda x: x[1], report[chain_id]["site_8"])]))
    pprint("Site 9:" +
           "".join([* map(lambda x: x[1], report[chain_id]["site_9"])]))

    with open(f'/home/rxz/dev/docker_ribxz/cli/scripts/PTC_COORDINATES/{rcsb_id.upper()}_PTC.json', 'w') as f:
        json.dump(report, f, indent=4)

    exit(1)

if args.fasta_profile:
    domain = 'bacteria'
    rcsb_id = argdict["target"]
    struct_profile = open_structure(rcsb_id, 'json')
    [chain_id, strand_target] = get_23SrRNA_strandseq(rcsb_id, custom_path=os.path.join(
        RIBETL_DATA, rcsb_id.upper(), f"{rcsb_id.upper()}_modified.cif"))

    fpath_23s = f'{rcsb_id.upper()}_{chain_id}_23SrRNA.fasta'
    domain_alignment: str = ''

    if domain == 'bacteria':
        domain_alignment = 'data/ribovision.bacteria.fasta'

    elif domain == 'eukarya':
        domain_alignment = 'data/ribovision.eukaryota.fasta'

    else:
        raise FileNotFoundError(
            "Domain misspecified. Must be either 'bacteria' or 'eukarya'.")
    seq_to_fasta(rcsb_id, strand_target, fpath_23s)
    muscle_combine_profile(domain_alignment, fpath_23s,'ribovision.bacteria.fasta')
    exit(1)

if args.generate:
    f    = open('rcsb_pdb_ids_20230106032038.txt', 'r')
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

# if not args.display_all:
    # print("\nTo display more residues per target structure, use additional --display_all flag.")



# 5afi | A| resi 2610-2611 |~ C3
# 3j7z | sele c. A and resi 2609 and name C4 | O4 | n4