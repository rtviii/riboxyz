
#!/usr/bin/env python3
# kdd@math.ubc.ca
# rtkushner@gmail.com
# This script is courtesy of the ribosome.xyz and its authors.
# This relies on the following packages to run
# - gemmi   : https://gemmi.readthedocs.io/en/latest/install.html
# - bipython: https://biopython.org/
# And additionally "requests" to download missing structures: https://pypi.org/project/requests/

# Distribute freely.
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
from extract_bsites import open_structure
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

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4574749/pdf/1719.pdf

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
parser.add_argument("--fuzzy", action='store_true')
parser.add_argument("--markers", action='store_true')
parser.add_argument("--fasta_profile", action='store_true')
parser.add_argument("--ptc", action='store_true')
parser.add_argument("--batch", type=int, required=False)

args = parser .parse_args()
argdict = vars(parser.parse_args())

bact_registry_file = open('rcsb_pdb_ids_20230106032038.txt', 'r')
line               = bact_registry_file.readline()
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

# ※ ---------------------------- 23SrRNA PTC residue locations ---------------------------- ※


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

def get_23SrRNA_strandseq(rcsb_id: str, custom_path=None) -> Tuple[str, str]:
    return get_one_letter_code_can_by_nomclass(rcsb_id, "23SrRNA", custom_path)


def get_one_letter_code_can_by_nomclass(rcsb_id: str, nomenclature_class: str, custom_path=None) -> Tuple[str, str]:

    default_path = f"{rcsb_id.upper()}_modified.cif" if custom_path == None else custom_path

    target = gemmi.cif.read_file(default_path)
    block = target.sole_block()
    model = gemmi.read_structure(default_path)[0]

    STRAND = None
    SEQ = None

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
        block.find_loop('_entity_poly.pdbx_seq_one_letter_code_can')
    ):
        # X-RAY structures have 'dual' chains. Split on comma to check both.
        if STRAND in chain_id.split(','):
            SEQ = str(one_letter_code).strip(";").strip("\n")

    if SEQ == None:
        print("Could not locate {} sequence in {} CIF file".format(
            nomenclature_class, rcsb_id))

    print("Located {} sequence in {} CIF file. (Chain {})".format(
        nomenclature_class, rcsb_id, STRAND))
    return (STRAND, SEQ)

class rRNA23S(SeqIO.SeqRecord):

    def __init__(self, seq):
        super().__init__(seq)

    def fuzzy_search_subseq(self, subseq: str):
        print("Fuzz seq", self._seq)
        print("Fuzz seq")


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
    _seq = _seq.replace("\n", "")
    seq_record = SeqRecord.SeqRecord(Seq(_seq).upper())
    seq_record.id = seq_record.description = rcsb_id
    SeqIO.write(seq_record, outfile, 'fasta',)


def muscle_combine_profile(msa_path1: str, msa_path2: str, out_filepath: str):
    """Combine two MSA-profiles into a single one. Used here to "append" a target sequence two the ribovision alignment. """
    cmd = ['/home/rxz/dev/docker_ribxz/cli/scripts/muscle3.8',  '-profile',
           '-in1', msa_path1, '-in2', msa_path2, '-out', out_filepath]
    subprocess.Popen(cmd, stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE, env=os.environ.copy()).wait()
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





if args.fuzzy:

    domain = 'b'
    rcsb_id = args.target.upper()

    [chain_id, strand_target] = get_23SrRNA_strandseq(
        rcsb_id,
        custom_path=os.path.join(
            RIBETL_DATA, rcsb_id.upper(), f"{rcsb_id.upper()}_modified.cif")
    )

    ptc_projected = {
        # "site_6": [],
        # "site_8": [],
        "site_9": []
    }

    def pick_match(ms, rna_length: int):
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

    # tgt_seq = retrieve_aligned_23s(rcsb_id, domain)
    rna23s = SeqIO.SeqRecord(strand_target)
    # found6 = find_near_matches(
    #     DORIS_ET_AL["subseq_6"], strand_target, max_l_dist=0)
    # found8 = find_near_matches(
    #     DORIS_ET_AL["subseq_8"], strand_target, max_l_dist=0)
    found9 = find_near_matches(
        DORIS_ET_AL["subseq_9"], strand_target, max_l_dist=0)

    # best_match6 = pick_match(found6, len(rna23s))
    # best_match8 = pick_match(found8, len(rna23s))
    best_match9 = pick_match(found9, len(rna23s))

    print("\tReturned matches ", )
    # print(best_match6)
    # print(best_match8)
    print(best_match9)

    # ptc_projected["site_6"] = [*range(best_match6.start, best_match6.end)]
    # ptc_projected["site_8"] = [*range(best_match8.start, best_match8.end)]
    ptc_projected["site_9"] = [*range(best_match9.start, best_match9.end)]


    report = {
            chain_id:{
                # "site_6": ptc_projected["site_6"],
                # "site_8": ptc_projected["site_8"],
                "site_9": ptc_projected["site_9"]
            }
        }
    pprint(ptc_projected)


    # pprint("Site 6:" +
    #        "".join([* map(lambda x: strand_target[x], report[chain_id]["site_6"])]))
    # pprint("Site 8:" +
    #        "".join([* map(lambda x: strand_target[x], report[chain_id]["site_8"])]))
    for res in [* map(lambda x: strand_target[x], report[chain_id]["site_9"])]:
        print("Res", res)


    ptcs_dir       = os.path.join(RIBETL_DATA, "PTC_COORDINATES")
    ptc_fuzzy_path = os.path.join(ptcs_dir, f"{rcsb_id.upper()}_FUZZY_PTC.json")

    with open(ptc_fuzzy_path, 'w') as f:
        json.dump(report, f, indent=4)
    print("[Saved {} successfully.]".format(ptc_fuzzy_path))
    
if args.markers:
    domain                    = 'bacteria'
    rcsb_id                   = argdict["target"].upper()
    struct_profile:Structure            = open_structure(rcsb_id, 'cif')
    print(struct_profile.child_dict[0])
    model:Model =struct_profile.child_dict[0]
    rna23s:Chain = model.child_dict['A']
    x:Residue;
    LANDMARK_IDS = [2610,2611,2612]
    SITE_DICT    = {}

    for res in rna23s.child_list:
        res:Residue
        seq_id_raw = res.id[1]
        if seq_id_raw in LANDMARK_IDS:

            if seq_id_raw not in SITE_DICT: 
                SITE_DICT[seq_id_raw] = {}

            for atom in res.child_dict.items():
                atom_name   = atom[0]
                atom_coords = atom[1].get_coord();
                SITE_DICT[seq_id_raw][atom_name] = list(map(lambda x : float(x),list(atom_coords)))

    pprint(SITE_DICT)
    markers_dir = os.path.join(RIBETL_DATA, "PTC_MARKERS")
    outfile     = os.path.join(markers_dir, f"{rcsb_id.upper()}_PTC_MARKERS.json")
    with open(outfile, 'w') as outf:
        json.dump(SITE_DICT, outf, indent=4)

        print("Saved {} successfully.".format(outfile))

    


    


if args.ptc:

    domain = 'b'
    rcsb_id = args.target.upper()

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
    domain                    = 'bacteria'
    rcsb_id                   = argdict["target"]
    struct_profile            = open_structure(rcsb_id, 'json')
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
    muscle_combine_profile(domain_alignment, fpath_23s,
                           'ribovision.bacteria.fasta')
    exit(1)

if args.generate:
    f = open('rcsb_pdb_ids_20230106032038.txt', 'r')
    line = f.readline()
    f.close()
    bacteria_structs = line.split(',')

    struct_ids = []
    parent_chain = []
    residue_name = []
    residue_seqid = []
    residue_x = []
    residue_y = []
    residue_z = []
    i = 0
    for struct in bacteria_structs[args.batch:args.batch+100]:
        struct = str.upper(struct)
        # get the ptc guess
        structpath = os.path.join(
            RIBETL_DATA,  struct, f"{struct}_modified.cif")

        # check whether structpath file exists

        if not os.path.isfile(structpath):
            print(f"Could not find {structpath} in RIBETL_DATA")
            continue
        try:
            ptc_guess = process_target_to_tuple(struct, args.mode, structpath)
            (res, posn, strand) = ptc_guess
            struct_ids.append(struct)
            residue_seqid.append(res.seqid)
            parent_chain.append(strand)
            residue_name.append(res.name)
            residue_x.append(posn[0])
            residue_y.append(posn[1])
            residue_z.append(posn[2])
            i = i+1
            print(f"Processed structs  : {i}")
        except:
            i = i+1
            print("err")
            print(f"Processed structs  : {i}")
            continue

    # add this to the dataframe
    df = pandas.DataFrame.from_dict({
        'struct_ids': struct_ids,
        'parent_chain': parent_chain,
        'residue_name': residue_name,
        'residue_seqid': residue_seqid,
        'residue_x': residue_x,
        'residue_y': residue_y,
        'residue_z': residue_z
    })

    # save the dataframe to a csv
    df.to_csv(f"ptc_100xbatch={args.batch}.csv")
    exit(1)

# if not args.display_all:
    # print("\nTo display more residues per target structure, use additional --display_all flag.")


# the task for today is to establish a robust pipeline for saving the visualized regions
# n regions of the 23S rRNA are conserved.
# - datasheet with the conserved regions: 6,8,9
# - robust way to visualize it with pymol



# 5afi | A| resi 2610-2611 |~ C3
# 3j7z | sele c. A and resi 2609 and name C4 | O4 | n4