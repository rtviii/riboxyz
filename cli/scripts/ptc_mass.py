
#!/usr/bin/env python3
# kdd@math.ubc.ca
# rtkushner@gmail.com
# This script is courtesy of the ribosome.xyz and its authors.
# This relies on the following packages to run
# - gemmi   : https://gemmi.readthedocs.io/en/latest/install.html
# - bipython: https://biopython.org/
# And additionally "requests" to download missing structures: https://pypi.org/project/requests/

# Distribute freely.
import os
import sys
from Bio import SeqRecord
from pprint import pprint
import subprocess
from typing import List, Tuple
from fuzzywuzzy import process
import re
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
DORIS_ET_AL = {
    "euk": {
        "site_6": [2400, 2401, 2402, 2403, 2404, 2405, 2406, 2407],
        "site_8": [2811, 2812, 2813, 2814, 2815, 2816, 2817, 2818],
        "site_9": [2941, 2942, 2943, 2944, 2945, 2946, 2947, 2948, 2949, 2950, 2951, 2952, 2953]
    },
    'bact': {
        "site_6": [2059, 2060, 2061, 2062, 2063, 2064, 2065, 2066],
        "site_8": [2446, 2447, 2448, 2449, 2450, 2451, 2452, 2453],
        "site_9": [2576, 2577, 2578, 2579, 2580, 2581, 2582, 2583, 2584, 2585, 2586, 2587, 2588]
    }
}

# -------------------------------------------------------------------------------------------|

parser = argparse.ArgumentParser(
    description='CLI for locating PTC residues of 23SrRNA in a given prokaryotic PDB file')

parser.add_argument("-t", "--target", type=str, required=False)
parser.add_argument("-d", "--domain", type=str, required=False, choices=['e', 'b'])
parser.add_argument("-m", "--mode",
                    type=str,  choices=['position', 'residue'], nargs='?', default='position', help="Display the position of the PTC residue, or the residue itself. Default: position")
parser.add_argument("--display_all"  ,      action          ='store_true')
parser.add_argument("--generate"     ,      action          ='store_true')
parser.add_argument("--fasta_profile",      action          ='store_true')
parser.add_argument("--ptc"        ,    action='store_true')
parser.add_argument("--batch"        ,    type=int,        required=False)

args               = parser .parse_args()
argdict            = vars(parser.parse_args())

bact_registry_file = open('rcsb_pdb_ids_20230106032038.txt', 'r')
line               = bact_registry_file.readline()
bact_registry_file.close()
bacteria_structs = line.split(',')


def util__backwards_match(alntgt: str, aln_resid: int, verbose:bool=False)->Tuple[int, str]:
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
                print("[ {} ] <-----> id.[aligned: {} | orgiginal: {} ]".format(alntgt[aln_resid],i, counter_proper))
            return ( counter_proper,  alntgt[i] )
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


def get_23SrRNA_strandseq(rcsb_id: str, custom_path=None) -> Tuple[str, str]:
    return get_one_letter_code_can_by_nomclass(rcsb_id, "23SrRNA", custom_path)


def get_one_letter_code_can_by_nomclass(rcsb_id: str, nomenclature_class: str, custom_path=None) -> Tuple[str, str]:

    default_path = f"{rcsb_id.upper()}_modified.cif" if custom_path == None else custom_path

    target = gemmi.cif.read_file(default_path)
    block  = target.sole_block()
    model  = gemmi.read_structure(default_path)[0]

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


class RibovisionAlignment:

    def __init__(self, mass_alignment_path:str) -> None:
        self.fasta_seqs: List[SeqIO.SeqRecord] = SeqIO.parse(open(mass_alignment_path), 'fasta')
        pass

    def find_aln_by_id(self, rcsb_id: str) -> SeqIO.SeqRecord:
        top_score = 0
        seq       = None
        for fasta in self.fasta_seqs:
            rat = max([process.fuzz.ratio(rcsb_id, fasta.description), process.fuzz.ratio(rcsb_id, fasta.id) ])
            if rat > top_score:
                top_score = rat
                seq       = fasta
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
        domain_alignment = 'ribovision.bacteria.fasta'
    elif domain == 'eukarya':
        domain_alignment = 'ribovision.eukaryota.fasta'
    else:
        raise FileNotFoundError(
            "Domain misspecified. Must be either 'bacteria' or 'eukarya'.")

    seq_to_fasta(rcsb_id, strand_target, fpath_23s)
    muscle_combine_profile(domain_alignment, fpath_23s,
                           f'combined_{rcsb_id.upper()}_ribovision_{domain}.fasta')


def seq_to_fasta(rcsb_id: str, _seq: str, outfile: str):
    from Bio.Seq import Seq
    _seq = _seq.replace("\n","")
    seq_record = SeqRecord.SeqRecord(Seq(_seq).upper())
    seq_record.id = seq_record.description = rcsb_id
    SeqIO.write(seq_record, outfile, 'fasta',)


def muscle_combine_profile(msa_path1: str, msa_path2: str, out_filepath: str):
    """Combine two MSA-profiles into a single one. Used here to "append" a target sequence two the ribovision alignment. """
    cmd = ['/home/rxz/dev/docker_ribxz/cli/scripts/muscle3.8',  '-profile', '-in1', msa_path1,'-in2', msa_path2, '-out', out_filepath]
    subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=os.environ.copy()).wait()
    sys.stdout.flush()

def get_seq_from_profiles(rcsb_id:str, domain:str):
    if domain == 'b':
        domain_alignment_path = 'ribovision.bacteria.fasta'
    elif domain == 'e':
        domain_alignment_path = 'ribovision.eukaryota.fasta'
    else:
        raise FileNotFoundError("Domain misspecified. Must be either 'bacteria' or 'eukarya'.")

    ribovision        = RibovisionAlignment(domain_alignment_path)
    target_seq_record = ribovision.find_aln_by_id(rcsb_id)

    return target_seq_record


if args.ptc:
    domain  = 'b'
    rcsb_id = args.target.upper()
    tgt_seq = get_seq_from_profiles(rcsb_id, domain)
    # print(tgt_seq.seq[2446:2456])

    for i in range(2576,2586):
        [ resid,res_projected ] = util__backwards_match(tgt_seq.seq, i, True)
    exit(1)

if args.fasta_profile:
    domain = 'bacteria'

    rcsb_id = argdict["target"]
    struct_profile = open_structure(rcsb_id, 'json')
    [chain_id, strand_target] = get_23SrRNA_strandseq(rcsb_id,custom_path=os.path.join(RIBETL_DATA, rcsb_id.upper(), f"{rcsb_id.upper()}_modified.cif"))

    fpath_23s             = f'{rcsb_id.upper()}_{chain_id}_23SrRNA.fasta'
    domain_alignment: str = ''

    if domain == 'bacteria':
        domain_alignment = 'ribovision.bacteria.fasta'

    elif domain == 'eukarya':
        domain_alignment = 'ribovision.eukaryota.fasta'

    else:
        raise FileNotFoundError("Domain misspecified. Must be either 'bacteria' or 'eukarya'.")
    seq_to_fasta(rcsb_id,strand_target,fpath_23s)
    muscle_combine_profile(domain_alignment, fpath_23s,'ribovision.bacteria.fasta')
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

if "targets" in argdict.keys():
    argdict["targets"] = [s.strip().upper()
                          for s in argdict["targets"].split(",")]

    if len(argdict) > 50:
        print("Please don't overload our servers. Paid out of pocket!:) \nInstead, get in touch for collaboration: rtkushner@gmail.com!")
        exit(1)

    for target in argdict["targets"]:

        if not args.display_all:
            target_ptc = process_target(target, args.mode)
            print("[\033[94m{}\033[0m] Approximate PTC position(1 of {} residues): \033[91m{}\033[0m".format(
                target, len(target_ptc), target_ptc[0]))

        else:
            print("[\033[94m{}\033[0m] PTC atom positions: ".format(target))
            for residue in process_target(target, args.mode):
                print(f"\t\033[91m{residue}\033[0m")

if not args.display_all:
    print("\nTo display more residues per target structure, use additional --display_all flag.")