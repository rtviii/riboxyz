
#!/usr/bin/env python3
# kdd@math.ubc.ca
# rtkushner@gmail.com
# This script is courtesy of the ribosome.xyz and its authors.
# This relies on the following packages to run
# - gemmi   : https://gemmi.readthedocs.io/en/latest/install.html
# - bipython: https://biopython.org/
# And additionally "requests" to download missing structures: https://pypi.org/project/requests/

# Distribute freely.
from genericpath import exists
import os
import re
import pandas
try:
    from Bio import pairwise2
    from Bio import SeqIO
except ImportError:
    print("Please install Biopython to use this script: pip install biopython")
    exit(1)
try:
    import gemmi
except:
    print("Please install gemmi to use this script: pip install gemmi")
    exit(1)
import pathlib
import argparse

RIBETL_DATA = os.environ.get('RIBETL_DATA')

# Change these two parameters to have a different "source" sequence to align *against*  ----|
#                                                                                           |
# List of conserved nucleotide sequences on 23s-28s can be found here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4574749/pdf/1719.pdf
# Some of the PTC residues in bacterial 23SrRNA                                             |
PTC_SEQ_IDS = [2445, 2446, 2447, 2448, 2449, 2450, 2451, 2452]  # |

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4574749/pdf/1719.pdf
DORIS_ET_AL = {
    "euk":{
        "site_6": [2400, 2401, 2402, 2403, 2404, 2405,2406,2407],
        "site_8": [2811, 2812, 2813, 2814, 2815, 2816, 2817, 2818],
        "site_9": [2941, 2942, 2943, 2944, 2945, 2946, 2947, 2948, 2949,2950,2951,2952,2953]

    },
    'bact':{
        "site_6": [2059, 2060,2061,2062,2063,2064,2065,2066],
        "site_8": [2446,2447,2448,2449,2450,2451,2452,2453],
        "site_9": [2576,2577,2578,2579,2580,2581,2582,2583,2584,2585,2586,2587,2588]
    }
}

# As per PDB 3J7Z ( https://www.rcsb.org/structure/3j7z )                                   |
# ECOLI23SRRNA = "GGUUAAGCGACUAAGCGUACACGGUGGAUGCCCUGGCAGUCAGAGGCGAUGAAGGACGUGCUAAUCUGCGAUAAGCGUCGGUAAGGUGAUAUGAACCGUUAUAACCGGCGAUUUCCGAAUGGGGAAACCCAGUGUGUUUCGACACACUAUCAUUAACUGAAUCCAUAGGUUAAUGAGGCGAACCGGGGGAACUGAAACAUCUAAGUACCCCGAGGAAAAGAAAUCAACCGAGAUUCCCCCAGUAGCGGCGAGCGAACGGGGAGCAGCCCAGAGCCUGAAUCAGUGUGUGUGUUAGUGGAAGCGUCUGGAAAGGCGCGCGAUACAGGGUGACAGCCCCGUACACAAAAAUGCACAUGCUGUGAGCUCGAUGAGUAGGGCGGGACACGUGGUAUCCUGUCUGAAUAUGGGGGGACCAUCCUCCAAGGCUAAAUACUCCUGACUGACCGAUAGUGAACCAGUACCGUGAGGGAAAGGCGAAAAGAACCCCGGCGAGGGGAGUGAAAAAGAACCUGAAACCGUGUACGUACAAGCAGUGGGAGCACGCUUAGGCGUGUGACUGCGUACCUUUUGUAUAAUGGGUCAGCGACUUAUAUUCUGUAGCAAGGUUAACCGAAUAGGGGAGCCGAAGGGAAACCGAGUCUUAACUGGGCGUUAAGUUGCAGGGUAUAGACCCGAAACCCGGUGAUCUAGCCAUGGGCAGGUUGAAGGUUGGGUAACACUAACUGGAGGACCGAACCGACUAAUGUUGAAAAAUUAGCGGAUGACUUGUGGCUGGGGGUGAAAGGCCAAUCAAACCGGGAGAUAGCUGGUUCUCCCCGAAAGCUAUUUAGGUAGCGCCUCGUGAAUUCAUCUCCGGGGGUAGAGCACUGUUUCGGCAAGGGGGUCAUCCCGACUUACCAACCCGAUGCAAACUGCGAAUACCGGAGAAUGUUAUCACGGGAGACACACGGCGGGUGCUAACGUCCGUCGUGAAGAGGGAAACAACCCAGACCGCCAGCUAAGGUCCCAAAGUCAUGGUUAAGUGGGAAACGAUGUGGGAAGGCCCAGACAGCCAGGAUGUUGGCUUAGAAGCAGCCAUCAUUUAAAGAAAGCGUAAUAGCUCACUGGUCGAGUCGGCCUGCGCGGAAGAUGUAACGGGGCUAAACCAUGCACCGAAGCUGCGGCAGCGACGCUUAUGCGUUGUUGGGUAGGGGAGCGUUCUGUAAGCCUGCGAAGGUGUGCUGUGAGGCAUGCUGGAGGUAUCAGAAGUGCGAAUGCUGACAUAAGUAACGAUAAAGCGGGUGAAAAGCCCGCUCGCCGGAAGACCAAGGGUUCCUGUCCAACGUUAAUCGGGGCAGGGUGAGUCGACCCCUAAGGCGAGGCCGAAAGGCGUAGUCGAUGGGAAACAGGUUAAUAUUCCUGUACUUGGUGUUACUGCGAAGGGGGGACGGAGAAGGCUAUGUUGGCCGGGCGACGGUUGUCCCGGUUUAAGCGUGUAGGCUGGUUUUCCAGGCAAAUCCGGAAAAUCAAGGCUGAGGCGUGAUGACGAGGCACUACGGUGCUGAAGCAACAAAUGCCCUGCUUCCAGGAAAAGCCUCUAAGCAUCAGGUAACAUCAAAUCGUACCCCAAACCGACACAGGUGGUCAGGUAGAGAAUACCAAGGCGCUUGAGAGAACUCGGGUGAAGGAACUAGGCAAAAUGGUGCCGUAACUUCGGGAGAAGGCACGCUGAUAUGUAGGUGAGGUCCCUCGCGGAUGGAGCUGAAAUCAGUCGAAGAUACCAGCUGGCUGCAACUGUUUAUUAAAAACACAGCACUGUGCAAACACGAAAGUGGACGUAUACGGUGUGACGCCUGCCCGGUGCCGGAAGGUUAAUUGAUGGGGUUAGCGCAAGCGAAGCUCUUGAUCGAAGCCCCGGUAAACGGCGGCCGUAACUAUAACGGUCCUAAGGUAGCGAAAUUCCUUGUCGGGUAAGUUCCGACCUGCACGAAUGGCGUAAUGAUGGCCAGGCUGUCUCCACCCGAGACUCAGUGAAAUUGAACUCGCUGUGAAGAUGCAGUGUACCCGCGGCAAGACGGAAAGACCCCGUGAACCUUUACUAUAGCUUGACACUGAACAUUGAGCCUUGAUGUGUAGGAUAGGUGGGAGGCUUUGAAGUGUGGACGCCAGUCUGCAUGGAGCCGACCUUGAAAUACCACCCUUUAAUGUUUGAUGUUCUAACGUUGACCCGUAAUCCGGGUUGCGGACAGUGUCUGGUGGGUAGUUUGACUGGGGCGGUCUCCUCCUAAAGAGUAACGGAGGAGCACGAAGGUUGGCUAAUCCUGGUCGGACAUCAGGAGGUUAGUGCAAUGGCAUAAGCCAGCUUGACUGCGAGCGUGACGGCGCGAGCAGGUGCGAAAGCAGGUCAUAGUGAUCCGGUGGUUCUGAAUGGAAGGGCCAUCGCUCAACGGAUAAAAGGUACUCCGGGGAUAACAGGCUGAUACCGCCCAAGAGUUCAUAUCGACGGCGGUGUUUGGCACCUCGAUGUCGGCUCAUCACAUCCUGGGGCUGAAGUAGGUCCCAAGGGUAUGGCUGUUCGCCAUUUAAAGUGGUACGCGAGCUGGGUUUAGAACGUCGUGAGACAGUUCGGUCCCUAUCUGCCGUGGGCGCUGGAGAACUGAGGGGGGCUGCUCCUAGUACGAGAGGACCGGAGUGGACGCAUCACUGGUGUUCGGGUUGUCAUGCCAAUGGCACUGCCCGGUAGCUAAAUGCGGAAGAGAUAAGUGCUGAAAGCAUCUAAGCACGAAACUUGCCCCGAGAUGAGUUCUCCCUGACCCUUUAAGGGUCCUGAAGGAACGUUGAAGACGACGACGUUGAUAGGCCGGGUGUGUAAGCGCAGCGAUGCGUUGAGCUAACCGGUACUAAUGAACCGUGAGGCUUAACCU"
RNA_3J7Z_23S = "GGUUAAGCGACUAAGCGUACACGGUGGAUGCCCUGGCAGUCAGAGGCGAUGAAGGACGUGCUAAUCUGCGAUAAGCGUCGGUAAGGUGAUAUGAACCGUUAUAACCGGCGAUUUCCGAAUGGGGAAACCCAGUGUGUUUCGACACACUAUCAUUAACUGAAUCCAUAGGUUAAUGAGGCGAACCGGGGGAACUGAAACAUCUAAGUACCCCGAGGAAAAGAAAUCAACCGAGAUUCCCCCAGUAGCGGCGAGCGAACGGGGAGCAGCCCAGAGCCUGAAUCAGUGUGUGUGUUAGUGGAAGCGUCUGGAAAGGCGCGCGAUACAGGGUGACAGCCCCGUACACAAAAAUGCACAUGCUGUGAGCUCGAUGAGUAGGGCGGGACACGUGGUAUCCUGUCUGAAUAUGGGGGGACCAUCCUCCAAGGCUAAAUACUCCUGACUGACCGAUAGUGAACCAGUACCGUGAGGGAAAGGCGAAAAGAACCCCGGCGAGGGGAGUGAAAAAGAACCUGAAACCGUGUACGUACAAGCAGUGGGAGCACGCUUAGGCGUGUGACUGCGUACCUUUUGUAUAAUGGGUCAGCGACUUAUAUUCUGUAGCAAGGUUAACCGAAUAGGGGAGCCGAAGGGAAACCGAGUCUUAACUGGGCGUUAAGUUGCAGGGUAUAGACCCGAAACCCGGUGAUCUAGCCAUGGGCAGGUUGAAGGUUGGGUAACACUAACUGGAGGACCGAACCGACUAAUGUUGAAAAAUUAGCGGAUGACUUGUGGCUGGGGGUGAAAGGCCAAUCAAACCGGGAGAUAGCUGGUUCUCCCCGAAAGCUAUUUAGGUAGCGCCUCGUGAAUUCAUCUCCGGGGGUAGAGCACUGUUUCGGCAAGGGGGUCAUCCCGACUUACCAACCCGAUGCAAACUGCGAAUACCGGAGAAUGUUAUCACGGGAGACACACGGCGGGUGCUAACGUCCGUCGUGAAGAGGGAAACAACCCAGACCGCCAGCUAAGGUCCCAAAGUCAUGGUUAAGUGGGAAACGAUGUGGGAAGGCCCAGACAGCCAGGAUGUUGGCUUAGAAGCAGCCAUCAUUUAAAGAAAGCGUAAUAGCUCACUGGUCGAGUCGGCCUGCGCGGAAGAUGUAACGGGGCUAAACCAUGCACCGAAGCUGCGGCAGCGACGCUUAUGCGUUGUUGGGUAGGGGAGCGUUCUGUAAGCCUGCGAAGGUGUGCUGUGAGGCAUGCUGGAGGUAUCAGAAGUGCGAAUGCUGACAUAAGUAACGAUAAAGCGGGUGAAAAGCCCGCUCGCCGGAAGACCAAGGGUUCCUGUCCAACGUUAAUCGGGGCAGGGUGAGUCGACCCCUAAGGCGAGGCCGAAAGGCGUAGUCGAUGGGAAACAGGUUAAUAUUCCUGUACUUGGUGUUACUGCGAAGGGGGGACGGAGAAGGCUAUGUUGGCCGGGCGACGGUUGUCCCGGUUUAAGCGUGUAGGCUGGUUUUCCAGGCAAAUCCGGAAAAUCAAGGCUGAGGCGUGAUGACGAGGCACUACGGUGCUGAAGCAACAAAUGCCCUGCUUCCAGGAAAAGCCUCUAAGCAUCAGGUAACAUCAAAUCGUACCCCAAACCGACACAGGUGGUCAGGUAGAGAAUACCAAGGCGCUUGAGAGAACUCGGGUGAAGGAACUAGGCAAAAUGGUGCCGUAACUUCGGGAGAAGGCACGCUGAUAUGUAGGUGAGGUCCCUCGCGGAUGGAGCUGAAAUCAGUCGAAGAUACCAGCUGGCUGCAACUGUUUAUUAAAAACACAGCACUGUGCAAACACGAAAGUGGACGUAUACGGUGUGACGCCUGCCCGGUGCCGGAAGGUUAAUUGAUGGGGUUAGCGCAAGCGAAGCUCUUGAUCGAAGCCCCGGUAAACGGCGGCCGUAACUAUAACGGUCCUAAGGUAGCGAAAUUCCUUGUCGGGUAAGUUCCGACCUGCACGAAUGGCGUAAUGAUGGCCAGGCUGUCUCCACCCGAGACUCAGUGAAAUUGAACUCGCUGUGAAGAUGCAGUGUACCCGCGGCAAGACGGAAAGACCCCGUGAACCUUUACUAUAGCUUGACACUGAACAUUGAGCCUUGAUGUGUAGGAUAGGUGGGAGGCUUUGAAGUGUGGACGCCAGUCUGCAUGGAGCCGACCUUGAAAUACCACCCUUUAAUGUUUGAUGUUCUAACGUUGACCCGUAAUCCGGGUUGCGGACAGUGUCUGGUGGGUAGUUUGACUGGGGCGGUCUCCUCCUAAAGAGUAACGGAGGAGCACGAAGGUUGGCUAAUCCUGGUCGGACAUCAGGAGGUUAGUGCAAUGGCAUAAGCCAGCUUGACUGCGAGCGUGACGGCGCGAGCAGGUGCGAAAGCAGGUCAUAGUGAUCCGGUGGUUCUGAAUGGAAGGGCCAUCGCUCAACGGAUAAAAGGUACUCCGGGGAUAACAGGCUGAUACCGCCCAAGAGUUCAUAUCGACGGCGGUGUUUGGCACCUCGAUGUCGGCUCAUCACAUCCUGGGGCUGAAGUAGGUCCCAAGGGUAUGGCUGUUCGCCAUUUAAAGUGGUACGCGAGCUGGGUUUAGAACGUCGUGAGACAGUUCGGUCCCUAUCUGCCGUGGGCGCUGGAGAACUGAGGGGGGCUGCUCCUAGUACGAGAGGACCGGAGUGGACGCAUCACUGGUGUUCGGGUUGUCAUGCCAAUGGCACUGCCCGGUAGCUAAAUGCGGAAGAGAUAAGUGCUGAAAGCAUCUAAGCACGAAACUUGCCCCGAGAUGAGUUCUCCCUGACCCUUUAAGGGUCCUGAAGGAACGUUGAAGACGACGACGUUGAUAGGCCGGGUGUGUAAGCGCAGCGAUGCGUUGAGCUAACCGGUACUAAUGAACCGUGAGGCUUAACCU"
# -------------------------------------------------------------------------------------------|

parser = argparse.ArgumentParser(description='CLI for locating PTC residues of 23SrRNA in a given prokaryotic PDB file')
parser.add_argument("-t","--targets", type=str, required=True)
parser.add_argument("-m","--mode",type=str,  choices=['position', 'residue'], nargs='?', default='position',help="Display the position of the PTC residue, or the residue itself. Default: position")
parser.add_argument("--display_all",action='store_true')
parser.add_argument("--generate",action='store_true')
parser.add_argument("--batch",type=int, required=False)

args    = parser.parse_args()
argdict = vars(parser.parse_args())

if "targets" in argdict.keys():
    argdict["targets"] = [s.strip().upper()
                          for s in argdict["targets"].split(",")]
    if len(argdict) > 50:
        print("Please don't overload our servers. Paid out of pocket!:) \nInstead, get in touch for collaboration: rtkushner@gmail.com!")
        exit(1)

def backwards_match(alntgt: str, resid: int):
    """Returns the target-sequence index of a residue in the (aligned) target sequence"""
    if resid > len(alntgt):
        exit(IndexError(
            f"Passed residue with invalid index ({resid}) to back-match to target.Seqlen:{len(alntgt)}"))
    counter_proper = 0
    for i, char in enumerate(alntgt):
        if i == resid:
            return counter_proper
        if char == '-':
            continue
        else:
            counter_proper += 1

def forwards_match(alnsrc: str, resid: int):
    """Returns the index of a source-sequence residue in the (aligned) source sequence."""
    count_proper = 0
    for alignment_indx, char in enumerate(alnsrc):
        if count_proper == resid:
            return alignment_indx
        if char == '-':
            continue
        else:
            count_proper += 1

def process_target_to_tuple(rcsb_id: str, result_as: str, custom_path=None):
    default_path = f"{rcsb_id.upper()}.cif" if custom_path == None else custom_path
    if not pathlib.Path(default_path).is_file():
        print("Could not find the fiel!!!! exti ")
        exit(1)
        # print(
        #     f"Could not locate file {default_path} in the current directory. Downloading via {f'https://api.ribosome.xyz/static_files/download_structure?struct_id={rcsb_id}'}.")
        # import requests
        # with open(default_path, 'wb') as outfile:
            # outfile.write(requests.get(
            #     f'https://api.ribosome.xyz/static_files/download_structure?struct_id={rcsb_id}').content)

    target = gemmi.cif.read_file(default_path)
    block = target.sole_block()
    model = gemmi.read_structure(default_path)[0]

    STRAND = None
    SEQ = None

    # Locate the chain of 23SrRNA class
    for (strand, nomclass) in zip(
        block.find_loop('_ribosome_nomenclature.entity_poly.pdbx_strand_id'),
        block.find_loop('_ribosome_nomenclature.polymer_class')
    ):
        if nomclass == '23SrRNA':
            STRAND = strand
            break

    # Now find sequence of this 23SrRNA
    for (chain_id, one_letter_code) in zip(
        block.find_loop('_entity_poly.pdbx_strand_id'),
        block.find_loop('_entity_poly.pdbx_seq_one_letter_code')
    ):
        # X-RAY structures have 'dual' chains. Split on comma to check both.
        if STRAND in chain_id.split(','):
            SEQ = str(one_letter_code).strip(";").strip("\n")

    if SEQ == None:
        print("Could not locate 23SrRNA sequence in {} CIF file".format(rcsb_id))

    alignment = pairwise2.align.globalxx(
        RNA_3J7Z_23S, SEQ, one_alignment_only=True)
    src_aln = alignment[0].seqA
    tgt_aln = alignment[0].seqB

    aln_ids = []
    tgt_ids = []

    for src_resid in PTC_SEQ_IDS:
        aln_ids.append(forwards_match(src_aln, src_resid))
    aln_ids = list(filter(lambda x: x != None, aln_ids))

    for aln_resid in aln_ids:
        if tgt_aln[aln_resid] == '-':
            continue
        tgt_ids.append(backwards_match(tgt_aln, aln_resid))

    guesses_residues     = [model[STRAND][ix] for ix in tgt_ids]
    guesses_alphaC_posns = [list(model[STRAND][ix][0].pos) for ix in tgt_ids]
    print(guesses_residues)
    print(guesses_alphaC_posns)
    return (guesses_residues[0], guesses_alphaC_posns[0], STRAND)

def process_target(rcsb_id: str, result_as: str, custom_path=None):
    default_path = f"{rcsb_id.upper()}.cif" if custom_path == None else custom_path
    # if not pathlib.Path(default_path).is_file():
    #     print(
    #         f"Could not locate file {default_path} in the current directory. Downloading via {f'https://api.ribosome.xyz/static_files/download_structure?struct_id={rcsb_id}'}.")
    #     import requests
    #     with open(default_path, 'wb') as outfile:
    #         outfile.write(requests.get(
    #             f'https://api.ribosome.xyz/static_files/download_structure?struct_id={rcsb_id}').content)

    target = gemmi.cif.read_file(default_path)
    block = target.sole_block()
    model = gemmi.read_structure(default_path)[0]

    STRAND = None
    SEQ = None

    # Locate the chain of 23SrRNA class
    for (strand, nomclass) in zip(
        block.find_loop('_ribosome_nomenclature.entity_poly.pdbx_strand_id'),
        block.find_loop('_ribosome_nomenclature.polymer_class')
    ):
        if nomclass == '23SrRNA':
            STRAND = strand
            break

    # Now find sequence of this 23SrRNA
    for (chain_id, one_letter_code) in zip(
        block.find_loop('_entity_poly.pdbx_strand_id'),
        block.find_loop('_entity_poly.pdbx_seq_one_letter_code')
    ):
        # X-RAY structures have 'dual' chains. Split on comma to check both.
        if STRAND in chain_id.split(','):
            SEQ = str(one_letter_code).strip(";").strip("\n")

    if SEQ == None:
        print("Could not locate 23SrRNA sequence in {} CIF file".format(rcsb_id))

    alignment = pairwise2.align.globalxx(RNA_3J7Z_23S, SEQ, one_alignment_only=True)

    src_aln = alignment[0].seqA
    tgt_aln = alignment[0].seqB

    aln_ids = []
    tgt_ids = []

    for src_resid in PTC_SEQ_IDS:
        aln_ids.append(forwards_match(src_aln, src_resid))
    aln_ids = list(filter(lambda x: x != None, aln_ids))

    for aln_resid in aln_ids:
        if tgt_aln[aln_resid] == '-':
            continue
        tgt_ids.append(backwards_match(tgt_aln, aln_resid))

    if result_as == "residue":
        return [model[STRAND][ix] for ix in tgt_ids]
    elif result_as == "position":
        print(
            "The alpha-carbon of each residue is taken to be its position for simplicity.")
        return [list(model[STRAND][ix][0].pos) for ix in tgt_ids]

if args.generate:
    f    = open('rcsb_pdb_ids_20230106032038.txt', 'r')
    line = f.readline()
    f.close()
    bacteria_structs = line.split(',')
    struct_ids       = []
    parent_chain     = []
    residue_name     = []
    residue_seqid    = []
    residue_x        = []
    residue_y        = []
    residue_z        = []
    i                = 0

    for struct in bacteria_structs[args.batch:args.batch+100]:
        struct = str.upper(struct)
        # get the ptc guess
        structpath = os.path.join(RIBETL_DATA,  struct, f"{struct}_modified.cif")

        # check whether structpath file exists
        
        if not os.path.isfile(structpath): 
            print(f"Could not find {structpath} in RIBETL_DATA")
            continue
        try:
            ptc_guess = process_target_to_tuple(struct, args.mode, structpath)
            (res, posn, strand)=ptc_guess
            struct_ids.append(struct)
            residue_seqid.append(res.seqid)
            parent_chain.append(strand)
            residue_name.append(res.name)
            residue_x.append(posn[0])
            residue_y.append(posn[1])
            residue_z.append(posn[2])
            i=i+1
            print(f"Processed structs  : {i}")
        except:
            i=i+1
            print("err")
            print(f"Processed structs  : {i}")
            continue

    # add this to the dataframe
    df = pandas.DataFrame.from_dict({
        'struct_ids'   : struct_ids   ,
        'parent_chain' : parent_chain ,
        'residue_name' : residue_name ,
        'residue_seqid': residue_seqid,
        'residue_x'    : residue_x    ,
        'residue_y'    : residue_y    ,
        'residue_z'    : residue_z
    })

    # save the dataframe to a csv
    df.to_csv(f"ptc_100xbatch={args.batch}.csv")
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
