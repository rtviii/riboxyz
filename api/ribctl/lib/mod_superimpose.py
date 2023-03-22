import argparse
import itertools
import os
from typing import Tuple
from ribctl.lib.types.types_ribosome import RibosomeStructure
from pymol import cmd
from rbxz_bend.settings import RIBETL_DATA
from ribctl.lib.mod_transpose_bsites import SeqMatch
from ribctl.lib.utils import open_structure


""" The goal is to have a module that,
given two files
and a residue range:
1. retrieves them
2. seq-aligns them
3. maps the range on the alignment
4. backtracks the new ranges to original chains
5. clips the backtracked ranges out of the chains in pymol
6. superimpose the individual snippets
"""


def ranged_align_by_polyclass(
        src_struct: str,
        tgt_struct: str,
        rng: Tuple[int, int],
        poly_class: str) -> Tuple[str, str, Tuple[int, int],
                                  str, str, Tuple[int, int]]:
    """Return a bundle of path + mapped range for a source and a target structure
    for a given polymer class. Feed this into pymol
    to chop up on the ranges and superimpose resultant snippets."""

    rstart, rend = rng

    json_src = RibosomeStructure.parse_obj(
        open_structure(src_struct.upper(), 'json'))
    json_tgt = RibosomeStructure.parse_obj(
        open_structure(tgt_struct.upper(), 'json'))

    seq_ids = {
        "src_auth_asym_id": '',
        "tgt_auth_asym_id": '',
        "src_seq": '',
        "tgt_seq": ''}

    for chain in itertools.chain(json_src.proteins, json_src.rnas if json_src.rnas else []):
        if poly_class in chain.nomenclature:
            seq_ids["src_auth_asym_id"] = chain.auth_asym_id
            seq_ids["src_seq"] = chain.entity_poly_seq_one_letter_code

    for chain in itertools.chain(json_tgt.proteins, json_tgt.rnas if json_tgt.rnas else []):
        if poly_class in chain.nomenclature:
            seq_ids["tgt_auth_asym_id"] = chain.auth_asym_id
            seq_ids["tgt_seq"] = chain.entity_poly_seq_one_letter_code

    if len(seq_ids["src_seq"]) < 1 or len(seq_ids["tgt_seq"]) < 1:
        print("""Could not retrieve either of the arguments:
			{}
		 Exiting.""".format(seq_ids))
        raise ValueError("One of the sequences is empty")

    indices = [*range(rstart, rend)]
    sm = SeqMatch(seq_ids["src_seq"], seq_ids["tgt_seq"], indices)

    target_range = (sm.tgt_ids[0], sm.tgt_ids[-1])
    source_range = (sm.src_ids[0], sm.src_ids[-1])

    src_chain_path = os.path.join(RIBETL_DATA, src_struct.upper(
    ), "CHAINS", "{}_STRAND_{}.cif".format(src_struct.upper(), seq_ids["src_auth_asym_id"]))

    tgt_chain_path = os.path.join(RIBETL_DATA, tgt_struct.upper(
    ), "CHAINS", "{}_STRAND_{}.cif".format(tgt_struct.upper(), seq_ids["tgt_auth_asym_id"]))

    print(sm.hl_ixs(sm.src, sm.src_ids))
    print("\n")
    print(sm.hl_ixs(sm.src_aln, sm.aligned_ids))
    print("\n")
    print(sm.hl_ixs(sm.tgt_aln, sm.aligned_ids))
    print("\n")
    print(sm.hl_ixs(sm.tgt, sm.tgt_ids))

    return (seq_ids["src_auth_asym_id"],
            src_chain_path,
            source_range,
            seq_ids["tgt_auth_asym_id"],
            tgt_chain_path, target_range)


def ranged_super(
        src_struct: str,
        src_auth_asym_id: str,
        tgt_struct: str,
        tgt_auth_asym_id: str,
        rng: Tuple[int, int],

) -> Tuple[str, Tuple[int, int], str, Tuple[int, int]]:
    """Return a bundle of path + mapped range for a source and a target structure
    for a given polymer class. Feed this into pymol 
    to chop up on the ranges and superimpose resultant snippets."""

    rstart, rend = rng

    json_src = open_structure(src_struct.upper(), 'json')
    json_tgt = open_structure(tgt_struct.upper(), 'json')

    tgt_seq, src_seq = [None, None]

    for chain in [*json_src['proteins'], *json_src['rnas']]:
        if src_auth_asym_id == chain['auth_asym_id']:
            src_seq = chain['entity_poly_seq_one_letter_code']

    for chain in [*json_tgt['proteins'], *json_tgt['rnas']]:
        if tgt_auth_asym_id == chain['auth_asym_id']:
            tgt_seq = chain['entity_poly_seq_one_letter_code']

    if None in [tgt_seq, src_seq]:
        print("""Could not retrieve either of the arguments:
			src_auth_asym, src_seq = [ {}, {} ],
			tgt_auth_asym, tgt_seq = [ {}, {} ]
		 Exiting.""".format(src_auth_asym_id, src_seq, tgt_auth_asym_id, tgt_seq))
        exit(1)

    ixs = [*range(rstart, rend)]
    sm = SeqMatch(src_seq, tgt_seq, ixs)

    target_range = (sm.tgt_ids[0], sm.tgt_ids[-1])
    source_range = (sm.src_ids[0], sm.src_ids[-1])

    src_chain_path = os.path.join(RIBETL_DATA, src_struct.upper(
    ), "CHAINS", "{}_STRAND_{}.cif".format(src_struct.upper(), src_auth_asym_id))
    tgt_chain_path = os.path.join(RIBETL_DATA, tgt_struct.upper(
    ), "CHAINS", "{}_STRAND_{}.cif".format(tgt_struct.upper(), tgt_auth_asym_id))

    print(sm.hl_ixs(sm.src, sm.src_ids))
    print("\n")
    print(sm.hl_ixs(sm.src_aln, sm.aligned_ids))
    print("\n")
    print(sm.hl_ixs(sm.tgt_aln, sm.aligned_ids))
    print("\n")
    print(sm.hl_ixs(sm.tgt, sm.tgt_ids))

    print("Aligning:\n{}\nvs\n{}".format(src_chain_path, tgt_chain_path))
    return (src_chain_path, source_range, tgt_chain_path, target_range)


def pymol_super(
    source_rcsb_id: str,
    source_range: tuple[int, int],
    source_auth_asym_id: str,

    target_rcsb_id: str,
    target_range: tuple[int, int],
    target_auth_asym_id: str,
):

    source_chain_path = os.path.join(RIBETL_DATA,
                                     source_rcsb_id.upper(),
                                     "CHAINS",
                                     "{}_STRAND_{}.cif".format(source_rcsb_id.upper(),
                                                               source_auth_asym_id))

    target_chain_path = os.path.join(RIBETL_DATA,
                                     target_rcsb_id.upper(),
                                     "CHAINS",
                                     "{}_STRAND_{}.cif".format(target_rcsb_id.upper(),
                                                               target_auth_asym_id))

    # Clip chains with pymol, create snippet objects, align those and save.
    print("loading ", source_chain_path)
    cmd.load(source_chain_path)
    cmd.select("resi {}-{}".format(source_range[0], source_range[1]))
    # cmd.select("resi {}-{} and m. {} ".format(source_range[0], source_range[1], source_rcsb_id))
    cmd.create("{}_{}".format(source_rcsb_id, source_auth_asym_id), "sele")
    cmd.delete(source_rcsb_id)

    cmd.load(target_chain_path)
    # cmd.select("resi {}-{} and m. {} ".format(target_range[0], target_range[1], target_rcsb_id))
    cmd.select("resi {}-{}".format(target_range[0], target_range[1]))
    cmd.create("{}_{}".format(target_rcsb_id, target_auth_asym_id), "sele")
    cmd.delete(target_rcsb_id)

    cmd.super("{}_{}".format(source_rcsb_id, source_auth_asym_id),
              "{}_{}".format(target_rcsb_id, target_auth_asym_id),
            #   object="superpose_{}_{}_{}_{}".format(source_rcsb_id, source_auth_asym_id, target_rcsb_id, target_auth_asym_id)
              )

    #  used to be:
    # "/home/rxz/dev/riboxyzbackend/ribetl/static/_TEMP_CHAIN.pdb"
    # cmd.save(os.environ.get("TEMP_CHAIN"))

    return cmd.get_cifstr(
        # selection="superpose_{}_{}_{}_{}".format(source_rcsb_id, source_auth_asym_id, target_rcsb_id, target_auth_asym_id)
        )

# if __name__ == "__main__":

#     prs = argparse.ArgumentParser()

#     prs.add_argument('-s', '--source_struct', type=str, required=True)
#     prs.add_argument('-t', '--target_struct', type=str, required=True)
#     prs.add_argument('-cs', '--chain_source', type=str, required=True)
#     prs.add_argument('-ct', '--chain_target', type=str, required=True)
#     prs.add_argument('-r', '--residue_range', type=str, required=True)

#     args = prs.parse_args()

#     src_struct = args.source_struct.upper()
#     tgt_struct = args.target_struct.upper()
#     chain_source = args.chain_source
#     chain_target = args.chain_target
#     rstart, rend = [* map(int, args.residue_range.split("-"))]

#     print(ranged_super(src_struct, chain_source,
#           tgt_struct, chain_target, (rstart, rend)))
