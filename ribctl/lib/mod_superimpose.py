import argparse
import itertools
import os
from typing import Tuple
from ribctl.lib.types.types_ribosome import RibosomeStructure
from pymol import cmd
from ribctl.lib.mod_transpose_bsites import SeqMatch
from ribctl.lib.utils import open_structure

RIBETL_DATA = os.environ.get("RIBETL_DATA")

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
            tgt_chain_path, 
            target_range)


def ranged_align_by_auth_asym_id(
        src_struct      : str,
        src_auth_asym_id: str,

        tgt_struct      : str,
        tgt_auth_asym_id: str,

        rng             : Tuple[int, int],
) -> Tuple[str,str, Tuple[int, int], str,str, Tuple[int, int]]:
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


    
    for src_chain in itertools.chain(json_src.proteins, json_src.rnas if json_src.rnas else []):
        if src_auth_asym_id == src_chain.auth_asym_id:
            seq_ids["src_seq"] = src_chain.entity_poly_seq_one_letter_code

    for chain in itertools.chain(json_tgt.proteins, json_tgt.rnas if json_tgt.rnas else []):
        if tgt_auth_asym_id == chain.auth_asym_id:
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
    return (seq_ids["src_auth_asym_id"],
            src_chain_path,
            source_range,

            seq_ids["tgt_auth_asym_id"],
            tgt_chain_path, 
            target_range)


def pymol_super(

    source_rcsb_id     : str,
    source_range       : tuple[int, int],
    source_auth_asym_id: str,

    target_rcsb_id     : str,
    target_range       : tuple[int, int],
    target_auth_asym_id: str,

):

    cmd.delete("all")
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


    # names
    clipped_src = "src_clipped_{}_{}".format(source_rcsb_id, source_auth_asym_id)
    clipped_tgt = "tgt_clipped_{}_{}".format(target_rcsb_id, target_auth_asym_id)



    cmd.load(source_chain_path, "strand1")
    cmd.select("resi {}-{} and m. strand1".format(source_range[0], source_range[1]))
    cmd.create(clipped_src, "sele")
    cmd.delete(source_rcsb_id)



    cmd.load(target_chain_path, "strand2")
    cmd.select("resi {}-{} and m. strand2".format(target_range[0], target_range[1]))
    cmd.create(clipped_tgt, "sele")
    cmd.delete(target_rcsb_id)

    superimpose_name = "superpose_{}_{}_{}_{}".format(source_rcsb_id, source_auth_asym_id, target_rcsb_id, target_auth_asym_id)
    cmd.super(clipped_src,clipped_tgt, object=superimpose_name)
    cmd.save("together.cif", selection=superimpose_name)
   


    #  used to be:
    # "/home/rxz/dev/riboxyzbackend/ribetl/static/_TEMP_CHAIN.pdb"
    # cmd.save(os.environ.get("TEMP_CHAIN"))
    x = cmd.get_cifstr(selection=superimpose_name)
    print(x)
    return  x
