from Bio import pairwise2
from Bio.pairwise2 import format_alignment


alignment = pairwise2.align.globalxx("AABBCC","AABBCC")
print(format_alignment(*alignment[0]))


# !!! ※ TEST FORWRADS AND BACKWARDS MATCH FUNCTIONS
def util__backwards_match(alntgt: str, resid: int):
    """Returns the target-sequence index of a residue in the (aligned) target sequence"""
    if resid > len(alntgt):
        raise IndexError(f"Passed residue with invalid index ({resid}) to back-match to target.Seqlen:{len(alntgt)}")

    counter_proper = 0
    for i, char in enumerate(alntgt):
        if i == resid:
            return counter_proper
        if char == '-':
            continue
        else:
            counter_proper += 1

def util__forwards_match(string: str, resid: int):
    """Returns the index of a source-sequence residue in the (aligned) source sequence."""
    if resid >= len(string):
        raise IndexError("Requested residue index({resid}) exceeds aligned(likely already gaps-extended) sequence. Something went wrong.")

    count_proper = 0
    for alignment_indx, char in enumerate(string):
        if count_proper == resid:
            return alignment_indx
        if char == '-':
            continue
        else:
            count_proper += 1


# ? Case Same-id-to-same (unperturbed seqs)

# ? Case SRC extended forward by 2 (after id)
# ? Case SRC extended forward by 2 (before id)

# ? Case TGT forward by 2 (after id)
# ? Case TGT extended forward by 2 (before id)