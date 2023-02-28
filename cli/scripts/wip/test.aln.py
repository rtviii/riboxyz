from Bio import pairwise2
from Bio.pairwise2 import format_alignment


alignment = pairwise2.align.globalxx("AABBCC","AABBCC")
# print(format_alignment(*alignment[0]))


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
SEQ1_1    = "AAABrfBCCC"
SEQ1_2    = "AAABrfBCCC"
IDX       = 3 #expected 3, B
alignment = pairwise2.align.globalxx(SEQ1_1,SEQ1_2)
src_aln   = alignment[0].seqA
tgt_aln   = alignment[0].seqB

src_aln_resid = util__forwards_match(src_aln, IDX)
tgt_aln_resid = util__backwards_match(tgt_aln, src_aln_resid)
assert(tgt_aln_resid == IDX)


# ? Case SRC extended forward by 2 (after id)
SEQ2_1    = "AAABrfBCCCii"
SEQ2_2    = "AAABrfBCCC"
IDX       = 3 #expected 3, B
alignment = pairwise2.align.globalxx(SEQ2_1,SEQ2_2)
src_aln   = alignment[0].seqA
tgt_aln   = alignment[0].seqB

src_aln_resid = util__forwards_match(src_aln, IDX)
tgt_aln_resid = util__backwards_match(tgt_aln, src_aln_resid)
assert(tgt_aln_resid == IDX)

# ? Case SRC extended forward by 2 (before id)

SEQ3_1    = "AAABrfBCCC"
SEQ3_2    = "AiiAABrfBCCC"
IDX       = 3 #expected 3, B
alignment = pairwise2.align.globalxx(SEQ3_1,SEQ3_2)
src_aln   = alignment[0].seqA
tgt_aln   = alignment[0].seqB

src_aln_resid = util__forwards_match(src_aln, IDX)
tgt_aln_resid = util__backwards_match(tgt_aln, src_aln_resid)
assert(tgt_aln_resid == IDX+2)

# ? Case TGT forward by 2 (after id)

SEQ4_1    = "AAABrfBCCC"
SEQ4_2    = "AAAiiBrfBCCC"
IDX       = 3 #expected 3, B
alignment = pairwise2.align.globalxx(SEQ4_1,SEQ4_2)
src_aln   = alignment[0].seqA
tgt_aln   = alignment[0].seqB

src_aln_resid = util__forwards_match(src_aln, IDX)
tgt_aln_resid = util__backwards_match(tgt_aln, src_aln_resid)
assert(tgt_aln_resid == IDX)
# ? Case TGT extended forward by 2 (before id)

SEQ5_1    = "AAABrfBCCC"
SEQ5_2    = "dAfAABrfBCCC"
IDX       = 3 #expected 3, B
alignment = pairwise2.align.globalxx(SEQ5_1,SEQ5_2)
src_aln   = alignment[0].seqA
tgt_aln   = alignment[0].seqB

src_aln_resid = util__forwards_match(src_aln, IDX)
tgt_aln_resid = util__backwards_match(tgt_aln, src_aln_resid)

print("original")
print(f"{SEQ5_1} \n{src_aln}    | [{SEQ5_1[IDX]}]")

print(f"computed :{tgt_aln_resid}")
print(f"{SEQ5_2} \n{tgt_aln}    | [{SEQ5_2[tgt_aln_resid]}]")