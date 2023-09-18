import os
from pyhmmer import hmmsearch,utils
import pyhmmer
from pyhmmer.plan7 import HMM, Pipeline
from pyhmmer.easel import SequenceFile, Sequence,DigitalSequence, Alphabet, SequenceBlock,DigitalSequenceBlock
from ribctl import ASSETS
from ribctl.ribosome_assets import RibosomeAssets
from ribctl.lib.types.types_poly_nonpoly_ligand import list_ProteinClass
from ribctl.lib.types.types_ribosome import ProteinClass


AMINO        = Alphabet.amino()

# def get_prot_hmm_dict() ->dict[ProteinClass, HMM]: 
#     "Get a dict from protein class to a corrsponding HMM"
#     prot_hmms_dict = {}

#     for hmm_class in list_ProteinClass:

#         class_hmm = os.path.join(ASSETS['hmm_ribosomal_proteins'], f"{hmm_class}.hmm")

#         with pyhmmer.plan7.HMMFile(class_hmm) as hmm_file:
#             hmm                       = hmm_file.read()
#             prot_hmms_dict[hmm_class] = hmm

#     return prot_hmms_dict

# def seq_prot_hmm_evalue(seq:str, hmm:HMM):
#     "Evaluate a sequence against a hmm"
#     seq_      = pyhmmer.easel.TextSequence(name=b"template", sequence=seq)
#     dsb       = DigitalSequenceBlock(AMINO, [seq_.digitize(AMINO)])
#     return pyhmmer.plan7.Pipeline(alphabet=AMINO).search_hmm(hmm,dsb)

# def seq_prot_against_protclasses(seq:str, hmm_dict:dict)->dict[ProteinClass, list[float]]:
#     "Evaluate a sequence against a dict of hmms and return a dict of scores"
#     _ = {}
#     for (prot_class, hmm) in hmm_dict.items():
#         result = seq_prot_hmm_evalue(seq, hmm)
#         _.update({prot_class: [] if len(result) == 0 else list(map(lambda x: x.evalue, result))})
#     return _

# --------------
# def hmmer_process_struct(rcsb_id:str):
#     """Given a ribosome structure (by rcsb_id), """
#     HMM_CLASS_DICT = get_prot_hmm_dict()
#     _, profile = RibosomeAssets(rcsb_id).get_struct_and_profile()

#     if len(profile.polymeric_factors) != 0:
#         for chain in profile.polymeric_factors:
#             _seq     = chain.entity_poly_seq_one_letter_code_can
#             res_dict = seq_prot_against_protclasses(_seq, HMM_CLASS_DICT)

#             nonzero_hits = {}
#             for (class_, hits) in res_dict.items():
#                 if len(hits) > 0:
#                     nonzero_hits.update({class_: hits})

#             if len(list(nonzero_hits.items()))>0:
#                 CRED = '\033[{}m'.format(94);CEND = '\033[0m'
#                 print(CRED + "Found hits for polymeric factor {} : {}".format(prot.entity_poly_strand_id, nonzero_hits)  + CEND)


#     for prot in profile.proteins:
#         _seq     = prot.entity_poly_seq_one_letter_code_can
#         res_dict = seq_prot_against_protclasses(_seq, HMM_CLASS_DICT)

#         nonzero_hits = {}
#         for (class_, hits) in res_dict.items():
#             if len(hits) > 0:
#                 nonzero_hits.update({class_: hits})

#         print("\nProtein {}.{} (old nomenclature {}) | {} :".format(prot.parent_rcsb_id,prot.entity_poly_strand_id, prot.nomenclature, prot.rcsb_pdbx_description))
#         print(nonzero_hits)

#         if len(list(nonzero_hits.items())) == 0:
#             raise LookupError("No hits found for protein {}".format(prot.entity_poly_strand_id))

#         if len(list(nonzero_hits.items()))>1:
#             CRED = '\033[{}m'.format(91);CEND = '\033[0m'
#             print(CRED + "Found hits for protein {} in multiple classes: {}".format(prot.entity_poly_strand_id, nonzero_hits)  + CEND)




for rcsb_id in ['7QV3']:
    hmmer_process_struct(rcsb_id)



#     _, profile = RibosomeAssets(rcsb_id).get_struct_and_profile()
#     for prot in profile.proteins:
#         assigned = ""
#         for prot_class in list_ProteinClass:
#             result = seq_hmm_evalue(prot.entity_poly_seq_one_letter_code_can, prot_class)
#             if result is not None:
#                 print("Identified {} as {} with evalue {}".format(prot.entity_poly_strand_id, prot_class, result))
#                 assigned = prot_class

    # sss   = "MSITKDQIIEAVAAMSVMDVVELISAMEEKFGVSAAAAVAVAAGPVEAAEEKTEFDVILKAAGANKVAVIKAVRGATGLGLKEAKDLVESAPAALKEGVSKDDAEALKKALEEAGAEVEVK"
    # seq   = pyhmmer.easel.TextSequence(name=str("template").encode(), sequence=sss)
    # dsb = DigitalSequenceBlock(AMINO, [seq.digitize(AMINO)])

    # class_hmm = os.path.join(HMM_PROFILES, f"{i}.hmm")

    # with pyhmmer.plan7.HMMFile(class_hmm) as hmm_file:
    #     hmm      = hmm_file.read()
    #     pipeline = pyhmmer.plan7.Pipeline(alphabet=AMINO)
    #     return pipeline.search_hmm(hmm,dsb)
    # # compare_seq_to_hmm(bl12path, class_hmm)




# seq          = SequenceFile("/home/rxz/dev/docker_ribxz/api/3J7Z_bL12.fasta", digital=True)

# sequence_file = SequenceFile(sequence_file_path)
# sequence = sequence_file.read_one_sequence()

# for hits in search:
#     print(f"HMM {hits.E} found {len(hits)} hits in the target sequences")

# alphabet      = pyhmmer.easel.Alphabet.amino()
# builder       = pyhmmer.plan7.Builder(alphabet)
# background    = pyhmmer.plan7.Background(alphabet)
# anonymous_msa = pyhmmer.easel.TextMSA(sequences=[])

# hmm, _, _ = builder.build_msa(msa, background)


# pipeline = pyhmmer.plan7.Pipeline(alphabet, background=background)
# with pyhmmer.easel.SequenceFile("data/seqs/LuxC.faa", digital=True, alphabet=alphabet) as seq_file:
#     hits = pipeline.search_hmm(hmm, seq_file)