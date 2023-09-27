import re
from Bio import SeqIO
from ete3 import NCBITaxa
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq



from ribctl.etl.ribosome_assets import RibosomeAssets
structs = RibosomeAssets.list_all_structs()



def rna_classify(poly_pdbx_description:str|None):
    rna_reg = {
        "5SrRNA"  : r"\b(5s)",
        "5.8SrRNA": r"\b(5\.8s)",
        "12SrRNA" : r"\b(12s)",
        "16SrRNA" : r"\b(16s)",
        "21SrRNA" : r"\b(21s)",
        "23SrRNA" : r"\b(23s)",
        "25SrRNA" : r"\b(25s)",
        "28SrRNA" : r"\b(28s)",
        "35SrRNA" : r"\b(35s)",
    }

    rnatypes = rna_reg.items()
    for i in rnatypes:
        matches = re.search(i[1], poly_pdbx_description if poly_pdbx_description !=None else '', flags=re.IGNORECASE | re.MULTILINE)
        if matches != None:
            return [i[0]]
    return []

to_keep =[]
stats ={}
for struct in structs:
    print(struct)
    try:
        profile   = RibosomeAssets(struct).profile()
        org_taxid = profile.src_organism_ids[0]
        if profile.rnas !=None:
            for rna in profile.rnas:
                if "5.8SrRNA" in rna.nomenclature and rna_classify(rna.rcsb_pdbx_description) == ["5.8SrRNA"]:
                    sequence    = rna.entity_poly_seq_one_letter_code_can
                    record_id   = str( org_taxid )
                    description = "{}.{}|{}".format(profile.rcsb_id, rna.auth_asym_id, rna.rcsb_pdbx_description)

                    sequence_obj = Seq(sequence)
                    record = SeqRecord(sequence_obj, id=record_id, description=description)
                    to_keep.append(record)

            
    except Exception as e:
        print(e)

dest = "5.8SrRNA.fasta"

with open(dest, "w") as output_handle:
    SeqIO.write(to_keep, output_handle, "fasta")





exit()

# ncbi    = NCBITaxa()
# dest    = "./m_12SrRNA.fasta"
# to_keep = []
# i       = 0
# specs   = {}

# with open(fasta_file, "r") as handle:
#     for record in SeqIO.parse(handle, "fasta"):
#         if "12S" in record.description and ( len(record.seq) < 1200 and len(record.seq) > 600 ):
#             try:
#                 taxname_huh = " ".join(record.description.split(" ")[1:3])
#                 taxid       = str(ncbi.get_name_translator([taxname_huh])[taxname_huh][0])
#                 record.id = str(taxid)

#                 if taxid not in specs:
#                     specs[taxid] = 1
#                 else:
#                     specs[taxid] += 1


#                 record.description = "rnacentral_id:{}|{}".format(record.id, record.description)
#                 if specs[taxid] < 5:
#                     to_keep.append(record)
#             except:
#                 ...

# with open(dest, "w") as output_handle:
#     SeqIO.write(to_keep, output_handle, "fasta")