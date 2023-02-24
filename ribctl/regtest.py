


import re


s5   = '5S ribosomal RNA'
trna = "HUMAN INITIATOR MET-TRNA-I"

def get_rna_nomenclature(text):
    # if re.search(r"\b(5s)", text, flags=re.I | re.M):
    #     return ["5SrRNA"]

    # matches  = re.match("trna", text, flags=re.I | re.M)
    # if re.search(r"(trna)|\b(transfer)\b", text, flags=re.I | re.M):
    #     return ["tRNA"]



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
        "mRNA"    : r"(mrna)|\b(messenger)\b",
        "tRNA"    : r"(trna)|\b(transfer)\b",
    }

    rnatypes = rna_reg.items()

    for i in rnatypes:
        matches = re.search(i[1], text, flags=re.I | re.M)
        if matches!=None:
            return [i[0]]
    return []

print(get_rna_nomenclature(trna))