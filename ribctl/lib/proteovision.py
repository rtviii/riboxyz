

# def PROTEOVISION_URL(SU, _class, taxids): return "https://ribovision3.chemistry.gatech.edu/showAlignment/{}/{}/{}".format(SU, _class, taxids)
# PROTEOVISION_MSA_FOLDER = '/home/rxz/dev/docker_ribxz/api/ribctl/__wip/data/msa_classes_proteovision/'
# def msa_class_proteovision_path(prot_class: ProteinClass):
#     def infer_subunit(protein_class: ProteinClass):
#         if protein_class in list_LSUProteinClass:
#             return "LSU"
#         elif protein_class in list_SSUProteinClass:
#             return "SSU"
#         else:
#             raise ValueError("Unknown protein class: {}".format(protein_class))
#     path = os.path.join(
#         '/home/rxz/dev/docker_ribxz/api',
#         'ribctl',
#         'assets',
#         'msa_profiles/{}/{}_ribovision.fasta'.format(
#             infer_subunit(prot_class), prot_class)
#     )
#     assert os.path.exists(path), "File not found: {}".format(path)
#     return path


# """download + tax-identify a protein class sequence from proteovision"""


# # RMPRD
# def process_proteovision_alignment(nomclass: ProteinClass):

#     def msa_add_taxonomic_ids(msa_path: str):

#         msa_main = parseMSA(msa_path)

#         def url(
#             protein_id): return f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=protein&id={protein_id}&retmode=json"

#         _sequences_new = []
#         _labels_new = []

#         if msa_main == None:
#             raise FileNotFoundError("File not found: {}".format(msa_path))

#         for seq in msa_main:
#             original_label = seq.getLabel()
#             nih_protein_id = original_label.split('|')[1]
#             sequence_proper = barr2str(seq.getArray())

#             response = requests.get(url(nih_protein_id))
#             if response.status_code == 200:
#                 protein_summary = response.json()
#                 # pprint(protein_summary)
#                 res = protein_summary["result"]
#                 try:
#                     uid = list(res.keys())[1]
#                     taxid = protein_summary["result"][uid]['taxid']
#                     new_label = f"{original_label}|{taxid}"

#                     _sequences_new.append(sequence_proper)
#                     _labels_new.append(new_label)
#                 except Exception as e:
#                     print("Failed to identify taxid in {}: {}".format(
#                         original_label, e))
#                     continue

#             else:
#                 print("Omtitting {} Error: ".format(
#                     original_label), response.status_code)

#         prd.writeMSA(msa_path, MSA(
#             np.array(_sequences_new), ' ', labels=_labels_new))

#     def save_aln(protein: ProteinClass):

#         # correct for the fact that the proteovision API requires two digits for the single-digit classes (i.e uL1 is uL01)
#         spec_letters, class_digit = protein.split(
#             "S" if infer_subunit(protein) == "SSU" else "L")

#         if int(class_digit) < 10:
#             class_digit = "0{}".format(class_digit)

#         def tax_string(spec_letters):
#             tax_string = ""
#             if 'e' in spec_letters:
#                 return TAXID_EUKARYA
#             if 'b' in spec_letters:
#                 return TAXID_BACTERIA
#             if 'u' in spec_letters:
#                 return ",".join([str(TAXID_EUKARYA), str(TAXID_ARCHEA), str(TAXID_BACTERIA)])

#             return tax_string[:-1]

#         URI = PROTEOVISION_URL(infer_subunit(protein), spec_letters + ("S" if infer_subunit(
#             protein) == "SSU" else "L") + class_digit, tax_string(spec_letters))
#         resp = requests.get(URI)
#         start, end = resp.text.find("<pre>"), resp.text.find("</pre>")
#         fasta_lines = resp.text[start+len("<pre>"):end].replace("&gt;", ">")
#         filename = "{}_ribovision.fasta".format(protein)
#         outpath = os.path.join(PROTEOVISION_MSA_FOLDER,
#                                infer_subunit(protein), filename)

#         print("{}\t| {} \t| URI: {}".format(protein, resp.status_code, URI))

#         if resp.status_code == 200:
#             with open(outpath, "w") as f:
#                 f.write(fasta_lines)

#     save_aln(nomclass)
#     msa_add_taxonomic_ids(msa_class_proteovision_path(nomclass))


# # RMPRD
# def display_msa_class(nomclass: ProteinClass):
#     msa_main = parseMSA(msa_class_proteovision_path(nomclass))
#     for seq in msa_main:
#         print(seq)