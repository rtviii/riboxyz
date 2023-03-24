import os
from pprint import pprint
import typing
import requests
from ribctl.lib.types.types_ribosome import ProteinClass
from  ribctl.lib.types.types_polymer import  list_LSU_Proteins,list_SSU_Proteins
import argparse
import subprocess
import sys
import numpy as np
from prody import MSA, Sequence, MSAFile,parseMSA
import prody as prd
import requests
prd.confProDy(verbosity='none')


TAXID_BACTERIA          = 2
TAXID_EUKARYA           = 2759
TAXID_ARCHEA            = 2157

PROTEOVISION_URL        = lambda SU, _class, taxids: "https://ribovision3.chemistry.gatech.edu/showAlignment/{}/{}/{}".format(SU,_class,taxids)
PROTEOVISION_MSA_FOLDER = '/home/rxz/dev/docker_ribxz/api/ribctl/__wip/data/msa_classes_proteovision/'

parser = argparse.ArgumentParser(description='Argument Parser for my ms code.')
parser.add_argument('-p', '--proteovision', required=False, help='Type of proteovision analysis to perform (lsu, ssu, single)')
parser.add_argument('-a','--addtaxid', type=PolymerClass, required=False, help='Additional taxonomic ID to include in the analysis')
parser.add_argument('-t','--lineage', type=str, required=False, help='Comma-separated list of taxonomic levels to include in the query')

args = parser.parse_args()

# Access the argument values
proteovision_type = args.proteovision
additional_taxid  = args.addtaxid
lineage           = args.lineage


def muscle_combine_profiles(msa_path1: str, msa_path2: str, out_filepath: str):
    """Combine two MSA-profiles into a single one. Used here to "append" a target sequence two the ribovision alignment. """
    cmd = ['/home/rxz/dev/docker_ribxz/cli/scripts/muscle3.8', '-profile','-in1', msa_path1, '-in2', msa_path2, '-out', out_filepath]

    subprocess.Popen(cmd,
                      stdout=subprocess.PIPE,
                      stderr=subprocess.PIPE,
                      env=os.environ.copy()
                      ).wait()

    sys.stdout.flush()

def barr2str (bArr):
    return ''.join([ x.decode("utf-8") for x in bArr])

def infer_subunit(protein_class:ProteinClass):
    if protein_class in list_LSU_Proteins:
        return "LSU"
    elif protein_class in list_SSU_Proteins:
        return "SSU"
    else:
        raise ValueError("Unknown protein class: {}".format(protein_class))

def msa_class_proteovision_path(_:ProteinClass):
    return '/home/rxz/dev/docker_ribxz/api/ribctl/__wip/data/msa_classes_proteovision/{}/{}_ribovision.fasta'.format(infer_subunit(_),_)

def msa_add_taxonomic_ids(msa_path:str):
    print("Adding taxonomic IDs to MSA: {}".format(msa_path))
    msa_main = parseMSA(msa_path)
    url = lambda protein_id: f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=protein&id={protein_id}&retmode=json"

    _sequences_new = []
    _labels_new   = []
    if msa_main == None:
        raise FileNotFoundError("File not found: {}".format(msa_path))
    for seq in msa_main:
        original_label  = seq.getLabel()
        nih_protein_id  = original_label.split('|')[1]
        sequence_proper = barr2str(seq.getArray())

        response        = requests.get(url(nih_protein_id))
        if response.status_code == 200:
            protein_summary = response.json()
            # pprint(protein_summary)
            res       = protein_summary["result"]
            try:
                uid       = list(res.keys())[1]
                taxid     = protein_summary["result"][uid]['taxid']
                new_label = f"{original_label}|{taxid}"

                _sequences_new.append(sequence_proper)
                _labels_new.append(new_label)
            except Exception as e:
                print("Failed to identify taxid in {}: {}".format(original_label, e))
                continue

        else:
            print("Omtitting {} Error: ".format(original_label), response.status_code)

    prd.writeMSA(msa_path, MSA(np.array(_sequences_new),' ',labels=_labels_new))

def save_aln(protein:ProteinClass):
        spec_letters, class_digit = protein.split("S" if infer_subunit(protein) == "SSU" else "L")

        # correct for the fact that the proteovision API requires two digits for the single-digit classes (i.e uL1 is uL01)
        if int(class_digit) < 10:
            class_digit = "0{}".format(class_digit)

        def tax_string(spec_letters):
            tax_string = ""
            if 'e' in spec_letters:
                return TAXID_EUKARYA
            if 'b' in spec_letters:
                return TAXID_BACTERIA
            if 'u' in spec_letters:
                return ",".join([str(TAXID_EUKARYA),str(TAXID_ARCHEA),str(TAXID_BACTERIA)])

            return tax_string[:-1]

        URI         = PROTEOVISION_URL(infer_subunit(protein),spec_letters + ( "S" if infer_subunit(protein) == "SSU" else "L" ) + class_digit, tax_string(spec_letters))
        resp        = requests.get(URI)
        start,end   = resp.text.find("<pre>"), resp.text.find("</pre>")
        fasta_lines = resp.text[start+len("<pre>"):end].replace("&gt;",">")
        filename    = "{}_ribovision.fasta".format(protein)
        outpath = os.path.join(PROTEOVISION_MSA_FOLDER,infer_subunit(protein), filename)

        print("{}\t| {} \t| URI: {}".format(protein,resp.status_code, URI) )

        if resp.status_code == 200:
            with open(outpath,"w") as f:
                f.write(fasta_lines)
                print("Wrote {} to {}".format(protein,filename))

                save_aln(nomclass)
                msa_add_taxonomic_ids(msa_class_proteovision_path(nomclass))

def process_proteovision_alignment(nomclass:ProteinClass):
    save_aln(nomclass)
    msa_add_taxonomic_ids(msa_class_proteovision_path(nomclass))

if proteovision_type:
        for nomclass in [*list_LSU_Proteins, *list_SSU_Proteins]:
            try:
                process_proteovision_alignment(nomclass)
            except Exception as e:
                print(e)
                ...


# star -----------------------------
# from ete3 import NcbiTaxa

# # create an instance of the NcbiTaxa class
# ncbi = NcbiTaxa()

# # define the tax ID of interest
# tax_id = 9606  # for human

# # retrieve the lineage
# lineage = ncbi.get_lineage(tax_id)

# # retrieve the scientific names of the lineage
# lineage_names = ncbi.get_taxid_translator(lineage)

# # print the lineage
# for taxid, name in lineage_names.items():
#     print(f"{taxid}\t{name}")