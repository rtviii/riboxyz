import os
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
parser.add_argument('-p', '--proteovision', choices=['lsu', 'ssu', 'single'], required=False,
                    help='Type of proteovision analysis to perform (lsu, ssu, single)')
parser.add_argument('-a','--addtaxid', type=int, required=False,
                    help='Additional taxonomic ID to include in the analysis')
parser.add_argument('-t','--lineage', type=str, required=False,
                    help='Comma-separated list of taxonomic levels to include in the analysis')

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
    
infer_subunit = lambda polymer_class : 'LSU' if 'L' in polymer_class else 'SSU'
datapath      = lambda subunit, polymer_class : '/home/rxz/dev/docker_ribxz/api/ribctl/__wip/data/msa_classes_proteovision/{}/{}_ribovision.fasta'.format(subunit, polymer_class)


nomclass      = 'uL23'
path          = datapath(infer_subunit(nomclass),nomclass)



def msa_add_taxonomic_ids(msa_path:str):
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
            taxid           = protein_summary["result"][nih_protein_id]["taxid"]
            new_label = f"{original_label}|{taxid}"

            _sequences_new.append(sequence_proper)
            _labels_new.append(new_label)
        else:
            print("Omtitting {} Error: ".format(original_label), response.status_code)

    prd.writeMSA(msa_path, MSA(np.array(_sequences_new),nomclass,labels=_labels_new))


def save_aln(protein:ProteinClass, subunit: typing.Literal["SSU","LSU"]):
        spec_letters, class_digit = protein.split("S" if subunit == "SSU" else "L")

        # correct for the fact that the proteovision API requires two digits for the single-digit classes (i.e uL1 is uL01)
        if int(class_digit) < 10:
            class_digit = "0{}".format(class_digit)

        def tax_string(spec_letters):
            tax_string = ""
            if 'e' in spec_letters:
                return E
            if 'b' in spec_letters:
                return B
            if 'u' in spec_letters:
                return ",".join([B,E,A])

            return tax_string[:-1]

        URI         = PROTEOVISION_URL(subunit,spec_letters + ( "S" if subunit == "SSU" else "L" ) + class_digit, tax_string(spec_letters))
        resp        = requests.get(URI)
        start,end   = resp.text.find("<pre>"), resp.text.find("</pre>")
        fasta_lines = resp.text[start+len("<pre>"):end].replace("&gt;",">")
        filename    = "{}_ribovision.fasta".format(protein)
        outpath = os.path.join(PROTEOVISION_MSA_FOLDER,subunit, filename)

        print("{}\t| {} \t| URI: {}".format(protein,resp.status_code, URI) )

        if resp.status_code == 200:
            with open(outpath,"w") as f:
                f.write(fasta_lines)
                print("Wrote {} to {}".format(protein,filename))

def save_alns(subunit_protein_names: list[str], subunit: typing.Literal["SSU","LSU"]):
    protnames = subunit_protein_names
    for protein in protnames:

        print("Fetching {}...".format(protein), end="")
        spec_letters, class_digit = protein.split("S" if subunit == "SSU" else "L")

        # correct for the fact that the proteovision API requires two digits for the single-digit classes (i.e uL1 is uL01)
        if int(class_digit) < 10:
            class_digit = "0{}".format(class_digit)

        def tax_string(spec_letters):
            tax_string = ""
            if 'e' in spec_letters:
                return E
            if 'b' in spec_letters:
                return B
            if 'u' in spec_letters:
                return ",".join([B,E,A])

            return tax_string[:-1]

        URI         = PROTEOVISION_URL(subunit,spec_letters + ( "S" if subunit == "SSU" else "L" ) + class_digit, tax_string(spec_letters))
        resp        = requests.get(URI)
        start,end   = resp.text.find("<pre>"), resp.text.find("</pre>")
        fasta_lines = resp.text[start+len("<pre>"):end].replace("&gt;",">")
        filename    = "{}_ribovision.fasta".format(protein)
        outpath = os.path.join(PROTEOVISION_MSA_FOLDER,subunit, filename)

        print("{}\t| {} \t| URI: {}".format(protein,resp.status_code, URI) )

        if resp.status_code == 200:
            with open(outpath,"w") as f:
                f.write(fasta_lines)
                print("Wrote {} to {}".format(protein,filename))
         
save_aln('uL23', "LSU")