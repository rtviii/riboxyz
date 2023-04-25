from io import StringIO
import io
import os
from pprint import pprint
import typing
import prody
import requests
from api.ribctl.etl.ribosome_assets import RibosomeAssets
from api.ribctl.lib.utils import open_structure
from ribctl.lib.types.types_ribosome import RNA, PolymericFactor, Protein, ProteinClass
from Bio import AlignIO
from  api.ribctl.lib.types.types_poly_nonpoly_ligand import  list_LSU_Proteins,list_SSU_Proteins
import argparse
import subprocess
import sys
import numpy as np
import gemmi
from Bio import pairwise2

from Bio.Seq import Seq

from Bio import SeqIO
from Bio import SeqRecord
from prody import MSA, Sequence, MSAFile, buildMSA,parseMSA
import prody as prd
import requests



prd.confProDy(verbosity='none')
RIBETL_DATA = os.environ.get('RIBETL_DATA')

TAXID_BACTERIA          = 2
TAXID_EUKARYA           = 2759
TAXID_ARCHEA            = 2157

PROTEOVISION_URL        = lambda SU, _class, taxids: "https://ribovision3.chemistry.gatech.edu/showAlignment/{}/{}/{}".format(SU,_class,taxids)
PROTEOVISION_MSA_FOLDER = '/home/rxz/dev/docker_ribxz/api/ribctl/__wip/data/msa_classes_proteovision/'

parser = argparse.ArgumentParser(description='Argument Parser for my ms code.')
parser.add_argument('-p', '--proteovision', required=False, help='Type of proteovision analysis to perform (lsu, ssu, single)')
parser.add_argument('-t','--lineage', type=str, required=False, help='Comma-separated list of taxonomic levels to include in the query')

args = parser.parse_args()

# Access the argument values
proteovision_type = args.proteovision
lineage           = args.lineage

#! util
def util__backwards_match(alntgt: str, aln_resid: int, verbose: bool = False) -> typing.Tuple[int, str, int]:
    """
    returns (projected i.e. "ungapped" residue id, the residue itself residue)
    """
    if aln_resid > len(alntgt):
        raise IndexError(f"Passed residue with invalid index ({aln_resid}) to backwards-match to target. Seqlen:{len(alntgt)}")

    counter_proper = 0
    for i, char in enumerate(alntgt):
        if i == aln_resid:
            if verbose:
                print("[ {} ] <-----> id.[aligned: {} | orgiginal: {} ]".format(
                    alntgt[aln_resid], i, counter_proper))
            return (counter_proper,  alntgt[i], aln_resid)
        if char == '-':
            continue
        else:
            counter_proper += 1

    raise LookupError("Could not find residue in aligned sequence.")

#! util
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

    raise LookupError("Could not find residue in aligned sequence.")

#! util
def barr2str (bArr):
    return ''.join([ x.decode("utf-8") for x in bArr])


#! util
def seq_to_fasta(rcsb_id: str, _seq: str, outfile: str):
    from Bio.Seq import Seq
    _seq          = _seq.replace("\n", "")
    seq_record    = SeqRecord.SeqRecord(Seq(_seq).upper())
    seq_record.id = seq_record.description = rcsb_id
    SeqIO.write(seq_record, outfile, 'fasta')
    

show = True

def msa_profiles_dict_prd()->dict[str,MSA ]:
    MSA_PROFILES_PATH  = '/home/rxz/dev/docker_ribxz/api/ribctl/assets/msa_profiles/'
    msa_dict  = {}

    # LSU
    LSU_path = os.path.join(MSA_PROFILES_PATH,"LSU")
    for msafile in os.listdir(LSU_path):
        classname = msafile.split("_")[0]
        class_msa       = prody.parseMSA(os.path.join(LSU_path, msafile))
        msa_dict = {f"{classname}": class_msa, **msa_dict}

    #SSU
    SSU_path = os.path.join(MSA_PROFILES_PATH,"SSU")
    for msafile in os.listdir(SSU_path):
        classname = msafile.split("_")[0]
        class_msa       = prody.parseMSA(os.path.join(SSU_path, msafile))
        msa_dict = {f"{classname}": class_msa, **msa_dict}
    return msa_dict
def msa_profiles_dict()->dict[str, AlignIO.MultipleSeqAlignment]:
    MSA_PROFILES_PATH  = '/home/rxz/dev/docker_ribxz/api/ribctl/assets/msa_profiles/'
    msa_dict  = {}

    # LSU
    LSU_path = os.path.join(MSA_PROFILES_PATH,"LSU")
    for msafile in os.listdir(LSU_path):
        classname = msafile.split("_")[0]
        class_msa       = AlignIO.read(os.path.join(LSU_path, msafile), "fasta")
        msa_dict = {f"{classname}": class_msa, **msa_dict}

    #SSU
    SSU_path = os.path.join(MSA_PROFILES_PATH,"SSU")
    for msafile in os.listdir(SSU_path):
        classname = msafile.split("_")[0]
        class_msa       = AlignIO.read(os.path.join(SSU_path, msafile), "fasta")
        msa_dict = {f"{classname}": class_msa, **msa_dict}
    return msa_dict

def get_fasta_string(chain:RNA|Protein| PolymericFactor)->str:
    fasta_description      = "[{}.{}] {} |{}| {}".format(chain.parent_rcsb_id,chain.auth_asym_id, chain.src_organism_names[0], "",  chain.src_organism_ids[0])
    _seq                   = chain.entity_poly_seq_one_letter_code_can.replace("\n","")
    seq_record             = SeqRecord.SeqRecord(Seq(_seq).upper())
    seq_record.id          = fasta_description
    seq_record.description = ""
    return seq_record.format('fasta')

def prot_class_msa_extend(rcsb_id:str, poly_class:ProteinClass)->AlignIO.MultipleSeqAlignment:
    R                 = RibosomeAssets(rcsb_id)
    chain             = R.get_chain_by_polymer_class(poly_class)
    if chain is None:
        raise LookupError("Could not find chain in {} for protein class: {}".format(rcsb_id,poly_class))

    fasta_target  = get_fasta_string(chain)
    class_profile = msa_class_proteovision_path(poly_class)

    cmd = [
        '/home/rxz/dev/docker_ribxz/api/ribctl/muscle3.8',
        '-profile',
        '-in1',
        class_profile,
        '-in2',
        '-',
        '-quiet']

    process = subprocess.Popen(cmd,
                               stdout=subprocess.PIPE,
                               stdin=subprocess.PIPE,
                               stderr=subprocess.PIPE, env=os.environ.copy())

    stdout, stderr = process.communicate(input=fasta_target.encode())
    out   ,err     = stdout.decode(), stderr.decode()

    process.wait()

    msa_file = StringIO(out)
    msa      = AlignIO.read(msa_file, "fasta")

    return msa


def msa_class_proteovision_path(prot_class:ProteinClass):
    def infer_subunit(protein_class:ProteinClass):
        if protein_class in list_LSU_Proteins:
            return "LSU"
        elif protein_class in list_SSU_Proteins:
            return "SSU"
        else:
            raise ValueError("Unknown protein class: {}".format(protein_class))
    path = os.path.join(
        '/home/rxz/dev/docker_ribxz/api',
        'ribctl',
        'assets',
        'msa_profiles/{}/{}_ribovision.fasta'.format(infer_subunit(prot_class),prot_class)
        )
    assert os.path.exists(path), "File not found: {}".format(path)
    return path

"""download + tax-identify a protein class sequence from proteovision"""
def process_proteovision_alignment(nomclass:ProteinClass):

    def msa_add_taxonomic_ids(msa_path:str):

        msa_main = parseMSA(msa_path)
        url      = lambda protein_id: f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=protein&id={protein_id}&retmode=json"

        _sequences_new = []
        _labels_new    = []

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

            # correct for the fact that the proteovision API requires two digits for the single-digit classes (i.e uL1 is uL01)
            spec_letters, class_digit = protein.split("S" if infer_subunit(protein) == "SSU" else "L")

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

    save_aln(nomclass)
    msa_add_taxonomic_ids(msa_class_proteovision_path(nomclass))

def display_msa_class(nomclass:ProteinClass):
    nomclass = 'uL23'
    msa_main = parseMSA(msa_class_proteovision_path(nomclass))

    for seq in msa_main:
        print(seq)

