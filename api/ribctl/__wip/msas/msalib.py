import os
from pprint import pprint
import typing
import requests
from api.ribctl.__wip.tunnel.ptc_mass import get_sequence_by_nomclass
from ribctl.lib.types.types_ribosome import ProteinClass
from  api.ribctl.lib.types.types_poly_nonpoly_ligand import  list_LSU_Proteins,list_SSU_Proteins
import argparse
import subprocess
import sys
import numpy as np
import gemmi
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
parser.add_argument('-t','--lineage', type=str, required=False, help='Comma-separated list of taxonomic levels to include in the query')

args = parser.parse_args()

# Access the argument values
proteovision_type = args.proteovision
lineage           = args.lineage

def muscle_combine_profile(msa_path1: str, msa_path2: str, out_filepath: str):
    """Combine two MSA-profiles into a single one. Used here to "append" a target sequence two the ribovision alignment. """
    cmd = ['/home/rxz/dev/docker_ribxz/cli/scripts/muscle3.8', '-profile','-in1', msa_path1, '-in2', msa_path2, '-out', out_filepath]
    subprocess.Popen(cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE, env=os.environ.copy()).wait()
    sys.stdout.flush()

def add_target_to_domain_alignment(rcsb_id: str, domain: str):
    """
    @param rcsb_id: PDB ID of target structure
    @param domain:  Domain of target structure (euk or bac)
    """

    rcsb_id = argdict["target"]
    struct_profile = open_structure(rcsb_id, 'json')
    [chain_id, strand_target] = get_23SrRNA_strandseq(
        rcsb_id,
        custom_path=os.path.join(
            RIBETL_DATA, rcsb_id.upper(), f"{rcsb_id.upper()}_modified.cif")
    )

    fpath_23s = f'{rcsb_id.upper()}_{chain_id}_23SrRNA.fasta'
    domain_alignment: str = ''
    if domain == 'bacteria':
        domain_alignment = 'data/ribovision.bacteria.fasta'
    elif domain == 'eukarya':
        domain_alignment = 'data/ribovision.eukaryota.fasta'
    else:
        raise FileNotFoundError(
            "Domain misspecified. Must be either 'bacteria' or 'eukarya'.")

    seq_to_fasta(rcsb_id, strand_target, fpath_23s)
    muscle_combine_profile(domain_alignment, fpath_23s,f'combined_{rcsb_id.upper()}_ribovision_{domain}.fasta')
def seq_to_fasta(rcsb_id: str, _seq: str, outfile: str):
    from Bio.Seq import Seq
    _seq          = _seq.replace("\n", "")
    seq_record    = SeqRecord.SeqRecord(Seq(_seq).upper())
    seq_record.id = seq_record.description = rcsb_id
    SeqIO.write(seq_record, outfile, 'fasta',)

def util__backwards_match(alntgt: str, aln_resid: int, verbose: bool = False) -> Tuple[int, str, int]:
    """
    returns (projected i.e. "ungapped" residue id, the residue itself residue)
    """
    if aln_resid > len(alntgt):
        raise IndexError(
            f"Passed residue with invalid index ({aln_resid}) to back-match to target. Seqlen:{len(alntgt)}")

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

    raise LookupError()

def util__forwards_match(string: str, resid: int):
    """Returns the index of a source-sequence residue in the (aligned) source sequence."""
    if resid >= len(string):
        raise IndexError(
            "Requested residue index({resid}) exceeds aligned(likely already gaps-extended) sequence. Something went wrong.")

    count_proper = 0
    for alignment_indx, char in enumerate(string):
        if count_proper == resid:
            return alignment_indx
        if char == '-':
            continue
        else:
            count_proper += 1

def get_sequence_by_nomclass(rcsb_id: str, nomenclature_class: str, canonical:bool=True,path:str=None) -> Tuple[str, str]:

    target       = gemmi.cif.read_file(path)
    block        = target.sole_block()
    model        = gemmi.read_structure(path)[0]

    STRAND = None
    SEQ    = None

    # Locate the chain of given nom. class
    for (strand, nomclass) in zip(
        block.find_loop('_ribosome_nomenclature.entity_poly.pdbx_strand_id'),
        block.find_loop('_ribosome_nomenclature.polymer_class')
    ):
        if nomclass == nomenclature_class:
            STRAND = strand
            break

    # Now find sequence of this class
    for (chain_id, one_letter_code) in zip(
        block.find_loop('_entity_poly.pdbx_strand_id'),
        block.find_loop('_entity_poly.pdbx_seq_one_letter_code_can') if canonical else block.find_loop('_entity_poly.pdbx_seq_one_letter_code')
    ):
        # X-RAY structures have 'dual' chains. Split on comma to check both.
        if STRAND in chain_id.split(','):
            SEQ = str(one_letter_code).strip(";").strip("\n")

    if SEQ == None:
        print("Could not locate {} sequence in {} CIF file".format(
            nomenclature_class, rcsb_id))
    return (STRAND, SEQ)

def retrieve_LSU_rRNA(rcsb_id, canonical:bool=True):
    # TODO: move method to RibosomeStructure please
    annotated_cifpath = os.path.join(RIBETL_DATA, rcsb_id.upper(), f"{rcsb_id.upper()}_modified.cif")
    rna_type          = ""
    #--------------
    [chain_id, strand_target] = get_sequence_by_nomclass(
        rcsb_id,
        "23SrRNA",
        canonical,
        path=os.path.join(
            RIBETL_DATA, rcsb_id.upper(), f"{rcsb_id.upper()}_modified.cif")
    )
    rna_type = "23SrRNA"

    if chain_id == None or strand_target == None:
        [chain_id, strand_target] = get_sequence_by_nomclass(
            rcsb_id,
            "25SrRNA",
            canonical,
            path=annotated_cifpath
        )
        rna_type = "25SrRNA"

    if chain_id == None or strand_target == None:
        [chain_id, strand_target] = get_sequence_by_nomclass(
            rcsb_id,
            "28SrRNA",
            canonical,
            path=annotated_cifpath)
        rna_type = "28SrRNA"

    if chain_id == None or strand_target == None:
        print("Failed to locate either 23S or 25S or 28 rRNA in {}".format(rcsb_id))
        exit(1)

    return [chain_id, strand_target, rna_type]

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

show = True

if show:

    nomclass = 'uL23'
    msa_main = parseMSA(msa_class_proteovision_path(nomclass))

    for seq in msa_main:
        print(seq)

if proteovision_type is not None: process_proteovision_alignment(proteovision_type)

if lineage is not None:
    
    # ※ The goal is to minimize evolutionary distance between the set of sequences that connect the polymer being considered to the the MSA *
    """given a substrand in a polymer class (say, uL23):
       - get taxid of its parent structure 
       - search the corresponding proteovision alignment for sequences with the following tax id priority:

               {the closest level in the lineage (ideally the aln has the same exact tax id) }--|
               {if not, the next closest populated level up in the lineage and its paralels  }  |->  "bridge" sequence

       - do the regular substrand lookup and backtracking into the 
            """


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