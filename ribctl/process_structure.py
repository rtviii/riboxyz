import functools
import itertools
from pprint import pprint
import requests
from urllib.parse import urlencode
from gql_querystrings import monolithic
from ribosome_types import Protein, RNA, Ligand, RibosomeStructure
import re
import json

# open both subunit maps:

LSU_map = {k.upper(): v for k, v in json.load(
    open('subunit_map_LSU.json', 'r')).items()}
SSU_map = {k.upper(): v for k, v in json.load(
    open('subunit_map_SSU.json', 'r')).items()}


def gql_monolith(rcsb_id): return monolithic.replace(
    "$RCSB_ID", rcsb_id.upper())
# gql_structs             = lambda rcsb_id: structure_string.replace("$RCSB_ID", rcsb_id.upper())
# gql_polymer_entities    = lambda rcsb_id: polymer_entities_string.replace("$RCSB_ID", rcsb_id.upper())
# gql_nonpolymer_entities = lambda rcsb_id: nonpolymer_entities_string.replace("$RCSB_ID", rcsb_id.upper())


def get_protein_nomenclature(protein):
    banregex = r"/\b([ueb][ls]\d{1,2})\b/gi"
    #     check authors's annotations. if classes are present --> use that.
    finds = re.search(
        banregex, protein['rcsb_polymer_entity']['pdbx_description'])
    print(finds)

    if (finds != None):
        firstcap: str = finds[0]
        clasname = firstcap[0].lower() + firstcap[1].upper() + firstcap[2:]
        return [clasname]

    elif protein['pfams'] == None or protein['pfams'] == []:
        return []
    else:
        pfamids = [pfam['rcsb_pfam_accession'] for pfam in protein['pfams']]

        nomenclature = []
        for pfam_id in pfamids:
            [nomenclature.append(kv[0]) if pfam_id in kv[1]
             ['pfamDomainAccession'] else ... for kv in LSU_map.items()]
            [nomenclature.append(kv[0]) if pfam_id in kv[1]
             ['pfamDomainAccession'] else ... for kv in SSU_map.items()]

    return list(set(nomenclature))


def get_rna_nomenclature(polymer):  

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
        matches = re.search(i[1], polymer['rcsb_polymer_entity'][ 'pdbx_description' ], flags=re.IGNORECASE | re.MULTILINE)
        if matches!=None:
            return [i[0]]
    return []


def inferOrganismsFromPolymers(polymers: list):

    host_organism_names: list[str] = []
    src_organism_names : list[str] = []
    host_organism_ids  : list[int] = []
    src_organism_ids   : list[int] = []

    for polymer in polymers:
        src_organism_names = [*src_organism_names, *polymer['src_organism_names']] if polymer[ 'src_organism_names' ]     != None else src_organism_names
        src_organism_ids   = [*src_organism_ids  , *polymer['src_organism_ids']] if polymer['src_organism_ids']           != None else src_organism_ids
        src_organism_names = [*src_organism_names, *polymer[ 'host_organism_names' ]] if polymer[ 'host_organism_names' ] != None else src_organism_names
        src_organism_ids   = [*src_organism_ids  , *polymer[ 'host_organism_ids' ]] if polymer[ 'host_organism_ids' ]     != None else src_organism_ids

    return {
        "src_organism_ids   ": list(set(src_organism_ids)),
        "src_organism_names ": list(set(src_organism_names)),
        "host_organism_ids  ": list(set(host_organism_ids)),
        "host_organism_names": list(set(host_organism_names))}


def extract_external_refs(external_refs):
    """
    external_refs: list[{ link: string; type: string; id: string }]
    """

    externalRefIds: list[str] = []
    externalRefTypes: list[str] = []
    externalRefLinks: list[str] = []

    if external_refs == None:
        ...
    else:
        for ref in external_refs:
            externalRefIds.append(ref['id'])
            externalRefTypes.append(ref['type'])
            externalRefLinks.append(ref['link'])

    return [externalRefIds, externalRefTypes, externalRefLinks]


def reshape_to_ligand(nonpoly):
    return {
        "pdbx_description   ": nonpoly['rcsb_nonpolymer_entity']['pdbx_description'],
        "formula_weight     ": nonpoly['rcsb_nonpolymer_entity']['formula_weight'],
        "chemicalId         ": nonpoly['pdbx_entity_nonpoly']['comp_id'],
        "chemicalName       ": nonpoly['pdbx_entity_nonpoly']['name'],
        "number_of_instances": nonpoly['rcsb_nonpolymer_entity']['pdbx_number_of_molecules']
    }


def is_ligand_like(polymer, nomenclature: list[str]):
    if 'tRNA' in nomenclature or 'mRNA' in nomenclature:
        return True
    #   // ? Look for enzymes, factors and antibiotics
    reg = r"/(\w*(?<!(cha|pro|dom|str))in\b)|(\b\w*zyme\b)|(factor)/gi"
    matches = re.search(
        reg, polymer['rcsb_polymer_entity']['pdbx_description'])

    if matches != None  \
            and not True in ['protein' not in _.lower() for _ in matches]  \
            and not ('protein' in polymer['rcsb_polymer_entity']['pdbx_description'].lower()) \
            and not ('rna' in polymer['rcsb_polymer_entity']['pdbx_description'].lower()):
        return True
    else:
        return False


def reshape_poly_to_rna(plm) -> list:

    src_organism_ids    = list(set([org['ncbi_taxonomy_id'] for org in plm['rcsb_entity_source_organism']])) if plm['rcsb_entity_source_organism'] != None else []
    host_organism_ids   = list(set([org['ncbi_taxonomy_id'] for org in plm['rcsb_entity_host_organism'  ]])) if plm['rcsb_entity_host_organism'  ] != None else []
    src_organism_names  = list(set([org['ncbi_taxonomy_id'] for org in plm['rcsb_entity_source_organism']])) if plm['rcsb_entity_source_organism'] != None else []
    host_organism_names = list(set([org['ncbi_taxonomy_id'] for org in plm['rcsb_entity_host_organism'  ]])) if plm['rcsb_entity_host_organism'  ] != None else []

    nomenclature = get_rna_nomenclature(plm)

    return [
        {
            "nomenclature": nomenclature,
            "ligand_like": is_ligand_like(plm, nomenclature),

            "asym_ids": plm['rcsb_polymer_entity_container_identifiers']['asym_ids'],
            "auth_asym_id": auth_asym_id,
            "parent_rcsb_id": plm['entry']['rcsb_id'],

            "host_organism_ids": host_organism_ids,
            "host_organism_names": host_organism_names,
            "src_organism_ids": src_organism_ids,
            "src_organism_names": src_organism_names,

            "rcsb_pdbx_description": "" if plm['rcsb_polymer_entity']['pdbx_description'] == None else plm['rcsb_polymer_entity']['pdbx_description'],

            "entity_poly_strand_id": plm['entity_poly']['pdbx_strand_id'],
            "entity_poly_seq_one_letter_code": plm['entity_poly']['pdbx_seq_one_letter_code'],
            "entity_poly_seq_one_letter_code_can": plm['entity_poly']['pdbx_seq_one_letter_code_can'],
            "entity_poly_seq_length": plm['entity_poly']['rcsb_sample_sequence_length'],
            "entity_poly_entity_type": plm['entity_poly']['type'],
            "entity_poly_polymer_type": plm['entity_poly']['rcsb_entity_polymer_type']
        }

        for auth_asym_id in plm['rcsb_polymer_entity_container_identifiers']['auth_asym_ids']]


def reshape_poly_to_protein(plm):
    if plm['pfams'] != None and len(plm['pfams']) > 0:

        pfam_comments = list(set([pfam['rcsb_pfam_comment']
                             for pfam in plm['pfams']]))
        pfam_descriptions = list(
            set([pfam['rcsb_pfam_description'] for pfam in plm['pfams']]))
        pfam_accessions = list(
            set([pfam['rcsb_pfam_accession'] for pfam in plm['pfams']]))

    else:
        pfam_comments = []
        pfam_descriptions = []
        pfam_accessions = []

    src_organism_ids = list(set([org['ncbi_taxonomy_id'] for org in plm['rcsb_entity_source_organism']])
                            ) if plm['rcsb_entity_host_organism'] != None else []
    src_organism_names = list(set([org['scientific_name'] for org in plm['rcsb_entity_source_organism']])
                              ) if plm['rcsb_entity_host_organism'] != None else []
    host_organism_ids = list(set([org['ncbi_taxnonomy_id'] for org in plm['rcsb_entity_host_organism']])
                             ) if plm['rcsb_entity_host_organism'] != None else []
    host_organism_names = list(set([org['scientific_name'] for org in plm['rcsb_entity_host_organism']])
                               ) if plm['rcsb_entity_host_organism'] != None else []

    nomenclature = get_protein_nomenclature(plm)

    return [
        {
            "nomenclature": nomenclature,
            "asym_ids": plm['rcsb_polymer_entity_container_identifiers']['asym_ids'],
            "parent_rcsb_id": plm['entry']['rcsb_id'],
            "auth_asym_id": auth_asym_id,
            "pfam_accessions": pfam_accessions,
            "pfam_comments": pfam_comments,
            "pfam_descriptions": pfam_descriptions,
            "ligand_like": is_ligand_like(plm, nomenclature),
            "host_organism_ids": host_organism_ids,
            "host_organism_names": host_organism_names,
            "src_organism_ids": src_organism_ids,
            "src_organism_names": src_organism_names,
            "uniprot_accession": [entry['rcsb_id'] for entry in plm['uniprots']] if plm['uniprots'] != None and len(plm['uniprots']) > 0 else [],
            "rcsb_pdbx_description": plm['rcsb_polymer_entity']['pdbx_description'],
            "entity_poly_strand_id": plm['entity_poly']['pdbx_strand_id'],
            "entity_poly_seq_one_letter_code": plm['entity_poly']['pdbx_seq_one_letter_code'],
            "entity_poly_seq_one_letter_code_can": plm['entity_poly']['pdbx_seq_one_letter_code_can'],
            "entity_poly_seq_length": plm['entity_poly']['rcsb_sample_sequence_length'],
            "entity_poly_entity_type": plm['entity_poly']['type'],
            "entity_poly_polymer_type": plm['entity_poly']['rcsb_entity_polymer_type']
        }


        for auth_asym_id in plm['rcsb_polymer_entity_container_identifiers']['auth_asym_ids']]


def process_pdb_record(pdb_api_response):
    """at the level of @entry, so response['data']['entry'] is the pdb record"""

    response = pdb_api_response
    polys = response['polymer_entities']
    ligands = response['nonpolymer_entities']

    proteins = polys
    def is_protein(
        poly): return poly['entity_poly']['rcsb_entity_polymer_type'] == 'Protein'

    proteins, rnas = [], []
    for poly in polys:
        proteins.append(poly) if is_protein(poly) else rnas.append(poly)

    reshaped_proteins = []
    reshaped_rnas     = []

    [reshaped_proteins.append(*reshape_poly_to_protein(poly))
     for poly in proteins]
    [reshaped_rnas.append(*reshape_poly_to_rna(poly)) for poly in rnas]

    reshaped_nonpoly = [reshape_to_ligand(nonpoly) for nonpoly in ligands] if ligands != None and len(ligands) > 0 else []
    organisms        = inferOrganismsFromPolymers(reshaped_proteins)
    externalRefs     = extract_external_refs(response['rcsb_external_references'])

    pub              = response['citation'][0]

    kwords_text = response['struct_keywords']['text'] if response['struct_keywords'] != None else None
    kwords = response['struct_keywords']['pdbx_keywords'] if response['struct_keywords'] != None else None

    reshaped = {
        "rcsb_id"               : response['rcsb_id'],
        "expMethod"             : response['exptl'][0]['method'],
        "resolution"            : response['rcsb_entry_info']['resolution_combined'][0],
        "rcsb_external_ref_id"  : externalRefs[0],
        "rcsb_external_ref_type": externalRefs[1],
        "rcsb_external_ref_link": externalRefs[2],
        "citation_year"         : pub['year'],
        "citation_rcsb_authors" : pub['rcsb_authors'],
        "citation_title"        : pub['title'],
        "citation_pdbx_doi"     : pub['pdbx_database_id_DOI'],
        "pdbx_keywords_text"    : kwords_text,
        "pdbx_keywords"         : kwords,
        "proteins"              : reshaped_proteins,
        "rnas"                  : reshaped_rnas,
        "ligands"               : reshaped_nonpoly,
        **organisms,
    }

    return reshaped


def query_rcsb_api(gql_string: str):
    reqstring = "https://data.rcsb.org/graphql?query={}".format(gql_string)

    try:
        resp = requests.get(reqstring)
        return resp.json()['data']['entry']

    except Exception as e:
        print("Could not land request to RCSB API. {}".format(e))


RCSB_ID = "4ug0"
mono    = query_rcsb_api(gql_monolith(RCSB_ID))

pprint(process_pdb_record(mono))
