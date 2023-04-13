import json
import os
from pathlib import Path
from typing import Any
import requests
from ribctl.lib.types.types_ribosome import RNA, NonpolymericLigand, Protein, RibosomeStructure
from ribctl.lib.gql_querystrings import monolithic
import re

p        = Path(__file__).parents[1]
lsu_path = os.path.join(p,'assets','subunit_map_LSU.json')
ssu_path = os.path.join(p,'assets','subunit_map_SSU.json')

LSU_map = {k: v for k, v in json.load(open(lsu_path, 'r')).items()}
SSU_map = {k: v for k, v in json.load(open(ssu_path, 'r')).items()}


def __get_protein_nomenclature(protein):
    banregex = r"/\b([ueb][ls]\d{1,2})\b/gi"
    #     check authors's annotations. if classes are present --> use that.
    finds = re.search(
        banregex, protein['rcsb_polymer_entity']['pdbx_description'])

    if (finds != None):
        firstcap: str = finds[0]
        classname: str = firstcap[0].lower(
        ) + firstcap[1].upper() + firstcap[2:]
        return [classname]

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

def __get_rna_nomenclature(polymer):

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
        matches = re.search(i[1], polymer['rcsb_polymer_entity']
                            ['pdbx_description'], flags=re.IGNORECASE | re.MULTILINE)
        if matches != None:
            return [i[0]]
    return []

def __infer_organisms_from_polymers(polymers: list[RNA|Protein]):

    host_organism_names: list[str] = []
    src_organism_names : list[str] = []
    host_organism_ids  : list[int] = []
    src_organism_ids   : list[int] = []

    for polymer in polymers:
        src_organism_names = [*src_organism_names, *polymer.src_organism_names
                              ] if polymer.src_organism_names != None else src_organism_names
        src_organism_ids = [*src_organism_ids, *polymer.src_organism_ids
                            ] if polymer.src_organism_ids != None else src_organism_ids
        src_organism_names = [*src_organism_names, *polymer.host_organism_names
                              ] if polymer.host_organism_names != None else src_organism_names
        src_organism_ids = [*src_organism_ids, *polymer.host_organism_ids
                            ] if polymer.host_organism_ids != None else src_organism_ids

    return {
        "src_organism_ids": list(map(int, set(src_organism_ids))),
        "src_organism_names": list(map(str, set(src_organism_names))),
        "host_organism_ids": list(map(int, set(host_organism_ids))),
        "host_organism_names": list(map(str, set(host_organism_names)))
    }

def __extract_external_refs(external_refs):
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

def __reshape_to_nonpolymericligand(nonpoly)->NonpolymericLigand:
    return NonpolymericLigand(
        chemicalId          = nonpoly['pdbx_entity_nonpoly']['comp_id'],
        chemicalName        = nonpoly['pdbx_entity_nonpoly']['name'],
        pdbx_description    = nonpoly['rcsb_nonpolymer_entity']['pdbx_description'],
        formula_weight      = nonpoly['rcsb_nonpolymer_entity']['formula_weight'],
        number_of_instances = nonpoly['rcsb_nonpolymer_entity']['pdbx_number_of_molecules'],
        #TODO
        nomenclature        = ["Analog"]
    )
    return NonpolymericLigand.parse_obj({
        "pdbx_description"   : nonpoly['rcsb_nonpolymer_entity']['pdbx_description'],
        "formula_weight"     : nonpoly['rcsb_nonpolymer_entity']['formula_weight'],
        "chemicalId"         : nonpoly['pdbx_entity_nonpoly']['comp_id'],
        "chemicalName"       : nonpoly['pdbx_entity_nonpoly']['name'],
        "number_of_instances": nonpoly['rcsb_nonpolymer_entity']['pdbx_number_of_molecules']
    })

def __is_ligand_like(polymer, nomenclature: list[str]):
    if 'tRNA' in nomenclature or 'mRNA' in nomenclature:
        return True
    #   // ? Look for enzymes, factors and antibiotics
    reg = r"/(\w*(?<!(cha|pro|dom|str))in\b)|(\b\w*zyme\b)|(factor)/gi"
    match = re.search(reg, polymer['rcsb_polymer_entity']['pdbx_description'])
    if match != None  \
            and not 'protein' not in match.string.lower()  \
            and not ('protein' in polymer['rcsb_polymer_entity']['pdbx_description'].lower()) \
            and not ('rna' in polymer['rcsb_polymer_entity']['pdbx_description'].lower()):
        return True
    else:
        return False

def __reshape_poly_to_rna(plm) -> list[RNA]:
    """this returns a list because certain polymers accounts for multiple RNA molecules"""

    # class RCSBOrganism:
    #     ncbi_taxonomy_id:int | None
    #     scientific_name :str | None


    host_organisms  :list[Any] | None = plm['rcsb_entity_host_organism']
    source_organisms:list[Any] | None = plm['rcsb_entity_source_organism']

    host_organism_ids   = []
    host_organism_names = []
    src_organism_ids   = []
    src_organism_names = []


    if host_organisms != None:
        for ho in host_organisms:
            if ho['ncbi_taxonomy_id'] != None:
                host_organism_ids.append(ho['ncbi_taxonomy_id'])
            if ho['scientific_name'] != None:
                host_organism_names.append(ho['scientific_name'])

    if source_organisms != None:
        for so in source_organisms:
            if so[ 'ncbi_taxonomy_id' ] != None:
                src_organism_ids.append(so['ncbi_taxonomy_id'])
            if so['scientific_name'] != None:
                src_organism_names.append(so['scientific_name'])

    host_organism_ids   = list(map(int, set(host_organism_ids)))
    host_organism_names = list(map(str, set(host_organism_names)))
    src_organism_ids    = list(map(int, set(src_organism_ids)))
    src_organism_names  = list(map(str, set(src_organism_names)))

    # # ------------
    # host_organism_ids   = list(map(int, set([org['ncbi_taxonomy_id'] for org in plm['rcsb_entity_host_organism'  ]]))) if plm['rcsb_entity_host_organism'  ] != None else []
    # host_organism_names = list(map(str, set([org['scientific_name'] for org in plm['rcsb_entity_host_organism'  ]]))) if plm['rcsb_entity_host_organism'  ] != None else []

    # src_organism_ids   = list(map(int, set([org['ncbi_taxonomy_id'] for org in plm['rcsb_entity_source_organism']]))) if plm['rcsb_entity_source_organism'] != None else []
    # src_organism_names = list(map(str, set([org['scientific_name'] for org in plm['rcsb_entity_source_organism']]))) if plm['rcsb_entity_source_organism']  != None else []

    nomenclature = __get_rna_nomenclature(plm)

    return [
        RNA(**{
            "nomenclature": nomenclature,
            "ligand_like": __is_ligand_like(plm, nomenclature),

            "asym_ids": plm['rcsb_polymer_entity_container_identifiers']['asym_ids'],
            "auth_asym_id": auth_asym_id,
            "parent_rcsb_id": plm['entry']['rcsb_id'],

            "host_organism_ids"  : host_organism_ids,
            "host_organism_names": host_organism_names,
            "src_organism_ids"   : src_organism_ids,
            "src_organism_names" : src_organism_names,

            "rcsb_pdbx_description": "" if plm['rcsb_polymer_entity']['pdbx_description'] == None else plm['rcsb_polymer_entity']['pdbx_description'],

            "entity_poly_strand_id": plm['entity_poly']['pdbx_strand_id'],
            "entity_poly_seq_one_letter_code": plm['entity_poly']['pdbx_seq_one_letter_code'],
            "entity_poly_seq_one_letter_code_can": plm['entity_poly']['pdbx_seq_one_letter_code_can'],
            "entity_poly_seq_length": plm['entity_poly']['rcsb_sample_sequence_length'],
            "entity_poly_entity_type": plm['entity_poly']['type'],
            "entity_poly_polymer_type": plm['entity_poly']['rcsb_entity_polymer_type']
        })

        for auth_asym_id in plm['rcsb_polymer_entity_container_identifiers']['auth_asym_ids']]

def __reshape_poly_to_protein(plm)->list[Protein]:
    if plm['pfams'] != None and len(plm['pfams']) > 0:

        pfam_comments = list(set([pfam['rcsb_pfam_comment']
                             for pfam in plm['pfams']]))
        pfam_descriptions = list(
            set([pfam['rcsb_pfam_description'] for pfam in plm['pfams']]))
        pfam_accessions = list(
            set([pfam['rcsb_pfam_accession'] for pfam in plm['pfams']]))

    else:
        pfam_comments     = []
        pfam_descriptions = []
        pfam_accessions   = []

    host_organisms  :list[Any] | None = plm['rcsb_entity_host_organism']
    source_organisms:list[Any] | None = plm['rcsb_entity_source_organism']

    host_organism_ids   = []
    host_organism_names = []

    src_organism_ids    = []
    src_organism_names  = []


    if host_organisms != None:
        for ho in host_organisms:
            if ho['ncbi_taxonomy_id'] != None:
                host_organism_ids.append(ho['ncbi_taxonomy_id'])
            if ho['scientific_name'] != None:
                host_organism_names.append(ho['scientific_name'])

    if source_organisms != None:
        for so in source_organisms:
            if so[ 'ncbi_taxonomy_id' ] != None:
                src_organism_ids.append(so['ncbi_taxonomy_id'])
            if so['scientific_name'] != None:
                src_organism_names.append(so['scientific_name'])

    host_organism_ids   = list(map(int, set(host_organism_ids)))
    host_organism_names = list(map(str, set(host_organism_names)))
    src_organism_ids    = list(map(int, set(src_organism_ids)))
    src_organism_names  = list(map(str, set(src_organism_names)))

    nomenclature = __get_protein_nomenclature(plm)

    return  [
        Protein(**{
            "nomenclature": nomenclature,
            "asym_ids": plm['rcsb_polymer_entity_container_identifiers']['asym_ids'],
            "parent_rcsb_id": plm['entry']['rcsb_id'],
            "auth_asym_id": auth_asym_id,
            "pfam_accessions": pfam_accessions,
            "pfam_comments": pfam_comments,
            "pfam_descriptions": pfam_descriptions,
            "ligand_like": __is_ligand_like(plm, nomenclature),
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
        }) for auth_asym_id in plm['rcsb_polymer_entity_container_identifiers']['auth_asym_ids']
    ]

def gql_monolith(rcsb_id): return monolithic.replace("$RCSB_ID", rcsb_id.upper())

def current_rcsb_structs() -> list[str]:
    """Return all structures in the rcsb that contain the phrase RIBOSOME and have more than 25 protein entities"""

    rcsb_search_api = "https://search.rcsb.org/rcsbsearch/v2/query"
    params = {
        "query":
        {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "group",
                    "logical_operator": "and",
                    "nodes": [
                        {"type": "group", "logical_operator": "and", "nodes": [{"type": "terminal", "service": "text", "parameters": {
                            "operator": "contains_phrase", "negation": False, "value": "RIBOSOME", "attribute": "struct_keywords.pdbx_keywords"}}]},
                        {"type": "group", "logical_operator": "and", "nodes": [{"type": "terminal", "service": "text", "parameters": {
                            "operator": "greater", "negation": False, "value": 25, "attribute": "rcsb_entry_info.polymer_entity_count_protein"}}]},
                    ],
                    "label": "text"
                }
            ],
            "label": "query-builder"
        },
        "return_type": "entry",
        "request_options": {
            "return_all_hits": True,
            "results_verbosity": "compact"
        }
    }

    query = rcsb_search_api + "?json=" + json.dumps(params)
    return requests.get(query).json()['result_set']

def query_rcsb_api(gql_string: str) -> dict:

    reqstring = "https://data.rcsb.org/graphql?query={}".format(gql_string)
    _resp = requests.get(reqstring)
    resp = _resp.json()

    if 'data' in resp and 'entry' in resp['data']:
        return resp['data']['entry']
    else:
        raise Exception("No data found for query: {}".format(gql_string))

def process_pdb_record(rcsb_id: str) -> RibosomeStructure:
    """
    returns dict of the shape types_RibosomeStructure 
    """

    response         = query_rcsb_api(gql_monolith(rcsb_id.upper()))
    poly_entities    = response['polymer_entities']
    nonpoly_entities = response['nonpolymer_entities']

    def is_protein(poly): return poly['entity_poly']['rcsb_entity_polymer_type'] == 'Protein'

    proteins, rnas = [], []
    for poly in poly_entities:
        proteins.append(poly) if is_protein(poly) else rnas.append(poly)

    reshaped_proteins:list[Protein] = []
    reshaped_rnas    :list[RNA]     = []

    [reshaped_rnas.extend(__reshape_poly_to_rna(poly)) for poly in rnas]
    [reshaped_proteins.extend(__reshape_poly_to_protein(poly)) for poly in proteins]

    #TODO: Polymeric Factors extraction

    reshaped_nonpoly:list[NonpolymericLigand] = [__reshape_to_nonpolymericligand(nonpoly) for nonpoly in nonpoly_entities] if nonpoly_entities != None and len(nonpoly_entities) > 0 else []
    organisms        = __infer_organisms_from_polymers(reshaped_proteins)  # type: ignore (only accessing commong fields)
    externalRefs     = __extract_external_refs(response['rcsb_external_references'])

    if response['citation'] != None and len(response['citation']) > 0:
        pub         = response['citation'][0]
    else:
        pub = {
            "year"                   : None,
            "rcsb_authors"           : None,
            "title"                  : None,
            "pdbx_database_id_DOI"   : None,
            "pdbx_database_id_PubMed": None
        }
    
    kwords_text = response['struct_keywords']['text'] if response['struct_keywords']          != None else None
    kwords      = response['struct_keywords']['pdbx_keywords'] if response['struct_keywords'] != None else None


    reshaped = RibosomeStructure(**{
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
    })

    reshaped = RibosomeStructure(
        rcsb_id                 = response['rcsb_id'],

        expMethod               = response['exptl'][0]['method'],
        resolution              = response['rcsb_entry_info']['resolution_combined'][0],
        rcsb_external_ref_id    = externalRefs[0],
        rcsb_external_ref_type  = externalRefs[1],
        rcsb_external_ref_link  = externalRefs[2],
        citation_year           = pub['year'],
        citation_rcsb_authors   = pub['rcsb_authors'],
        citation_title          = pub['title'],
        citation_pdbx_doi       = pub['pdbx_database_id_DOI'],
        pdbx_keywords_text      = kwords_text,
        pdbx_keywords           = kwords,
        src_organism_ids= organisms['src_organism_ids'],
        src_organism_names= organisms['src_organism_names'],
        
        host_organism_ids= organisms['host_organism_ids'],
        host_organism_names= organisms['host_organism_names'],

        proteins                = reshaped_proteins,
        rnas                    = reshaped_rnas,
        nonpolymeric_ligands = reshaped_nonpoly,

        polymeric_factors= ...,
        assembly_map=...

    )

    return reshaped