import json
from pprint import pprint
from typing import Any
import requests
from api.ribctl.lib.types.types_ribosome import RNA, Ligand, Protein, RibosomeStructure
from ribctl.lib import LSU_map, SSU_map
from ribctl.lib.gql_querystrings import monolithic
import re


def gql_monolith(rcsb_id): return monolithic.replace(
    "$RCSB_ID", rcsb_id.upper())

# gql_structs             = lambda rcsb_id: structure_string.replace("$RCSB_ID", rcsb_id.upper())
# gql_polymer_entities    = lambda rcsb_id: polymer_entities_string.replace("$RCSB_ID", rcsb_id.upper())
# gql_nonpolymer_entities = lambda rcsb_id: nonpolymer_entities_string.replace("$RCSB_ID", rcsb_id.upper())


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

def __reshape_to_ligand(nonpoly)->Ligand:
    return Ligand(**{
        "pdbx_description": nonpoly['rcsb_nonpolymer_entity']['pdbx_description'],
        "formula_weight": nonpoly['rcsb_nonpolymer_entity']['formula_weight'],
        "chemicalId": nonpoly['pdbx_entity_nonpoly']['comp_id'],
        "chemicalName": nonpoly['pdbx_entity_nonpoly']['name'],
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

    nomenclature = __get_protein_nomenclature(plm)

    prots =  [
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
    print( "returning reshaped:",  len(prots) )
    return prots

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
        print("Processing poly:", poly['entry']['rcsb_id'])
        proteins.append(poly) if is_protein(poly) else rnas.append(poly)

    reshaped_proteins:list[Protein] = []
    reshaped_rnas    :list[RNA]     = []

    [reshaped_rnas.extend(__reshape_poly_to_rna(poly)) for poly in rnas]

    print("Proteins is :", proteins)
    for poly in proteins:
        print("invoking reshe")
        protein  = __reshape_poly_to_protein(poly)
        print("GOt least with proteins:", protein)
        reshaped_proteins.extend(__reshape_poly_to_protein(poly)) 

    print("Reshaped proteins:", reshaped_proteins)
    reshaped_nonpoly:list[Ligand] = [__reshape_to_ligand(nonpoly) for nonpoly in nonpoly_entities] if nonpoly_entities != None and len(nonpoly_entities) > 0 else []
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

    print("returning 6s0k", reshaped.dict())
    print("reshp rptoeins", reshaped_proteins)

    return reshaped

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



# Functionality to migrate from older cli:

# 4 Scripts
# Assets verification

        # await x.initialize_assets(flags.commit || flags.repair)
        # await x.initialize_ligands(flags.commit || flags.repair, x.structure as RibosomeStructure)

        # if (flags.commit) {
        #     commit_struct_to_Db(rcsb_id)
# export class RibosomeAssets {
#     rcsb_id: string

#     constructor(
    #     rcsb_id: string
    # ) {
    #     this.rcsb_id = rcsb_id
    # }
    # envcheck = () => {
    #     if (!process.env["RIBETL_DATA"]) {
    #         throw Error("RIBETL_DATA environment variable not set. Cannot access assets.")
    #     }
    # }

    # folder_path = () => {
    #     this.envcheck()
    #     return `${process.env["RIBETL_DATA"]}/${this.rcsb_id}`
    # }

    # cif_filepath = () => {
    #     this.envcheck()
    #     return `${this.folder_path()}/${this.rcsb_id}.cif`
    # }

    # cif_modified_filepath = () => {
    #     this.envcheck()
    #     return `${this.folder_path()}/${this.rcsb_id}_modified.cif`
    # }

    # json_profile_filepath = () => {
    #     this.envcheck()
    #     return `${this.folder_path()}/${this.rcsb_id}.json`
    # }

    # chains_folder = () => {
    #     this.envcheck()
    #     return `${this.folder_path()}/CHAINS`
    # }

    # png_thumbnail_filepath = () => {
    #     this.envcheck()
    #     return `${this.folder_path()}/_ray_${this.rcsb_id}.png`
    # }


#     async __verify_cif(obtain: boolean = false):Promise<boolean> {
#         if (existsSync(this.cif_filepath())) {
#             return true
#         } else {

#             if (obtain) {
#                 await download_unpack_place(this.rcsb_id); return true
#             } else return false
#         }
#     }
#     async __verify_cif_modified(obtain: boolean = false) {
#         if (existsSync(this.cif_modified_filepath())) { return true } else {
#             if (obtain) {
#                 this.__verify_cif(true)
#                 let y = exec(`${process.env["PYTHONBIN"]} ${process.env["SPLIT_RENAME_PY"]} -s ${this.rcsb_id}`,
#                     (err, stdout, stderr) => {
#                         console.log(err);
#                         console.log(stdout.toString());
#                         console.log(stderr.toString());
#                     })
#                 console.log(y.stdout)
#             } else return false
#         }
#     }
#     async __verify_json_profile(obtain: boolean = false) {
#         if (existsSync(this.json_profile_filepath())) { return true } else {
#             if (obtain) {
#                 let ribosome = await processPDBRecord(this.rcsb_id)
#                 let filename = await save_struct_profile(ribosome)
#                 process.stdout.write(`Saved structure profile:\t${filename}`);
#             } else return false
#         }
#     }
#     async __verify_png_thumbnail(obtain: boolean = false) {
#         if (existsSync(this.png_thumbnail_filepath())) { return true } else {
#             if (obtain) {
#                await new Promise<void>((rs,rj)=>{
#                exec(`${process.env["PYTHONBIN"]} ${process.env["RENDER_THUMBNAIL_PY"]} -s ${this.rcsb_id}`,
#                     (err, stdout, stderr) => {
#                         if (err !==0 ){
#                             rj(stderr)
#                         }
#                         else{rs()}
#                     })
#                })
#             } else return false
#         }
#     }
#     async __verify_chains_folder(obtain: boolean = false) {
#         if (existsSync(this.chains_folder())) { return true } else {
#             if (obtain) {
#                 this.__verify_cif(true)
#                 exec(`${process.env["PYTHONBIN"]} ${process.env["SPLIT_RENAME_PY"]} -s ${this.rcsb_id}`, (err, stdout, stderr) => {
#                     console.log(err);
#                     console.log(stdout.toString());
#                     console.log(stderr.toString());
#                 })
#             } else return false
#         }
#     }
# #     //verify that each chain file exists
#     async __verify_chain_files(parent_structure: RibosomeStructure): Promise<boolean> {
#         if (!this.__verify_chains_folder()) {
#             return false
#         }
#         const chain_files_all = [
#             parent_structure.proteins?.map((c) => { return c.auth_asym_id }),

#             (parent_structure.rnas || []).map((c) => { return c.auth_asym_id })
#         ].map(
#             (chain_id) => {
#                 if (!existsSync(`${this.chains_folder()}/${this.rcsb_id}_STRAND_${chain_id}.cif`)) {
#                     console.log(`[${this.rcsb_id}]: NOT FOUND ${this.chains_folder()}/${this.rcsb_id}_STRAND_${chain_id}.cif`)
#                     return false
#                 } else return true
#             }
#         ).reduce((prev, cur) => { return prev && cur }, true)
#         return chain_files_all
#     }

#     async __verify_ligands_and_polymers(struct: RibosomeStructure) {
#         let ligs = struct.ligands && struct.ligands.map((lig_chem_id) => {
#             if (!existsSync(`${this.folder_path()}/LIGAND_${lig_chem_id}.json`)) {
#                 console.log(`[${this.rcsb_id}]: NOT FOUND ${this.folder_path()}/LIGAND_${lig_chem_id}.json (Is it an ION? Expected.)`)
#                 return false
#             } else return true
#         }).reduce((prev, cur) => { return prev && cur }, true)


#         let ligandlike: string[] = []
#         let ligands: string[] = []

#         for (var chain of [...struct.proteins, ...(struct.rnas || [])]) {
#             if (chain.ligand_like) {
#                 ligandlike = [...ligandlike, chain.auth_asym_id]
#             }
#         }
#         if (!struct.ligands) {
#         } else {
#             for (var lig of struct.ligands) {
#                 ligands = [...ligands, lig.chemicalId]
#             }
#         }


#         let polys = ligandlike.map((polymer_id) => {
#             if (!existsSync(`${this.folder_path()}/POLYMER_${polymer_id}.json`)) {
#                 console.log(`[${this.rcsb_id}]: NOT FOUND ${this.folder_path()}/POLYMER_${polymer_id}.json`)
#                 return false
#             } else return true
#         })

#         if ((!polys || !ligs)) {
#             console.log("Some ligands are missing. Calling script:", process.env["PYTHONBIN"], process.env["EXTRACT_BSITES_PY"])
#             exec(`${process.env["PYTHONBIN"]} ${process.env["EXTRACT_BSITES_PY"]} -s ${this.rcsb_id} --save`, (err, stdout, stderr) => {
#                 console.log(err);
#                 console.log(stdout);
#                 console.log(stderr);
#             })
#         }
#     }


#     async init_assets(obtain: boolean = false) {
#     }
# }

# export class StructureFolder {

#     rcsb_id   : string;
#     assets    : RibosomeAssets;
#     structure?: RibosomeStructure;

#     constructor(rcsb_id: string) {
#         this.rcsb_id = rcsb_id.toUpperCase()
#         this.assets  = new RibosomeAssets(this.rcsb_id)
#     }

#     async initialize_assets(obtain: boolean) {
#         console.log("Initializing assets");

#         await this.assets.__verify_json_profile(true)
#         await this.assets.__verify_cif(true)
#         await this.assets.__verify_cif_modified(true)
#         await this.assets.__verify_chains_folder(true)
#         // await this.assets.__verify_png_thumbnail(true)

#         console.log("Initializing assets");
#         if (!this.assets.__verify_json_profile(obtain)) {
#             throw Error(`Structure ${this.rcsb_id} assets not found. Cannot initiate resource.`)
#         }

#         console.log("accessing filepath");
#         this.structure = JSON.parse(readFileSync(this.assets.json_profile_filepath(), 'utf-8'))
#     }


#     async initialize_ligands(obtain: boolean, ribosome: RibosomeStructure) {
#         let ligs = ribosome.ligands && ribosome.ligands.map((lig) => {
#             console.log("Looking through chemids" , lig)
#             if (!existsSync(`${this.assets.folder_path()}/LIGAND_${lig}.json`)) {
#                 console.log(`[${this.rcsb_id}]: NOT FOUND ${this.assets.folder_path()}/LIGAND_${lig.chemicalId}.json (Is it an ION? Expected.)`)
#                 return false
#             } else return true
#         }).reduce((prev, cur) => { return prev && cur }, true)

#         let ligandlike: string[] = []
#         for (var chain of [...ribosome.proteins, ...(ribosome.rnas || [])]) {
#             if (chain.ligand_like) {
#                 ligandlike = [...ligandlike, chain.auth_asym_id]
#             }
#         }

#         let polys = ligandlike && ligandlike.map((polymer_id) => {
#             if (!existsSync(`${this.assets.folder_path()}/POLYMER_${polymer_id}.json`)) {
#                 console.log(`[${this.rcsb_id}]: NOT FOUND ${this.assets.folder_path()}/POLYMER_${polymer_id}.json`)
#                 return false
#             } else return true
#         })

#         if (obtain) {
#             console.log("Some ligands are missing. Calling script:", `${process.env["PYTHONBIN"]} ${process.env["EXTRACT_BSITES_PY"]} -s ${this.rcsb_id} --save`)
#             exec(`${process.env["PYTHONBIN"]} ${process.env["EXTRACT_BSITES_PY"]} -s ${this.rcsb_id} --save`, (err, stdout, stderr) => {
#                 console.log(err);
#                 console.log(stdout);
#                 console.log(stderr);
#             })
#         }

#     }
# }

# /**
#  * Request and display the state of the given rcsb_id structure in the database instance.
#  */
# export const queryStructDb = (rcsb_id: string) => {
#     return new Promise<string[]>((resolve, reject) => {
#         let y = cp.exec(`echo \"match (struct:RibosomeStructure {rcsb_id:\\"${rcsb_id.toUpperCase()}\\"}) return struct.rcsb_id\" | cypher-shell -a \"${process.env["NEO4J_URI"]}\" --format plain -u ${process.env["NEO4J_USER"]} -p ${process.env["NEO4J_PASSWORD"]} --database ${process.env["NEO4J_CURRENTDB"]}`,
#             { env: process.env },
#             (err, stdout, stderr) => {
#                 if (err && err?.code != 0) {
#                     process.stdout.write("Got shell error " + stderr + stdout)
#                     console.log("Got Error code:", err?.code)
#                     reject(err)
#                 }

#                 const dbstructs = stdout != null ? (stdout as string).replace(/"/g, '').split("\n").filter(r => r.length === 4) : []
#                 console.log(dbstructs)
#                 resolve(dbstructs)
#             })
#     })

# }
# /**
#  * Download a .cif model of the structure.
#  * @param struct_id
#  */
# export const download_unpack_place = async (struct_id: string) => {
#     const BASE_URL = "http://files.rcsb.org/download/"
#     const FORMAT = ".cif.gz"

#     const structid = struct_id.toUpperCase()
#     let url = BASE_URL + structid + FORMAT
#     let compressed: Buffer = await axios.get(url, { responseType: 'arraybuffer' }).then(r => { return r.data })
#         .catch(e => { process.stdout.write(`Structure ${structid} failed: `, e); return []; })
#     let decompressed = await ungzip(compressed);

#     // let destination_chains = path.join(
#     //   process.env["RIBETL_DATA"] as string,
#     //   `${structid}`,
#     //   `CHAINS`)

#     // if (!existsSync(destination_chains)) {
#     //   mkdirSync(destination_chains)
#     //   process.stdout.write(`Created directory ${destination_chains}.`);
#     // }
#     let structfile = path.join(
#         process.env["RIBETL_DATA"] as string,
#         `${structid}`,
#         `${structid}.cif`)
#     writeFileSync(structfile, decompressed)
# }

# export const save_struct_profile = (r: RibosomeStructure): string => {
#     var rcsb_id = r.rcsb_id;
#     var target_filename = path.join(
#         process.env["RIBETL_DATA"] as string,
#         rcsb_id.toUpperCase(),
#         rcsb_id.toUpperCase() + ".json"
#     );

#     if (!existsSync(path.dirname(target_filename))) {
#         mkdirSync(path.dirname(target_filename));
#     }
#     writeFileSync(target_filename, JSON.stringify(r, null, 4));
#     return target_filename
# }

# export const commit_struct_to_Db = (rcsb_id: string) => {
#     console.log(`Commiting ${rcsb_id} to the database`)
#     const commit_script = process.env["COMMIT_STRUCTURE_SH"]
#     let   current_db    = process.env["NEO4J_CURRENTDB"]
#     let   uri           = process.env["NEO4J_URI"]
#     let   invocation    = `${commit_script} -s ${rcsb_id} -d ${current_db} -a "${uri}"`
#     console.log("Invoking:", invocation);
#     let   proc          = cp.exec(invocation)

#     if (proc.stderr !== null) {
#         proc.stderr.on("data",
#             (data) => {
#                 console.log(data)
#             }
#         )
#     }

#     proc.stdout?.on("data", (data) => { console.log(data) })
#     proc.on("exit", () => { process.exit() })
# }
# export const commit_struct_to_db_sync = (rcsb_id: string): Promise<void> => {
#     console.log(`Commiting ${rcsb_id} to the database`)
#     const commit_script = process.env["COMMIT_STRUCTURE_SH"]
#     let current_db = process.env["NEO4J_CURRENTDB"]
#     let uri = process.env["NEO4J_URI"]
#     let invocation = `${commit_script} -s ${rcsb_id} -d ${current_db} -a "${uri}"`
#     return new Promise<void>((resolve, reject) => {
#         let proc = cp.exec(invocation)
#         if (proc.stderr !== null) {
#             proc.stderr.on("data",
#                 (data) => {
#                     console.log(data)
#                 }
#             )
#         }

#         proc.stdout?.on("data", (data) => { console.log(data) })
#         proc.stderr?.on("data", (data) => { console.log(data); reject(data) })
#         proc.on("exit", () => { process.exit(); resolve() })

#     })

# }
