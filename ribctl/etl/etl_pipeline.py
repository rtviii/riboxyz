import functools
import json
from pprint import pprint
from typing import Any, Optional
from pyhmmer.plan7 import HMM
import requests
from ribctl.etl.etl_polypeptides import (
    factor_classify,
    protein_classify,
    rna_classify,
    rp_hmm_dict_init,
    seq_prot_against_protclasses,
)
from ribctl.lib.types.types_ribosome import (
    RNA,
    AssemblyInstancesMap,
    NonpolymericLigand,
    Polymer,
    PolymericFactor,
    Protein,
    ProteinClass,
    RibosomeStructure,
)
from ribctl.etl.gql_querystrings import single_structure_graphql_template


"""These methods take a single rcsb_id and observes that the corresponding structure is:\
     - retrieved from RCSB
     - processed according to the pipeline of methods 
"""


def current_rcsb_structs() -> list[str]:
    """Return all structures in the rcsb that contain the phrase RIBOSOME and have more than 25 protein entities"""

    rcsb_search_api = "https://search.rcsb.org/rcsbsearch/v2/query"
    params = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "group",
                    "logical_operator": "and",
                    "nodes": [
                        {
                            "type": "group",
                            "logical_operator": "and",
                            "nodes": [
                                {
                                    "type": "terminal",
                                    "service": "text",
                                    "parameters": {
                                        "operator": "contains_phrase",
                                        "negation": False,
                                        "value": "RIBOSOME",
                                        "attribute": "struct_keywords.pdbx_keywords",
                                    },
                                }
                            ],
                        },
                        {
                            "type": "group",
                            "logical_operator": "and",
                            "nodes": [
                                {
                                    "type": "terminal",
                                    "service": "text",
                                    "parameters": {
                                        "operator": "greater",
                                        "negation": False,
                                        "value": 25,
                                        "attribute": "rcsb_entry_info.polymer_entity_count_protein",
                                    },
                                }
                            ],
                        },
                    ],
                    "label": "text",
                }
            ],
            "label": "query-builder",
        },
        "return_type": "entry",
        "request_options": {"return_all_hits": True, "results_verbosity": "compact"},
    }

    query = rcsb_search_api + "?json=" + json.dumps(params)
    return requests.get(query).json()["result_set"]

def query_rcsb_api(gql_string: str) -> dict:
    """This defines a query in the RCSB search language that identifies the structures we view as 'current' i.e. 40+ proteins, smaller than 4A resolution etc."""

    reqstring = "https://data.rcsb.org/graphql?query={}".format(gql_string)
    _resp = requests.get(reqstring)
    resp = _resp.json()

    if "data" in resp and "entry" in resp["data"]:
        return resp["data"]["entry"]
    else:
        raise Exception("No data found for query: {}".format(gql_string))

def rcsb_single_structure_graphql(rcsb_id):
    return single_structure_graphql_template.replace("$RCSB_ID", rcsb_id.upper())

class ReannotationPipeline:
    """
    ETL Pipeline as it currently stands takes care of injesting a graphql profile from RCSB and reshaping it into a RibosomeStructure.
    This is the class that oversees all of the annotations that added/edited on the _semantic profile_(as opposed to the structural files) by us, including:

    - metadat reshaping
    - rRNA classification
    - rProtein classification
    - factors annotation
    - ligand annotation



    """

    # ? Input data:
    rcsb_data_dict: dict

    # ? Initialized classification resources:
    hmm_ribosomal_proteins: dict[ProteinClass, HMM]
    # hmm_ribosomal_rnas:dict[RNAClass, HMM]
    # hmm_ribosomal_factors:dict[PolymericFactorClass, HMM]
    # hmm_ribosomal_ligands:dict[NonpolymericLigandClass, HMM]

    # ? Housekeeping
    asm_maps: list[AssemblyInstancesMap]
    rcsb_polymers: int
    rcsb_nonpolymers: int

    # rProteins: list[Protein] | None
    # rRNA     : list[RNA] | None

    def __init__(self, response: dict):
        self.rcsb_data_dict = response
        self.hmm_ribosomal_proteins = rp_hmm_dict_init()
        self.asm_maps = self.asm_parse(response["assemblies"])

        self.rcsb_polymers = len(self.rcsb_data_dict["polymer_entities"])
        self.rcsb_nonpolymers = (
            len(self.rcsb_data_dict["nonpolymer_entities"])
            if self.rcsb_data_dict["nonpolymer_entities"] != None
            else 0
        )

        # What is this garbage, you ask? Let me tell you.
        # The rcsb_data_dict contains a list of polymer entities, each of which has 4 properties(hidden in ..._container_identifiers):

        # - auth_asym_ids: this is the closest you have to an id of the polymer in the structure, but if there are multiple assemblies, this is an array of the two chains that are the same polymer in each assembly.
        # - asym_ids: this is the ids of the this AND all other chains that are the same assymteric unit as the given polymer. So, useless for purposes of identification:
        # - entity_id: so far as I can tell, this only serves to connect a polymer to its assembly container ( why this is not done via auth_asym_id directly is beyond me, legacy reasons this and that probably)
        # - pdbx_strand_id: confusing piece of garbage that sometimes overlaps with auth_asym_id

        # Ex. 4V8E.DD
        # assembly_id=1 asym_ids=['BB', 'CD', 'ED', 'ZA'] auth_asym_id='DD' parent_rcsb_id='4V8E' src_organism_names=[] host_organism_names=[] src_organism_ids=[] host_organism_ids=[] rcsb_pdbx_description='TRNA-TYR' entity_poly_strand_id='BB,BD,DB,DD' entity_poly_seq_one_letter_code='GGUGGGGUUCCCGAGCGGCCAAAGGGAGCAGACUGUA(MIA)AUCUGCCGUCAUCGACUUCGAAGGUUCGAAUCCUUCCCCCACCACCA' entity_poly_seq_one_letter_code_can='GGUGGGGUUCCCGAGCGGCCAAAGGGAGCAGACUGUAAAUCUGCCGUCAUCGACUUCGAAGGUUCGAAUCCUUCCCCCACCACCA' entity_poly_seq_length=85 entity_poly_polymer_type='RNA' entity_poly_entity_type='polyribonucleotide' nomenclature=['tRNA']

        # Given this state of affairs, if we want to be able to uniquely (coordinate-wise) identify each chain, we need to account for cases where PDB has provided two and sometimes four polymer collapsed into one representation.
        # We do this by counting assymetric_ids of each polymer and hold that as a

        self.flattened_polymers_target = functools.reduce(
            lambda count, poly: count
            + len(poly["rcsb_polymer_entity_container_identifiers"]["asym_ids"]),
            self.rcsb_data_dict["polymer_entities"],
            0,
        )

    #! Reshaping

    def infer_organisms_from_polymers(self, polymers: list[RNA | Protein]):
        """Grabbing taxid from every polymer in the structure to see which taxid prevails proportionally. Only needed because rcsb does not provide unequivocal taxid for structures (sometimes it's host+source)"""

        host_organism_names: list[str] = []
        src_organism_names : list[str] = []
        host_organism_ids  : list[int] = []
        src_organism_ids   : list[int] = []

        for polymer in polymers:
            src_organism_names = (
                [*src_organism_names, *polymer.src_organism_names]
                if polymer.src_organism_names != None
                else src_organism_names
            )
            src_organism_ids = (
                [*src_organism_ids, *polymer.src_organism_ids]
                if polymer.src_organism_ids != None
                else src_organism_ids
            )
            src_organism_names = (
                [*src_organism_names, *polymer.host_organism_names]
                if polymer.host_organism_names != None
                else src_organism_names
            )
            src_organism_ids = (
                [*src_organism_ids, *polymer.host_organism_ids]
                if polymer.host_organism_ids != None
                else src_organism_ids
            )

        return {
            "src_organism_ids"   : list(map(int, set(src_organism_ids))),
            "src_organism_names" : list(map(str, set(src_organism_names))),
            "host_organism_ids"  : list(map(int, set(host_organism_ids))),
            "host_organism_names": list(map(str, set(host_organism_names))),
        }

    def extract_external_refs(self, external_refs):
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
                externalRefIds.append(ref["id"])
                externalRefTypes.append(ref["type"])
                externalRefLinks.append(ref["link"])

        return [externalRefIds, externalRefTypes, externalRefLinks]

    def nonpoly_reshape_to_ligand(self, nonpoly: dict) -> NonpolymericLigand:
        # print("got a nonpoly elem")
        # pprint(nonpoly)

        # nonpoly["nonpolymer_comp"]["chem_comp"]["id"]
        # nonpoly["nonpolymer_comp"]["chem_comp"]["name"]
        # nonpoly["nonpolymer_comp"]["chem_comp"]["three_letter_code"]

        # nonpoly["nonpolymer_comp"]["drugbank"]["drugbank_container_identifiers"][
        #     "drugbank_id"
        # ]
        # nonpoly["nonpolymer_comp"]["drugbank"]["drugbank_info"]["cas_number"]
        # nonpoly["nonpolymer_comp"]["drugbank"]["drugbank_info"]["description"]

        # nonpoly["nonpolymer_comp"]["rcsb_chem_comp_target"]

        # {
        #     "interaction_type": "target",
        #     "name": "Spermine synthase",
        #     "provenance_source": "DrugBank",
        #     "reference_database_accession_code": "P52788",
        #     "reference_database_name": "UniProt",
        # },

        return NonpolymericLigand(
            nonpolymer_comp=nonpoly["nonpolymer_comp"],
            chemicalId=nonpoly["pdbx_entity_nonpoly"]["comp_id"],
            chemicalName=nonpoly["pdbx_entity_nonpoly"]["name"],
            pdbx_description=nonpoly["rcsb_nonpolymer_entity"]["pdbx_description"],
            formula_weight=nonpoly["rcsb_nonpolymer_entity"]["formula_weight"],
            number_of_instances=nonpoly["rcsb_nonpolymer_entity"][
                "pdbx_number_of_molecules"
            ],
        )

    def asm_parse(self, dictionaries: list[dict]) -> list[AssemblyInstancesMap]:
        return list(map(AssemblyInstancesMap.parse_obj, dictionaries))

    def poly_assign_to_asm(self, auth_asym_id: str) -> int:
        """
        Every structure in PDB might have 1 or more "assemblies", i.e. *almost* identical models sitting next to each other in space.
        The purpose of this method is to assign a polymer to the correct assembly, given its auth_asym_id.
        """

        if len(self.asm_maps) == 1:
            return 0
        else:
            for assembly_map in self.asm_maps:
                for polymer_instance in assembly_map.polymer_entity_instances:
                    if (
                        polymer_instance.rcsb_polymer_entity_instance_container_identifiers.auth_asym_id
                        == auth_asym_id
                    ):
                        return int(assembly_map.rcsb_id.split("-")[1]) - 1
            else:
                raise LookupError()

    def process_polypeptides(self) -> tuple[list[Protein], list[PolymericFactor]]:
        poly_entities = self.rcsb_data_dict["polymer_entities"]

        reshaped_proteins         : list[Protein]         = []
        reshaped_polymeric_factors: list[PolymericFactor] = []

        def is_protein(poly):
            #* According to RCSB schema, polymer_entites include Proteins, RNA but also DNA, NA-Hybrids and "Other". 
            #* We only make the distinction between Proteins and RNA and Other for purposes of simplicity
            return poly["entity_poly"]["rcsb_entity_polymer_type"] == "Protein"

        for poly in poly_entities:
            # A polymer can be either rna or protein
            if is_protein(poly):
                # a protein can be either a ribosomal protein, a factor or a nascent chain
                if (factor_classify(poly["rcsb_polymer_entity"]["pdbx_description"])!= None):

                    # TODO: MOVE TO HMM BASED CLASSIFICATION METHOD
                    reshaped_polymeric_factors.extend(self.poly_reshape_to_rFactor(poly))

                else:
                    # TODO: MOVE TO HMM BASED CLASSIFICATION METHOD
                    reshaped_proteins.extend(self.poly_reshape_to_rProtein(poly))
            else:
                # print("Filtered out a protein")
                ...

        self.rProteins = reshaped_proteins
        return (reshaped_proteins, reshaped_polymeric_factors)

    def process_polynucleotides(self) -> tuple[list[RNA], list[PolymericFactor]]:
        poly_entities = self.rcsb_data_dict["polymer_entities"]
        rnas          = []


        def is_rna(poly):
            #* According to RCSB schema, polymer_entites include Proteins, RNA but also DNA, NA-Hybrids and "Other". 
            #* We only make the distinction between Proteins and RNA and Other for purposes of simplicity
            return poly["entity_poly"]["rcsb_entity_polymer_type"] == "RNA"

        for poly in poly_entities:
            rnas.append(poly) if is_rna(poly) else ...

        print("==================================================================================================")

        reshaped_rnas             : list[RNA]             = []
        reshaped_polymeric_factors: list[PolymericFactor] = []

        for j, poly_rna in enumerate(rnas):
            if (
                factor_classify(poly_rna["rcsb_polymer_entity"]["pdbx_description"])
                != None
            ):
                # TODO: HMM WORKFLOW
                reshaped_polymeric_factors.extend(
                    self.poly_reshape_to_rFactor(poly_rna)
                )
            else:
                # TODO: DIFFERENTIATE BETWEEN MRNA (regex) AND TRNA (HMM workflow )
                reshaped_rnas.extend(self.poly_reshape_to_rRNA(poly_rna))

        self.rRNA = reshaped_rnas
        return (reshaped_rnas, reshaped_polymeric_factors)

    def process_other_polymers(self)->list[Polymer]:

        poly_entities = self.rcsb_data_dict["polymer_entities"]
        other = []

        # print("Processing polypeptides :")
        # print(len(poly_entities))

        def is_not_rna_protein_polymer(poly):
            #* According to RCSB schema, polymer_entites include Proteins, RNA but also DNA, NA-Hybrids and "Other". 
            #* We only make the distinction between Proteins and RNA and Other for purposes of simplicity
            return poly["entity_poly"]["rcsb_entity_polymer_type"] not in ["RNA", "Protein"]

        for poly in poly_entities:
            other.append(poly) if is_not_rna_protein_polymer(poly) else ...

        flat_other = []        
        for poly in other:
            flat_other = [*flat_other, *self.poly_reshape_to_Other(poly)]

        return flat_other

    def process_nonpolymers(self) -> list[NonpolymericLigand]:
        nonpoly_entities = self.rcsb_data_dict["nonpolymer_entities"]
        reshaped_nonpoly = (
            [self.nonpoly_reshape_to_ligand(nonpoly) for nonpoly in nonpoly_entities]
            if nonpoly_entities != None and len(nonpoly_entities) > 0
            else []
        )

        return reshaped_nonpoly

    def process_metadata(self):
        organisms = self.infer_organisms_from_polymers([*self.rProteins, *self.rRNA])
        externalRefs = self.extract_external_refs(
            self.rcsb_data_dict["rcsb_external_references"]
        )

        if (
            self.rcsb_data_dict["citation"] != None
            and len(self.rcsb_data_dict["citation"]) > 0
        ):
            pub = self.rcsb_data_dict["citation"][0]
        else:
            pub = {
                "year": None,
                "rcsb_authors": None,
                "title": None,
                "pdbx_database_id_DOI": None,
                "pdbx_database_id_PubMed": None,
            }

        kwords_text = (
            self.rcsb_data_dict["struct_keywords"]["text"]
            if self.rcsb_data_dict["struct_keywords"] != None
            else None
        )
        kwords = (
            self.rcsb_data_dict["struct_keywords"]["pdbx_keywords"]
            if self.rcsb_data_dict["struct_keywords"] != None
            else None
        )

        return [organisms, externalRefs, pub, kwords_text, kwords]

    def poly_reshape_to_rProtein(self, rpotein_polymer_obj) -> list[Protein]:
        if rpotein_polymer_obj["pfams"] != None and len(rpotein_polymer_obj["pfams"]) > 0:
            pfam_comments = list(
                set([pfam["rcsb_pfam_comment"] for pfam in rpotein_polymer_obj["pfams"]])
            )
            pfam_descriptions = list(
                set([pfam["rcsb_pfam_description"] for pfam in rpotein_polymer_obj["pfams"]])
            )
            pfam_accessions = list(
                set([pfam["rcsb_pfam_accession"] for pfam in rpotein_polymer_obj["pfams"]])
            )

        else:
            pfam_comments     = []
            pfam_descriptions = []
            pfam_accessions   = []

        host_organisms: list[Any] | None = rpotein_polymer_obj["rcsb_entity_host_organism"]
        source_organisms: list[Any] | None = rpotein_polymer_obj["rcsb_entity_source_organism"]

        host_organism_ids = []
        host_organism_names = []

        src_organism_ids = []
        src_organism_names = []

        if host_organisms != None:
            for ho in host_organisms:
                if ho["ncbi_taxonomy_id"] != None:
                    host_organism_ids.append(ho["ncbi_taxonomy_id"])
                if ho["scientific_name"] != None:
                    host_organism_names.append(ho["scientific_name"])

        if source_organisms != None:
            for so in source_organisms:
                if so["ncbi_taxonomy_id"] != None:
                    src_organism_ids.append(so["ncbi_taxonomy_id"])
                if so["scientific_name"] != None:
                    src_organism_names.append(so["scientific_name"])

        host_organism_ids   = list(map(int, set(host_organism_ids)))
        host_organism_names = list(map(str, set(host_organism_names)))
        src_organism_ids    = list(map(int, set(src_organism_ids)))
        src_organism_names  = list(map(str, set(src_organism_names)))

        # ? Compare prot sequence against all HMMs (returns a dict), pick the class with the lowest e-value
        hmm_resulsts = seq_prot_against_protclasses(
            rpotein_polymer_obj["entity_poly"]["pdbx_seq_one_letter_code_can"],
            self.hmm_ribosomal_proteins,
        )
        results_hmm = dict(
            sorted(
                {key: value for key, value in hmm_resulsts.items() if value}.items(),
                key=lambda item: min(item[1]),
            )
        )
        nomenclature = (
            [list(results_hmm.keys())[0]] if len(list(results_hmm.keys())) > 0 else []
        )

        return [
            Protein(
                assembly_id=self.poly_assign_to_asm(auth_asym_id),
                nomenclature=nomenclature,
                asym_ids=rpotein_polymer_obj["rcsb_polymer_entity_container_identifiers"]["asym_ids"],
                parent_rcsb_id=rpotein_polymer_obj["entry"]["rcsb_id"],
                auth_asym_id=auth_asym_id,
                pfam_accessions=pfam_accessions,
                pfam_comments=pfam_comments,
                pfam_descriptions=pfam_descriptions,
                host_organism_ids=host_organism_ids,
                host_organism_names=host_organism_names,
                src_organism_ids=src_organism_ids,
                src_organism_names=src_organism_names,
                uniprot_accession=[entry["rcsb_id"] for entry in rpotein_polymer_obj["uniprots"]]
                if rpotein_polymer_obj["uniprots"] != None and len(rpotein_polymer_obj["uniprots"]) > 0
                else [],
                rcsb_pdbx_description=rpotein_polymer_obj["rcsb_polymer_entity"]["pdbx_description"],
                entity_poly_strand_id=rpotein_polymer_obj["entity_poly"]["pdbx_strand_id"],
                entity_poly_seq_one_letter_code=rpotein_polymer_obj["entity_poly"][
                    "pdbx_seq_one_letter_code"
                ],
                entity_poly_seq_one_letter_code_can=rpotein_polymer_obj["entity_poly"][
                    "pdbx_seq_one_letter_code_can"
                ],
                entity_poly_seq_length=rpotein_polymer_obj["entity_poly"][
                    "rcsb_sample_sequence_length"
                ],
                entity_poly_entity_type=rpotein_polymer_obj["entity_poly"]["type"],
                entity_poly_polymer_type=rpotein_polymer_obj["entity_poly"]["rcsb_entity_polymer_type"],
            )
            for auth_asym_id in rpotein_polymer_obj["rcsb_polymer_entity_container_identifiers"][
                "auth_asym_ids"
            ]
        ]

    def poly_reshape_to_rRNA(self, rrna_polymer_obj) -> list[RNA]:
        """this returns a list because certain polymers accounts for multiple RNA molecules"""

        host_organisms: list[Any] | None = rrna_polymer_obj["rcsb_entity_host_organism"]
        source_organisms: list[Any] | None = rrna_polymer_obj["rcsb_entity_source_organism"]

        host_organism_ids = []
        host_organism_names = []
        src_organism_ids = []
        src_organism_names = []

        if host_organisms != None:
            for ho in host_organisms:
                if ho["ncbi_taxonomy_id"] != None:
                    host_organism_ids.append(ho["ncbi_taxonomy_id"])
                if ho["scientific_name"] != None:
                    host_organism_names.append(ho["scientific_name"])

        if source_organisms != None:
            for so in source_organisms:
                if so["ncbi_taxonomy_id"] != None:
                    src_organism_ids.append(so["ncbi_taxonomy_id"])
                if so["scientific_name"] != None:
                    src_organism_names.append(so["scientific_name"])

        host_organism_ids = list(map(int, set(host_organism_ids)))
        host_organism_names = list(map(str, set(host_organism_names)))
        src_organism_ids = list(map(int, set(src_organism_ids)))
        src_organism_names = list(map(str, set(src_organism_names)))

        # # ------------
        # host_organism_ids   = list(map(int, set([org['ncbi_taxonomy_id'] for org in plm['rcsb_entity_host_organism'  ]]))) if plm['rcsb_entity_host_organism'  ] != None else []
        # host_organism_names = list(map(str, set([org['scientific_name'] for org in plm['rcsb_entity_host_organism'  ]]))) if plm['rcsb_entity_host_organism'  ] != None else []

        # src_organism_ids   = list(map(int, set([org['ncbi_taxonomy_id'] for org in plm['rcsb_entity_source_organism']]))) if plm['rcsb_entity_source_organism'] != None else []
        # src_organism_names = list(map(str, set([org['scientific_name'] for org in plm['rcsb_entity_source_organism']]))) if plm['rcsb_entity_source_organism']  != None else []

        nomenclature = rna_classify(rrna_polymer_obj["rcsb_polymer_entity"]["pdbx_description"])

        return [
            RNA(
                assembly_id=self.poly_assign_to_asm(auth_asym_id),
                nomenclature=nomenclature,
                asym_ids=rrna_polymer_obj["rcsb_polymer_entity_container_identifiers"]["asym_ids"],
                auth_asym_id=auth_asym_id,
                parent_rcsb_id=rrna_polymer_obj["entry"]["rcsb_id"],
                host_organism_ids=host_organism_ids,
                host_organism_names=host_organism_names,
                src_organism_ids=src_organism_ids,
                src_organism_names=src_organism_names,
                rcsb_pdbx_description=""
                if rrna_polymer_obj["rcsb_polymer_entity"]["pdbx_description"] == None
                else rrna_polymer_obj["rcsb_polymer_entity"]["pdbx_description"],
                entity_poly_strand_id=rrna_polymer_obj["entity_poly"]["pdbx_strand_id"],
                entity_poly_seq_one_letter_code=rrna_polymer_obj["entity_poly"][
                    "pdbx_seq_one_letter_code"
                ],
                entity_poly_seq_one_letter_code_can=rrna_polymer_obj["entity_poly"][
                    "pdbx_seq_one_letter_code_can"
                ],
                entity_poly_seq_length=rrna_polymer_obj["entity_poly"][
                    "rcsb_sample_sequence_length"
                ],
                entity_poly_entity_type=rrna_polymer_obj["entity_poly"]["type"],
                entity_poly_polymer_type=rrna_polymer_obj["entity_poly"]["rcsb_entity_polymer_type"],
            )
            for auth_asym_id in rrna_polymer_obj["rcsb_polymer_entity_container_identifiers"][
                "auth_asym_ids"
            ]
        ]

    def poly_reshape_to_rFactor(self, factor_polymer_obj) -> list[PolymericFactor]:
        host_organisms: list[Any] | None = factor_polymer_obj["rcsb_entity_host_organism"]
        source_organisms: list[Any] | None = factor_polymer_obj["rcsb_entity_source_organism"]

        host_organism_ids = []
        host_organism_names = []
        src_organism_ids = []
        src_organism_names = []

        if host_organisms != None:
            for ho in host_organisms:
                if ho["ncbi_taxonomy_id"] != None:
                    host_organism_ids.append(ho["ncbi_taxonomy_id"])
                if ho["scientific_name"] != None:
                    host_organism_names.append(ho["scientific_name"])

        if source_organisms != None:
            for so in source_organisms:
                if so["ncbi_taxonomy_id"] != None:
                    src_organism_ids.append(so["ncbi_taxonomy_id"])
                if so["scientific_name"] != None:
                    src_organism_names.append(so["scientific_name"])

        host_organism_ids   = list(map(int, set(host_organism_ids)))
        host_organism_names = list(map(str, set(host_organism_names)))
        src_organism_ids    = list(map(int, set(src_organism_ids)))
        src_organism_names  = list(map(str, set(src_organism_names)))

        nomenclature = [
            *filter(
                lambda x: x is not None,
                [factor_classify(factor_polymer_obj["rcsb_polymer_entity"]["pdbx_description"])],
            )
        ]

        return [
            PolymericFactor(
                assembly_id=self.poly_assign_to_asm(auth_asym_id),
                nomenclature=nomenclature,
                asym_ids=factor_polymer_obj["rcsb_polymer_entity_container_identifiers"]["asym_ids"],
                auth_asym_id=auth_asym_id,
                parent_rcsb_id=factor_polymer_obj["entry"]["rcsb_id"],
                host_organism_ids=host_organism_ids,
                host_organism_names=host_organism_names,
                src_organism_ids=src_organism_ids,
                src_organism_names=src_organism_names,
                rcsb_pdbx_description=""
                if factor_polymer_obj["rcsb_polymer_entity"]["pdbx_description"] == None
                else factor_polymer_obj["rcsb_polymer_entity"]["pdbx_description"],
                entity_poly_strand_id=factor_polymer_obj["entity_poly"]["pdbx_strand_id"],
                entity_poly_seq_one_letter_code=factor_polymer_obj["entity_poly"][
                    "pdbx_seq_one_letter_code"
                ],
                entity_poly_seq_one_letter_code_can=factor_polymer_obj["entity_poly"][
                    "pdbx_seq_one_letter_code_can"
                ],
                entity_poly_seq_length=factor_polymer_obj["entity_poly"][
                    "rcsb_sample_sequence_length"
                ],
                entity_poly_entity_type=factor_polymer_obj["entity_poly"]["type"],
                entity_poly_polymer_type=factor_polymer_obj["entity_poly"]["rcsb_entity_polymer_type"],
            )
            for auth_asym_id in factor_polymer_obj["rcsb_polymer_entity_container_identifiers"][
                "auth_asym_ids"
            ]
        ]

    def poly_reshape_to_Other(self, other_polymer_obj) -> list[Polymer]:

        host_organisms  : Optional[ list[Any] ] = other_polymer_obj["rcsb_entity_host_organism"]
        source_organisms: Optional[ list[Any] ] = other_polymer_obj["rcsb_entity_source_organism"]

        host_organism_ids   = []
        host_organism_names = []
        src_organism_ids    = []
        src_organism_names  = []

        if host_organisms != None:
            for ho in host_organisms:
                if ho["ncbi_taxonomy_id"] != None:
                    host_organism_ids.append(ho["ncbi_taxonomy_id"])
                if ho["scientific_name"] != None:
                    host_organism_names.append(ho["scientific_name"])

        if source_organisms != None:
            for so in source_organisms:
                if so["ncbi_taxonomy_id"] != None:
                    src_organism_ids.append(so["ncbi_taxonomy_id"])
                if so["scientific_name"] != None:
                    src_organism_names.append(so["scientific_name"])

        host_organism_ids   = list(map(int, set(host_organism_ids)))
        host_organism_names = list(map(str, set(host_organism_names)))
        src_organism_ids    = list(map(int, set(src_organism_ids)))
        src_organism_names  = list(map(str, set(src_organism_names)))


        return [
            Polymer(
                assembly_id=self.poly_assign_to_asm(auth_asym_id),
                nomenclature=[], # Nomenclature does not exist for arbitrary polymers
                asym_ids=other_polymer_obj["rcsb_polymer_entity_container_identifiers"]["asym_ids"],
                auth_asym_id=auth_asym_id,
                parent_rcsb_id=other_polymer_obj["entry"]["rcsb_id"],
                host_organism_ids=host_organism_ids,
                host_organism_names=host_organism_names,
                src_organism_ids=src_organism_ids,
                src_organism_names=src_organism_names,
                rcsb_pdbx_description=""
                if other_polymer_obj["rcsb_polymer_entity"]["pdbx_description"] == None
                else other_polymer_obj["rcsb_polymer_entity"]["pdbx_description"],
                entity_poly_strand_id=other_polymer_obj["entity_poly"]["pdbx_strand_id"],
                entity_poly_seq_one_letter_code=other_polymer_obj["entity_poly"][
                    "pdbx_seq_one_letter_code"
                ],
                entity_poly_seq_one_letter_code_can=other_polymer_obj["entity_poly"][
                    "pdbx_seq_one_letter_code_can"
                ],
                entity_poly_seq_length=other_polymer_obj["entity_poly"][
                    "rcsb_sample_sequence_length"
                ],
                entity_poly_entity_type=other_polymer_obj["entity_poly"]["type"],
                entity_poly_polymer_type=other_polymer_obj["entity_poly"]["rcsb_entity_polymer_type"],
            )
            for auth_asym_id in other_polymer_obj["rcsb_polymer_entity_container_identifiers"][
                "auth_asym_ids"
            ]
        ]

    def process_structure(self):
        [reshaped_proteins,reshaped_polymeric_factors_prot] = self.process_polypeptides()
        [reshaped_rnas, reshaped_polymeric_factors_rna] = self.process_polynucleotides()
        other_polymers = self.process_other_polymers()

        # for i in reshaped_rnas:
        #     pprint(i)
        # for i in reshaped_proteins:
        #     pprint(i)
        # for i in reshaped_polymeric_factors_prot:
        #     pprint(i)
        # for i in reshaped_polymeric_factors_rna:
        #     pprint(i)

        # print("target falt poly", self.flattened_polymers_target)
        # print(len(reshaped_proteins))
        # print(len(reshaped_rnas))
        # print(len(reshaped_polymeric_factors_rna))
        # print(len(reshaped_polymeric_factors_prot))

        # print("poly identifiers len: ",len( self.rcsb_data_dict["polymer_entities"] ), )

        assert (
            len(reshaped_proteins)
            + len(reshaped_rnas)
            + len(reshaped_polymeric_factors_rna)
            + len(reshaped_polymeric_factors_prot)
            + len(other_polymers)
        ) == self.flattened_polymers_target

        reshaped_nonpolymers = self.process_nonpolymers()
        [organisms, externalRefs, pub, kwords_text, kwords] = self.process_metadata()

        reshaped = RibosomeStructure(
            rcsb_id=self.rcsb_data_dict["rcsb_id"],
            expMethod=self.rcsb_data_dict["exptl"][0]["method"],
            resolution=self.rcsb_data_dict["rcsb_entry_info"]["resolution_combined"][0],
            rcsb_external_ref_id=externalRefs[0],
            rcsb_external_ref_type=externalRefs[1],
            rcsb_external_ref_link=externalRefs[2],
            citation_year=pub["year"],
            citation_rcsb_authors=pub["rcsb_authors"],
            citation_title=pub["title"],
            citation_pdbx_doi=pub["pdbx_database_id_DOI"],
            pdbx_keywords_text=kwords_text,
            pdbx_keywords=kwords,
            src_organism_ids=organisms["src_organism_ids"],
            src_organism_names=organisms["src_organism_names"],
            host_organism_ids=organisms["host_organism_ids"],
            host_organism_names=organisms["host_organism_names"],
            proteins=reshaped_proteins,
            rnas=reshaped_rnas,
            polymeric_factors=[
                *reshaped_polymeric_factors_prot,
                *reshaped_polymeric_factors_rna,
            ],
            nonpolymeric_ligands=reshaped_nonpolymers,
            other_polymers=other_polymers,
            assembly_map=self.asm_maps,
        )

        return reshaped
