import asyncio
from enum import   auto
import enum
import json
import os
from pprint import pprint
import subprocess
import typing
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from loguru import logger
import requests
from ribctl import AMINO_ACIDS_3_TO_1_CODE, ASSETS_PATH, CHAINSPLITTER_PATH, CLASSIFICATION_REPORTS
from ribctl.etl.assets_structure import AssetClass, StructureAssets
from ribctl.lib.libtax import PhylogenyNode, PhylogenyRank, Taxid
from Bio.PDB.Structure import Structure
from Bio.PDB.MMCIFParser import FastMMCIFParser
from ribctl.lib.landmarks.ptc_via_doris import ptc_resdiues_get, ptc_residues_calculate_midpoint
from ribctl.lib.utils import download_unpack_place
from ribctl.lib.schema.types_ribosome import ( RNA, PTCInfo, Polymer, PolymerClass, PolynucleotideClass, PolynucleotideClass, PolypeptideClass, RibosomeStructure, RibosomeStructureMetadata, )
from ribctl import RIBETL_DATA
from ribctl.logs.loggers import get_etl_logger
from ribctl.ribosome_ops import RibosomeOps




class GlobalAssets:

    @staticmethod
    def status_vs_rcsb() -> list[str]:
        """Return a list of structures that are in the RCSB but not in the local database."""
        return list(set(GlobalAssets.current_rcsb_structs()) - set(GlobalAssets.list_all_structs()))

    @staticmethod
    def current_rcsb_structs() -> list[str]:
        """Return all structures in the rcsb that contain the phrase RIBOSOME and have more than 25 protein entities"""

        rcsb_search_api = "https://search.rcsb.org/rcsbsearch/v2/query"
        
        q2 = {
            "query": {
                "type"            : "group",
                "logical_operator": "and",
                "nodes"           : [
                    {
                        "type"      : "terminal",
                        "service"   : "text",
                        "parameters": {
                            "operator" : "contains_phrase",
                            "negation" : False,
                            "value"    : "RIBOSOME",
                            "attribute": "struct_keywords.pdbx_keywords",
                        },
                    },
                    {
                        "type"      : "terminal",
                        "service"   : "text",
                        "parameters": {
                            "operator" : "greater",
                            "negation" : False,
                            "value"    : 12,
                            "attribute": "rcsb_entry_info.polymer_entity_count_protein",
                        },
                    },
                ],
                "label": "query-builder",
            },
            "return_type": "entry",
            "request_options": {"return_all_hits": True, "results_verbosity": "compact"},
        }
        query = rcsb_search_api + "?json=" + json.dumps(q2)
        return sorted(requests.get(query).json()["result_set"])
    @staticmethod
    def list_all_structs()->list[str]:
        
        profiles_exist = [ rcsb_id if os.path.exists(RibosomeOps(rcsb_id).assets.paths.profile)  else None for rcsb_id in os.listdir(RIBETL_DATA)]
        return list(filter(lambda x: x!=None, profiles_exist))

    @staticmethod
    def ptc_references( ribosome_type:typing.Literal['arch','bact','euk','mito']):
        filename = "ptc_reference_residues_{}.pickle".format(ribosome_type.upper())
        return os.path.join(ASSETS_PATH, "landmarks_cache", filename )

    @staticmethod
    def collect_all_taxa() -> set[PhylogenyNode]:
        _ = set()
        for struct in GlobalAssets.list_all_structs():
            rp = RibosomeOps(struct).profile
            for org in [*rp.src_organism_ids, *rp.host_organism_ids]:
                # if Taxid.rank(org) not in list(typing.get_args(PhylogenyRank)):
                #     org = Taxid.coerce_to_rank(org, "species")
                assert org is not None
                try:
                    pn = PhylogenyNode(
                        ncbi_tax_id     = org,
                        scientific_name = Taxid.get_name(org),
                        rank            = Taxid.rank(org),
                    )
                except Exception as e:
                    print( struct, Taxid.get_name(org), "|\t", Taxid.rank(org), "->", Taxid.get_lineage(org), )
                    print("Error with", org, struct)
                    print(e)
                _.add(pn)
        return _

    @staticmethod
    def global_status()->dict[str, dict[AssetClass, bool]]:
        _ = {}
        for struct in GlobalAssets.list_all_structs():
            _[struct] = StructureAssets(struct).assets_status()
        return _