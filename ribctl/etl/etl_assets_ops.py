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
from ribctl import AMINO_ACIDS_3_TO_1_CODE, CHAINSPLITTER_PATH, CLASSIFICATION_REPORTS
from ribctl.lib.libtax import PhylogenyNode, PhylogenyRank, Taxid
from Bio.PDB.Structure import Structure
from Bio.PDB.MMCIFParser import FastMMCIFParser
from ribctl.lib.tunnel import ptc_resdiues_get, ptc_residues_calculate_midpoint
from ribctl.lib.utils import download_unpack_place
from ribctl.lib.schema.types_ribosome import ( RNA, PTCInfo, Polymer, PolymerClass, PolynucleotideClass, PolynucleotideClass, PolypeptideClass, RibosomeStructure, )
from ribctl import RIBETL_DATA
from ribctl.logs.loggers import get_etl_logger


RCSB_ID = typing.NewType("RCSB_ID", str)


class AssetClass(enum.StrEnum):
    profile   = auto()
    cif       = auto()
    ptc       = auto()
    # chains    = auto()
    thumbnail = auto()

    @staticmethod
    def from_str(_:str):
        assets = [status for status in dir( AssetClass) if not status.startswith('_')]
        if _ in assets:
            return getattr(AssetClass, _)
        return None

class AssetPath:
    rcsb_id:str
    def __init__(self, rcsb_id) -> None:
        self.rcsb_id = rcsb_id
        pass

    def binding_site(self, chemId:str):
        return f"{self.dir}/{self.rcsb_id.upper()}_LIG_{chemId.upper()}.json"

    def binding_site_prediction(self, chemId:str, source_struct:str):
        return f"{self.dir}/{self.rcsb_id.upper()}_LIG_{chemId.upper()}_PREDICTION_VIA_{source_struct.upper()}.json"
    
    @property
    def dir(self):
        return os.path.join(RIBETL_DATA, self.rcsb_id)

    @property
    def cif(self):
        return f"{self.dir}/{self.rcsb_id}.cif"

    @property
    def ptc(self):
        return os.path.join( self.dir, "{}_PTC_COORDINATES.json".format(self.rcsb_id) )

    @property
    def profile(self):
        return os.path.join(self.dir, f"{self.rcsb_id}.json")

    # @property
    # def chains_dir(self):
    #     return f"{self.dir}/CHAINS"

    @property
    def classification_report(self):
        return os.path.join(self.dir, f"classification_report_{self.rcsb_id}.json")

    @property
    def thumbnail(self):
        return f"{self.dir}/{self.rcsb_id}.png"

    @property
    def tunnel_dir(self):
        return os.path.join(self.dir, "TUNNELS")

class RibosomeOps:
    rcsb_id: str
    paths: AssetPath

    def __init__(self, rcsb_id: str) -> None:
        if not RIBETL_DATA:
            raise Exception("RIBETL_DATA environment variable not set. Cannot access assets." )
        self.rcsb_id = rcsb_id.upper()
        self.paths   = AssetPath(self.rcsb_id)


    def nomenclature_table(self, verbose: bool = False) -> dict[str, dict]:
        prof = self.profile()
        m    = {}

        for p in prof.other_polymers:
            m[p.auth_asym_id] = {
                "nomenclature": list(map(lambda x: x.value, p.nomenclature)),
            }

            if verbose:
                m[p.auth_asym_id].update(
                    {
                        "entity_poly_strand_id": p.entity_poly_strand_id,
                        "rcsb_pdbx_description": p.rcsb_pdbx_description,
                    }
                )
        for prot in prof.proteins:
            m[prot.auth_asym_id] = {
                "nomenclature": list(map(lambda x: x.value, prot.nomenclature)),
            }
            if verbose:
                m[prot.auth_asym_id].update(
                    {
                        "entity_poly_strand_id": prot.entity_poly_strand_id,
                        "rcsb_pdbx_description": prot.rcsb_pdbx_description,
                    }
                )
        if prof.rnas != None:
            for rna in prof.rnas:
                m[rna.auth_asym_id] = {
                    "nomenclature": list(map(lambda x: x.value, rna.nomenclature)),
                }
                if verbose:
                    m[rna.auth_asym_id].update(
                        {
                            "entity_poly_strand_id": rna.entity_poly_strand_id,
                            "rcsb_pdbx_description": rna.rcsb_pdbx_description,
                        }
                    )

        return m

    def ptc(self) -> PTCInfo:
        with open(self.paths.ptc, "r") as infile:
            _ = json.load(infile)
            return PTCInfo.model_validate(_)


    def profile(self) -> RibosomeStructure:
        with open(self.paths.profile, "r") as f:
            return RibosomeStructure.model_validate(json.load(f))

    def biopython_structure(self)-> Structure:
        cifpath = RibosomeOps(self.rcsb_id).paths.cif
        return FastMMCIFParser(QUIET=True).get_structure(self.rcsb_id, cifpath)

    def write_own_json_profile(self, new_profile: dict, overwrite: bool = False):
        """Update self, basically."""
        if os.path.exists(self.paths.profile) and not overwrite:
            print( "You are about to overwrite {}. Specify `overwrite=True` explicitly.".format( self.paths.profile )
            )
        elif overwrite:
            with open(self.paths.profile, "w") as f:
                json.dump(new_profile, f)
                logger.debug(f"Updated profile for {self.rcsb_id}")

    def get_taxids(self) -> tuple[list[int], list[int]]:
        p = self.profile()
        return (p.src_organism_ids, p.host_organism_ids)

    def get_chain_by_polymer_class(
        self, poly_class: PolynucleotideClass, assembly: int = 0
    ) -> Polymer | None:

        profile = self.profile()

        if poly_class in [v.value for v in PolynucleotideClass]:
            if profile.rnas is not None:
                for rna in profile.rnas:
                    if (
                        poly_class in [r.value for r in rna.nomenclature]
                        and rna.assembly_id == assembly
                    ):
                        return rna
            else:
                return None

        elif poly_class in [v.value for v in PolypeptideClass]:
            for prot in profile.proteins:
                if (
                    poly_class in [v.value for v in prot.nomenclature]
                    and prot.assembly_id == assembly
                ):
                    return prot
            else:
                return None

        else:
            for polyf in profile.other_polymers:
                if (
                    poly_class in [p.value for p in polyf.nomenclature]
                    and polyf.assembly_id == assembly
                ):
                    return polyf
        return None

    def get_poly_by_auth_asym_id( self, auth_asym_id: str ) -> Polymer :

        profile = self.profile()
        for chain in [ *profile.proteins, *profile.rnas, *profile.other_polymers]:
            if chain.auth_asym_id == auth_asym_id:
                return chain
        raise KeyError("No chain found with auth_asym_id: {}".format(auth_asym_id))

    def get_poly_by_polyclass(
        self, class_: PolymerClass, assembly: int = 0
    ) ->Polymer | None:
        """@assembly here stands to specify which of the two or more models the rna comes from
        in the case that a structure contains multiple models (ex. 4V4Q XRAY)"""

        profile = self.profile()

        # if profile.rnas == None:
        #     return None

        # for rna in profile.rnas:
        #     if class_ in rna.nomenclature and rna.assembly_id == assembly:
        #         return rna

        for polymer in [*profile.rnas, *profile.other_polymers, *profile.proteins]: 
            if class_ in  polymer.nomenclature and polymer.assembly_id == assembly:
                return polymer
        

    def get_LSU_rRNA(self, assembly: int = 0) -> RNA:
        """retrieve the largest rRNA sequence in the structure
        @returns (seq, auth_asym_id, rna_type)
        """

        rna = self.get_poly_by_polyclass(PolynucleotideClass("23SrRNA"), assembly)
        if rna == None:
            rna = self.get_poly_by_polyclass(PolynucleotideClass("25SrRNA"), assembly)
        if rna == None:
            rna = self.get_poly_by_polyclass(PolynucleotideClass("28SrRNA"), assembly)
        if rna == None:
            rna = self.get_poly_by_polyclass(PolynucleotideClass("mt16SrRNA"), assembly)
        if rna == None:
            raise Exception("No LSU rRNA found in structure")
        else:
            return rna

    def biopython_get_chain(self, auth_asym_id: str) -> Chain:
        return self.biopython_structure().child_dict[0].child_dict[auth_asym_id]

    #TODO : Utilities, various
    @staticmethod
    def biopython_chain_get_seq(
        struct: Structure,
        auth_asym_id: str,
        protein_rna: typing.Literal["protein", "rna"],
        sanitized: bool = False,
    ) -> str:

        chain3d = struct.child_dict[0].child_dict[auth_asym_id]
        ress = chain3d.child_list

        seq = ""
        for i in ress:
            if i.resname not in AMINO_ACIDS_3_TO_1_CODE.keys():
                print("Unknown residue:", i.resname)
                continue

            if protein_rna == "rna":
                seq += i.resname
            else:
                seq += AMINO_ACIDS_3_TO_1_CODE[i.resname]

        return seq

    # ※ -=-=-=-=-=-=-=-=-=-=-=-=-=-=-= Verification =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    def _verify_dir_exists(self):
        if not os.path.exists(self.paths.dir):
            os.umask(0)
            os.makedirs(self.paths.dir, 0o755)

class Assets:

    rcsb_id: str
    ro     : RibosomeOps
    paths  : AssetPath
    
    def __init__(self, rcsb_id: str) -> None:
        self.rcsb_id = rcsb_id
        self.ro      = RibosomeOps(rcsb_id)
        self.paths   = AssetPath(rcsb_id)
    
    
    @staticmethod
    def global_status()->dict[RCSB_ID, dict[AssetClass, bool]]:
        _ = {}
        for struct in Assets.list_all_structs():
            _[struct] = Assets(struct).assets_status()
        return _

    @staticmethod
    def collect_all_taxa() -> set[PhylogenyNode]:
        _ = set()
        for struct in Assets.list_all_structs():
            rp = RibosomeOps(struct).profile()
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
    def list_all_structs()->list[str]:
        
        profiles_exist = [ rcsb_id if os.path.exists(RibosomeOps(rcsb_id).paths.profile)  else None for rcsb_id in os.listdir(RIBETL_DATA)]
        return list(filter(lambda x: x!=None, profiles_exist))


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
    def status_vs_rcsb() -> list[str]:
        """Return a list of structures that are in the RCSB but not in the local database."""
        return list(set(Assets.current_rcsb_structs()) - set(Assets.list_all_structs()))

    def assets_status(self)->dict[AssetClass, bool]:
        _ = {}
        for a in AssetClass:
            _[a]= self.verify_exists(a)
        return _

    def verify_exists(self, asset: AssetClass) -> bool:
        match asset:
            case AssetClass.profile:
                return os.path.exists(self.paths.profile)
            case AssetClass.cif:
                return os.path.exists(self.paths.cif)
            case AssetClass.ptc:
                return os.path.exists(self.paths.ptc)
            # case AssetClass.chains:
            #     return os.path.exists(self.paths.chains_dir) and len(os.listdir(self.paths.chains_dir))> 0
            case AssetClass.thumbnail:
                return os.path.exists(self.paths.thumbnail)
            case _:
                raise KeyError("AssetClass {} does not exist.".format(asset))

    async def upsert_cif(self, overwrite: bool = False):
        if not os.path.exists(self.paths.cif):
            await download_unpack_place(self.rcsb_id)
            print("Saved cif file: \t", self.paths.cif)
        else:
            if overwrite:
                await download_unpack_place(self.rcsb_id)
                print("Overwrote cif file: \t", self.paths.cif)

    #TODO: ChimeraX image gen
    def upsert_thumbnail(self, overwrite: bool = False) -> bool:
        ...

    #TODO: ChimeraX splitchain
    async def upsert_chains(self):
        try:
            #TODO: chimerax should be a env variable pointint to the binary, so too should the script
            # Command to run ChimeraX with the script
            command = ['chimerax',   '--nogui', '--offscreen' ,'--cmd' ,'open /home/rtviii/dev/riboxyz/chimerax/chainsplitter.py; open {}; chainsplitter {}; exit'.format(self.rcsb_id, self.rcsb_id)]
            process = subprocess.Popen(
                command,
                stdout = subprocess.PIPE,
                stderr = subprocess.PIPE,
                text   = True
            )
            # Capture output in real-time
            while True:
                output = process.stdout.readline()
                if output == '' and process.poll() is not None:
                    break
                if output:
                    print("Output:", output.strip())

            # Capture any remaining output and errors
            stdout, stderr = process.communicate()

            if process.returncode == 0:
                print("ChimeraX script executed successfully.")
            else:
                print("ChimeraX script execution failed.")
                print("Error:", stderr)

            return process.returncode, stdout, stderr

        except Exception as e:
            print(f"An error occurred: {str(e)}")
            return None, None, str(e)

        
    async def upsert_ptc(self, overwrite: bool = False):

        etllogger = get_etl_logger()
        asset_ptc_coords_path = os.path.join(self.paths.ptc)

        if os.path.exists(asset_ptc_coords_path) and not overwrite:
            logger.debug(f"PTC coordinates already exist for {self.rcsb_id} and overwrite is set to False" )

        ress, auth_asym_id = ptc_resdiues_get( self.ro.biopython_structure(), self.ro.profile().rnas )
        midpoint_coords    = ptc_residues_calculate_midpoint(ress, auth_asym_id)

        writeout = {
            "site_9_residues"      : [(res.get_resname(), res.id[1]) for res in ress],
            "LSU_rRNA_auth_asym_id": auth_asym_id,
            "midpoint_coordinates" : midpoint_coords,
            "nomenclature_table"   : self.ro.nomenclature_table(),
        }

        with open(asset_ptc_coords_path, "w") as f:
            json.dump(writeout, f)
            etllogger.info( f"Saved PTC coordinates for {self.rcsb_id} to {asset_ptc_coords_path}" )