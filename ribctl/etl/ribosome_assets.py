import json
import os
from pprint import pprint
import typing
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from typing import Optional
from ribctl import AMINO_ACIDS_3_TO_1_CODE
from ribctl.lib.ribosome_types.types_binding_site import BindingSite
from ribctl.lib.mod_extract_bsites import  struct_ligand_ids, bsite_ligand
from ribctl.lib.mod_split_rename import split_rename
from ribctl.etl.etl_pipeline import current_rcsb_structs, ReannotationPipeline, rcsb_single_structure_graphql, query_rcsb_api
from ribctl.lib.tunnel import ptc_resdiues_get, ptc_residues_calculate_midpoint
# from ribctl.lib.mod_render_thumbnail import render_thumbnail
from ribctl.lib.utils import download_unpack_place, open_structure
from ribctl.lib.ribosome_types.types_ribosome import RNA, LifecycleFactorClass, PolymerClass, PolynucleotideClass, Protein, CytosolicProteinClass, RibosomeStructure
from ribctl import RIBETL_DATA
from pydantic import BaseModel
from concurrent.futures import ALL_COMPLETED, Future, ProcessPoolExecutor, ThreadPoolExecutor, wait
from ribctl.logs.loggers import get_etl_logger


logger = get_etl_logger()

class Assetlist(BaseModel)   : 
      profile                : Optional[bool]
      ptc_coords             : Optional[bool]
      cif                    : Optional[bool]
      cif_modified_and_chains: Optional[bool]
      ligands                : Optional[bool]
      png_thumbnail          : Optional[bool]

class RibosomeAssets():
    rcsb_id: str

    def __init__(self, rcsb_id: str) -> None:
        self.rcsb_id = rcsb_id.upper()

    def _envcheck(self):
        if not RIBETL_DATA:
            raise Exception(
                "RIBETL_DATA environment variable not set. Cannot access assets.")

    def _dir_path(self):
        self._envcheck()
        return os.path.join(RIBETL_DATA, self.rcsb_id)

    def _emdb_path(self):
        self._envcheck()
        emdb_path = os.path.join(RIBETL_DATA, self.rcsb_id.upper(), "{}.map.gz")
        return os.path.join(RIBETL_DATA, self.rcsb_id)

    def _cif_filepath(self):
        self._envcheck()
        return f"{self._dir_path()}/{self.rcsb_id}.cif"

    def _cif_modified_filepath(self):
        self._envcheck()
        return f"{self._dir_path()}/{self.rcsb_id}_modified.cif"

    def _ptc_residues(self) -> dict[str, dict[str, list[float]]]:
        PTC_RESIDUES_PATH = os.path.join(RIBETL_DATA, self.rcsb_id, "{}_PTC_COORDINATES.json".format(self.rcsb_id))
        with open(PTC_RESIDUES_PATH, 'r') as infile:
            return json.load(infile)

    def _json_profile_filepath(self):
        self._envcheck()
        return os.path.join(self._dir_path(),f"{self.rcsb_id}.json")

    def profile(self) -> RibosomeStructure:
        with open(self._json_profile_filepath(), "r") as f:
            return RibosomeStructure.model_validate(json.load(f))

    def _nomenclature_table(self) -> dict[str, dict]:
        #TODO: update getter
        prof = self.profile()
        m    = {}
        
        for p in prof.other_polymers:
            m[p.auth_asym_id]=  {
                   "nomenclature"         : list(map(lambda x: x.name, p.nomenclature)),
               "entity_poly_strand_id": p.entity_poly_strand_id,
               "rcsb_pdbx_description": p.rcsb_pdbx_description
               }

        for prot in prof.proteins:
            m[prot.auth_asym_id] = {
                   "nomenclature"         : list(map(lambda x: x.name, prot.nomenclature)),
                   "entity_poly_strand_id": prot.entity_poly_strand_id,
                   "rcsb_pdbx_description":prot.rcsb_pdbx_description
                   }
        if prof.rnas!=None:
            for rna in prof.rnas:
                m[rna.auth_asym_id] = {
                   "nomenclature"         : list(map(lambda x: x.name, rna.nomenclature)),
                   "entity_poly_strand_id": rna.entity_poly_strand_id,
                   "rcsb_pdbx_description": rna.rcsb_pdbx_description
                   }


        return m
                
    def biopython_structure(self):
        return open_structure(self.rcsb_id, 'cif')

    def chains_dir(self):
        self._envcheck()
        return f"{self._dir_path()}/CHAINS"

    def _png_thumbnail_filepath(self):
        self._envcheck()
        return f"{self._dir_path()}/_ray_{self.rcsb_id}.png"

    def write_own_json_profile(self, new_profile: dict, overwrite: bool = False):
        """Update self, basically."""
        if os.path.exists(self._json_profile_filepath()) and not overwrite:
            print("You are about to overwrite {}. Specify `overwrite=True` explicitly.".format(self._json_profile_filepath()))
        elif overwrite:
            with open(self._json_profile_filepath(), "w") as f:
                json.dump(new_profile, f)
                logger.debug(f"Updated profile for {self.rcsb_id}")

    @staticmethod
    def list_all_structs():
        return os.listdir(RIBETL_DATA)

    # ※ -=-=-=-=-=-=-=-=-=-=-=-=-=-=-= Getters =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    def get_taxids(self) -> tuple[list[int], list[int]]:
        p = self.profile()
        return (p.src_organism_ids, p.host_organism_ids)

    def get_struct_and_profile(self) -> tuple[Structure, RibosomeStructure]:
        return self.biopython_structure(), self.profile()

    def get_chain_by_polymer_class(self, poly_class: PolymerClass , assembly: int = 0) -> PolymerClass | None:
        profile = self.profile()
        # print("Got provfile for ", profile.rcsb_id)
        # print("searching for ", poly_class)
        for prot in profile.proteins:
            # print(prot.nomenclature)
            if poly_class in [v.value for v in prot.nomenclature] and prot.assembly_id == assembly:
                return prot

        if profile.rnas is not None:
            for rna in profile.rnas:
                if poly_class in [r.value for r in rna.nomenclature] and rna.assembly_id == assembly:
                    return rna

        if profile.other_polymers:
            for polyf in profile.other_polymers:
                if poly_class in [ p.value for p in polyf.nomenclature ] and polyf.assembly_id == assembly:
                    return polyf
        return None

    def get_chain_by_auth_asym_id(self, auth_asym_id: str) -> tuple[
            RNA | Protein  | None,
            typing.Literal["RNA", "Protein", "PolymericFactor"] | None]:

        profile = self.profile()
        for chain in profile.proteins:
            if chain.auth_asym_id == auth_asym_id:
                return (chain, "Protein")

        if profile.rnas is not None:
            for chain in profile.rnas:
                if chain.auth_asym_id == auth_asym_id:
                    return (chain, "RNA")

        #TODO: update getter
        if profile.polymeric_factors is not None:
            for chain in profile.polymeric_factors:
                if chain.auth_asym_id == auth_asym_id:
                    return (chain, "PolymericFactor")

        return (None, None)

    def get_rna_by_nomclass(self, class_: PolynucleotideClass, assembly: int = 0) -> RNA | None:
        """@assembly here stands to specify which of the two or more models the rna comes from
        in the case that a structure contains multiple models (ex. 4V4Q XRAY)"""

        profile = self.profile()

        if profile.rnas == None:
            return None

        for rna in profile.rnas:
            if class_ in rna.nomenclature and rna.assembly_id == assembly:
                return rna

    def get_prot_by_nomclass(self, class_: CytosolicProteinClass, assembly: int = 0) -> Protein | None:
        _auth_asym_id = {
            v: k for k, v in self.____nomenclature_v2().items()}.get(class_, None)
        if _auth_asym_id != None:
            return self.get_chain_by_auth_asym_id(_auth_asym_id)[0]

        for prot in self.profile().proteins:
            if class_ in prot.nomenclature and prot.assembly_id == assembly:
                return prot
        else:
            return None

    def get_LSU_rRNA(self, assembly: int = 0) -> RNA:
        """retrieve the largest rRNA sequence in the structure
        @returns (seq, auth_asym_id, rna_type)
        """

        rna = self.get_rna_by_nomclass(RNAClass("23SrRNA"), assembly);
        if rna == None:
            rna = self.get_rna_by_nomclass(RNAClass("25SrRNA"), assembly)
        if rna == None:
            rna = self.get_rna_by_nomclass(RNAClass("28SrRNA"), assembly)
        if rna == None:
            raise Exception("No LSU rRNA found in structure")
        else:
            return rna

    def biopython_get_chain(self, auth_asym_id: str) -> Chain:
        return self.biopython_structure().child_dict[0].child_dict[auth_asym_id]

    @staticmethod
    def biopython_chain_get_seq(struct: Structure, auth_asym_id: str, protein_rna: typing.Literal["protein", "rna"], sanitized: bool = False) -> str:

        chain3d = struct.child_dict[0].child_dict[auth_asym_id]
        ress = chain3d.child_list

        # if sanitized == True:
        #     print("sanitized")
        #     ress = [*filter(lambda r: r.get_resname()
        #                     in ["A", "C", "G", "U", "PSU"], ress)]

        #     for _r in ress:
        #         if _r.get_resname() == "PSU":
        #             _r.resname = "U"
        seq = ''
        for i in ress:
            if i.resname not in AMINO_ACIDS_3_TO_1_CODE.keys():
                print("Unknown residue:", i.resname)
                continue

            if protein_rna == 'rna':
                seq += i.resname
            else:
                seq += AMINO_ACIDS_3_TO_1_CODE[i.resname]

        # return reduce(lambda x, y: x + y.resname if protein_rna == 'rna' else AMINO_ACIDS_3_TO_1_CODE[y.resname], ress, '')
        return seq

    # ※ -=-=-=-=-=-=-=-=-=-=-=-=-=-=-= Verification =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    def _verify_dir_exists(self):
        if not os.path.exists(self._dir_path()):
            os.umask(0)
            os.makedirs(self._dir_path(), 0o777)

    async def _update_cif(self, overwrite: bool = False) -> bool:
        if not os.path.exists(self._cif_filepath()):
            await download_unpack_place(self.rcsb_id)
            print("Saved structure file:\t", self._cif_filepath())
            return True
        else:
            if overwrite:
                await download_unpack_place(self.rcsb_id)
                print("Saved structure file:\t", self._cif_filepath())
                return True
            else:
                return False

    async def _update_cif_modified_and_chains(self, overwrite: bool = False) -> bool:
        if not os.path.isdir(self.chains_dir()):
            os.makedirs(self.chains_dir())
            await split_rename(self.rcsb_id)
            return True
        else:
            if overwrite:
                await split_rename(self.rcsb_id)
                return True
            else:
                return False

    async def _update_json_profile(self, overwrite: bool = False) -> bool:
        self._verify_dir_exists()
        if not os.path.isfile(self._json_profile_filepath()):
            ribosome = ReannotationPipeline(query_rcsb_api(rcsb_single_structure_graphql(self.rcsb_id.upper()))).process_structure()
            if not RibosomeStructure.model_validate(ribosome):
                raise Exception("Created invalid ribosome profile (Schema validation failed). Not writing")
            self.write_own_json_profile( ribosome.model_dump(), overwrite=True)
            logger.debug("Processing {} profile.".format(self.rcsb_id))
        else:
            if overwrite:
                ribosome = ReannotationPipeline(query_rcsb_api(rcsb_single_structure_graphql(self.rcsb_id.upper()))).process_structure()
                self.write_own_json_profile(ribosome.model_dump(),overwrite)
                logger.debug("Processing {} profile: already exists for {}. Overwriting.".format(self.rcsb_id, self.rcsb_id))
            else:
                logger.debug("Processing {} profile: already exists for {}.".format(self.rcsb_id, self.rcsb_id))
        return True
        
    def _update_png_thumbnail(self, overwrite: bool = False) -> bool:
        if overwrite:
            print("Obtaning thumbnail...")
            render_thumbnail(self.rcsb_id)
            return True
        else:
            if os.path.exists(self._png_thumbnail_filepath()):
                return True
            else:
                return False

    async def _update_chains_dir(self):
        split_rename(self.rcsb_id)

    async def _update_ptc_coordinates(self, overwrite:bool=False):

        etllogger = get_etl_logger()
        asset_ptc_coords_path = os.path.join(self._dir_path(),f'{self.rcsb_id}_PTC_COORDINATES.json')

        if os.path.exists(asset_ptc_coords_path) and not overwrite:
            raise Exception(f'PTC coordinates already exist for {self.rcsb_id} and overwrite is set to False')

        ress, auth_asym_id = ptc_resdiues_get(self.biopython_structure(),  self.profile().rnas)
        midpoint_coords = ptc_residues_calculate_midpoint(ress, auth_asym_id)

        writeout = {
            "site_9_residues"      : [(res.get_resname(), res.get_segid(), res.full_id) for res in ress],
            "LSU_rRNA_auth_asym_id": auth_asym_id,
            "midpoint_coordinates" : midpoint_coords,
            'nomenclature_table'   : self._nomenclature_table()
        }

        with open(asset_ptc_coords_path, 'w') as f:
            json.dump(writeout, f)
            etllogger.info(f"Saved PTC coordinates for {self.rcsb_id} to {asset_ptc_coords_path}")


    async def _update_ligands(self, overwrite:bool=False):

        ligands           = struct_ligand_ids(self.rcsb_id, self.profile())
        
        for ligand in ligands:
            ligand_chemid  = ligand.chemicalId
            if not os.path.exists(BindingSite.path_nonpoly_ligand(self.rcsb_id, ligand_chemid)):
                bsite = bsite_ligand( ligand_chemid, self.biopython_structure())
                bsite.save(bsite.path_nonpoly_ligand( self.rcsb_id, ligand_chemid))
            else:
                if overwrite:
                    bsite = bsite_ligand( ligand_chemid, self.biopython_structure())
                    bsite.save(bsite.path_nonpoly_ligand( self.rcsb_id, ligand_chemid))
                else:
                    ...
    # async def _verify_ligads_and_ligandlike_polys(self, overwrite: bool = False):

    #     # def ligand_path(chem_id): return os.path.join(self._dir_path(), f"polymer_{chem_id.upper()}.json")
    #     # def poly_factor_path(auth_asym_id): return os.path.join(self._dir_path(), f"polymer_{auth_asym_id.upper()}.json")

    #     ligands           = struct_ligand_ids(self.rcsb_id, self.profile())
        # polymeric_factors = struct_polymeric_factor_ids(self.profile())
    #     all_verified_flag = True

    #     for ligand_chemid in ligands:
    #         if not os.path.exists(BindingSite.path_nonpoly_ligand(self.rcsb_id, ligand_chemid)):
    #             all_verified_flag = False
    #             bsite = bsite_ligand( ligand_chemid, self.biopython_structure())

    #             bsite.save(bsite.path_nonpoly_ligand( self.rcsb_id, ligand_chemid))
    #         else:
    #             if overwrite:
    #                 bsite = bsite_ligand( ligand_chemid, self.biopython_structure())
    #                 bsite.save(bsite.path_nonpoly_ligand( self.rcsb_id, ligand_chemid))
    #             else:
    #                 ...

        # if polymeric_factors is not None:
        #     for poly in polymeric_factors:
        #         if not os.path.exists(BindingSite.path_poly_factor(self.rcsb_id, poly.nomenclature[0], poly.auth_asym_id)):
        #             all_verified_flag = False
        #             bsite = bsite_extrarbx_polymer(
        #                 poly.auth_asym_id, self.biopython_structure())
        #             bsite.save(bsite.path_poly_factor(
        #                 self.rcsb_id, poly.nomenclature[0], poly.auth_asym_id))
        #         else:
        #             if overwrite:
        #                 bsite = bsite_extrarbx_polymer(
        #                     poly.auth_asym_id, self.biopython_structure())
        #                 bsite.save(bsite.path_poly_factor(
        #                     self.rcsb_id, poly.nomenclature[0], poly.auth_asym_id))
        #             else:
        #                 ...

        # return all_verified_flag


def classify_struct_by_proportions(ribosome: RibosomeStructure) -> int:
    ids = []
    if ribosome.rnas is not None:
        for rna in ribosome.rnas:
            ids = [*rna.src_organism_ids, *ids]

    for protein in ribosome.proteins:
        ids = [*protein.src_organism_ids, *ids]

    proportions = {}
    for i in set(ids):
        proportions[i] = ids.count(i) / len(ids)

    return max(proportions, key=proportions.get)