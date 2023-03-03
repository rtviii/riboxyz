import os
import typing
from pydantic import BaseModel
from ribctl.lib.types.types_polymer import LSU_Proteins, RNAClass, SSU_Proteins
import json
from pydantic import parse_obj_as
from ribctl.lib import RIBETL_DATA
from ribctl.lib.utils import download_unpack_place, open_structure
from ribctl.lib.struct_render_thumbnail import render_thumbnail
from ribctl.lib.struct_rcsb_api import process_pdb_record
from ribctl.lib.struct_split_rename import split_rename
from ribctl.lib.struct_extract_bsites import get_ligands, get_liglike_polymers, render_ligand, render_liglike_polymer
from ninja import Schema

ProteinClass = LSU_Proteins | SSU_Proteins




class LastUpdate(BaseModel):
    date: str
    added_structure: str

class Protein(BaseModel):
    asym_ids: list[str]
    auth_asym_id: str

    parent_rcsb_id: str
    pfam_accessions: list[str]
    pfam_comments: list[str]
    pfam_descriptions: list[str]

    src_organism_names: list[str]
    host_organism_names: list[str]
    src_organism_ids: list[int]
    host_organism_ids: list[int]

    ligand_like: bool

    uniprot_accession: list[str]

    rcsb_pdbx_description: str | None

    entity_poly_strand_id: str
    entity_poly_seq_one_letter_code: str
    entity_poly_seq_one_letter_code_can: str
    entity_poly_seq_length: int
    entity_poly_polymer_type: str
    entity_poly_entity_type: str

    nomenclature:  list[ProteinClass]


# ※ TESTING NINJA
class RibosomeResponse(Schema):
    id: int
    name: str
    proteins: list[Protein]

class RNA(BaseModel):

    asym_ids: list[str]

    auth_asym_id: str
    nomenclature: list[RNAClass]
    parent_rcsb_id: str

    src_organism_names: list[str]
    host_organism_names: list[str]
    src_organism_ids: list[int]
    host_organism_ids: list[int]

    rcsb_pdbx_description: str | None
    # entity_polymer
    entity_poly_strand_id: str
    entity_poly_seq_one_letter_code: str
    entity_poly_seq_one_letter_code_can: str
    entity_poly_seq_length: int
    entity_poly_polymer_type: str
    entity_poly_entity_type: str

    ligand_like: bool



class Ligand(BaseModel):
    chemicalId: str
    chemicalName: str
    formula_weight: float
    pdbx_description: str
    number_of_instances: int


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
        return f"{RIBETL_DATA}/{self.rcsb_id}"

    def _cif_filepath(self):
        self._envcheck()
        return f"{self._dir_path()}/{self.rcsb_id}.cif"

    def _cif_modified_filepath(self):
        self._envcheck()
        return f"{self._dir_path()}/{self.rcsb_id}_modified.cif"

    def _json_profile_filepath(self):
        self._envcheck()
        return f"{self._dir_path()}/{self.rcsb_id}.json"

    def json_profile(self):
        with open(self._json_profile_filepath(), "r") as f:
            return json.load(f)

    def biopython_sturcture(self):
        return open_structure(self.rcsb_id, 'cif')

    def chains_dir(self):
        self._envcheck()
        return f"{self._dir_path()}/CHAINS"

    def _png_thumbnail_filepath(self):
        self._envcheck()
        return f"{self._dir_path()}/_ray_{self.rcsb_id}.png"

    def save_json_profile(self, filepath: str, profile: dict):
        with open(filepath, "w") as f:
            json.dump(profile, f)

    def _verify_cif(self, obtain: bool = False) -> bool:
        if os.path.exists(self._cif_filepath()):
            return True
        else:
            if obtain:
                download_unpack_place(self.rcsb_id)
                return True
            else:
                return False

    def _verify_cif_modified(self, obtain: bool = False) -> bool:
        if os.path.exists(self._cif_modified_filepath()):
            return True
        else:
            if obtain:
                split_rename(self.rcsb_id)
                return True
            else:
                return False

    def _verify_json_profile(self, obtain: bool = False) -> bool:
        if os.path.exists(self._json_profile_filepath()):
            return True
        else:
            if obtain:
                ribosome = process_pdb_record(self.rcsb_id)
                if not parse_obj_as(RibosomeStructure, ribosome):
                    raise Exception("Invalid ribosome structure profile.")

                self.save_json_profile(self._json_profile_filepath(), ribosome)
                print(
                    f"Saved structure profile:\t{self._json_profile_filepath()}")
                return True
            else:
                return False

    def _verify_png_thumbnail(self, obtain: bool = False) -> bool:
        if os.path.exists(self._png_thumbnail_filepath()):
            return True
        else:
            if obtain:
                print("Obtaning thumbnail...")
                render_thumbnail(self.rcsb_id)

                return True
            else:
                return False

    def _verify_chains_dir(self):
        split_rename(self.rcsb_id)

    def _verify_ligads_and_ligandlike_polys(self, obtain: bool = False):

        def ligand_path(chem_id): return os.path.join(
            self._dir_path(), f"LIGAND_{chem_id.upper()}.json")
        def liglike_poly_path(auth_asym_id): return os.path.join(
            self._dir_path(), f"POLYMER_{auth_asym_id.upper()}.json")

        ligands = get_ligands(self.rcsb_id, self.json_profile())
        ligandlike_polymers = get_liglike_polymers(self.json_profile())

        _flag = True

        for ligand in ligands:
            if not os.path.exists(ligand_path(ligand[0])):
                _flag = False
                render_ligand(
                    self.rcsb_id, ligand[0], self.biopython_sturcture(), obtain)

        for ligandlike_poly in ligandlike_polymers:
            if not os.path.exists(liglike_poly_path(ligandlike_poly.auth_asym_id)):
                _flag = False
                render_liglike_polymer(
                    self.rcsb_id, ligandlike_poly.auth_asym_id, self.biopython_sturcture(), obtain)

        return _flag


class RibosomeStructure(BaseModel):

    rcsb_id:    str
    expMethod:  str
    resolution: float

    pdbx_keywords:      str | None
    pdbx_keywords_text: str | None

    rcsb_external_ref_id: list[str]
    rcsb_external_ref_type: list[str]
    rcsb_external_ref_link: list[str]

    citation_year: int
    citation_rcsb_authors: list[str]
    citation_title: str
    citation_pdbx_doi: str

    src_organism_ids: list[int]
    src_organism_names: list[str]

    host_organism_ids: list[int]
    host_organism_names: list[str]

    proteins: list[Protein]
    rnas: list[RNA] | None
    ligands: list[Ligand] | None



    @staticmethod
    def from_json_profile(rcsb_id: str):

        _rib = RibosomeAssets(rcsb_id.upper()).json_profile()
        rib  = RibosomeStructure(**_rib)

        return rib

    # def _ingres_split_rename():
    #     ...

    # def _ingres_extract_bsites():
    #     ...

    # def _ingres_render_thumbnail():
    #     ...

    # def _ingres_commit_structure():
    #     ...

# ※--------------------------------------------------------


class InterProFamily(BaseModel):
    family_id: str
    type: str
    description: str


class GOClass(BaseModel):
    class_id: str
    annotation: str


class PFAMFamily(BaseModel):
    family_id: str
    annotation: str
    family_type: str


class NomeclatureClass(BaseModel):
    class_id:  ProteinClass | RNAClass
