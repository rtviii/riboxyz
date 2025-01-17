import typing
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from ribctl import AMINO_ACIDS_3_TO_1_CODE
from ribctl.asset_manager.assets_structure import StructureAssets
from Bio.PDB.Structure import Structure
from ribctl.lib.types.polymer.base import CytosolicRNAClass, MitochondrialRNAClass
from ribctl.lib.schema.types_ribosome import ( RNA, PTCInfo, Polymer, PolymerClass,  RibosomeStructure, RibosomeStructureMetadata, )
from ribctl import RIBETL_DATA

class RibosomeOps:
    rcsb_id: str
    assets: StructureAssets

    def __init__(self, rcsb_id: str) -> None:
        if not RIBETL_DATA:
            raise Exception("RIBETL_DATA environment variable not set. Cannot access assets." )
        self.rcsb_id = rcsb_id.upper()
        self.assets  = StructureAssets(self.rcsb_id)

    @property
    def taxid(self)->int:
        return self.assets.profile().src_organism_ids[0]

    @property
    def profile(self) -> RibosomeStructure:
        return self.assets.profile()

    def get_biopython_chain_by_polymer_class(self, polymer_class: PolymerClass) -> Chain:
        model = self.assets.biopython_structure()[0]
        poly = self.get_poly_by_polyclass(polymer_class)
        if poly == None:
            raise KeyError("No polymer found with class: {}".format(polymer_class))
        return model.child_dict[poly.auth_asym_id]

    def get_biopython_chain_by_auth_asym_id(self, auth_asym_id: str) -> Chain:
        model = self.assets.biopython_structure()[0]
        return model.child_dict[auth_asym_id]

    def nomenclature_table(self, verbose: bool = False) -> dict[str, dict]:
        prof = self.assets.profile()
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

    def get_taxids(self) -> tuple[list[int], list[int]]:
        p = self.assets.profile()
        return (p.src_organism_ids, p.host_organism_ids)

    def get_poly_by_auth_asym_id( self, auth_asym_id: str ) -> Polymer :
        profile = self.assets.profile()
        for chain in [ *profile.proteins, *profile.rnas, *profile.other_polymers]:
            if chain.auth_asym_id == auth_asym_id:
                return chain
        raise KeyError("No chain found with auth_asym_id: {}".format(auth_asym_id))

    def get_poly_by_polyclass( self, class_: PolymerClass, assembly: int = 0 ) ->Polymer | None:
        """@assembly here stands to specify which of the two or more models the rna comes from
        in the case that a structure contains multiple models (ex. 4V4Q XRAY)"""

        profile = self.assets.profile()
        polymer:Polymer
        for polymer in [*profile.rnas, *profile.other_polymers, *profile.proteins]: 

            if class_ in  polymer.nomenclature and polymer.assembly_id == assembly and polymer.entity_poly_seq_length > 30:
                return polymer
    def get_LSU_rRNA(self, assembly: int = 0) -> RNA:
        """retrieve the largest rRNA sequence in the structure
        @returns (seq, auth_asym_id, rna_type)
        """

        rna = self.get_poly_by_polyclass(CytosolicRNAClass.rRNA_23S, assembly)
        if rna == None:
            rna = self.get_poly_by_polyclass(CytosolicRNAClass.rRNA_25S, assembly)
        if rna == None:
            rna = self.get_poly_by_polyclass(CytosolicRNAClass.rRNA_28S, assembly)
        if rna == None:
            rna = self.get_poly_by_polyclass(MitochondrialRNAClass.mtrRNA16S, assembly)
        if rna == None:
            raise Exception("No LSU rRNA found in structure")
        else:
            return rna

    @staticmethod
    def biopython_chain_get_seq(
        struct: Structure,
        auth_asym_id: str,
        protein_rna: typing.Literal["protein", "rna"],
        sanitized: bool = False,
    ) -> str:

        chain3d = struct.child_dict[0].child_dict[auth_asym_id]
        ress    = chain3d.child_list

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