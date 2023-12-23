from enum import Enum
import json
from typing import List, Literal, Optional, Union
from pydantic import BaseModel, field_serializer

from ribctl.lib.ribosome_types.types_ribosome import CytosolicProteinClass, MitochondrialProteinClass, PolymerClass

raw_chain = {
            "assembly_id": 1,
            "asym_ids": [
                "B",
                "CC"
            ],
            "auth_asym_id": "S0",
            "parent_rcsb_id": "5DAT",
            "src_organism_names": [
                "Saccharomyces cerevisiae (strain ATCC 204508 / S288c)"
            ],
            "host_organism_names": [],
            "src_organism_ids": [
                559292
            ],
            "host_organism_ids": [],
            "rcsb_pdbx_description": "40S ribosomal protein S0-A",
            "entity_poly_strand_id": "S0,s0",
            "entity_poly_seq_one_letter_code": "SLPATFDLTPEDAQLLLAANTHLGARNVQVHQEPYVFNARPDGVHVINVGKTWEKLVLAARIIAAIPNPEDVVAISSRTFGQRAVLKFAAHTGATPIAGRFTPGSFTNYITRSFKEPRLVIVTDPRSDAQAIKEASYVNIPVIALTDLDSPSEFVDVAIPCNNRGKHSIGLIWYLLAREVLRLRGALVDRTQPWSIMPDLYFYRDPEEVEQQVAEEATTEEAGEEEAKEEVTEEQAEATEWAEENADNVEW",
            "entity_poly_seq_one_letter_code_can": "SLPATFDLTPEDAQLLLAANTHLGARNVQVHQEPYVFNARPDGVHVINVGKTWEKLVLAARIIAAIPNPEDVVAISSRTFGQRAVLKFAAHTGATPIAGRFTPGSFTNYITRSFKEPRLVIVTDPRSDAQAIKEASYVNIPVIALTDLDSPSEFVDVAIPCNNRGKHSIGLIWYLLAREVLRLRGALVDRTQPWSIMPDLYFYRDPEEVEQQVAEEATTEEAGEEEAKEEVTEEQAEATEWAEENADNVEW",
            "entity_poly_seq_length": 251,
            "entity_poly_polymer_type": "Protein",
            "entity_poly_entity_type": "polypeptide(L)",
            "nomenclature": [
                "uS2"
            ],
            "pfam_accessions": [
                "PF00318"
            ],
            "pfam_comments": [
                "NULL"
            ],
            "pfam_descriptions": [
                "Ribosomal protein S2"
            ],
            "uniprot_accession": [
                "P32905"
            ]
        }

class Polymer(BaseModel):

    def __hash__(self):
        return hash(self.auth_asym_id + self.parent_rcsb_id)

    # def to_dict(self):
    #     """A hack for enum.union to work with pydantic BaseModel. Otherwise EnumUnion instances are represented as <MitochondrialProteinClass.mL64: 'mL64'> etc.(Correct is "mL64")"""
    #     return json.loads(self.json())


    # FIXME:
    # def to_SeqRecord(self)->SeqRecord:
    #     return SeqRecord(self.entity_poly_seq_one_letter_code_can, id=f"{self.parent_rcsb_id}.{self.auth_asym_id}", 
    #                      description="{}|{}".format(self.src_organism_ids[0],self.rcsb_pdbx_description))
    def __str__(self) -> str:
        return super().model_dump_json()

    assembly_id: int

    asym_ids    : list[str]
    auth_asym_id: str

    parent_rcsb_id   : str

    src_organism_names : list[str]
    host_organism_names: list[str]

    src_organism_ids   : list[int]
    host_organism_ids  : list[int]

    rcsb_pdbx_description: Optional[str] 

    entity_poly_strand_id              : str
    entity_poly_seq_one_letter_code    : str
    entity_poly_seq_one_letter_code_can: str
    entity_poly_seq_length             : int
    entity_poly_polymer_type           : str
    entity_poly_entity_type            : str

    @field_serializer('nomenclature')
    def serialize_dt(self, nomenclature:list[PolymerClass], _info):
        return [x.name for x in nomenclature]
    nomenclature:  list[PolymerClass]


class Protein(Polymer):


    pfam_accessions  : list[str]
    pfam_comments    : list[str]
    pfam_descriptions: list[str]
    uniprot_accession: list[str]

    def __hash__(self):
        return hash(self.auth_asym_id + self.parent_rcsb_id)

    @staticmethod
    def from_polymer(p: Polymer,
                      **kwargs): 

        if ( kwargs["pfams"] != None and len(kwargs["pfams"]) > 0 ):
            pfam_comments     = list( set( [ pfam["rcsb_pfam_comment"    ] for pfam in kwargs["pfams"] ] ) )
            pfam_descriptions = list( set( [ pfam["rcsb_pfam_description"] for pfam in kwargs["pfams"] ] ) )
            pfam_accessions   = list( set( [ pfam["rcsb_pfam_accession"  ] for pfam in kwargs["pfams"] ] ) )

        else:
            pfam_comments     = []
            pfam_descriptions = []
            pfam_accessions   = []

        
        return Protein(**{
              **p.model_dump(),
              "pfam_accessions"   : pfam_accessions,
              "pfam_comments"     : pfam_comments,
              "pfam_descriptions" : pfam_descriptions,
              "uniprot_accession" : [ entry["rcsb_id"] for entry in kwargs["uniprots"] ] if kwargs["uniprots"] != None and len(kwargs["uniprots"]) > 0 else []
        })


    def to_polymer(self)->Polymer:
        return Polymer(**self.model_dump())



with open('test.json', 'r') as f:
    poly = Polymer.model_validate(json.load(f))
    print("Uploaded poly", poly)

with open('test1.json', 'w') as f:
    json.dump(poly.model_dump(), f)
