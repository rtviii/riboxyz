/* tslint:disable */
/* eslint-disable */
/**
/* This file was automatically generated from pydantic models by running pydantic2ts.
/* Do not modify it by hand - just update the pydantic models and then re-run the script
*/

export interface GOClass {
  class_id: string;
  annotation: string;
}
export interface InterProFamily {
  family_id: string;
  type: string;
  description: string;
}
export interface LastUpdate {
  date: string;
  added_structure: string;
}
export interface Ligand {
  chemicalId: string;
  chemicalName: string;
  formula_weight: number;
  pdbx_description: string;
  number_of_instances: number;
}
export interface NomeclatureClass {
  class_id:
    | (
        | "uL1"
        | "uL2"
        | "uL3"
        | "uL4"
        | "uL5"
        | "uL6"
        | "eL6"
        | "eL8"
        | "bL9"
        | "uL10"
        | "uL11"
        | "bL12"
        | "uL13"
        | "eL13"
        | "uL14"
        | "eL14"
        | "uL15"
        | "eL15"
        | "uL16"
        | "bL17"
        | "uL18"
        | "eL18"
        | "bL19"
        | "eL19"
        | "bL20"
        | "eL20"
        | "bL21"
        | "eL21"
        | "uL22"
        | "eL22"
        | "uL23"
        | "uL24"
        | "eL24"
        | "bL25"
        | "bL27"
        | "eL27"
        | "bL28"
        | "eL28"
        | "uL29"
        | "eL29"
        | "uL30"
        | "eL30"
        | "bL31"
        | "eL31"
        | "bL32"
        | "eL32"
        | "bL33"
        | "eL33"
        | "bL34"
        | "eL34"
        | "bL35"
        | "bL36"
        | "eL36"
        | "eL37"
        | "eL38"
        | "eL39"
        | "eL40"
        | "eL41"
        | "eL42"
        | "eL43"
        | "P1/P2"
      )
    | (
        | "bS1"
        | "eS1"
        | "uS2"
        | "uS3"
        | "uS4"
        | "eS4"
        | "uS5"
        | "bS6"
        | "eS6"
        | "uS7"
        | "eS7"
        | "uS8"
        | "eS8"
        | "uS9"
        | "uS10"
        | "eS10"
        | "uS11"
        | "uS12"
        | "eS12"
        | "uS13"
        | "uS14"
        | "uS15"
        | "bS16"
        | "uS17"
        | "eS17"
        | "bS18"
        | "uS19"
        | "eS19"
        | "bS20"
        | "bS21"
        | "bTHX"
        | "eS21"
        | "eS24"
        | "eS25"
        | "eS26"
        | "eS27"
        | "eS28"
        | "eS30"
        | "eS31"
        | "RACK1"
      )
    | (
        | "5SrRNA"
        | "5.8SrRNA"
        | "12SrRNA"
        | "16SrRNA"
        | "21SrRNA"
        | "23SrRNA"
        | "25SrRNA"
        | "28SrRNA"
        | "35SrRNA"
        | "mRNA"
        | "tRNA"
      );
}
export interface PFAMFamily {
  family_id: string;
  annotation: string;
  family_type: string;
}
export interface Protein {
  asym_ids: string[];
  auth_asym_id: string;
  parent_rcsb_id: string;
  pfam_accessions: string[];
  pfam_comments: string[];
  pfam_descriptions: string[];
  src_organism_names: string[];
  host_organism_names: string[];
  src_organism_ids: number[];
  host_organism_ids: number[];
  ligand_like: boolean;
  uniprot_accession: string[];
  rcsb_pdbx_description?: string;
  entity_poly_strand_id: string;
  entity_poly_seq_one_letter_code: string;
  entity_poly_seq_one_letter_code_can: string;
  entity_poly_seq_length: number;
  entity_poly_polymer_type: string;
  entity_poly_entity_type: string;
  nomenclature: (
    | (
        | "uL1"
        | "uL2"
        | "uL3"
        | "uL4"
        | "uL5"
        | "uL6"
        | "eL6"
        | "eL8"
        | "bL9"
        | "uL10"
        | "uL11"
        | "bL12"
        | "uL13"
        | "eL13"
        | "uL14"
        | "eL14"
        | "uL15"
        | "eL15"
        | "uL16"
        | "bL17"
        | "uL18"
        | "eL18"
        | "bL19"
        | "eL19"
        | "bL20"
        | "eL20"
        | "bL21"
        | "eL21"
        | "uL22"
        | "eL22"
        | "uL23"
        | "uL24"
        | "eL24"
        | "bL25"
        | "bL27"
        | "eL27"
        | "bL28"
        | "eL28"
        | "uL29"
        | "eL29"
        | "uL30"
        | "eL30"
        | "bL31"
        | "eL31"
        | "bL32"
        | "eL32"
        | "bL33"
        | "eL33"
        | "bL34"
        | "eL34"
        | "bL35"
        | "bL36"
        | "eL36"
        | "eL37"
        | "eL38"
        | "eL39"
        | "eL40"
        | "eL41"
        | "eL42"
        | "eL43"
        | "P1/P2"
      )
    | (
        | "bS1"
        | "eS1"
        | "uS2"
        | "uS3"
        | "uS4"
        | "eS4"
        | "uS5"
        | "bS6"
        | "eS6"
        | "uS7"
        | "eS7"
        | "uS8"
        | "eS8"
        | "uS9"
        | "uS10"
        | "eS10"
        | "uS11"
        | "uS12"
        | "eS12"
        | "uS13"
        | "uS14"
        | "uS15"
        | "bS16"
        | "uS17"
        | "eS17"
        | "bS18"
        | "uS19"
        | "eS19"
        | "bS20"
        | "bS21"
        | "bTHX"
        | "eS21"
        | "eS24"
        | "eS25"
        | "eS26"
        | "eS27"
        | "eS28"
        | "eS30"
        | "eS31"
        | "RACK1"
      )
  )[];
}
export interface RNA {
  asym_ids: string[];
  auth_asym_id: string;
  nomenclature: (
    | "5SrRNA"
    | "5.8SrRNA"
    | "12SrRNA"
    | "16SrRNA"
    | "21SrRNA"
    | "23SrRNA"
    | "25SrRNA"
    | "28SrRNA"
    | "35SrRNA"
    | "mRNA"
    | "tRNA"
  )[];
  parent_rcsb_id: string;
  src_organism_names: string[];
  host_organism_names: string[];
  src_organism_ids: number[];
  host_organism_ids: number[];
  rcsb_pdbx_description?: string;
  entity_poly_strand_id: string;
  entity_poly_seq_one_letter_code: string;
  entity_poly_seq_one_letter_code_can: string;
  entity_poly_seq_length: number;
  entity_poly_polymer_type: string;
  entity_poly_entity_type: string;
  ligand_like: boolean;
}
export interface RibosomeResponse {
  id: number;
  name: string;
  proteins: Protein[];
}
export interface RibosomeStructure {
  rcsb_id: string;
  expMethod: string;
  resolution: number;
  pdbx_keywords?: string;
  pdbx_keywords_text?: string;
  rcsb_external_ref_id: string[];
  rcsb_external_ref_type: string[];
  rcsb_external_ref_link: string[];
  citation_year: number;
  citation_rcsb_authors: string[];
  citation_title: string;
  citation_pdbx_doi: string;
  src_organism_ids: number[];
  src_organism_names: string[];
  host_organism_ids: number[];
  host_organism_names: string[];
  proteins: Protein[];
  rnas?: RNA[];
  ligands?: Ligand[];
}
export interface Schema {}
