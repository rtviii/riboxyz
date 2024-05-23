assembly_identification_string = """
assemblies{
    rcsb_id 
   	nonpolymer_entity_instances{
  
      rcsb_nonpolymer_entity_instance_container_identifiers{
        
        comp_id
        auth_asym_id
        rcsb_id
        auth_seq_id
      }
    }
    polymer_entity_instances{
       rcsb_polymer_entity_instance_container_identifiers {        
        asym_id
        auth_asym_id
        entry_id
        entity_id
      }
    }
  }"""

entry_info_string = """
    struct_keywords {
      pdbx_keywords
      text
    }
    rcsb_entry_info {
      resolution_combined
    }
    rcsb_external_references {
      link
      type
      id
    }
    exptl {
      method
    }
    citation {
      rcsb_authors
      year
      title
      pdbx_database_id_DOI
    }
    struct_keywords {
      pdbx_keywords
      text
    }

  rcsb_accession_info{
    deposit_date
  }

"""

single_structure_graphql_template = """{
  entry(entry_id: "$RCSB_ID") {
    rcsb_id

    assemblies{
        rcsb_id 
       	nonpolymer_entity_instances{
      
          rcsb_nonpolymer_entity_instance_container_identifiers{
            auth_asym_id
            auth_seq_id
            entity_id
          }
        }
        polymer_entity_instances{
           rcsb_polymer_entity_instance_container_identifiers {        
            auth_asym_id
            entity_id
          }
        }
      }



    struct_keywords {
        pdbx_keywords
        text
      }
      rcsb_entry_info {
        resolution_combined
      }
      rcsb_external_references {
        link
        type
        id
      }
      exptl {
        method
      }
      citation {
        rcsb_authors
        year
        title
        pdbx_database_id_DOI
      }
      struct_keywords {
        pdbx_keywords
        text
      }



    polymer_entities {
      entry{
        rcsb_id
      }
      rcsb_polymer_entity_container_identifiers {
        asym_ids
        auth_asym_ids
        entry_id
        entity_id
      }
      pfams {
        rcsb_pfam_accession
        rcsb_pfam_comment
        rcsb_pfam_description
      }
      rcsb_entity_source_organism {
        ncbi_taxonomy_id
        scientific_name
      }

      rcsb_entity_host_organism {
        ncbi_taxonomy_id
        scientific_name
      }

      uniprots {
        rcsb_id
      }
      rcsb_polymer_entity {
        pdbx_description
      }
      entity_poly {
        pdbx_seq_one_letter_code
        pdbx_seq_one_letter_code_can
        pdbx_strand_id
        rcsb_entity_polymer_type
        rcsb_sample_sequence_length
        type
      }
    }


    nonpolymer_entities {
      pdbx_entity_nonpoly {
        comp_id
        name
        entity_id
      }
      rcsb_nonpolymer_entity {
        details
        formula_weight
        pdbx_description
        pdbx_number_of_molecules
      }
      nonpolymer_comp {
        chem_comp {
          id
          name
          three_letter_code
        }
        rcsb_chem_comp_target {
          interaction_type
          name
          provenance_source
          reference_database_accession_code
          reference_database_name
        }

        drugbank {
          drugbank_container_identifiers {
            drugbank_id
          }
          drugbank_info {
            cas_number
            description
            indication
            mechanism_of_action
            name
            pharmacology
          }
        }
      }
      
    }


  }
}"""

polymer_entities_string = """{
  entry(entry_id: "$RCSB_ID") {
   assemblies {
      polymer_entity_instances {
          rcsb_id
        }
      }


    polymer_entities {
      rcsb_polymer_entity_container_identifiers {
        asym_ids
        auth_asym_ids
        entry_id
        entity_id
      }
      pfams {
        rcsb_pfam_accession
        rcsb_pfam_comment
        rcsb_pfam_description
      }
      rcsb_entity_source_organism {
        ncbi_taxonomy_id
        scientific_name
      }
      rcsb_entity_host_organism {
        ncbi_taxonomy_id
        scientific_name
      }

      uniprots {
        rcsb_id
      }
      rcsb_polymer_entity {
        pdbx_description
      }
      entity_poly {
        pdbx_seq_one_letter_code
        pdbx_seq_one_letter_code_can
        pdbx_strand_id
        rcsb_entity_polymer_type
        rcsb_sample_sequence_length
        type
      }
       rcsb_polymer_entity_annotation {
        annotation_id
        assignment_version
        description
        name
        provenance_source
        type
      }
    }
  }
}"""

nonpolymer_entities_string = """{
  entry(entry_id: "$RCSB_ID") {
     nonpolymer_entities {
      pdbx_entity_nonpoly {
        comp_id
        name
        entity_id
      }
      rcsb_nonpolymer_entity {
        details
        formula_weight
        pdbx_description
        pdbx_number_of_molecules
      }
      nonpolymer_comp {
        chem_comp {
          id
          name
          three_letter_code
        }
        rcsb_chem_comp_target {
          interaction_type
          name
          provenance_source
          reference_database_accession_code
          reference_database_name
        }

        drugbank {
          drugbank_container_identifiers {
            drugbank_id
          }
          drugbank_info {
            cas_number
            description
            indication
            mechanism_of_action
            name
            pharmacology
          }
        }
      }
      
    }
  }
}"""

structure_string = """{
  entry(entry_id:"$RCSB_ID"){
    rcsb_id
    struct_keywords {
      pdbx_keywords
      text
    }
    rcsb_entry_info {
      resolution_combined
    }
    rcsb_external_references {
      link
      type
      id
    }
    exptl {
      method
    }
    citation {
      rcsb_authors
      year
      title
      pdbx_database_id_DOI
    }
    struct_keywords {
      pdbx_keywords
      text
    }
    assemblies{


      rcsb_assembly_info {
        assembly_id
        atom_count
        branched_atom_count
        branched_entity_count
        branched_entity_instance_count
        hydrogen_atom_count
        modeled_polymer_monomer_count
        na_polymer_entity_types
        nonpolymer_atom_count
        nonpolymer_entity_count
        nonpolymer_entity_instance_count
        num_heterologous_interface_entities
        num_heteromeric_interface_entities
        num_homomeric_interface_entities
        num_interface_entities
        num_interfaces
        num_isologous_interface_entities
        num_na_interface_entities
        num_prot_na_interface_entities
        num_protein_interface_entities
        polymer_atom_count
        polymer_composition
        polymer_entity_count
        polymer_entity_count_DNA
        polymer_entity_count_RNA
        polymer_entity_count_nucleic_acid
        polymer_entity_count_nucleic_acid_hybrid
        polymer_entity_count_protein
        polymer_entity_instance_count
        polymer_entity_instance_count_DNA
        polymer_entity_instance_count_RNA
        polymer_entity_instance_count_nucleic_acid
        polymer_entity_instance_count_nucleic_acid_hybrid
        polymer_entity_instance_count_protein
        polymer_monomer_count
        selected_polymer_entity_types
        solvent_atom_count
        solvent_entity_count
        solvent_entity_instance_count
        total_assembly_buried_surface_area
        total_number_interface_residues
        unmodeled_polymer_monomer_count
      }
    }
  }
}"""
