import axios from "axios";
import { pull, uniq } from "lodash";
import {
  ProteinClass,
  Ligand,
  Protein,
  RibosomeStructure,
  RNAClass,
  RNA,
} from "./ribosome_types";

// both should be served by django
import { small_subunit_map } from "./subunit_maps/small-subunit-map";
import { large_subunit_map } from "./subunit_maps/large-subunit-map";

import {
  Nonpolymer_Entity,
  PDBGQLResponse,
  Polymer_Entity,
} from "./rcsb_graphql_schema";




// Renamed, added so far:  host organisms, host ids, srcids srcnames
const query_template = (pdbid: string) => {
  const GQL_QUERY_SHAPE = `{
  entry(entry_id: "${pdbid.toUpperCase()}") {
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
        formula_weight
        pdbx_description
        pdbx_number_of_molecules
      }
    }
  }
}`
  return `https://data.rcsb.org/graphql?query=${GQL_QUERY_SHAPE}`;
};

const matchRPNomenclature = (
  polymer: Polymer_Entity
): ProteinClass[] => {

  var banregex = /\b([ueb][ls]\d{1,2})\b/gi

  // * check authors's annotations. if classes are present --> use that.
  var finds = banregex.exec(polymer.rcsb_polymer_entity.pdbx_description)
  if (finds !== null) {
    var firstcap = finds[0]
    var name = (
      firstcap[0].toLowerCase()
      + firstcap[1].toUpperCase()
      + firstcap.slice(2)
    ) as ProteinClass
    return [name]
  }
  // * otherwise resort to pfams
  if (!polymer.pfams) {
    return []
  } else {
    var pfamIds = polymer.pfams.reduce<string[]>((acc, item) => {
      acc.push(item.rcsb_pfam_accession);
      return acc;
    }, []);
    var nomenclature: ProteinClass[] = [];
    pfamIds!.map(id => {
      Object.entries(large_subunit_map).map(entry => {
        if (entry[1].pfamDomainAccession.includes(id)) { nomenclature.push(entry[0] as ProteinClass); }
      });
      Object.entries(small_subunit_map).map(entry => {
        if (entry[1].pfamDomainAccession.includes(id)) {
          nomenclature.push(entry[0] as ProteinClass);
        }
      });
    });

    return uniq(nomenclature)
  }
};

const matchRNANomenclature = (polymer: Polymer_Entity): RNAClass[] => {

  const rna_reg: Record<RNAClass, RegExp> = {
    "5SrRNA": /\b(5s)/gi,
    "5.8SrRNA": /\b(5\.8s)/gi,
    "12SrRNA": /\b(12s)/gi,
    "16SrRNA": /\b(16s)/gi,
    "21SrRNA": /\b(21s)/gi,
    "23SrRNA": /\b(23s)/gi,
    "25SrRNA": /\b(25s)/gi,
    "28SrRNA": /\b(28s)/gi,
    "35SrRNA": /\b(35s)/gi,
    "mRNA": /(mrna)|\b(messenger)\b/gi,
    "tRNA": /(trna)|\b(transfer)\b/gi,
  }

  var rnatypes = Object.keys(rna_reg) as RNAClass[]
  for (var c of rnatypes) {
    var matches = rna_reg[c].exec(polymer.rcsb_polymer_entity.pdbx_description)
    if (matches) {
      return [c as RNAClass]
    }
  }
  return []

}

interface Organisms {
  src_organism_ids: number[]
  src_organism_names: string[]
  host_organism_ids: number[]
  host_organism_names: string[]

}
const inferOrganismsFromPolymers = (proteins: Protein[]): Organisms => {

  var host_organism_names: string[] = [];
  var src_organism_names: string[] = [];

  var host_organism_ids: number[] = [];
  var src_organism_ids: number[] = [];

  proteins.map(protein => {
    protein.src_organism_names ? src_organism_names.push(...protein.src_organism_names) : null;
    protein.src_organism_ids ? src_organism_ids.push(...protein.src_organism_ids) : null;
    protein.host_organism_names ? src_organism_names.push(...protein.host_organism_names) : null;
    protein.host_organism_ids ? src_organism_ids.push(...protein.host_organism_ids) : null;
  });

  return {
    src_organism_ids: uniq(src_organism_ids),
    src_organism_names: uniq(src_organism_names),
    host_organism_ids: uniq(host_organism_ids),
    host_organism_names: uniq(host_organism_names),
  };
};

const extractRefs = (
  external_refs: { link: string; type: string; id: string }[]
) => {

  var externalRefIds: string[] = [];
  var externalRefTypes: string[] = [];
  var externalRefLinks: string[] = [];

  external_refs === null ? (() => { }) : external_refs.map(ref => {
    externalRefIds.push(ref.id);
    externalRefTypes.push(ref.type);
    externalRefLinks.push(ref.link);
  });

  return [externalRefIds, externalRefTypes, externalRefLinks];
};

const reshape_ToLigand = (nonpoly: Nonpolymer_Entity): Ligand => {
  return {
    pdbx_description: nonpoly.rcsb_nonpolymer_entity.pdbx_description,
    formula_weight: nonpoly.rcsb_nonpolymer_entity.formula_weight,
    chemicalId: nonpoly.pdbx_entity_nonpoly.comp_id,
    chemicalName: nonpoly.pdbx_entity_nonpoly.name,
    number_of_instances: nonpoly.rcsb_nonpolymer_entity.pdbx_number_of_molecules
  };
};

const is_ligand_like = (plm: Polymer_Entity, nomenclature: Array<RNAClass | ProteinClass>) => {
  if (nomenclature.includes('tRNA') || nomenclature.includes('mRNA')) {
    return true
  }

  // ? Look for enzymes, factors and antibiotics
  var reg = /(\w*(?<!(cha|pro|dom|stra|pl\w*))in\b)|(\b\w*zyme\b)|(factor)/gi;
  const matches = plm.rcsb_polymer_entity.pdbx_description.match(reg)
  if (matches !== null &&
    !(matches.map(_ => _.toLowerCase()).includes('protein')) &&
    !(plm.rcsb_polymer_entity.pdbx_description.toLowerCase().includes('protein')) &&
    !(plm.rcsb_polymer_entity.pdbx_description.toLowerCase().includes('rrna'))
  ) {
    return true
  }
  return false
}

const reshape_PolyEntity_to_rRNA = (plm: Polymer_Entity): RNA[] => {

  var src_organism_ids: number[] = plm.rcsb_entity_source_organism ? uniq(plm.rcsb_entity_source_organism.map(org => org.ncbi_taxonomy_id)) : [];
  var src_organism_names: string[] = plm.rcsb_entity_source_organism ? uniq(plm.rcsb_entity_source_organism.map(org => org.scientific_name)) : [];

  var host_organism_ids: number[] = plm.rcsb_entity_host_organism ? uniq(plm.rcsb_entity_host_organism.map(org => org.ncbi_taxonomy_id)) : [];
  var host_organism_names: string[] = plm.rcsb_entity_host_organism ? uniq(plm.rcsb_entity_host_organism.map(org => org.scientific_name)) : [];
  var nomenclature = matchRNANomenclature(plm)


  return plm.rcsb_polymer_entity_container_identifiers.auth_asym_ids.map(auth_asym_id => (
    {
      nomenclature: nomenclature,
      ligand_like: is_ligand_like(plm, nomenclature),
      asym_ids: plm.rcsb_polymer_entity_container_identifiers.asym_ids,
      auth_asym_id: auth_asym_id,
      parent_rcsb_id: plm.entry.rcsb_id,

      host_organism_ids: host_organism_ids,
      host_organism_names: host_organism_names,
      src_organism_ids: src_organism_ids,
      src_organism_names: src_organism_names,

      rcsb_pdbx_description: plm.rcsb_polymer_entity.pdbx_description === null ? "" : plm.rcsb_polymer_entity.pdbx_description,
      entity_poly_strand_id: plm.entity_poly.pdbx_strand_id,
      entity_poly_seq_one_letter_code: plm.entity_poly.pdbx_seq_one_letter_code,
      entity_poly_seq_one_letter_code_can: plm.entity_poly.pdbx_seq_one_letter_code_can,
      entity_poly_seq_length: plm.entity_poly.rcsb_sample_sequence_length,
      entity_poly_entity_type: plm.entity_poly.type,
      entity_poly_polymer_type: plm.entity_poly.rcsb_entity_polymer_type
    }
  ))




}
const reshape_PolyEntity_to_RibosomalProtein = (plm: Polymer_Entity): Protein[] => {
  if (plm.pfams) {
    var pfam_comments = uniq(
      plm.pfams.map(pfam => pfam.rcsb_pfam_comment)
    );
    var pfam_descriptions = uniq(
      plm.pfams.map(pfam => pfam.rcsb_pfam_description)
    );
    var pfam_accessions = uniq(
      plm.pfams.map(pfam => pfam.rcsb_pfam_accession)
    );
  }
  else {
    var pfam_comments: string[] = [];
    var pfam_descriptions: string[] = [];
    var pfam_accessions: string[] = [];
  }

  var src_organism_ids: number[] = uniq(plm.rcsb_entity_source_organism.map(org => org.ncbi_taxonomy_id));
  var src_organism_names: string[] = uniq(plm.rcsb_entity_source_organism.map(org => org.scientific_name));
  var host_organism_ids: number[] = plm.rcsb_entity_host_organism ? uniq(plm.rcsb_entity_host_organism.map(org => org.ncbi_taxonomy_id)) : [];
  var host_organism_names: string[] = plm.rcsb_entity_host_organism ? uniq(plm.rcsb_entity_host_organism.map(org => org.scientific_name)) : [];
  var nomenclature = matchRPNomenclature(plm)
  return plm.rcsb_polymer_entity_container_identifiers.auth_asym_ids.map(auth_asym_id => (
    {
      nomenclature     : nomenclature                                                        ,
      asym_ids         : plm              .rcsb_polymer_entity_container_identifiers.asym_ids,
      auth_asym_id     : auth_asym_id                                                        ,
      parent_rcsb_id   : plm              .entry.rcsb_id                                     ,
      pfam_accessions  : pfam_accessions                                                     ,
      pfam_comments    : pfam_comments                                                       ,
      pfam_descriptions: pfam_descriptions                                                   ,

      ligand_like: is_ligand_like(plm, nomenclature),

      host_organism_ids: host_organism_ids,
      host_organism_names: host_organism_names,
      src_organism_ids: src_organism_ids,
      src_organism_names: src_organism_names,

      uniprot_accession: plm.uniprots ? plm.uniprots.map(entry => entry.rcsb_id) : [],
      rcsb_pdbx_description: plm.rcsb_polymer_entity.pdbx_description,
      entity_poly_strand_id: plm.entity_poly.pdbx_strand_id,
      entity_poly_seq_one_letter_code: plm.entity_poly.pdbx_seq_one_letter_code,
      entity_poly_seq_one_letter_code_can: plm.entity_poly.pdbx_seq_one_letter_code_can,
      entity_poly_seq_length: plm.entity_poly.rcsb_sample_sequence_length,
      entity_poly_entity_type: plm.entity_poly.type,
      entity_poly_polymer_type: plm.entity_poly.rcsb_entity_polymer_type
    }
  ))

}

export const processPDBRecord = async (
  pdbid: string
): Promise<RibosomeStructure> => {
  return await axios.get(query_template(pdbid))
    .then(response => {


      var pdbRecord: PDBGQLResponse = response.data.data;
      var proteins = pdbRecord.entry.polymer_entities.filter(poly => poly.entity_poly.rcsb_entity_polymer_type === "Protein");
      var rnas = pdbRecord.entry.polymer_entities.filter(poly => poly.entity_poly.rcsb_entity_polymer_type === "RNA");
      var ligands = pdbRecord.entry.nonpolymer_entities;
      var reshaped_proteins: Protein[] = [];
      var reshaped_rnas: RNA[] = [];


      proteins.map(protein => { reshaped_proteins.push(...reshape_PolyEntity_to_RibosomalProtein(protein)) });
      rnas.map(rna => { reshaped_rnas.push(...reshape_PolyEntity_to_rRNA(rna)) });
      var reshaped_nonpoly: Ligand[] | null = ligands == null ? null : ligands.map(r => reshape_ToLigand(r));

      var organisms = inferOrganismsFromPolymers(reshaped_proteins)
      var externalRefs = extractRefs(pdbRecord.entry.rcsb_external_references);
      var pub = pdbRecord.entry.citation[0];

      var kwords_text = pdbRecord.entry.struct_keywords
        ? pdbRecord.entry.struct_keywords.text
        : null;
      var kwords = pdbRecord.entry.struct_keywords
        ? pdbRecord.entry.struct_keywords.pdbx_keywords
        : null;

      var reshaped: RibosomeStructure = {
        rcsb_id: pdbRecord.entry.rcsb_id,
        expMethod: pdbRecord.entry.exptl[0].method,
        resolution: pdbRecord.entry.rcsb_entry_info.resolution_combined[0],

        rcsb_external_ref_id: externalRefs[0],
        rcsb_external_ref_type: externalRefs[1],
        rcsb_external_ref_link: externalRefs[2],

        citation_year: pub.year,
        citation_rcsb_authors: pub.rcsb_authors,
        citation_title: pub.title,
        citation_pdbx_doi: pub.pdbx_database_id_DOI,

        ...organisms,

        pdbx_keywords_text: kwords_text,
        pdbx_keywords: kwords,
        proteins: reshaped_proteins,
        rnas: reshaped_rnas,
        ligands: reshaped_nonpoly,
      };

      return reshaped;
    });

};