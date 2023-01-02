"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.processPDBRecord = void 0;
const tslib_1 = require("tslib");
const axios_1 = tslib_1.__importDefault(require("axios"));
const lodash_1 = require("lodash");
// both should be served by django
const small_subunit_map_1 = require("../subunit_maps/small-subunit-map");
const large_subunit_map_1 = require("../subunit_maps/large-subunit-map");
// Renamed, added so far:  host organisms, host ids, srcids srcnames
const query_template = (pdbid) => {
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
}`;
    return `https://data.rcsb.org/graphql?query=${GQL_QUERY_SHAPE}`;
};
const matchRPNomenclature = (polymer) => {
    var banregex = /\b([ueb][ls]\d{1,2})\b/gi;
    // * check authors's annotations. if classes are present --> use that.
    var finds = banregex.exec(polymer.rcsb_polymer_entity.pdbx_description);
    if (finds !== null) {
        var firstcap = finds[0];
        var name = (firstcap[0].toLowerCase()
            + firstcap[1].toUpperCase()
            + firstcap.slice(2));
        return [name];
    }
    // * otherwise resort to pfams
    if (!polymer.pfams) {
        return [];
    }
    else {
        var pfamIds = polymer.pfams.reduce((acc, item) => {
            acc.push(item.rcsb_pfam_accession);
            return acc;
        }, []);
        var nomenclature = [];
        pfamIds.map(id => {
            Object.entries(large_subunit_map_1.large_subunit_map).map(entry => {
                if (entry[1].pfamDomainAccession.includes(id)) {
                    nomenclature.push(entry[0]);
                }
            });
            Object.entries(small_subunit_map_1.small_subunit_map).map(entry => {
                if (entry[1].pfamDomainAccession.includes(id)) {
                    nomenclature.push(entry[0]);
                }
            });
        });
        return (0, lodash_1.uniq)(nomenclature);
    }
};
const matchRNANomenclature = (polymer) => {
    const rna_reg = {
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
    };
    var rnatypes = Object.keys(rna_reg);
    for (var c of rnatypes) {
        var matches = rna_reg[c].exec(polymer.rcsb_polymer_entity.pdbx_description);
        if (matches) {
            return [c];
        }
    }
    return [];
};
const inferOrganismsFromPolymers = (proteins) => {
    var host_organism_names = [];
    var src_organism_names = [];
    var host_organism_ids = [];
    var src_organism_ids = [];
    proteins.map(protein => {
        protein.src_organism_names ? src_organism_names.push(...protein.src_organism_names) : null;
        protein.src_organism_ids ? src_organism_ids.push(...protein.src_organism_ids) : null;
        protein.host_organism_names ? src_organism_names.push(...protein.host_organism_names) : null;
        protein.host_organism_ids ? src_organism_ids.push(...protein.host_organism_ids) : null;
    });
    return {
        src_organism_ids: (0, lodash_1.uniq)(src_organism_ids),
        src_organism_names: (0, lodash_1.uniq)(src_organism_names),
        host_organism_ids: (0, lodash_1.uniq)(host_organism_ids),
        host_organism_names: (0, lodash_1.uniq)(host_organism_names),
    };
};
const extractRefs = (external_refs) => {
    var externalRefIds = [];
    var externalRefTypes = [];
    var externalRefLinks = [];
    external_refs === null ? (() => { }) : external_refs.map(ref => {
        externalRefIds.push(ref.id);
        externalRefTypes.push(ref.type);
        externalRefLinks.push(ref.link);
    });
    return [externalRefIds, externalRefTypes, externalRefLinks];
};
const reshape_ToLigand = (nonpoly) => {
    return {
        pdbx_description: nonpoly.rcsb_nonpolymer_entity.pdbx_description,
        formula_weight: nonpoly.rcsb_nonpolymer_entity.formula_weight,
        chemicalId: nonpoly.pdbx_entity_nonpoly.comp_id,
        chemicalName: nonpoly.pdbx_entity_nonpoly.name,
        number_of_instances: nonpoly.rcsb_nonpolymer_entity.pdbx_number_of_molecules
    };
};
const is_ligand_like = (plm, nomenclature) => {
    if (nomenclature.includes('tRNA') || nomenclature.includes('mRNA')) {
        return true;
    }
    // ? Look for enzymes, factors and antibiotics
    var reg = /(\w*(?<!(cha|pro|dom|stra|pl\w*))in\b)|(\b\w*zyme\b)|(factor)/gi;
    const matches = plm.rcsb_polymer_entity.pdbx_description.match(reg);
    if (matches !== null &&
        !(matches.map(_ => _.toLowerCase()).includes('protein')) &&
        !(plm.rcsb_polymer_entity.pdbx_description.toLowerCase().includes('protein')) &&
        !(plm.rcsb_polymer_entity.pdbx_description.toLowerCase().includes('rrna'))) {
        return true;
    }
    return false;
};
const reshape_PolyEntity_to_rRNA = (plm) => {
    var src_organism_ids = plm.rcsb_entity_source_organism ? (0, lodash_1.uniq)(plm.rcsb_entity_source_organism.map(org => org.ncbi_taxonomy_id)) : [];
    var src_organism_names = plm.rcsb_entity_source_organism ? (0, lodash_1.uniq)(plm.rcsb_entity_source_organism.map(org => org.scientific_name)) : [];
    var host_organism_ids = plm.rcsb_entity_host_organism ? (0, lodash_1.uniq)(plm.rcsb_entity_host_organism.map(org => org.ncbi_taxonomy_id)) : [];
    var host_organism_names = plm.rcsb_entity_host_organism ? (0, lodash_1.uniq)(plm.rcsb_entity_host_organism.map(org => org.scientific_name)) : [];
    var nomenclature = matchRNANomenclature(plm);
    return plm.rcsb_polymer_entity_container_identifiers.auth_asym_ids.map(auth_asym_id => ({
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
    }));
};
const reshape_PolyEntity_to_RibosomalProtein = (plm) => {
    if (plm.pfams) {
        var pfam_comments = (0, lodash_1.uniq)(plm.pfams.map(pfam => pfam.rcsb_pfam_comment));
        var pfam_descriptions = (0, lodash_1.uniq)(plm.pfams.map(pfam => pfam.rcsb_pfam_description));
        var pfam_accessions = (0, lodash_1.uniq)(plm.pfams.map(pfam => pfam.rcsb_pfam_accession));
    }
    else {
        var pfam_comments = [];
        var pfam_descriptions = [];
        var pfam_accessions = [];
    }
    var src_organism_ids = (0, lodash_1.uniq)(plm.rcsb_entity_source_organism.map(org => org.ncbi_taxonomy_id));
    var src_organism_names = (0, lodash_1.uniq)(plm.rcsb_entity_source_organism.map(org => org.scientific_name));
    var host_organism_ids = plm.rcsb_entity_host_organism ? (0, lodash_1.uniq)(plm.rcsb_entity_host_organism.map(org => org.ncbi_taxonomy_id)) : [];
    var host_organism_names = plm.rcsb_entity_host_organism ? (0, lodash_1.uniq)(plm.rcsb_entity_host_organism.map(org => org.scientific_name)) : [];
    var nomenclature = matchRPNomenclature(plm);
    return plm.rcsb_polymer_entity_container_identifiers.auth_asym_ids.map(auth_asym_id => ({
        nomenclature: nomenclature,
        asym_ids: plm.rcsb_polymer_entity_container_identifiers.asym_ids,
        auth_asym_id: auth_asym_id,
        parent_rcsb_id: plm.entry.rcsb_id,
        pfam_accessions: pfam_accessions,
        pfam_comments: pfam_comments,
        pfam_descriptions: pfam_descriptions,
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
    }));
};
const processPDBRecord = async (pdbid) => {
    return await axios_1.default.get(query_template(pdbid))
        .then(response => {
        var pdbRecord = response.data.data;
        var proteins = pdbRecord.entry.polymer_entities.filter(poly => poly.entity_poly.rcsb_entity_polymer_type === "Protein");
        var rnas = pdbRecord.entry.polymer_entities.filter(poly => poly.entity_poly.rcsb_entity_polymer_type === "RNA");
        var ligands = pdbRecord.entry.nonpolymer_entities;
        var reshaped_proteins = [];
        var reshaped_rnas = [];
        proteins.map(protein => { reshaped_proteins.push(...reshape_PolyEntity_to_RibosomalProtein(protein)); });
        rnas.map(rna => { reshaped_rnas.push(...reshape_PolyEntity_to_rRNA(rna)); });
        var reshaped_nonpoly = ligands == null ? null : ligands.map(r => reshape_ToLigand(r));
        var organisms = inferOrganismsFromPolymers(reshaped_proteins);
        var externalRefs = extractRefs(pdbRecord.entry.rcsb_external_references);
        var pub = pdbRecord.entry.citation[0];
        var kwords_text = pdbRecord.entry.struct_keywords
            ? pdbRecord.entry.struct_keywords.text
            : null;
        var kwords = pdbRecord.entry.struct_keywords
            ? pdbRecord.entry.struct_keywords.pdbx_keywords
            : null;
        var reshaped = {
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
exports.processPDBRecord = processPDBRecord;
