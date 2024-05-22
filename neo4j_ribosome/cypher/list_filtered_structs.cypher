// match (r:RibosomeStructure) 
// with r
// order by r.rcsb_id desc
// where toLower(r.citation_title) 
// + toLower(r.pdbx_keywords_text) 
// + apoc.text.join(r.citation_rcsb_authors, "")  contains "complex" 
// and rib.citation_year > 2020 
// and rib.resolution < 3
// and ALL(x in ["uL4", "uL22"] where x in apoc.coll.flatten(collect{match (rib)-[]-(p:Polymer) return p.nomenclature }) )
// and ANY(tax in [9606] where tax in apoc.coll.flatten(collect{ match (rib)-[:source]-(p:PhylogenyNode)-[:descendant_of*]-(s:PhylogenyNode) return [p.ncbi_tax_id, s.ncbi_tax_id]}) )

with collect(rib)[..10] as rib, count(rib) as total_count 
unwind rib as ribosomes
optional match (l:Ligand)-[]-(ribosomes) 
with collect(PROPERTIES(l)) as ligands, ribosomes, total_count

match (rps:Protein)-[]-(ribosomes) 
with collect(PROPERTIES(rps)) as proteins, ligands, ribosomes, total_count

optional match (rna:RNA)-[]-(ribosomes) 
with collect(PROPERTIES(rna)) as rnas, proteins, ligands, ribosomes, total_count

with apoc.map.mergeList([{proteins:proteins},{nonpolymeric_ligands:ligands},{rnas:rnas},{other_polymers:[]}]) as rest, ribosomes, total_count
return collect(apoc.map.merge(ribosomes, rest)),  collect(distinct total_count)[0]