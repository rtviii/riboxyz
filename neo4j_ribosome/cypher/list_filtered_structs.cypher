match (rib:RibosomeStructure) 
with rib
order by rib.rcsb_id desc
where toLower(rib.citation_title) 
+ toLower(rib.pdbx_keywords_text) 
+ apoc.text.join(rib.citation_rcsb_authors, "")  contains "complex" 
and rib.citation_year > 2020 
and rib.resolution < 3
and ALL(x in ["uL4", "uL22"] where x in apoc.coll.flatten(collect{match (rib)-[]-(p:Polymer) return p.nomenclature }) )
and ANY(tax in [9606] where tax in apoc.coll.flatten(collect{ match (rib)-[:source]-(p:PhylogenyNode)-[:descendant_of*]-(s:PhylogenyNode) return [p.ncbi_tax_id, s.ncbi_tax_id]}) )
with collect(rib)[..10] as rib, count(rib) as total_count 

unwind rib as ribs
optional match (l:Ligand)-[]-(ribs) 
with collect(PROPERTIES(l)) as ligands, ribs, total_count

match (rps:Protein)-[]-(ribs) 
with collect(PROPERTIES(rps)) as proteins, ligands, ribs, total_count

optional match (rna:RNA)-[]-(ribs) 
with collect(PROPERTIES(rna)) as rnas, proteins, ligands, ribs, total_count

with apoc.map.mergeList([{proteins:proteins},{nonpolymeric_ligands:ligands},{rnas:rnas},{other_polymers:[]}]) as rest, ribs, total_count
return collect(apoc.map.merge(ribs, rest)),  collect(distinct total_count)[0]