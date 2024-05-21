match (rib:RibosomeStructure) 
with rib 
order by rib.rcsb_id desc 
limit 20

return rib.rcsb_id


// optional match (l:Ligand)-[]-(rib) 
// with collect(PROPERTIES(l)) as ligands, rib

// match (rps:Protein)-[]-(rib) 
// with collect(PROPERTIES(rps)) as proteins, ligands, rib

// optional match (rna:RNA)-[]-(rib) 
// with collect(PROPERTIES(rna)) as rnas, proteins, ligands, rib

// with  apoc.map.mergeList([{proteins:proteins},{nonpolymeric_ligands:ligands},{rnas:rnas},{other_polymers:[]}]) as rest, rib
// return apoc.map.merge(rib, rest)