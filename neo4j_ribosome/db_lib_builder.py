import sys

from ribctl.global_ops import GlobalOps

sys.dont_write_bytecode = True

from neo4j_ribosome.node_ligand import link__ligand_to_struct, node__ligand
from ribctl.lib.libtax import PhylogenyNode, Taxid
from neo4j_ribosome.node_phylogeny import link__phylogeny, node__phylogeny
from neo4j.graph import Node
from neo4j import Driver, GraphDatabase
from neo4j_ribosome.node_polymer import (
    link__polymer_to_polymer_class,
    link__polymer_to_structure,
    node__polymer,
    upsert_polymer_to_protein,
    upsert_polymer_to_rna,
    node__polymer_class,
)
from neo4j_ribosome.node_structure import (
    link__structure_to_lineage_member,
    link__structure_to_organism,
    node__structure,
    struct_exists,
)
from ribctl.lib.schema.types_ribosome import (
    RibosomeStructure,
    RibosomeStructureMetadata,
)
from ribctl.lib.types.polymer.hierarchies import Polynucleotides, Polypeptides, Polymers
from ribctl.lib.types.polymer import     MitochondrialProteinClass, PolynucleotideClass, PolynucleotideClass
from ribctl.ribosome_ops import StructureAssets, RibosomeOps, Structure
from neo4j import GraphDatabase, Driver, ManagedTransaction, Transaction
from ribctl.lib.schema.types_ribosome import (
    NonpolymericLigand,
)

NODE_CONSTRAINTS = [
    """CREATE CONSTRAINT rcsb_id_unique IF NOT EXISTS FOR (ribosome:RibosomeStructure) REQUIRE ribosome.rcsb_id IS UNIQUE;""",
    """CREATE CONSTRAINT polymer_class_unique IF NOT EXISTS FOR (poly_class:PolymerClass) REQUIRE poly_class.class_id IS UNIQUE;""",
    """CREATE CONSTRAINT taxid_unique IF NOT EXISTS FOR (phylonode:PhylogenyNode) REQUIRE phylonode.ncbi_tax_id IS UNIQUE;""",
    """CREATE CONSTRAINT chemicalId IF NOT EXISTS FOR (ligand:Ligand) REQUIRE ligand.chemicalId IS UNIQUE;""",
]

# If you are connecting via a shell or programmatically via a driver,
# just issue a `ALTER CURRENT USER SET PASSWORD FROM 'current password' TO 'new password'` statement against
# the system database in the current session, and then restart your driver with the new password configured.

class Neo4jAdapter:

    driver: Driver
    uri: str
    user: str
    databases: list[str]

    def __init__(
        self,
        uri: str,
        user: str,
        current_db: str,
        password: str | None = None,
    ) -> None:

        self.uri = uri
        self.user = user
        self.password = password

        try:
            self.driver = GraphDatabase.driver(
                uri, auth=(user, password), database=current_db
            )
        except Exception as ae:
            print(ae)

    def initialize_new_instance(self):
        self.init_constraints()
        self.init_polymer_classes()
        self.init_phylogenies()

    def init_constraints(self) -> None:
        with self.driver.session() as session:
            for c in NODE_CONSTRAINTS:
                session.execute_write(lambda tx: tx.run(c))
                print("\nAdded constraint: ", c)

    def init_polymer_classes(self):
        with self.driver.session() as session:
            for polymer_class in [*list(Polymers)]:
                session.execute_write(node__polymer_class(polymer_class.value))
            print("\nAdded polymer classes: ", [*list(Polymers)])

    def add_phylogeny_node(self, taxid: int) -> Node:
        with self.driver.session() as session:
            node = session.execute_write(
                node__phylogeny(PhylogenyNode.from_taxid(taxid))
            )
            return node

    def init_phylogenies(self):
        taxa = GlobalOps.collect_all_taxa()
        for taxon in taxa:
            self._create_lineage(taxon.ncbi_tax_id)

    def _create_lineage(self, taxid: int) -> None:
        lin = Taxid.get_lineage(taxid)
        lin.reverse()
        previous_id: int | None = None

        with self.driver.session() as session:
            for taxid in lin:
                node = session.execute_write(
                    node__phylogeny(PhylogenyNode.from_taxid(taxid))
                )
                if previous_id == None:  # initial (superkingdom has no parent node)
                    previous_id = taxid
                    continue

                session.execute_write(
                    link__phylogeny(taxid, previous_id)
                )  # link current tax to parent
                previous_id = taxid
        return

    def link_structure_to_phylogeny(
        self, rcsb_id: str, profile: RibosomeStructure | None = None, verbose: bool = False
    ):
        rcsb_id = rcsb_id.upper()

        if profile is None:
            profile = RibosomeOps(rcsb_id).profile

        _ = []
        with self.driver.session() as s:
            for organism_host in profile.host_organism_ids:

                self._create_lineage(organism_host)
                s.execute_write(
                    link__structure_to_organism(rcsb_id, organism_host, "host")
                )
                lineage_memebers_host = Taxid.get_lineage(organism_host)
                for org in lineage_memebers_host:
                    s.execute_write(
                        link__structure_to_lineage_member(
                            rcsb_id, org, "belongs_to_lineage_host"
                        )
                    )
                _.extend(lineage_memebers_host)

            for organism_src in profile.src_organism_ids:
                self._create_lineage(organism_src)
                s.execute_write(
                    link__structure_to_organism(rcsb_id, organism_src, "source")
                )
                lineage_memebers_source = Taxid.get_lineage(organism_src)
                for org in lineage_memebers_source:
                    s.execute_write(
                        link__structure_to_lineage_member(
                            rcsb_id, org, "belongs_to_lineage_source"
                        )
                    )
                _.extend(lineage_memebers_source)

    def check_structure_exists(self, rcsb_id: str) -> bool:
        rcsb_id = rcsb_id.upper()
        with self.driver.session() as session:
            return session.execute_read(struct_exists(rcsb_id))

    def upsert_ligand_node(
        self, ligand: NonpolymericLigand, parent_rcsb_id: str | None = None
    ):
        with self.driver.session() as s:
            ligand_node = s.execute_write(node__ligand(ligand))
            if parent_rcsb_id is not None:
                s.execute_write(link__ligand_to_struct(ligand_node, parent_rcsb_id))

    def add_total_structure(self, rcsb_id: str, disable_exists_check: bool = False):
        rcsb_id = rcsb_id.upper()

        if not disable_exists_check:
            if self.check_structure_exists(rcsb_id):
                print("\nStruct node {} already exists.".format(rcsb_id))
                return
    
        R: RibosomeStructure = RibosomeOps(rcsb_id).profile

        with self.driver.session() as s:

            structure_node = s.execute_write(node__structure(R))
            self.link_structure_to_phylogeny(rcsb_id, R)

            for protein in R.proteins:
                protein_node = s.execute_write(node__polymer(protein))

                s.execute_write(
                    link__polymer_to_structure(protein_node, protein.parent_rcsb_id)
                )
                s.execute_write(link__polymer_to_polymer_class(protein_node))
                s.execute_write(upsert_polymer_to_protein(protein_node, protein))

            if R.rnas is not None:
                for rna in R.rnas:
                    rna_node = s.execute_write(node__polymer(rna))
                    s.execute_write(
                        link__polymer_to_structure(rna_node, rna.parent_rcsb_id)
                    )
                    s.execute_write(link__polymer_to_polymer_class(rna_node))
                    s.execute_write(upsert_polymer_to_rna(rna_node, rna))

            if R.other_polymers is not None:
                for polymer in R.other_polymers:
                    other_poly_node = s.execute_write(node__polymer(polymer))
                    s.execute_write(
                        link__polymer_to_structure(
                            other_poly_node, polymer.parent_rcsb_id
                        )
                    )
                    s.execute_write(link__polymer_to_polymer_class(other_poly_node))

            if R.nonpolymeric_ligands is not None:
                for ligand in R.nonpolymeric_ligands:
                    ligand_node = s.execute_write(node__ligand(ligand))
                    s.execute_write(link__ligand_to_struct(ligand_node, R.rcsb_id))

        print("Successfully initialized structure {}.".format(rcsb_id))
        return structure_node

    def upsert_structure_node(self, rcsb_id: str):
        rcsb_id = rcsb_id.upper()
        R: RibosomeStructure = RibosomeOps(rcsb_id).profile
        with self.driver.session() as s:
            structure_node = s.execute_write(node__structure(R))
        print("Successfully merged structure {}.".format(rcsb_id))
        return structure_node

    def delete_structure(self, rcsb_id: str, dry_run: bool = False) -> dict[str, int]:
        """
        Delete a structure and its polymers, preserving ligand nodes while removing their relationships.
        
        Args:
            rcsb_id: The RCSB ID of the structure to delete
            dry_run: If True, only count nodes/relationships that would be deleted without actually deleting
            
        Returns:
            Dictionary with counts of deleted/affected nodes and relationships
        """
        rcsb_id = rcsb_id.upper()
        
        # Query to count and optionally delete the structure and polymers
        query = """
        MATCH (s:RibosomeStructure {rcsb_id: $rcsb_id})
        
        // Count all connected nodes and relationships
        OPTIONAL MATCH (s)-[r1]-(p:Polymer)
        OPTIONAL MATCH (s)-[r2]-(l:Ligand)
        OPTIONAL MATCH (s)-[r3]-(t:PhylogenyNode)
        
        WITH s,
             collect(DISTINCT p) as polymers,
             collect(DISTINCT r1) as polymerRels,
             collect(DISTINCT r2) as ligandRels,
             collect(DISTINCT r3) as taxaRels,
             count(DISTINCT p) as polymerCount,
             count(DISTINCT r1) + count(DISTINCT r2) + count(DISTINCT r3) as relCount
        
        // Delete relationships to ligands (preserve nodes)
        FOREACH (r IN ligandRels | DELETE r)
        
        // Delete relationships to taxa (preserve nodes)
        FOREACH (r IN taxaRels | DELETE r)
        
        // Delete polymers and their relationships
        FOREACH (p IN polymers | DETACH DELETE p)
        
        // Finally delete the structure
        DETACH DELETE s
        
        RETURN {
            structure: 1,
            polymers: polymerCount,
            relationships: relCount
        } as counts
        """
        
        # Only count query for dry run
        dry_run_query = """
        MATCH (s:RibosomeStructure {rcsb_id: $rcsb_id})
        
        OPTIONAL MATCH (s)-[r1]-(p:Polymer)
        OPTIONAL MATCH (s)-[r2]-(l:Ligand)
        OPTIONAL MATCH (s)-[r3]-(t:PhylogenyNode)
        
        RETURN {
            structure: 1,
            polymers: count(DISTINCT p),
            relationships: count(DISTINCT r1) + count(DISTINCT r2) + count(DISTINCT r3)
        } as counts
        """
        
        with self.driver.session() as session:
            # First check if structure exists
            if not self.check_structure_exists(rcsb_id):
                raise ValueError(f"Structure {rcsb_id} does not exist in database")
                
            try:
                # Execute appropriate query based on dry_run flag
                result = session.execute_write(
                    lambda tx: tx.run(
                        dry_run_query if dry_run else query,
                        rcsb_id=rcsb_id
                    ).single()
                )
                
                counts = result['counts'] if result else {
                    'structure': 0,
                    'polymers': 0,
                    'relationships': 0
                }
                
                action = "Would delete" if dry_run else "Deleted"
                print(f"\n{action} for {rcsb_id}:")
                print(f"- Structure node: {counts['structure']}")
                print(f"- Polymer nodes: {counts['polymers']}")
                print(f"- Relationships: {counts['relationships']}")
                print("(Ligand nodes preserved)")
                
                return counts
                
            except Exception as e:
                print(f"Error deleting structure {rcsb_id}: {str(e)}")
                raise