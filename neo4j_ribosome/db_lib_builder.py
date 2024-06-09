import sys
sys.dont_write_bytecode = True
from ribctl.lib.libtax import PhylogenyNode, Taxid
from neo4j_ribosome.node_phylogeny import link__phylogeny, node__phylogeny
from neo4j.graph import Node
from neo4j import Driver, GraphDatabase
from neo4j_ribosome.node_polymer import  link__polymer_to_polymer_class, link__polymer_to_structure, node__polymer, upsert_polymer_to_protein, upsert_polymer_to_rna,node__polymer_class
from neo4j_ribosome.node_structure import   link__ligand_to_struct, link__structure_to_phylogeny, node__ligand, node__structure, struct_exists
from ribctl.etl.etl_pipeline import current_rcsb_structs
from ribctl.lib.schema.types_ribosome import MitochondrialProteinClass, PolymerClass, PolynucleotideClass, RibosomeStructure
from ribctl.etl.ribosome_assets import RibosomeAssets
from neo4j import GraphDatabase, Driver, ManagedTransaction, Transaction
from ribctl.lib.schema.types_ribosome import  NonpolymericLigand,  CytosolicProteinClass, RibosomeStructure

NODE_CONSTRAINTS = [
    """CREATE CONSTRAINT rcsb_id_unique IF NOT EXISTS FOR (ribosome:RibosomeStructure) REQUIRE ribosome.rcsb_id IS UNIQUE;""",
    """CREATE CONSTRAINT polymer_class_unique IF NOT EXISTS FOR (poly_class:PolymerClass) REQUIRE poly_class.class_id IS UNIQUE;""",
    """CREATE CONSTRAINT taxid_unique IF NOT EXISTS FOR (phylonode:PhylogenyNode) REQUIRE phylonode.ncbi_tax_id IS UNIQUE;""",
]

# If you are connecting via a shell or programmatically via a driver,
# just issue a `ALTER CURRENT USER SET PASSWORD FROM 'current password' TO 'new password'` statement against
# the system database in the current session, and then restart your driver with the new password configured.

class Neo4jBuilder():

    driver   : Driver
    uri      : str
    user     : str
    databases: list[str]


    def __init__(self, uri: str, user: str, password: str|None=None) -> None:
        self.uri      = uri
        self.user     = user
        self.password = password 

        try:
            self.driver = GraphDatabase.driver(uri, auth=(user, password))
            print("Established connection to ", self.uri)
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
                print("Added constraint: ", c)

    def init_polymer_classes(self):
        with self.driver.session() as session:
            for polymer_class in [*list(PolymerClass)]:
                session.execute_write(node__polymer_class(polymer_class.value))
                print("Added polymer class: ", polymer_class.value)

    def add_phylogeny_node(self, taxid:int)->Node:
        with self.driver.session() as session:
            node = session.execute_write(node__phylogeny(PhylogenyNode.from_taxid(taxid)))
            return node

    def init_phylogenies(self):
        taxa = RibosomeAssets.collect_all_taxa()
        for taxon in taxa:
            self._create_lineage(taxon.ncbi_tax_id)

    def _create_lineage(self,taxid:int)->None:
        lin = Taxid.get_lineage(taxid)
        print("lineage: " ,lin)
        lin.reverse()
        previous_id: int|None = None

        with self.driver.session() as session:
            for taxid in lin:
                print("creating node for " , taxid)
                node = session.execute_write(node__phylogeny(PhylogenyNode.from_taxid(taxid)))
                if previous_id == None: # initial (superkingdom has no parent node)
                    previous_id = taxid
                    continue

                session.execute_write(link__phylogeny( taxid , previous_id)) # link current tax to parent
                previous_id = taxid
        return
        
    def link_structure_to_phylogeny(self,rcsb_id:str):

        rcsb_id = rcsb_id.upper()
        R:RibosomeStructure = RibosomeAssets(rcsb_id).profile()

        with self.driver.session() as s:
            if self.check_structure_exists(rcsb_id):
                print("Struct node {} already exists.".format(rcsb_id))
                ...
            else:
                s.execute_write(node__structure(R))

            for organism_host in R.host_organism_ids:
                s.execute_write(link__structure_to_phylogeny(rcsb_id, organism_host, 'host_organism'))
            for organism_src in R.src_organism_ids:
                s.execute_write(link__structure_to_phylogeny(rcsb_id, organism_src, 'source_organism'))
        print("Linked structure {} to phylogeny".format(rcsb_id))

    def check_structure_exists(self, rcsb_id:str)->bool:
        rcsb_id = rcsb_id.upper()
        with self.driver.session() as session:
            return session.execute_read(struct_exists(rcsb_id))

    def upsert_ligand_node(self, ligand:NonpolymericLigand, parent_rcsb_id:str):

        with self.driver.session() as s:
            ligand_node = s.execute_write(node__ligand(ligand))
            s.execute_write(link__ligand_to_struct(ligand_node, parent_rcsb_id))

    def add_structure(self, rcsb_id:str, disable_exists_check:bool=False):
        rcsb_id = rcsb_id.upper()
        if not disable_exists_check and self.check_structure_exists(rcsb_id):
            print("Struct node {} already exists.".format(rcsb_id))
            return

        R:RibosomeStructure = RibosomeAssets(rcsb_id).profile()

        with self.driver.session() as s:
            structure_node = s.execute_write(node__structure(R))

            for organism in R.host_organism_ids:
                self._create_lineage(organism)
                s.execute_write(link__structure_to_phylogeny(rcsb_id, organism, 'host_organism'))
            for organism in R.src_organism_ids:
                self._create_lineage(organism)
                s.execute_write(link__structure_to_phylogeny(rcsb_id, organism, 'source_organism'))

            for protein in R.proteins:
                   protein_node = s.execute_write(node__polymer(protein))
                   s.execute_write(link__polymer_to_structure    (protein_node   , protein.parent_rcsb_id))
                   s.execute_write(link__polymer_to_polymer_class(protein_node                           ))
                   s.execute_write(upsert_polymer_to_protein     (protein_node   , protein               ))

            if R.rnas is not None:
                for rna in R.rnas:
                    rna_node = s.execute_write(node__polymer(rna))
                    s.execute_write(link__polymer_to_structure(rna_node, rna.parent_rcsb_id))
                    s.execute_write(link__polymer_to_polymer_class(rna_node))
                    s.execute_write(upsert_polymer_to_rna(rna_node, rna))

            if R.other_polymers is not None:
                for polymer in R.other_polymers:
                    other_poly_node = s.execute_write(node__polymer(polymer))
                    s.execute_write(link__polymer_to_structure(other_poly_node, polymer.parent_rcsb_id))
                    s.execute_write(link__polymer_to_polymer_class(other_poly_node))

            if R.nonpolymeric_ligands is not None:
                for ligand in R.nonpolymeric_ligands:
                    ligand_node = s.execute_write(node__ligand(ligand))
                    s.execute_write(link__ligand_to_struct(ligand_node, R.rcsb_id))
    
        print("Successfully initialized structure {}.".format(rcsb_id))
        return structure_node



    # ------------------- OLD SHIT --------------------



