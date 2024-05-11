import sys

from ribctl.lib.libtax import PhylogenyNode, Taxid


sys.dont_write_bytecode = True

from neo4j_adapter.node_phylogeny import link__phylogeny, node__phylogeny
from neo4j.graph import Node
from concurrent.futures import ALL_COMPLETED, Future, ThreadPoolExecutor, wait
from neo4j import Driver, GraphDatabase
from neo4j_adapter.node_polymer import  link__polymer_to_polymer_class, link__polymer_to_structure, node__polymer, upsert_polymer_to_protein, upsert_polymer_to_rna
from neo4j_adapter.node_protein import  node__polymer_class
from neo4j_adapter.node_structure import   link__ligand_to_struct, node__ligand, node__structure, struct_exists
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

class Neo4jAdapter():

    driver   : Driver
    uri      : str
    user     : str
    databases: list[str]

    def init_constraints(self) -> None:
        with self.driver.session() as session:
            for c in NODE_CONSTRAINTS:
                session.execute_write(lambda tx: tx.run(c))

    def init_polymer_classes(self):
        with self.driver.session() as session:
            for polymer_class in [*list(PolymerClass)]:
                session.execute_write(node__polymer_class(polymer_class.value))

    def initialize_new_instance(self):

        self.init_constraints()
        self.init_polymer_classes()
        self.init_phylogenies()

    def __init__(self, uri: str, user: str, password: str|None=None) -> None:
        self.uri      = uri
        self.user     = user
        self.password = password 

        try:
            self.driver = GraphDatabase.driver(uri, auth=(user, password))
            print("Established connection to ", self.uri)
        except Exception as ae:
            print(ae)

    def init_phylogenies(self):
        taxa = RibosomeAssets.collect_all_taxa()
        for taxon in taxa:
            self.create_lineage(taxon.ncbi_tax_id)
        

    def add_structure(self, rcsb_id:str):

        rcsb_id = rcsb_id.upper()
        if self.check_structure_exists(rcsb_id):
            print("Struct node {} already exists.".format(rcsb_id))
            return

        R:RibosomeStructure = RibosomeAssets(rcsb_id).profile()

        with self.driver.session() as s:
            structure_node = s.execute_write(node__structure(R))

            for protein in R.proteins:
                   protein_node = s.execute_write(node__polymer                 (protein                        ))
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
    
        return structure_node

    def add_phylogeny_node(self, taxid:int)->Node:
        with self.driver.session() as session:
            node = session.execute_write(node__phylogeny(PhylogenyNode.from_taxid(taxid)))
            return node

    def create_lineage(self,taxid:int)->None:
        lin = Taxid.get_lineage(taxid)
        lin.reverse()
        previous_id: int|None = None
        with self.driver.session() as session:
            for taxid in lin:
                node = session.execute_write(node__phylogeny(PhylogenyNode.from_taxid(taxid)))
                print("Created node {} with taxid {}".format(node, taxid))
                if previous_id == None: # initial (superkingdom has no parent node)
                    previous_id = taxid
                    continue

                session.execute_write(link__phylogeny( taxid , previous_id)) # link current tax to parent
                previous_id = taxid
        print("Created lineage: ", lin)
        return



    # ------------------- OLD SHIT --------------------

    def see_current_auth(self):
        with self.driver.session(database='system') as s:
            r = s.run("""show users""")
            users_array = r.data()
            return users_array

    def sync_with_rcsb(self, workers:int)->None:

        synced   = self.get_all_structs()
        unsynced = sorted(current_rcsb_structs())
        futures:list[Future] =  []
        print("Syncing over the following structs:", unsynced)

        with ThreadPoolExecutor(max_workers=workers) as executor:
            for rcsb_id in list(set(unsynced ) - set(synced)):
                fut    = executor.submit(self.add_structure, rcsb_id)
                futures.append(fut)

        print(wait(futures, return_when=ALL_COMPLETED))

    def get_all_structs(self):
        with self.driver.session() as s:
            struct_ids = []
            [struct_ids.extend(struct) for struct in s.execute_read(lambda tx: tx.run("""//
            match (n:RibosomeStructure) return n.rcsb_id;
            """).values('n.rcsb_id'))]
            return struct_ids

    def get_individual_ligand(self, chemId: str) -> NonpolymericLigand:
        with self.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                return tx.run("""match (ligand:Ligand{chemicalId: $CHEM_ID}) return ligand""", {"CHEM_ID": chemId}).data()[0]['ligand']
            return session.execute_read(_)

    def check_structure_exists(self, rcsb_id:str)->bool:
        rcsb_id = rcsb_id.upper()
        with self.driver.session() as session:
            return session.execute_read(struct_exists(rcsb_id))

    def get_RibosomeStructure(self, rcsb_id: str) -> RibosomeStructure:
        with self.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                return tx.run("""//
            match (n:RibosomeStructure {rcsb_id:$RCSB_ID})

            optional match (rr:RNA)-[]-(n)
            with n, collect(rr) as rrna

            optional match (rp:Protein)-[]-(n)
            with n, rrna,  collect(rp) as rps

            optional match (l:Ligand)-[]-(n)
            with n, rrna, rps, collect(l) as ligs

            with {
                        rcsb_id               : n.rcsb_id               ,
                        expMethod             : n.expMethod             ,
                        resolution            : n.resolution            ,

                        pdbx_keywords         : n.pdbx_keywords         ,
                        pdbx_keywords_text    : n.pdbx_keywords_text    ,

                        rcsb_external_ref_id  : n.rcsb_external_ref_id  ,
                        rcsb_external_ref_type: n.rcsb_external_ref_type,
                        rcsb_external_ref_link: n.rcsb_external_ref_link,

                        citation_year         : n.citation_year         ,
                        citation_rcsb_authors : n.citation_rcsb_authors ,
                        citation_title        : n.citation_title        ,
                        citation_pdbx_doi     : n.citation_pdbx_doi     ,

                        src_organism_ids      : n.src_organism_ids      ,
                        src_organism_names    : n.src_organism_names    ,

                        host_organism_ids     : n.host_organism_ids     ,
                        host_organism_names   : n.host_organism_names   ,

                        proteins : rps,
                        rnas     : rrna,
                        ligands  : ligs
                    } as structure

                    return structure
                        """, {"RCSB_ID": rcsb_id.upper()}).data()[0]['structure']
            return session.execute_read(_)

    def match_structs_w_proteins(self, targets: list[CytosolicProteinClass]) -> list[str]:
        with self.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                return [s[0] for s in tx.run("""//
                    match (n:RibosomeStructure)-[]-(rp:Protein)
                    with n, rp, [] as strnoms 
                    unwind rp.nomenclature as unwound
                    with collect(unwound) as unwound, n, $T as tgts
                    where all(x in tgts where x in unwound)
                    with n.rcsb_id as rcsb_id
                    return rcsb_id
                        """, {"T": targets}).values('rcsb_id')]

            return session.execute_read(_)

    def proteins_number(self) -> int:
        with self.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                return tx.run("""//
                match (n:Protein) return count(n)
                    """).value()[0]
            return session.execute_read(_)

    def number_of_structures(self) -> int:
        with self.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                return tx.run("""//
                match (n:RibosomeStructure) return count(n)
                    """).value()[0]
            return session.execute_read(_)



    def get_banclass_for_chain(self, rcsb_id: str, auth_asym_id) -> list[CytosolicProteinClass]:
        # TODO: This method should handle both protein and RNA chains(atm only proteins)
        with self.driver.session() as session:

            def _(tx: Transaction | ManagedTransaction):
                return [pc[0] for pc in tx.run("""//
                match (n:RibosomeStructure {rcsb_id:$RCSB_ID})-[]-(c:Protein{auth_asym_id:$AUTH_ASYM_ID})-[]-(pc:ProteinClass) return pc.class_id
                    """, {"RCSB_ID": rcsb_id, "AUTH_ASYM_ID": auth_asym_id}).values('pc.class_id')]
            return session.execute_read(_)



