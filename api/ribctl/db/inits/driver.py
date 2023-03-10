from ast import List
from pprint import pprint
from typing import  Self
from neo4j import Driver, GraphDatabase, Result, Transaction
from pyparsing import Any
from ribctl.db.inits.proteins import add_protein, node__protein_class
from ribctl.db.inits.rna import add_rna, node__rna_class
from ribctl.db.inits.structure import add_ligand, node__structure
from ribctl.lib.types.types_ribosome import RibosomeStructure
from ribctl.lib.types.types_ribosome_assets import RibosomeAssets
from ribctl.lib.types.types_polymer import list_LSU_Proteins, list_SSU_Proteins, list_RNAClass

"""Functions of the form create_node_xxxx return a closure over their target [xxxx] because the neo4j expects a 'unit-of-work'
function that takes a transaction as its argument."""

# ※ ----------------[ 0.Database  inits: constraints & nomenclature classes]
NODE_CONSTRAINTS = [
    """CREATE CONSTRAINT IF NOT EXISTS ON (ipro:InterProFamily) ASSERT ipro.family_id  IS UNIQUE;""",
    """CREATE CONSTRAINT IF NOT EXISTS ON (go:GOClass) ASSERT go.class_id IS UNIQUE;""",
    """CREATE CONSTRAINT IF NOT EXISTS ON (q:RibosomeStructure) Assert q.rcsb_id IS UNIQUE;""",
    """CREATE CONSTRAINT IF NOT EXISTS ON (pf:PFAMFamily) assert pf.family_id  is unique;""",
    """CREATE CONSTRAINT IF NOT EXISTS ON (lig:Ligand) assert lig.chemicalId is unique;""",
    """CREATE CONSTRAINT IF NOT EXISTS ON (nc:NomenclatureClass) assert nc.class_id is unique;""",
]

# ※ ----------------[ 1.Structure Nodes]
# ※ ----------------[ 2.RNA Nodes]
# ※ ----------------[ 3.Protein Nodes]
# ※ ----------------[ 4.Ligand Nodes]
# ※ ----------------[ 5.Ingress]


def init_driver(uri: str = "neo4j://localhost:7687", username: str = "neo4j", password="neo4j"):
    # TODO: Create an instance of the driver here
    api = GraphDatabase.driver(uri, auth=(username, password))
    return api


driver = init_driver("neo4j://localhost:7687", "neo4j", "neo4j")


class Neo4jDB():

    driver: Driver

    def __init__(self, **kwargs) -> None:
        self.driver = init_driver(**kwargs)

    def see_constraints(self) -> list[dict[str, Any]]:
        with self.driver.session() as s:
            r = s.run("""//
            CALL db.constraints;
                  """)
            return r.data()
    def get_all_structs(self) :
        with self.driver.session() as s:
            struct_ids = []
            [struct_ids.extend(struct) for struct in  s.read_transaction(lambda tx: tx.run("""//
            match (n:RibosomeStructure) return n.rcsb_id;
            """).values('n.rcsb_id'))]
            return struct_ids

    def get_any(self) -> list[dict[str, Any]]:
        with self.driver.session() as s:
            return s.read_transaction(lambda tx: tx.run("""//
            match (n) return n limit 10;
            """).data())

    def init_db(self):
        self.__init_constraints()
        self.__init_protein_classes()
        self.__init_rna_classes()

    def __init_constraints(self) -> None:
        with self.driver.session() as s:
            for c in NODE_CONSTRAINTS:
                s.write_transaction(lambda tx: tx.run(c))

    def __init_protein_classes(self):
        with self.driver.session() as s:
            for protein_class in [*list_LSU_Proteins, *  list_SSU_Proteins]:
                s.execute_write(node__protein_class(protein_class))

    def __init_rna_classes(self,):
        with self.driver.session() as s:
            for rna_class in list_RNAClass:
                s.execute_write(node__rna_class(rna_class))

    def add_structure(self, struct_assets: RibosomeAssets):

        R = RibosomeStructure(**struct_assets.json_profile())

        global struct_node_result;
        with self.driver.session() as s:
           struct_node_result = s.write_transaction(node__structure(R))


        for protein in R.proteins:
            add_protein(self.driver, protein)

        if R.rnas is not None:
            for rna in R.rnas:
                add_rna(self.driver, rna)

        if R.ligands is not None:
            for ligand in R.ligands:
                add_ligand(self.driver, ligand, R.rcsb_id)

        print("Added structure: ", R.rcsb_id)
        return struct_node_result.data()

