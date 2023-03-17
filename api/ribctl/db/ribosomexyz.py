from ast import List
from pprint import pprint
from venv import logger
from neo4j import Driver, GraphDatabase, Result, Transaction
from neo4j.exceptions import ServiceUnavailable, ClientError
from pyparsing import Any
from ribctl.lib.struct_rcsb_api import current_rcsb_structs
from ribctl.db.inits.proteins import add_protein, node__protein_class
from ribctl.db.inits.rna import add_rna, node__rna_class
from ribctl.db.inits.structure import add_ligand, node__structure
from ribctl.lib.types.types_ribosome import RibosomeStructure
from ribctl.lib.types.types_ribosome_assets import RibosomeAssets
from ribctl.lib.types.types_polymer import list_LSU_Proteins, list_SSU_Proteins, list_RNAClass
from concurrent.futures import Future, ProcessPoolExecutor

"""Functions of the form create_node_xxxx return a closure over their target [xxxx] because the neo4j expects a 'unit-of-work'
function that takes a transaction as its argument."""

# ※ ----------------[ 0.Database  inits: constraints & nomenclature classes]
NODE_CONSTRAINTS = [
    """CREATE CONSTRAINT IF NOT EXISTS ON (ipro:InterProFamily) ASSERT ipro.family_id  IS UNIQUE;""",
    """CREATE CONSTRAINT IF NOT EXISTS ON (go:GOClass) ASSERT go.class_id IS UNIQUE;""",
    """CREATE CONSTRAINT IF NOT EXISTS ON (q:RibosomeStructure) Assert q.rcsb_id IS UNIQUE;""",
    """CREATE CONSTRAINT IF NOT EXISTS ON (lig:Ligand) assert lig.chemicalId is unique;""",
    """CREATE CONSTRAINT IF NOT EXISTS ON (nc:NomenclatureClass) assert nc.class_id is unique;""",
]

# ※ ----------------[ 1.Structure Nodes]
# ※ ----------------[ 2.RNA Nodes]
# ※ ----------------[ 3.Protein Nodes]
# ※ ----------------[ 4.Ligand Nodes]
# ※ ----------------[ 5.Ingress]



# If you are connecting via a shell or programmatically via a driver,
# just issue a `ALTER CURRENT USER SET PASSWORD FROM 'current password' TO 'new password'` statement against
# the system database in the current session, and then restart your driver with the new password configured.

class riboxyzDB():

    driver: Driver

    # def sync_with_rcsb(self)->None:

    #     def _():
    #         D        = riboxyzDB()
    #         synced   = D.get_all_structs()
    #         unsynced = sorted(current_rcsb_structs())

    #         for rcsb_id in set(unsynced ) - set(synced):
    #             assets = RibosomeAssets(rcsb_id)
    #             try:
    #                 assets._verify_json_profile(True)
    #                 D.add_structure(assets)
    #             except Exception as e:
    #                 print(e)
    #                 logger.error("Exception occurred:", exc_info=True)

    #     with ProcessPoolExecutor() as executor:
    #         future:Future = executor.submit(sync_with_rcsb)
    #         print(future)

    def __init__(self,uri:str,username:str, password:str) -> None:
        self.driver = GraphDatabase.driver(uri, auth=(username, password))
    
    def change_default_neo4j_credentials(self):
        with self.driver.session(database='system') as system_s:
            r = system_s.run("""//
ALTER CURRENT USER SET PASSWORD FROM "neo4j" TO "ribosomexyz";
                  """)
            return r.data()

    def see_current_auth(self):
        with self.driver.session(database='system') as s:
            r = s.run("""show users""")
            # {
            #   "user": "neo4j",
            #   "roles": null,
            #   "passwordChangeRequired": false,
            #   "suspended": null,
            #   "home": null
            # }
            users_array = r.data()
            if users_array[0]["passwordChangeRequired"] == True:
                self.change_default_neo4j_credentials()
            return users_array



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

