from neo4j import  Driver, ManagedTransaction, Transaction
from api.ribctl.db.proteins import node__protein_class
from api.ribctl.db.rna import node__rna_class
from ribctl.lib.types.types_polymer import list_LSU_Proteins, list_SSU_Proteins, list_RNAClass


def init_protein_classes(driver:Driver):
    with driver.session() as s:
        for protein_class in [*list_LSU_Proteins,*  list_SSU_Proteins]:
            s.execute_write(node__protein_class(protein_class))

def init_rna_classes(driver:Driver):
    with driver.session() as s:
        for rna_class in list_RNAClass:
            s.execute_write(node__rna_class(rna_class))

def create_constraints(tx: Transaction | ManagedTransaction):
    tx.run("""//
    CREATE CONSTRAINT IF NOT EXISTS ON (ipro:InterProFamily) ASSERT ipro.family_id  IS UNIQUE;
    CREATE CONSTRAINT IF NOT EXISTS ON (go:GOClass) ASSERT go.class_id IS UNIQUE;
    CREATE CONSTRAINT IF NOT EXISTS ON (q:RibosomeStructure) Assert q.rcsb_id IS UNIQUE;
    CREATE CONSTRAINT IF NOT EXISTS ON (pf:PFAMFamily) assert pf.family_id  is unique;
    CREATE CONSTRAINT IF NOT EXISTS ON (lig:Ligand) assert lig.chemicalId is unique;
    CREATE CONSTRAINT IF NOT EXISTS ON (nc:NomenclatureClass) assert nc.class_id is unique;
    CREATE CONSTRAINT IF NOT EXISTS ON (nc:NomenclatureClass) assert nc.class_id is unique;
    """)

def get_all_constraints(tx: Transaction | ManagedTransaction):
    return tx.run("""//
    CALL db.constraints;
    """)