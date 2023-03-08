from typing import Callable
from neo4j import GraphDatabase, Driver, ManagedTransaction, Record, Result, Transaction
from neo4j.graph import Node, Relationship
from neo4j import ManagedTransaction, Transaction
from api.ribctl.db.inits.driver import Neo4jDB
from api.ribctl.lib.types.types_ribosome import Ligand, Protein
from api.schema.v0 import LigandInstance, NeoStruct
from ribctl.lib.types.types_polymer import list_LSU_Proteins, list_SSU_Proteins, list_RNAClass


class QueryOps(Neo4jDB):
    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)

    def get_all_ligands(self)->list[LigandInstance]:
        """this is used by "request all ligands" in the frontend, binding sites action types."""

        with self.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                return tx.run("""//
                match (l:Ligand)-[]-(r:RibosomeStructure)  where 
                not l.chemicalName  contains "ION" 
                and not l.chemicalName contains "CLUSTER"
                and not l.chemicalName contains "["
                and r.expMethod <> "X-RAY DIFFRACTION"
                with collect({
                    polymer:false,
                    description:l.chemicalName,
                    chemicalId:l.chemicalId,
                    presentIn:{
                        src_organism_ids: r.src_organism_ids,
                        description:l.chemicalName,
                        citation_title:r.citation_title,
                        expMethod:r.expMethod,
                        rcsb_id:r.rcsb_id,
                        resolution:r.resolution}
                        }) as lig_instance
                return lig_instance
                """).value()[0]
            return session.read_transaction(_)

    def get_individual_ligand(self, chemId: str) -> Ligand:
        with self.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                return tx.run("""match (ligand:Ligand{chemicalId: $CHEM_ID}) return ligand""", {"CHEM_ID": chemId}).data()[0]['ligand']
            return session.read_transaction(_)

    def get_struct(self, rcsb_id: str) -> NeoStruct:

        with self.driver.session() as session:

            def _(tx: Transaction | ManagedTransaction):
                return tx.run("""//
                    match (struct:RibosomeStructure{rcsb_id:$RCSB_ID})
                    optional match (rrna:RNA)-[]-(struct)
                    with struct, collect(rrna) as rnas
                    optional match (rp:Protein)-[]-(struct)
                    with struct, rnas,  collect(rp) as rps
                    optional match (l:Ligand)-[]-(struct)
                    with struct, rnas, rps, collect(l.chemicalId) as ligands
                    return struct, ligands,rnas, rps
                                        """, {"RCSB_ID": rcsb_id.upper()}).data()[0]

            return session.read_transaction(_)

    def get_all_structures(self) -> list[NeoStruct]:

        with self.driver.session() as session:

            def _(tx: Transaction | ManagedTransaction):
                return tx.run("""//
                        match (ribs:RibosomeStructure) 
                                unwind ribs as struct

                                optional match (l:Ligand)-[]-(struct)
                                with collect(l.chemicalId) as ligands, struct

                                optional match (rps:Protein)-[]-(struct)
                                with ligands, struct, collect({
                                    auth_asym_id                   : rps.auth_asym_id,
                                    nomenclature                   : rps.nomenclature,
                                    entity_poly_seq_one_letter_code: rps.entity_poly_seq_one_letter_code
                                    }) as rps

                                optional match (struct_rnas:RNA)-[]-(struct)
                                with ligands, struct, rps, collect({
                                    auth_asym_id                   : struct_rnas.auth_asym_id,
                                    nomenclature                   : struct_rnas.nomenclature,
                                    entity_poly_seq_one_letter_code: struct_rnas.entity_poly_seq_one_letter_code
                                    }) as rnas
                                with ligands, rps, rnas, keys(struct) as keys, struct 
                                return struct, ligands,rps,rnas limit 3
                                        """).data()

            return session.read_transaction(_)
