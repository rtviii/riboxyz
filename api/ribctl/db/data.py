from typing import Callable
from neo4j import GraphDatabase, Driver, ManagedTransaction, Record, Result, ResultSummary, Transaction
from neo4j.graph import Node, Relationship
from neo4j import ManagedTransaction, Transaction
from api.ribctl.db.inits.driver import Neo4jDB
from api.ribctl.lib.types.types_ribosome import Ligand, Protein, ProteinClass, RibosomeStructure
from api.schema.data_requests import LigandsByStruct
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

    def get_all_ligandlike(self)->list[LigandInstance]:
        with self.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                return tx.run("""//
                                match (l {ligand_like:true})-[]-(r:RibosomeStructure) 
                        where r.expMethod <> "X-RAY DIFFRACTION"
                        with collect ( {
                            polymer     : true,
                            description : l.rcsb_pdbx_description,
                            presentIn  : {
                                auth_asym_id    : l.auth_asym_id,
                                src_organism_ids: r.src_organism_ids,
                                description     : l.rcsb_pdbx_description,
                                citation_title  : r.citation_title,
                                expMethod       : r.expMethod,
                                rcsb_id         : r.rcsb_id,
                                resolution      : r.resolution
                            }
                        } ) as liglike
                        return liglike""").data()[0]['liglike']
            return session.read_transaction(_)
            

    def get_RibosomeStructure(self,rcsb_id:str)->RibosomeStructure:
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
                        """,{"RCSB_ID":rcsb_id.upper()}).data()[0]['structure']
            return session.read_transaction(_)

    def get_ligands_by_struct(self)->list[LigandsByStruct]:
        with self.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                return tx.run("""//
        match (n:RibosomeStructure)-[]-(l:Ligand)
           with { title: n.citation_title, struct:n.rcsb_id, organism:n.src_organism_names, taxid:n.src_organism_ids, 
           ligands: collect({ chemid: l.chemicalId, name:l.chemicalName, number:l.number_of_instances })} as struct_ligs
           return struct_ligs
                        """).data()[0]['struct_ligs']
            return session.read_transaction(_)

    def match_structs_w_proteins(self,targets:list[ProteinClass])->list[str]:
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
                        """,{"T":targets}).values('rcsb_id')]

            return session.read_transaction(_)


    def get_full_structure(self,rcsb_id:str)->NeoStruct:
        #TODO: This method is identical to "get_struct" and was merely queried in a different way. DEPRECATE ONE (AS WELL AS THE API ENDPOINT)
        with self.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                return tx.run("""//
    match (rib:RibosomeStructure {rcsb_id:$RCSB_ID}) 
        unwind rib as struct
        optional match (l:Ligand)-[]-(struct)
        with collect(l.chemicalId) as ligands, struct
        optional match (rps:Protein)-[]-(struct)
        with ligands, struct, collect({auth_asym_id:rps.auth_asym_id, nomenclature:rps.nomenclature, entity_poly_seq_one_letter_code: rps.entity_poly_seq_one_letter_code}) as rps
        optional match (rnas:RNA)-[]-(struct)
        with ligands, struct, rps, collect({auth_asym_id: rnas.auth_asym_id, nomenclature: rnas.nomenclature, entity_poly_seq_one_letter_code:rnas.entity_poly_seq_one_letter_code}) as rnas
        return struct, ligands,rps, rnas
                        """,{"RCSB_ID":rcsb_id.upper()}).data()[0]
            return session.read_transaction(_)