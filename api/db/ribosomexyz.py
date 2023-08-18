from concurrent.futures import ALL_COMPLETED, Future, ThreadPoolExecutor, wait
from api.logs.loggers import get_updates_logger
import typing
from neo4j.exceptions import AuthError
from pyparsing import Any
from neo4j import Driver, GraphDatabase
from ribctl.etl.etl_pipeline import current_rcsb_structs
from rbxz_bend.settings import NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER
from api.db.inits.proteins import add_protein, node__protein_class
from api.db.inits.rna import add_rna, node__rna_class
from api.db.inits.structure import add_ligand, node__structure
from ribctl.lib.types.types_ribosome import RibosomeStructure
from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.types.types_poly_nonpoly_ligand import list_LSUProteinClass, list_SSUProteinClass, list_RNAClass
from neo4j import GraphDatabase, Driver, ManagedTransaction, Transaction
from ribctl.lib.types.types_ribosome import  NonpolymericLigand,  ProteinClass, RibosomeStructure
from schema.data_requests import LigandsByStruct
from schema.v0 import ExogenousRNAByStruct,BanClassMetadata, LigandInstance, NeoStruct, NomenclatureClass, NomenclatureClassMember
from ribctl.lib.types.types_poly_nonpoly_ligand import RNAClass, list_LSUProteinClass, list_SSUProteinClass, list_RNAClass


# ※ ----------------[ 0.Database  inits: constraints & nomenclature classes]
NODE_CONSTRAINTS = [
    """CREATE CONSTRAINT IF NOT EXISTS ON (ipro:InterProFamily) ASSERT ipro.family_id  IS UNIQUE;""",
    """CREATE CONSTRAINT IF NOT EXISTS ON (go:GOClass) ASSERT go.class_id IS UNIQUE;""",
    """CREATE CONSTRAINT IF NOT EXISTS ON (ribosome:RibosomeStructure) Assert q.rcsb_id IS UNIQUE;""",
    # """CREATE CONSTRAINT IF NOT EXISTS ON (lig:Ligand) assert lig.chemicalId is unique;""",
    """CREATE CONSTRAINT IF NOT EXISTS ON (nc:NomenclatureClass) assert nc.class_id is unique;""",
]

# !>>>>>>> Enterprise edition needed
# The NODE KEY constraint allows you to specify multiple properties that must be unique across all nodes with a specific label. Here's an example of how to create a NODE KEY constraint on the "name" and "email" properties of nodes with the "Person" label:
# CREATE CONSTRAINT ON (p:Person) ASSERT (p.name, p.email) IS NODE KEY;
# !<<<<<<<

# ※ ----------------[ 1.Structure Nodes]
# ※ ----------------[ 2.RNA Nodes]
# ※ ----------------[ 3.Protein Nodes]
# ※ ----------------[ 4.Ligand Nodes]
# ※ ----------------[ 5.Ingress]

# If you are connecting via a shell or programmatically via a driver,
# just issue a `ALTER CURRENT USER SET PASSWORD FROM 'current password' TO 'new password'` statement against
# the system database in the current session, and then restart your driver with the new password configured.

class ribosomexyzDB():
    driver: Driver
    uri      : str
    password : str
    user     : str
    databases: list[str]


    def change_default_pass(self):
        with GraphDatabase.driver(self.uri, auth=("neo4j", "neo4j")).session(database='system') as s:
            s.run("""ALTER CURRENT USER SET PASSWORD FROM "neo4j" TO "ribosomexyz";""")
        print("[INIT]:Changed the default Neo4j password. Initializing new database instance")
   
    def initialize_new_instance(self):

        self.__init_constraints()
        print("[INIT]: Initialized constraints.")
   
        self.__init_protein_classes()
        print("[INIT]: Initialized protein classes.")

        self.__init_rna_classes()
        print("[INIT]: Initialized rna classes.")

        print("[INIT]: Done with constraints and classes. >>>Sync with RCSB PDB is manual<<<")

    def __init__(self, uri: str, user: str, password: str) -> None:
        self.uri      = uri
        self.user     = user
        self.password = password
        try:
            self.driver = GraphDatabase.driver(uri, auth=(user, password))
        except AuthError as ae:
            self.change_default_pass()
            self.initialize_new_instance()

    def see_current_auth(self):
        print("see_current_auth")
        print(f"NEO4J_VAR: {NEO4J_URI} {NEO4J_USER} {NEO4J_PASSWORD}")
        with self.driver.session(database='system') as s:
            # {
            #   "user": "neo4j",
            #   "roles": null,
            #   "passwordChangeRequired": false,
            #   "suspended": null,
            #   "home": null
            # }
            r = s.run("""show users""")
            users_array = r.data()
            return users_array

    def show_dbs(self):
        with self.driver.session(database='system') as s:
            r = s.run("""show databases""")
            return r.data()

    def write(self, cypher:str):
        with self.driver.session() as s:
            r = s.run(cypher)
            return r.data()

    def see_constraints(self) -> list[dict[str, Any]]:
        with self.driver.session() as s:
            r = s.run("""//
            CALL db.constraints;
                  """)
            return r.data()


    #※----------------------------------------------------------------------------------------


    def get_all_structs(self):
        with self.driver.session() as s:
            struct_ids = []
            [struct_ids.extend(struct) for struct in s.execute_read(lambda tx: tx.run("""//
            match (n:RibosomeStructure) return n.rcsb_id;
            """).values('n.rcsb_id'))]
            return struct_ids

    def get_any(self) -> list[dict[str, Any]]:
        with self.driver.session() as s:

            return s.execute_read(lambda tx: tx.run("""//
            match (n) return n limit 10;
            """).data())

    # def init_db(self,driver):
    #     self.__init_constraints(di)
    #     self.__init_protein_classes()
    #     self.__init_rna_classes()

        # TODO : Ligand classes

    def add_structure(self, struct_assets: RibosomeAssets):

        R = RibosomeStructure.parse_obj(struct_assets.profile())

        with self.driver.session() as s:
            struct_node_result = s.execute_write(node__structure(R))

        for protein in R.proteins:
            add_protein(self.driver, protein)

        if R.rnas is not None:
            for rna in R.rnas:
                add_rna(self.driver, rna)

        if R.ligands is not None:
            for ligand in R.ligands:
                add_ligand(self.driver, ligand, R.rcsb_id)

       
        return struct_node_result.data()

    def get_all_ligands(self) -> list[LigandInstance]:
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
            return session.execute_read(_)

    def get_individual_ligand(self, chemId: str) -> NonpolymericLigand:
        with self.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                return tx.run("""match (ligand:Ligand{chemicalId: $CHEM_ID}) return ligand""", {"CHEM_ID": chemId}).data()[0]['ligand']
            return session.execute_read(_)

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

            return session.execute_read(_)

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
                                return struct, ligands,rps,rnas
                                        """).data()

            return session.execute_read(_)

    def get_all_ligandlike(self) -> list[LigandInstance]:
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
            return session.execute_read(_)

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

    def get_ligands_by_struct(self) -> list[LigandsByStruct]:
        with self.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                return tx.run("""//
        match (n:RibosomeStructure)-[]-(l:Ligand)
           with { title: n.citation_title, struct:n.rcsb_id, organism:n.src_organism_names, taxid:n.src_organism_ids, 
           ligands: collect({ chemid: l.chemicalId, name:l.chemicalName, number:l.number_of_instances })} as struct_ligs
           return struct_ligs
                        """).data()[0]['struct_ligs']
            return session.execute_read(_)

    def match_structs_w_proteins(self, targets: list[ProteinClass]) -> list[str]:
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

    def get_full_structure(self, rcsb_id: str) -> NeoStruct:
        # TODO: This method is identical to "get_struct" and was merely queried in a different way. DEPRECATE ONE (AS WELL AS THE API ENDPOINT)
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
                        """, {"RCSB_ID": rcsb_id.upper()}).data()[0]
            return session.execute_read(_)

    def get_banclass_for_chain(self, rcsb_id: str, auth_asym_id) -> list[ProteinClass]:
        # TODO: This method should handle both protein and RNA chains(atm only proteins)
        with self.driver.session() as session:

            def _(tx: Transaction | ManagedTransaction):
                return [pc[0] for pc in tx.run("""//
                match (n:RibosomeStructure {rcsb_id:$RCSB_ID})-[]-(c:Protein{auth_asym_id:$AUTH_ASYM_ID})-[]-(pc:ProteinClass) return pc.class_id
                    """, {"RCSB_ID": rcsb_id, "AUTH_ASYM_ID": auth_asym_id}).values('pc.class_id')]
            return session.execute_read(_)

    def get_banclasses_metadata(self, family: typing.Literal['b', 'e', 'u'], subunit: typing.Literal['SSU', 'LSU']) -> list[BanClassMetadata]:
        # TODO: This method should handle both protein and RNA chains(atm only proteins)
        with self.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):

                def flag_into_filter(_subunit: str):
                    if _subunit == "SSU":
                        return 'toLower(n.class_id) contains "s" or toLower(n.class_id) contains "bthx" or toLower(n.class_id) contains "rack"'
                    elif _subunit == "LSU":
                        return 'toLower(n.class_id) contains "l"'
                    else:
                        raise ValueError(
                            "Subunit must be either 'ssu' or 'lsu'")

                fstring = flag_into_filter(subunit)
                return tx.run("""//
                        match (n:ProteinClass)-[]-(rp:Protein)-[]-(s:RibosomeStructure) where  toLower(n.class_id) contains "{FAMILY}" and {SUBUNIT} 
                        unwind s.`src_organism_ids` as orgid
                        with collect(distinct orgid) as organisms, n.class_id as banClass, collect(s.rcsb_id) as structs, collect(distinct rp.pfam_comments) as comments
                        return  banClass, organisms, comments, structs""".format_map({"FAMILY": family, "SUBUNIT": fstring})).data()  # type: ignore
            return session.execute_read(_)

    def list_nom_classes(self) -> list[NomenclatureClass]:
        # TODO: Merge into get_banclasses_metadata
        with self.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                return tx.run("""//
                    match (b:ProteinClass)-[]-(rp)-[]-(str:RibosomeStructure)
                    with collect(str.rcsb_id) as structs, b.class_id as banClass, collect({
                        organism_desc: rp.src_organism_names,
                        organism_id  : rp.src_organism_ids,
                        uniprot      : rp.uniprot_accession,
                        parent       : str.rcsb_id,
                        parent_reso  : str.resolution,
                        strand_id    : rp.entity_poly_strand_id
                        }) as rps
                    return structs, banClass, rps
                 """).data()
            return session.execute_read(_)

    def gmo_nom_class(self, class_id: ProteinClass) -> list[NomenclatureClassMember]:

        with self.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                return tx.run("""//
                    match (rib:RibosomeStructure)-[]-(n:Protein)-[]-(nc:ProteinClass{class_id:$BANCLASS})
                    with collect({  
                    parent_rcsb_id                     : rib.parent_rcsb_id,
                    parent_resolution                  : rib.resolution,
                    parent_citation                    : rib.citation_title,
                    parent_year                        : rib.citation_year,
                    parent_method                      : rib.expMethod,

                    pfam_accessions                    : n.pfam_accessions,
                    pfam_comments                      : n.pfam_comments,
                    pfam_descriptions                  : n.pfam_descriptions,

                    asym_ids                            : n.asym_ids                           ,
                    auth_asym_id                        : n.auth_asym_id                       ,

                    src_organism_ids                   : n.src_organism_ids,
                    src_organism_names                 : n.src_organism_names,

                    host_organism_names                 : n.host_organism_names                ,
                    host_organism_ids                   : n.host_organism_ids                  ,

                    uniprot_accession                  : n.uniprot_accession,
                    rcsb_pdbx_description              : n.rcsb_pdbx_description,

                    entity_poly_strand_id              : n.entity_poly_strand_id,
                    entity_poly_seq_one_letter_code    : n.entity_poly_seq_one_letter_code,
                    entity_poly_seq_one_letter_code_can: n.entity_poly_seq_one_letter_code_can,
                    entity_poly_seq_length             : n.entity_poly_seq_length,
                    entity_poly_polymer_type           : n.entity_poly_polymer_type,
                    entity_poly_entity_type            : n.entity_poly_entity_type,

                    nomenclature                       : n.nomenclature,
                    ligand_like                        : n.ligand_like
                        }) as member
                        return member
                 """, {"BANCLASS": class_id}).value('member')[0]
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

    def get_rnas_by_struct(self) -> list[ExogenousRNAByStruct]:
        with self.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                return tx.run("""//
match (n:RibosomeStructure)-[]-(r:RNA) 
where toLower(r.rcsb_pdbx_description) contains "mrna" or 
toLower(r.rcsb_pdbx_description) contains "trna" or
toLower(r.rcsb_pdbx_description)  contains "m-rna" or
toLower(r.rcsb_pdbx_description)  contains "t-rna" or
toLower(r.rcsb_pdbx_description)  contains "messenger"
or toLower(r.rcsb_pdbx_description) contains "transfer"
with n.rcsb_id as struct, collect(r.rcsb_pdbx_description) as rnas
        return struct,rnas
                    """).data()
            return session.execute_read(_)

    def get_rna_class(self, class_id: RNAClass) -> list[NomenclatureClassMember]:
        with self.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                return [ rna[0] for rna in tx.run("""//
                    match (c:RNAClass { class_id:$RNA_CLASS })-[]-(n)-[]-(rib:RibosomeStructure)
            with {
            parent_year                         : rib.citation_year                    ,
            parent_resolution                   : rib.resolution                       ,
            parent_citation                     : rib.citation_title                   ,
            parent_rcsb_id                      : rib.parent_rcsb_id                   ,
            parent_method                       : rib.expMethod                        ,

            pfam_accessions                    : n.pfam_accessions                     ,
            pfam_comments                      : n.pfam_comments                       ,
            pfam_descriptions                  : n.pfam_descriptions                   ,

            asym_ids                            : n.asym_ids                           ,
            auth_asym_id                        : n.auth_asym_id                       ,

            src_organism_names                  : n.src_organism_names                 ,
            src_organism_ids                    : n.src_organism_ids                   ,

            host_organism_names                 : n.host_organism_names                ,
            host_organism_ids                   : n.host_organism_ids                  ,
            uniprot_accesion                    : n.uniprot_accesion                   ,
            rcsb_pdbx_description               : n.rcsb_pdbx_description              ,

            entity_poly_strand_id               : n.entity_poly_strand_id              ,
            entity_poly_seq_one_letter_code     : n.entity_poly_seq_one_letter_code    ,
            entity_poly_seq_one_letter_code_can : n.entity_poly_seq_one_letter_code_can,
            entity_poly_seq_length              : n.entity_poly_seq_length             ,
            entity_poly_polymer_type            : n.entity_poly_polymer_type           ,
            entity_poly_entity_type             : n.entity_poly_entity_type            ,

            nomenclature                        : [c.class_id]                         ,
            ligand_like                         : n.ligand_like                        
        } as rna
        return rna""", {"RNA_CLASS": class_id}).values('rna')]

            return session.execute_read(_)

    def sync_with_rcsb(self, workers:int)->None:

        logger = get_updates_logger()
        synced   = self.get_all_structs()
        unsynced = sorted(current_rcsb_structs())
        futures:list[Future] =  []

        logger.info("Started syncing with RCSB") 

        def log_commit_result( rcsb_id:str):
            def _(f:Future):
                if not None == f.exception():
                    logger.error(rcsb_id + ":" +f.exception().__str__())
                else:
                    logger.debug(rcsb_id + ":" + f.result().__str__())
            return _


        with ThreadPoolExecutor(max_workers=workers) as executor:
            for rcsb_id in list(set(unsynced ) - set(synced)):

                assets = RibosomeAssets(rcsb_id)
                assets._verify_json_profile(True)
                fut = executor.submit(self.add_structure, assets)
                fut.add_done_callback(log_commit_result(rcsb_id))
                futures.append(fut)

        wait(futures, return_when=ALL_COMPLETED)
        logger.info("Finished syncing with RCSB")

    def __init_constraints(self) -> None:
        with self.driver.session() as session:
            for c in NODE_CONSTRAINTS:
                session.execute_write(lambda tx: tx.run(c))

    def __init_protein_classes(self):
        with self.driver.session() as session:
            for protein_class in [*list_LSUProteinClass, *  list_SSUProteinClass]:
                session.execute_write(node__protein_class(protein_class))

    def __init_rna_classes(self):
        with self.driver.session() as session:
            for rna_class in list_RNAClass:
                session.execute_write(node__rna_class(rna_class))