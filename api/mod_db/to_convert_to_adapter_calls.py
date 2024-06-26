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

    def get_rna_class(self, class_id: PolynucleotideClass) -> list[NomenclatureClassMember]:
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
                        match (n:ProteinClass)-[]-(rp:Protein)-[]-(s:RibosomeStructure) where toLower(n.class_id) contains "{FAMILY}" and {SUBUNIT} 
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

    def gmo_nom_class(self, class_id: CytosolicProteinClass) -> list[NomenclatureClassMember]:

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