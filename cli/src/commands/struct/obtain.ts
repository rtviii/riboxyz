import { Flags } from "@oclif/core";
import { exec } from "node:child_process";
import { BaseCommand } from "../..";
import { commit_struct_to_Db, queryStructDb, StructureFolder } from "../../structure_processing/structure";
import { RibosomeStructure } from "../../structure_processing/structure_types";

export default class Obtain extends BaseCommand {
    static description = 'Query structure in the database'
    static flags = {
        repair: Flags.boolean({ char: 'R' }),
        commit: Flags.boolean({ char: 'C' }),
    }

    static args = [
        { name: 'rcsb_id', required: true }
    ]

    public async run(): Promise<void> {
        const { args, flags } = await this.parse(Obtain)
        const rcsb_id = String(args.rcsb_id).toUpperCase();
        this.log(`Obtaining structure ${rcsb_id}`)
        let x = new StructureFolder(rcsb_id)

        this.log("----------11-------------------------------")
        await x.initialize_assets(flags.repair)
        this.log("---------12------------------------")
        await x.initialize_ligands(flags.repair, x.structure as RibosomeStructure)
        if (flags.commit) {
            commit_struct_to_Db(rcsb_id)
        }

    }
}

