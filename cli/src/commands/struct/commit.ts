import { exec } from "node:child_process";
import { BaseCommand } from "../..";
import { queryStructDb, StructureFolder } from "../../structure_processing/structure";
import { RibosomeStructure } from "../../structure_processing/structure_types";

export default class Commit extends BaseCommand {

    static description = 'Query structure in the database'
    static flags = {}

    static args = [
        { name: 'rcsb_id', required: true }
    ]

    public async run(): Promise<void> {
        const { args, flags } = await this.parse(Commit)
        const rcsb_id = String(args.rcsb_id).toUpperCase();
        this.log(`Comitting structure ${rcsb_id}`)
        let x = new StructureFolder(rcsb_id, true)
        await x.initialize_assets(true)
        await x.initialize_ligands(true, x.structure as RibosomeStructure )

    }
}
