import { exec } from "node:child_process";
import { BaseCommand } from "../..";
import { queryStructDb, StructureFolder } from "../../structure_processing/structure";

export default class Obtain extends BaseCommand {

    static description = 'Query structure in the database'
    static flags = {}

    static args = [
        { name: 'rcsb_id', required: true }
    ]

    public async run(): Promise<void> {
        const { args, flags } = await this.parse(Obtain)
        const rcsb_id = String(args.rcsb_id).toUpperCase();
        this.log(`Obtaining structure ${rcsb_id}`)
        let x = new StructureFolder(rcsb_id, true)
        await x.initialize_assets(true)
    }
}

