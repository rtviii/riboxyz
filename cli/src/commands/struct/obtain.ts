import { Flags } from "@oclif/core";
import { exec } from "node:child_process";
import { BaseCommand } from "../..";
import { queryStructDb, StructureFolder } from "../../structure_processing/structure";
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
        let x = new StructureFolder(rcsb_id, true)
            await x.initialize_assets(flags.repair)
            await x.initialize_ligands(flags.repair, x.structure as RibosomeStructure )

        if (flags.commit) {
          this.log(`Commiting ${rcsb_id} to the database`)
          const commit_script = process.env["COMMIT_STRUCTURE_SH"]
          let   current_db    = process.env["NEO4J_CURRENTDB"]
          let   uri           = process.env["NEO4J_URI"]
          let   invocation    = `${commit_script} -s ${rcsb_id} -d ${current_db} -a "${uri}"`
          let   cp            = exec(invocation)
          if (cp.stderr !== null) {
            cp.stderr.on("data",
              (data) => {
                console.log(data)
              }
            )
          }
          cp.stdout?.on("data", (data) => { console.log(data) })
        }

    }
}

