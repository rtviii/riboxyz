import { Command, Flags } from '@oclif/core'
import { ShellString } from 'shelljs'
import { exec } from 'child_process'
import { BaseCommand } from '../..'



export default class Show extends BaseCommand {


  static description = 'Query structure in the database'

  static flags = {
    assets: Flags.boolean({ char: 'a' }),
    db: Flags.boolean({ char: 'd' }),
  }

  static args = [{ name: 'rcsb_id' }]

  public async run(): Promise<void> {
    const { args, flags } = await this.parse(Show)
    const rcsb_id: string = args.rcsb_id;
    this.log(`got flags: `, flags)
    this.log(`got env: `, process.env["NEO4J_URI"])
    let dbquery = new Promise<string[]>((resolve, reject) => {
      let y = exec(`echo \"match (struct:RibosomeStructure {rcsb_id:\"${rcsb_id.toUpperCase()}\"}) return struct.rcsb_id\" | cypher-shell -a \"${process.env["NEO4J_URI"]}\" --format plain -u ${process.env["NEO4J_USER"]} -p ${process.env["NEO4J_PASSWORD"]} --database ${process.env["NEO4J_CURRENTDB"]}`,
        { env: process.env },
        (err, stdout, stderr) => {
          if (err?.code != 0) {
            process.stdout.write("Got shell error " + stderr + stdout)
            reject(err)
          }
          const dbstructs = (stdout as string).replace(/"/g, '').split("\n").filter(r => r.length === 4)
          this.log("got dbstructs: ", dbstructs)
          resolve(dbstructs)
        })

      // exec(`echo \"match (struct:RibosomeStructure {rcsb_id:${rcsb_id.toUpperCase()}}) return struct.rcsb_id\" | cypher-shell -a \"${process.env["NEO4J_URI"]}\" --format plain -u ${process.env["NEO4J_USER"]} -p ${process.env["NEO4J_PASSWORD"]} --database ${process.env[ "NEO4J_CURRENTDB" ]}`,
      //     function (err:number, stdout:string, stderr:any) {
      //         if (err != 0) {
      //             process.stdout.write("Got shell error " + stderr  + stdout)
      //             reject(err)
      //         }
      //         const dbstructs = (stdout as string).replace(/"/g, '').split("\n").filter(r => r.length === 4)
      //         resolve(dbstructs)
      //     })
    })

    console.log(await dbquery)

  }
}
