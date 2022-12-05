import { Command, Flags } from '@oclif/core'
import { ShellString } from 'shelljs'
import { exec } from 'child_process'
import { BaseCommand } from '../..'
import { flags } from '@oclif/core/lib/parser'



export default class Show extends BaseCommand {


  static description = 'Query structure in the database'
  static flags = {
    files: Flags.boolean({}),
    db: Flags.boolean({}),
  }

  static args = [
    { name: 'rcsb_id', required: true }
  ]

  public async run(): Promise<void> {
    const { args, flags } = await this.parse(Show)
    const rcsb_id = args.rcsb_id;

    if (flags.files) {
      queryStructAssets(rcsb_id)
    } else if (flags.db) {
      // console.log(await dbquery)
    }

  }
}


class StructureFolder {

  cif_filepath            : string
  cif_modified_filepath   : string
  json_profile_filepath   : string
  chains_folder           : string
  ligand_profile_filepath?: string[]

  constructor(rcsb_id: string) {
    if (!process.env["RIBETL_DATA"]) {
      throw Error("RIBETL_DATA environment variable not set. Cannot access assets.")
    }
    rcsb_id                    = rcsb_id.toUpperCase()
    this.cif_filepath          = `${process.env["RIBETL_DATA"]}/${rcsb_id}/${rcsb_id}.cif`
    this.cif_modified_filepath = `${process.env["RIBETL_DATA"]}/${rcsb_id}/${rcsb_id}_modified.cif`
    this.json_profile_filepath = `${process.env["RIBETL_DATA"]}/${rcsb_id}/${rcsb_id}.json`
    this.chains_folder         = `${process.env["RIBETL_DATA"]}/${rcsb_id}/CHAINS`
  }


  verify_assets() {

  }


}


interface StructureFolder {
}


/**
 * Request and display files associated with the given rcsb_id on disk
 */
const queryStructAssets = (rcsb_id: string) => {
  type StructureAsset = 'ChainsFolder' | 'CifModel' | 'CifModel_modified' | 'JsonProfile' | 'LigandProfile' | 'CifModelChain'
  console.log("hihihc assets call")
  console.log(process.env["RIBETL_DATA"])
  let x = new StructureFolder(rcsb_id);



}



/**
 * Request and display the state of the given rcsb_id structure in the database instance.
 */
const queryStructDb = (rcsb_id: string) => {
  let dbquery = new Promise<string[]>((resolve, reject) => {
    let y = exec(`echo \"match (struct:RibosomeStructure {rcsb_id:\"${rcsb_id.toUpperCase()}\"}) return struct.rcsb_id\" | cypher-shell -a \"${process.env["NEO4J_URI"]}\" --format plain -u ${process.env["NEO4J_USER"]} -p ${process.env["NEO4J_PASSWORD"]} --database ${process.env["NEO4J_CURRENTDB"]}`,
      { env: process.env },
      (err, stdout, stderr) => {
        if (err?.code != 0) {
          process.stdout.write("Got shell error " + stderr + stdout)
          reject(err)
        }
        const dbstructs = (stdout as string).replace(/"/g, '').split("\n").filter(r => r.length === 4)
        resolve(dbstructs)
      })

  })

}