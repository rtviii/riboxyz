import { Command, Flags } from '@oclif/core'
import { ShellString } from 'shelljs'
import { exec } from 'child_process'
import { BaseCommand } from '../..'
import { flags, flagUsages } from '@oclif/core/lib/parser'
import { existsSync, readFileSync } from 'fs'
import { string } from 'yargs'
import { RibosomeStructure } from '../../RibosomeTypes'


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
  png_thumbnail_filepath  : string
  rcsb_id                 : string
  __structure             : RibosomeStructure
  constructor(rcsb_id: string) {
    this.rcsb_id = rcsb_id.toUpperCase()
    if (!process.env["RIBETL_DATA"]) {
      throw Error("RIBETL_DATA environment variable not set. Cannot access assets.")
    }
    this.cif_filepath           = `${process.env["RIBETL_DATA"]}/${this.rcsb_id}/${this.rcsb_id}.cif`
    this.cif_modified_filepath  = `${process.env["RIBETL_DATA"]}/${this.rcsb_id}/${this.rcsb_id}_modified.cif`
    this.json_profile_filepath  = `${process.env["RIBETL_DATA"]}/${this.rcsb_id}/${this.rcsb_id}.json`
    this.chains_folder          = `${process.env["RIBETL_DATA"]}/${this.rcsb_id}/CHAINS`
    this.png_thumbnail_filepath = `${process.env["RIBETL_DATA"]}/${this.rcsb_id}/_ray_${this.rcsb_id}.png`
    this.__structure            = JSON.parse(readFileSync(this.json_profile_filepath, 'utf-8'))
  }



  private __verify_cif() {
    if (existsSync(this.cif_filepath)) {
      return true
    } else return false
  }
  private __verify_cif_modified() {
    if (existsSync(this.cif_modified_filepath)) {
      return true
    } else return false
  }
  private __verify_json_profile() {
    if (existsSync(this.cif_modified_filepath)) {
      return true
    } else return false
  }
  private __verify_png_thumbnail() {
    if (existsSync(this.png_thumbnail_filepath)) {
      return true
    } else return false
  }
  private __verify_chains_folder() {
    if (existsSync(this.chains_folder)) {
      return true
    } else return false
  }
  //verify that each chain file exists
  private __verify_chain_files() {
    const chain_files = [ 
      ...this.__structure.proteins?.map((c)=>{return c.auth_asym_id}),  
     ...( this.__structure.rnas || [] ).map((c)=>{return c.auth_asym_id}) 
    ].forEach(
      (chain_id) => {
        if (!existsSync(`${this.chains_folder}/${this.rcsb_id}_STRAND_${chain_id}.cif`)){
          console.log(`[${this.rcsb_id}]: NOT FOUND ${this.chains_folder}/${this.rcsb_id}_STRAND_${chain_id}.cif`)
          return false
        }
      }
    )
    
  }

  assets_verify() {
    if (!this.__verify_cif          ()) { console.log(`[${this.rcsb_id}]: NOT FOUND ${this.cif_filepath          }`) }
    if (!this.__verify_cif_modified ()) { console.log(`[${this.rcsb_id}]: NOT FOUND ${this.cif_modified_filepath }`) }
    if (!this.__verify_json_profile ()) { console.log(`[${this.rcsb_id}]: NOT FOUND ${this.json_profile_filepath }`) }
    if (!this.__verify_png_thumbnail()) { console.log(`[${this.rcsb_id}]: NOT FOUND ${this.png_thumbnail_filepath}`) }
    if (!this.__verify_chains_folder()) { console.log(`[${this.rcsb_id}]: NOT FOUND ${this.chains_folder          }`) }
    this.__verify_chain_files()
  }

  db_verify(){
    
  }
}

interface StructureFolder {
}

/**
 * Request and display files associated with the given rcsb_id on disk
 */
const queryStructAssets = (rcsb_id: string) => {
  type StructureAsset = 'ChainsFolder' | 'CifModel' | 'CifModel_modified' | 'JsonProfile' | 'LigandProfile' | 'CifModelChain' | 'PngThumbnail'
  console.log(process.env["RIBETL_DATA"])
  let x = new StructureFolder(rcsb_id);
  x.assets_verify()
  // var struct:RibosomeStructure = JSON.parse(readFileSync(x.json_profile_filepath, 'utf-8'))

  // console.log(struct)


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