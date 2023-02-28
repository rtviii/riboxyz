import { Flags } from '@oclif/core'
import { exec} from 'child_process'
import { BaseCommand } from '../..'
import { existsSync, readFileSync } from 'fs'
import { RibosomeStructure } from '../../structure_processing/structure_types'
// https://search.rcsb.org/#introduction
import axios from "axios";
import {  ungzip } from 'node-gzip'
import { mkdirSync, writeFileSync } from 'fs'
import { processPDBRecord } from "../../structure_processing/structure_json_profile";
import path = require("path");
import { StructureFolder } from '../../structure_processing/structure'

export default class Show extends BaseCommand {

  static description = 'Query structure in the database'
  static flags = {
    files : Flags.boolean({}),
    db    : Flags.boolean({}),
    dryrun: Flags.boolean(),
    force : Flags.boolean({ char: 'f' }),
  }

  static args = [
    { name: 'rcsb_id', required: true }
  ]

  public async run(): Promise<void> {
    const { args, flags } = await this.parse(Show)
    const rcsb_id = String(args.rcsb_id).toUpperCase();
    if (flags.files) {
      let struct = new StructureFolder(rcsb_id);
      console.log("\t\tFiles at ", struct.assets.folder_path() )
      struct.assets.init_assets()
    } 
    if (flags.db) {
      console.log("\t\tDB at ", process.env["NEO4J_URI"])
      console.log(`The following structs have been found for ${rcsb_id}`, await queryStructDb(rcsb_id))
    }
  }
}







/**
 * Request and display the state of the given rcsb_id structure in the database instance.
 */
const queryStructDb = (rcsb_id: string) => {
  return new Promise<string[]>((resolve, reject) => {
    let y = exec(`echo \"match (struct:RibosomeStructure {rcsb_id:\\"${rcsb_id.toUpperCase()}\\"}) return struct.rcsb_id\" | cypher-shell -a \"${process.env["NEO4J_URI"]}\" --format plain -u ${process.env["NEO4J_USER"]} -p ${process.env["NEO4J_PASSWORD"]} --database ${process.env["NEO4J_CURRENTDB"]}`,
      { env: process.env },
      (err, stdout, stderr) => {
        if (err && err?.code != 0) {
          process.stdout.write("Got shell error " + stderr + stdout)
          console.log("Got Error code:", err?.code)
          reject(err)
        }

        const dbstructs = stdout != null ? (stdout as string).replace(/"/g, '').split("\n").filter(r => r.length === 4) : []
        console.log(dbstructs)
        resolve(dbstructs)
      })
  })

}


/**
 * Download a .cif model of the structure.
 * @param struct_id 
 */
const download_unpack_place = async (struct_id: string) => {
  const BASE_URL = "http://files.rcsb.org/download/"
  const FORMAT = ".cif.gz"

  const structid = struct_id.toUpperCase()
  let url = BASE_URL + structid + FORMAT
  let compressed: Buffer = await axios.get(url, { responseType: 'arraybuffer' }).then(r => { return r.data })
    .catch(e => { process.stdout.write(`Structure ${structid} failed: `, e); return []; })
  let decompressed = await ungzip(compressed);

  // let destination_chains = path.join(
  //   process.env["RIBETL_DATA"] as string,
  //   `${structid}`,
  //   `CHAINS`)

  // if (!existsSync(destination_chains)) {
  //   mkdirSync(destination_chains)
  //   process.stdout.write(`Created directory ${destination_chains}.`);
  // }
  let structfile = path.join(
    process.env["RIBETL_DATA"] as string,
    `${structid}`,
    `${structid}.cif`)
  writeFileSync(structfile, decompressed)
}

export const save_struct_profile = (r: RibosomeStructure): string => {
  var rcsb_id = r.rcsb_id;
  var target_filename = path.join(
    process.env["RIBETL_DATA"] as string,
    rcsb_id.toUpperCase(),
    rcsb_id.toUpperCase() + ".json"
  );

  if (!existsSync(path.dirname(target_filename))) {
    mkdirSync(path.dirname(target_filename));
  }
  writeFileSync(target_filename, JSON.stringify(r, null, 4));
  return target_filename
};
