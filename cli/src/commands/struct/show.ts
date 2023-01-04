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

export default class Show extends BaseCommand {

  static description = 'Query structure in the database'
  static flags = {
    files : Flags.boolean({}),
    db    : Flags.boolean({}),
    repair: Flags.boolean({ char: 'R' }),
    commit: Flags.boolean({ char: 'C' }),
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
      let assets = new StructureFolder(rcsb_id);
      console.log("\t\tFiles at ", assets.folder_path )
      assets.assets_verify(flags.repair)
    } 
    if (flags.db) {
      console.log("\t\tDB at ", process.env["NEO4J_URI"])
      console.log(`The following structs have been found for ${rcsb_id}`, await queryStructDb(rcsb_id))
    }
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




export class StructureFolder {

  folder_path: string
  private cif_filepath: string
  private cif_modified_filepath: string
  private json_profile_filepath: string
  private chains_folder: string
  private png_thumbnail_filepath: string
  private rcsb_id: string
  __structure: RibosomeStructure
  private ligands: string[] = []
  private ligand_like_polymers: string[] = []

  constructor(rcsb_id: string) {
    this.rcsb_id = rcsb_id.toUpperCase()

    if (!process.env["RIBETL_DATA"]) {
      throw Error("RIBETL_DATA environment variable not set. Cannot access assets.")
    }
    this.folder_path = `${process.env["RIBETL_DATA"]}/${this.rcsb_id}`

    this.cif_filepath = `${this.folder_path}/${this.rcsb_id}.cif`
    this.cif_modified_filepath = `${this.folder_path}/${this.rcsb_id}_modified.cif`
    this.json_profile_filepath = `${this.folder_path}/${this.rcsb_id}.json`
    this.chains_folder = `${this.folder_path}/CHAINS`
    this.png_thumbnail_filepath = `${this.folder_path}/_ray_${this.rcsb_id}.png`
    this.__structure = JSON.parse(readFileSync(this.json_profile_filepath, 'utf-8'))

    for (var chain of [...this.__structure.proteins, ...(this.__structure.rnas || [])]) {
      if (chain.ligand_like) {
        this.ligand_like_polymers = [...this.ligand_like_polymers, chain.auth_asym_id]
      }
    }

    if (!this.__structure.ligands) {
      console.log("No ligands found");
    } else {
      for (var lig of this.__structure.ligands) {
        this.ligands = [...this.ligands, lig.chemicalId]
      }
    }
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
  private __verify_chain_files(): boolean {
    if (!this.__verify_chains_folder()) {
      return false
    }
    const chain_files_all = [
      ...this.__structure.proteins?.map((c) => { return c.auth_asym_id }),
      ...(this.__structure.rnas || []).map((c) => { return c.auth_asym_id })
    ].map(
      (chain_id) => {
        if (!existsSync(`${this.chains_folder}/${this.rcsb_id}_STRAND_${chain_id}.cif`)) {
          console.log(`[${this.rcsb_id}]: NOT FOUND ${this.chains_folder}/${this.rcsb_id}_STRAND_${chain_id}.cif`)
          return false
        } else return true
      }
    ).reduce((prev, cur) => { return prev && cur }, true)
    return chain_files_all

  }

  private __verify_ligands_and_polymers(try_fix: boolean = false) {
    let ligs = this.ligands && this.ligands.map((lig_chem_id) => {
      if (!existsSync(`${this.folder_path}/LIGAND_${lig_chem_id}.json`)) {
        console.log(`[${this.rcsb_id}]: NOT FOUND ${this.folder_path}/LIGAND_${lig_chem_id}.json (Is it an ION? Expected.)`)
        return false
      } else return true
    }).reduce((prev, cur) => { return prev && cur }, true)

    let polys = this.ligand_like_polymers && this.ligand_like_polymers.map((polymer_id) => {
      if (!existsSync(`${this.folder_path}/POLYMER_${polymer_id}.json`)) {
        console.log(`[${this.rcsb_id}]: NOT FOUND ${this.folder_path}/POLYMER_${polymer_id}.json`)
        return false
      } else return true
    })

    if ((!polys || !ligs) && try_fix) {
      console.log("Some ligands are missing. Calling script:", process.env["PYTHONBIN"], process.env["EXTRACT_BSITES_PY"])
      exec(`${process.env["PYTHONBIN"]} ${process.env["EXTRACT_BSITES_PY"]} -s ${this.rcsb_id} --save`, (err, stdout, stderr) => {
        console.log(err);
        console.log(stdout);
        console.log(stderr);
      })
    }
  }

  async assets_verify(try_fix: boolean = false) {
    if (!this.__verify_cif()) {
      console.log(`[${this.rcsb_id}]: NOT FOUND ${this.cif_filepath}`)
      if (try_fix) {
        download_unpack_place(this.rcsb_id)
      }
    }

    if (!this.__verify_cif_modified()) {
      console.log(`[${this.rcsb_id}]: NOT FOUND ${this.cif_modified_filepath}`)
      if (try_fix) {
        let y = exec(`${process.env["PYTHONBIN"]} ${process.env["SPLIT_RENAME_PY"]} -s ${this.rcsb_id}`,
          (err, stdout, stderr) => {
            console.log(err);
            console.log(stdout);
            console.log(stderr);
          })
        console.log(y.stdout)
      }
    }

    if (!this.__verify_json_profile()) {

      if (try_fix) {
        let ribosome = await processPDBRecord(this.rcsb_id)
        let filename = await save_struct_profile(ribosome)
        process.stdout.write(`Saved structure profile:\t${filename}`);
      }
      console.log(`[${this.rcsb_id}]: NOT FOUND ${this.json_profile_filepath}`)
    }
    if (!this.__verify_png_thumbnail()) {
      console.log(`[${this.rcsb_id}]: NOT FOUND ${this.png_thumbnail_filepath}`)
      if (try_fix) {
        let proc = exec(`${process.env["PYTHONBIN"]} ${process.env["RENDER_THUMBNAIL_PY"]} -s ${this.rcsb_id}`,
          (err, stdout, stderr) => {
            console.log(err);
            console.log(stdout);
            console.log(stderr);
          })
      }
    }

    if (!this.__verify_chain_files()) {
      if (try_fix) {
        exec(`${process.env["PYTHONBIN"]} ${process.env["SPLIT_RENAME_PY"]} -s ${this.rcsb_id}`, (err, stdout, stderr) => {
          console.log(err);
          console.log(stdout);
          console.log(stderr);
        })
      }
    }
    this.__verify_ligands_and_polymers(try_fix)
  }

  db_verify() {

    // TODO

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


// Options describing how to process a given structure
interface IngressOptions {
  rcsb_id?: string,
  acquirePDBRecord?: boolean;
  downloadCifFile?: boolean;
  splitRenameChains?: boolean;
  renderStructHero?: boolean;
  extractBindingSites?: boolean;
  commitToNeo4j?: boolean;
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

// const process_structure = async (opts: IngressOptions) => {
//   if (!opts.rcsb_id) throw new Error("No structure ID provided.")

//   let struct_id = opts.rcsb_id.toUpperCase()

//   if (opts.acquirePDBRecord) {
//     let ribosome = await processPDBRecord(struct_id)
//     let filename = await save_struct_profile(ribosome)
//     process.stdout.write(`Saved structure profile:\t${filename}`);
//   }

//   if (opts.downloadCifFile) {
//     await download_unpack_place(struct_id)
//     process.stdout.write(`Saved structure cif file ${struct_id}.cif`);
//   }
//   if (opts.renderStructHero) {
//     exec(`${process.env["PYTHONBIN"]} ${process.env["RENDER_THUMBNAIL_PY"]} -s ${struct_id}`)
//   }

//   if (opts.splitRenameChains) {
//     exec(`${process.env["PYTHONBIN"]} ${process.env["SPLIT_RENAME_PY"]} -s ${struct_id}`)
//   }

//   if (opts.extractBindingSites) {
//     exec(`${process.env["PYTHONBIN"]}   \
//          ${process.env["EXTRACT_BSITES_PY"]}       \
//          -s ${struct_id}                           \
//          --save`)
//   }

//   if (opts.commitToNeo4j) {
//     var shellstring = `${process.env["COMMIT_STRUCTURE_SH"]} -s ${opts.rcsb_id.toUpperCase()} \
//          -d ${process.env.NEO4J_CURRENTDB} \
//          -a "${process.env.NEO4J_URI}"`;
//     exec(shellstring)
//   }

// }