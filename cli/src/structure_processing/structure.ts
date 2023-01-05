import { existsSync, mkdirSync, readFileSync, writeFileSync } from "fs"
import { exec } from "shelljs"
import { processPDBRecord } from "./structure_json_profile"
import { RibosomeStructure } from "./structure_types"
import * as cp from 'child_process'
import { ungzip } from "node-gzip"
import axios from "axios"
import path = require("path")
import { Polymer_Entity } from "./rcsb_graphql_schema"



export class RibosomeAssets {
    rcsb_id: string

    constructor(
        rcsb_id: string
    ) {
        this.rcsb_id = rcsb_id
    }
    envcheck = () => {
        if (!process.env["RIBETL_DATA"]) {
            throw Error("RIBETL_DATA environment variable not set. Cannot access assets.")
        }
    }

    folder_path = () => {
        this.envcheck()
        return `${process.env["RIBETL_DATA"]}/${this.rcsb_id}`
    }

    cif_filepath = () => {
        this.envcheck()
        return `${this.folder_path()}/${this.rcsb_id}.cif`
    }

    cif_modified_filepath = () => {
        this.envcheck()
        return `${this.folder_path()}/${this.rcsb_id}_modified.cif`
    }

    json_profile_filepath = () => {
        this.envcheck()
        return `${this.folder_path()}/${this.rcsb_id}.json`
    }

    chains_folder = () => {
        this.envcheck()
        return `${this.folder_path()}/CHAINS`
    }

    png_thumbnail_filepath = () => {
        this.envcheck()
        return `${this.folder_path()}/_ray_${this.rcsb_id}.png`
    }



    async __verify_cif(obtain: boolean = false) {
        if (existsSync(this.cif_filepath())) { return true } else {
            if (obtain) {
                download_unpack_place(this.rcsb_id); return true
            } else return false
        }
    }
    async __verify_cif_modified(obtain: boolean = false) {
        if (existsSync(this.cif_modified_filepath())) { return true } else {
            if (obtain) {
                let y = exec(`${process.env["PYTHONBIN"]} ${process.env["SPLIT_RENAME_PY"]} -s ${this.rcsb_id}`,
                    (err, stdout, stderr) => {
                        console.log(err);
                        console.log(stdout);
                        console.log(stderr);
                    })
                console.log(y.stdout)
            } else return false
        }
    }
    async __verify_json_profile(obtain: boolean = false) {
        if (existsSync(this.cif_modified_filepath())) { return true } else {
            if (obtain) {
                console.log("Trying to call");
                
                let ribosome = await processPDBRecord(this.rcsb_id)
                console.log("Trying to call over");
                let filename = await save_struct_profile(ribosome)
                process.stdout.write(`Saved structure profile:\t${filename}`);
            } else return false
        }
    }
    async __verify_png_thumbnail(obtain: boolean = false) {
        if (existsSync(this.png_thumbnail_filepath())) { return true } else {
            if (obtain) {
                let proc = exec(`${process.env["PYTHONBIN"]} ${process.env["RENDER_THUMBNAIL_PY"]} -s ${this.rcsb_id}`,
                    (err, stdout, stderr) => {
                        console.log(err);
                        console.log(stdout);
                        console.log(stderr);
                    })
            } else return false
        }
    }
    async __verify_chains_folder(obtain: boolean = false) {
        if (existsSync(this.chains_folder())) { return true } else {
            if (obtain) {
                exec(`${process.env["PYTHONBIN"]} ${process.env["SPLIT_RENAME_PY"]} -s ${this.rcsb_id}`, (err, stdout, stderr) => {
                    console.log(err);
                    console.log(stdout);
                    console.log(stderr);
                })
            } else return false
        }
    }
    //verify that each chain file exists
    async __verify_chain_files(parent_structure: RibosomeStructure): Promise<boolean> {
        if (!this.__verify_chains_folder()) {
            return false
        }
        const chain_files_all = [
            parent_structure.proteins?.map((c) => { return c.auth_asym_id }),
            (parent_structure.rnas || []).map((c) => { return c.auth_asym_id })
        ].map(
            (chain_id) => {
                if (!existsSync(`${this.chains_folder()}/${this.rcsb_id}_STRAND_${chain_id}.cif`)) {
                    console.log(`[${this.rcsb_id}]: NOT FOUND ${this.chains_folder()}/${this.rcsb_id}_STRAND_${chain_id}.cif`)
                    return false
                } else return true
            }
        ).reduce((prev, cur) => { return prev && cur }, true)
        return chain_files_all

    }

    async __verify_ligands_and_polymers(struct: RibosomeStructure) {
        let ligs = struct.ligands && struct.ligands.map((lig_chem_id) => {
            if (!existsSync(`${this.folder_path()}/LIGAND_${lig_chem_id}.json`)) {
                console.log(`[${this.rcsb_id}]: NOT FOUND ${this.folder_path()}/LIGAND_${lig_chem_id}.json (Is it an ION? Expected.)`)
                return false
            } else return true
        }).reduce((prev, cur) => { return prev && cur }, true)


        let ligandlike: string[] = []
        let ligands: string[] = []

        for (var chain of [...struct.proteins, ...(struct.rnas || [])]) {
            if (chain.ligand_like) {
                ligandlike = [...ligandlike, chain.auth_asym_id]
            }
        }
        if (!struct.ligands) {
        } else {
            for (var lig of struct.ligands) {
                ligands = [...ligands, lig.chemicalId]
            }
        }


        let polys = ligandlike.map((polymer_id) => {
            if (!existsSync(`${this.folder_path()}/POLYMER_${polymer_id}.json`)) {
                console.log(`[${this.rcsb_id}]: NOT FOUND ${this.folder_path()}/POLYMER_${polymer_id}.json`)
                return false
            } else return true
        })

        if ((!polys || !ligs)) {
            console.log("Some ligands are missing. Calling script:", process.env["PYTHONBIN"], process.env["EXTRACT_BSITES_PY"])
            exec(`${process.env["PYTHONBIN"]} ${process.env["EXTRACT_BSITES_PY"]} -s ${this.rcsb_id} --save`, (err, stdout, stderr) => {
                console.log(err);
                console.log(stdout);
                console.log(stderr);
            })
        }
    }

    async init_assets(obtain: boolean) {
        await this.__verify_cif(obtain)
        await this.__verify_json_profile(obtain)
        await this.__verify_png_thumbnail(obtain)
        await this.__verify_chains_folder(obtain)
    }
}


export class StructureFolder {

    rcsb_id: string;
    assets: RibosomeAssets;
    structure?: RibosomeStructure;


    constructor(rcsb_id: string, obtain: boolean = false, ribosome_structure?: RibosomeStructure) {
        this.rcsb_id = rcsb_id.toUpperCase()
        console.log("before assets");
        this.assets = new RibosomeAssets(this.rcsb_id)
        console.log("after");


    }

    async initialize_assets(obtain: boolean) {
        await this.assets.init_assets(obtain)
        if (!this.assets.__verify_json_profile(obtain)) {
            throw Error(`Structure ${this.rcsb_id} assets not found. Cannot initiate resource.`)
        }
        this.structure = JSON.parse(readFileSync(this.assets.json_profile_filepath(), 'utf-8'))
    }

}


/**
 * Request and display the state of the given rcsb_id structure in the database instance.
 */
export const queryStructDb = (rcsb_id: string) => {
    return new Promise<string[]>((resolve, reject) => {
        let y = cp.exec(`echo \"match (struct:RibosomeStructure {rcsb_id:\\"${rcsb_id.toUpperCase()}\\"}) return struct.rcsb_id\" | cypher-shell -a \"${process.env["NEO4J_URI"]}\" --format plain -u ${process.env["NEO4J_USER"]} -p ${process.env["NEO4J_PASSWORD"]} --database ${process.env["NEO4J_CURRENTDB"]}`,
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
export const download_unpack_place = async (struct_id: string) => {
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
}