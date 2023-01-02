"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.save_struct_profile = exports.StructureFolder = void 0;
const tslib_1 = require("tslib");
const core_1 = require("@oclif/core");
const child_process_1 = require("child_process");
const __1 = require("../..");
const fs_1 = require("fs");
// https://search.rcsb.org/#introduction
const axios_1 = tslib_1.__importDefault(require("axios"));
const node_gzip_1 = require("node-gzip");
const fs_2 = require("fs");
const structure_json_profile_1 = require("../../structure_processing/structure_json_profile");
const path = require("path");
class Show extends __1.BaseCommand {
    async run() {
        var _a;
        const { args, flags } = await this.parse(Show);
        const rcsb_id = String(args.rcsb_id).toUpperCase();
        if (flags.files) {
            let assets = new StructureFolder(rcsb_id);
            console.log("\t\tFiles at ", assets.folder_path);
            assets.assets_verify(flags.repair);
        }
        if (flags.db) {
            console.log("\t\tDB at ", process.env["NEO4J_URI"]);
            console.log(`The following structs have been found for ${rcsb_id}`, await queryStructDb(rcsb_id));
        }
        if (flags.commit) {
            this.log(`Commiting ${rcsb_id} to the database`);
            const commit_script = process.env["COMMIT_STRUCTURE_SH"];
            let current_db = process.env["NEO4J_CURRENTDB"];
            let uri = process.env["NEO4J_URI"];
            let invocation = `${commit_script} -s ${rcsb_id} -d ${current_db} -a "${uri}"`;
            let cp = (0, child_process_1.exec)(invocation);
            if (cp.stderr !== null) {
                cp.stderr.on("data", (data) => {
                    console.log(data);
                });
            }
            (_a = cp.stdout) === null || _a === void 0 ? void 0 : _a.on("data", (data) => { console.log(data); });
        }
    }
}
exports.default = Show;
Show.description = 'Query structure in the database';
Show.flags = {
    files: core_1.Flags.boolean({}),
    db: core_1.Flags.boolean({}),
    repair: core_1.Flags.boolean({ char: 'R' }),
    commit: core_1.Flags.boolean({ char: 'C' }),
    dryrun: core_1.Flags.boolean(),
    force: core_1.Flags.boolean({ char: 'f' }),
};
Show.args = [
    { name: 'rcsb_id', required: true }
];
class StructureFolder {
    constructor(rcsb_id) {
        this.ligands = [];
        this.ligand_like_polymers = [];
        this.rcsb_id = rcsb_id.toUpperCase();
        if (!process.env["RIBETL_DATA"]) {
            throw Error("RIBETL_DATA environment variable not set. Cannot access assets.");
        }
        this.folder_path = `${process.env["RIBETL_DATA"]}/${this.rcsb_id}`;
        this.cif_filepath = `${this.folder_path}/${this.rcsb_id}.cif`;
        this.cif_modified_filepath = `${this.folder_path}/${this.rcsb_id}_modified.cif`;
        this.json_profile_filepath = `${this.folder_path}/${this.rcsb_id}.json`;
        this.chains_folder = `${this.folder_path}/CHAINS`;
        this.png_thumbnail_filepath = `${this.folder_path}/_ray_${this.rcsb_id}.png`;
        this.__structure = JSON.parse((0, fs_1.readFileSync)(this.json_profile_filepath, 'utf-8'));
        for (var chain of [...this.__structure.proteins, ...(this.__structure.rnas || [])]) {
            if (chain.ligand_like) {
                this.ligand_like_polymers = [...this.ligand_like_polymers, chain.auth_asym_id];
            }
        }
        if (!this.__structure.ligands) {
            console.log("No ligands found");
        }
        else {
            for (var lig of this.__structure.ligands) {
                this.ligands = [...this.ligands, lig.chemicalId];
            }
        }
    }
    __verify_cif() {
        if ((0, fs_1.existsSync)(this.cif_filepath)) {
            return true;
        }
        else
            return false;
    }
    __verify_cif_modified() {
        if ((0, fs_1.existsSync)(this.cif_modified_filepath)) {
            return true;
        }
        else
            return false;
    }
    __verify_json_profile() {
        if ((0, fs_1.existsSync)(this.cif_modified_filepath)) {
            return true;
        }
        else
            return false;
    }
    __verify_png_thumbnail() {
        if ((0, fs_1.existsSync)(this.png_thumbnail_filepath)) {
            return true;
        }
        else
            return false;
    }
    __verify_chains_folder() {
        if ((0, fs_1.existsSync)(this.chains_folder)) {
            return true;
        }
        else
            return false;
    }
    //verify that each chain file exists
    __verify_chain_files() {
        var _a;
        if (!this.__verify_chains_folder()) {
            return false;
        }
        const chain_files_all = [
            ...(_a = this.__structure.proteins) === null || _a === void 0 ? void 0 : _a.map((c) => { return c.auth_asym_id; }),
            ...(this.__structure.rnas || []).map((c) => { return c.auth_asym_id; })
        ].map((chain_id) => {
            if (!(0, fs_1.existsSync)(`${this.chains_folder}/${this.rcsb_id}_STRAND_${chain_id}.cif`)) {
                console.log(`[${this.rcsb_id}]: NOT FOUND ${this.chains_folder}/${this.rcsb_id}_STRAND_${chain_id}.cif`);
                return false;
            }
            else
                return true;
        }).reduce((prev, cur) => { return prev && cur; }, true);
        return chain_files_all;
    }
    __verify_ligands_and_polymers(try_fix = false) {
        let ligs = this.ligands && this.ligands.map((lig_chem_id) => {
            if (!(0, fs_1.existsSync)(`${this.folder_path}/LIGAND_${lig_chem_id}.json`)) {
                console.log(`[${this.rcsb_id}]: NOT FOUND ${this.folder_path}/LIGAND_${lig_chem_id}.json (Is it an ION? Expected.)`);
                return false;
            }
            else
                return true;
        }).reduce((prev, cur) => { return prev && cur; }, true);
        let polys = this.ligand_like_polymers && this.ligand_like_polymers.map((polymer_id) => {
            if (!(0, fs_1.existsSync)(`${this.folder_path}/POLYMER_${polymer_id}.json`)) {
                console.log(`[${this.rcsb_id}]: NOT FOUND ${this.folder_path}/POLYMER_${polymer_id}.json`);
                return false;
            }
            else
                return true;
        });
        if ((!polys || !ligs) && try_fix) {
            console.log("Some ligands are missing. Calling script:", process.env["PYTHONBIN"], process.env["EXTRACT_BSITES_PY"]);
            (0, child_process_1.exec)(`${process.env["PYTHONBIN"]} ${process.env["EXTRACT_BSITES_PY"]} -s ${this.rcsb_id} --save`, (err, stdout, stderr) => {
                console.log(err);
                console.log(stdout);
                console.log(stderr);
            });
        }
    }
    async assets_verify(try_fix = false) {
        if (!this.__verify_cif()) {
            console.log(`[${this.rcsb_id}]: NOT FOUND ${this.cif_filepath}`);
            if (try_fix) {
                download_unpack_place(this.rcsb_id);
            }
        }
        if (!this.__verify_cif_modified()) {
            console.log(`[${this.rcsb_id}]: NOT FOUND ${this.cif_modified_filepath}`);
            if (try_fix) {
                let y = (0, child_process_1.exec)(`${process.env["PYTHONBIN"]} ${process.env["SPLIT_RENAME_PY"]} -s ${this.rcsb_id}`, (err, stdout, stderr) => {
                    console.log(err);
                    console.log(stdout);
                    console.log(stderr);
                });
                console.log(y.stdout);
            }
        }
        if (!this.__verify_json_profile()) {
            if (try_fix) {
                let ribosome = await (0, structure_json_profile_1.processPDBRecord)(this.rcsb_id);
                let filename = await (0, exports.save_struct_profile)(ribosome);
                process.stdout.write(`Saved structure profile:\t${filename}`);
            }
            console.log(`[${this.rcsb_id}]: NOT FOUND ${this.json_profile_filepath}`);
        }
        if (!this.__verify_png_thumbnail()) {
            console.log(`[${this.rcsb_id}]: NOT FOUND ${this.png_thumbnail_filepath}`);
            if (try_fix) {
                let proc = (0, child_process_1.exec)(`${process.env["PYTHONBIN"]} ${process.env["RENDER_THUMBNAIL_PY"]} -s ${this.rcsb_id}`, (err, stdout, stderr) => {
                    console.log(err);
                    console.log(stdout);
                    console.log(stderr);
                });
            }
        }
        if (!this.__verify_chain_files()) {
            if (try_fix) {
                (0, child_process_1.exec)(`${process.env["PYTHONBIN"]} ${process.env["SPLIT_RENAME_PY"]} -s ${this.rcsb_id}`, (err, stdout, stderr) => {
                    console.log(err);
                    console.log(stdout);
                    console.log(stderr);
                });
            }
        }
        this.__verify_ligands_and_polymers(try_fix);
    }
    db_verify() {
        // TODO
    }
}
exports.StructureFolder = StructureFolder;
/**
 * Request and display the state of the given rcsb_id structure in the database instance.
 */
const queryStructDb = (rcsb_id) => {
    return new Promise((resolve, reject) => {
        let y = (0, child_process_1.exec)(`echo \"match (struct:RibosomeStructure {rcsb_id:\\"${rcsb_id.toUpperCase()}\\"}) return struct.rcsb_id\" | cypher-shell -a \"${process.env["NEO4J_URI"]}\" --format plain -u ${process.env["NEO4J_USER"]} -p ${process.env["NEO4J_PASSWORD"]} --database ${process.env["NEO4J_CURRENTDB"]}`, { env: process.env }, (err, stdout, stderr) => {
            if (err && (err === null || err === void 0 ? void 0 : err.code) != 0) {
                process.stdout.write("Got shell error " + stderr + stdout);
                console.log("Got Error code:", err === null || err === void 0 ? void 0 : err.code);
                reject(err);
            }
            const dbstructs = stdout != null ? stdout.replace(/"/g, '').split("\n").filter(r => r.length === 4) : [];
            console.log(dbstructs);
            resolve(dbstructs);
        });
    });
};
/**
 * Download a .cif model of the structure.
 * @param struct_id
 */
const download_unpack_place = async (struct_id) => {
    const BASE_URL = "http://files.rcsb.org/download/";
    const FORMAT = ".cif.gz";
    const structid = struct_id.toUpperCase();
    let url = BASE_URL + structid + FORMAT;
    let compressed = await axios_1.default.get(url, { responseType: 'arraybuffer' }).then(r => { return r.data; })
        .catch(e => { process.stdout.write(`Structure ${structid} failed: `, e); return []; });
    let decompressed = await (0, node_gzip_1.ungzip)(compressed);
    // let destination_chains = path.join(
    //   process.env["RIBETL_DATA"] as string,
    //   `${structid}`,
    //   `CHAINS`)
    // if (!existsSync(destination_chains)) {
    //   mkdirSync(destination_chains)
    //   process.stdout.write(`Created directory ${destination_chains}.`);
    // }
    let structfile = path.join(process.env["RIBETL_DATA"], `${structid}`, `${structid}.cif`);
    (0, fs_2.writeFileSync)(structfile, decompressed);
};
const save_struct_profile = (r) => {
    var rcsb_id = r.rcsb_id;
    var target_filename = path.join(process.env["RIBETL_DATA"], rcsb_id.toUpperCase(), rcsb_id.toUpperCase() + ".json");
    if (!(0, fs_1.existsSync)(path.dirname(target_filename))) {
        (0, fs_2.mkdirSync)(path.dirname(target_filename));
    }
    (0, fs_2.writeFileSync)(target_filename, JSON.stringify(r, null, 4));
    return target_filename;
};
exports.save_struct_profile = save_struct_profile;
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
