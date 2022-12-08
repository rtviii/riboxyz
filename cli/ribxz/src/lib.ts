#!/usr/bin/env ts-node-script
// https://search.rcsb.org/#introduction
import axios from "axios";
import { gzip, ungzip } from 'node-gzip'
import fs from 'fs'
import yargs from 'yargs'
import path from "path";
import shell from "shelljs";
import { processPDBRecord } from "./rcsb_data_api/requestGqlProfile";
import { RibosomeStructure } from "./RibosomeTypes";
import _ from "lodash";




// Options describing how to process a given structure
interface IngressOptions {
    rcsb_id?            : string,
    acquirePDBRecord?   : boolean;
    downloadCifFile?    : boolean;
    splitRenameChains?  : boolean;
    renderStructHero?   : boolean;
    extractBindingSites?: boolean;
    commitToNeo4j?      : boolean;
}

const check_env_var = () => {
    let MUST_SET = [
        "RIBETL_DATA",

        "NEO4J_URI",
        "NEO4J_CURRENTDB",
        "NEO4J_USER",
        "NEO4J_PASSWORD",

        // Custom scripts involved in a the processing: should be found in the __ingress/scripts directory
        "EXTRACT_BSITES_PY", // Parse out the ligands and other non-ribosomal molecules, currently:"bsite_mixed.py")
        "RENDER_THUMBNAIL_PY", // Render a thumbnail of the structure, currently:"render_thumbnail.py")
        "SPLIT_RENAME_PY", // Extract individual chains of the structre into separate files in $RIBETL_DATA/$PDBID/CHAINS and a ban-nomenclature table to the .cif file.
        "COMMIT_STRUCTURE_SH", // Update the database (shell piping cypher into the db, currently "neo4j_commit_structure.sh")
    ]

    for (var envvar of MUST_SET) {
        if (!process.env[envvar]) {
            process.stdout.write(`Environment variable is not set: ${envvar}. (Use "--env ${envvar}=XXX")`);
        }
    }
}


// class DbUpdate {
//     state         : "scheduled" | "running" | "done";
//     completed     : boolean;
//     new_structures: string[];

//     constructor(struct_list: string[]) {
//     }
//     commit_structures(): void {

//     };
// }


const main = async () => {


    // Development defaults: 
    process.env["NEO4J_USER"]      = "neo4j"
    process.env["NEO4J_PASSWORD"]  = "rrr"
    process.env["NEO4J_CURRENTDB"] = "neo4j"
    process.env["NEO4J_URI"]       = "bolt://0.0.0.0:7687"

    process.env["RIBETL_DATA"] = "/home/rxz/dev/static"

    process.env["EXTRACT_BSITES_PY"]   = "/home/backend/ingress/scripts/extract_bsites.py"
    process.env["RENDER_THUMBNAIL_PY"] = "/home/backend/ingress/scripts/render_thumbnail.py"
    process.env["COMMIT_STRUCTURE_SH"] = "/home/backend/ingress/scripts/commit_structure.sh"
    process.env["SPLIT_RENAME_PY"]     = "/home/backend/ingress/scripts/split_rename.py"


    // add a help message if command received less than 1 subcommand
    const args = yargs(process.argv.slice(2))
        // .demandCommand(1, 'You need at least one command before moving on.')
        .command(require('./cmd_structure'))
        .command('check', 'Verify the health of the database', () => { }, (argv) => {
            console.log("Health check")
        })
        .options({
            structure: { type: "string", demandOption: false },
            pythonbin: { type: "string", demandOption: false, default: "/usr/bin/python3" },
            dryrun   : { type: "boolean", demandOption: false },
        })
        .option('env', {
            type: 'array',
            alias: 'e',
            description: `Add "env" argument to yargs that takes any number of key-value pairs: e.g. --env NEO4J_URI="bolt://localhost:7687"`,
            coerce: (arg: string[]) => {
                arg.forEach(kvp => {
                    const [key, value] = kvp.split('=');
                    process.env[key] = value;
                    return { [key]: value };
                })
            },
        })
        // .boolean('update_all_missing')
        .parseAsync()

    check_env_var()

    // if (args.dryrun) {
    //     process.stdout.write("Dry run. No changes will be made to the database.")
    //     let missing = await missing_structures()
    //     missing.forEach(m => process.stdout.write(m))
    //     process.stdout.write("Finished parser")
    // }

    // if (args.update_all_missing) {
    //     let missing = await missing_structures()
    //     process.stdout.write(`Attempting to process ${missing.length} structures`)
    //     missing.forEach(async (struct_id) => {
    //         let processing_opts = Object.assign({}, IngressOptions, { rcsb_id: struct_id })
    //         await process_structure(processing_opts)
    //     })

    //     process.exit(1)
    // }


    // if (args.structure) {

    //     if (args.structure != null) {
    //         process.stdout.write("Enter a valid RCSB PDB ID to build from RCSB's gql api.");
    //         process.exit(2);
    //     }

    //     try {
    //         await process_structure({});
    //     } catch (e: any) {
    //         process.stdout.write("Failed: \n\n");
    //         process.stdout.write(e);
    //     }
    // }


}

const missing_structures = async () => {
    process.stdout.write("Getting missing structures");
    var rcsb_search_api = "https://search.rcsb.org/rcsbsearch/v2/query"
    const params = {
        "query":
        {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "group",
                    "logical_operator": "and",
                    "nodes": [
                        { "type": "group", "logical_operator": "and", "nodes": [{ "type": "terminal", "service": "text", "parameters": { "operator": "contains_phrase", "negation": false, "value": "RIBOSOME", "attribute": "struct_keywords.pdbx_keywords" } }] },
                        { "type": "group", "logical_operator": "and", "nodes": [{ "type": "terminal", "service": "text", "parameters": { "operator": "greater", "negation": false, "value": 25, "attribute": "rcsb_entry_info.polymer_entity_count_protein" } }] },
                        // { "type": "group", "logical_operator": "and", "nodes": [{ "type": "terminal", "service": "text", "parameters": { "operator": "less"           , "negation": false, "value": 20         , "attribute": "rcsb_entry_info.resolution_combined"          } }] }
                    ],
                    "label": "text"
                }
            ],
            "label": "query-builder"
        },
        "return_type": "entry",
        "request_options": {
            "return_all_hits": true,
            "results_verbosity": "compact"
        }
    };
    let query = rcsb_search_api + "?json=" + encodeURIComponent(JSON.stringify(params))

    // let cypherstring = "match (struct:RibosomeStructure) return struct.rcsb_id"
    // cypherstring = encodeURIComponent(cypherstring);
    // let ribxz_query = `http://localhost:8000/neo4j/cypher/?cypher=${cypherstring}`

    let dbquery = new Promise<string[]>((resolve, reject) => {
        shell.config.silent = true
        shell.exec(`echo \"match (struct:RibosomeStructure) return struct.rcsb_id\" | cypher-shell -a \"${process.env["NEO4J_URI"]}\" --format plain -u ${process.env.NEO4J_USER} -p ${process.env.NEO4J_PASSWORD} --database ${process.env.NEO4J_CURRENTDB}`,
            function (err, stdout, stderr) {
                if (err != 0) {
                    process.stdout.write("Got shell error")
                    reject(err)
                }
                const dbstructs = (stdout as string).replace(/"/g, '').split("\n").filter(r => r.length === 4)
                resolve(dbstructs)
            })
    })


    return Promise.all([dbquery, axios.get(query)]).then(r => {
        var ribxz_structs: string[] = r[0]
        var rcsb_structs: string[] = r[1].data.result_set

        var missing_from_ribxz = rcsb_structs.filter(struct => {
            if (!ribxz_structs.includes(struct)) {
                return true
            } else return false
        })
        process.stdout.write(`riboxyz contains ${ribxz_structs.length} structures. Up-to-date RCSB API contains ${rcsb_structs.length} structures.`)
        process.stdout.write(`Structs absent from ribosome.xyz: ${missing_from_ribxz.length} `,)
        return missing_from_ribxz
    })
        .catch(e => { process.stdout.write(`Got error: ${e}`); return [] })
}

const download_unpack_place = async (struct_id: string) => {
    const BASE_URL = "http://files.rcsb.org/download/"
    const FORMAT = ".cif.gz"

    const structid = struct_id.toUpperCase()
    let url = BASE_URL + structid + FORMAT
    let compressed: Buffer = await axios.get(url, { responseType: 'arraybuffer' }).then(r => { return r.data })
        .catch(e => { process.stdout.write(`Structure ${structid} failed: `, e); return []; })
    let decompressed = await ungzip(compressed);

    let destination_chains = path.join(
        process.env["RIBETL_DATA"] as string,
        `${structid}`,
        `CHAINS`)

    let structfile = path.join(
        process.env["RIBETL_DATA"] as string,
        `${structid}`,
        `${structid}.cif`)
    if (!fs.existsSync(destination_chains)) {
        fs.mkdirSync(destination_chains)
        process.stdout.write(`Created directory ${destination_chains}.`);
    }
    fs.writeFileSync(structfile, decompressed)
}

const save_struct_profile = (r: RibosomeStructure): string => {
    var rcsb_id = r.rcsb_id;
    var target_filename = path.join(
        process.env["RIBETL_DATA"] as string,
        rcsb_id.toUpperCase(),
        rcsb_id.toUpperCase() + ".json"
    );

    if (!fs.existsSync(path.dirname(target_filename))) {
        shell.mkdir("-p", path.dirname(target_filename));
    }
    fs.writeFileSync(target_filename, JSON.stringify(r, null, 4));
    return target_filename
};

const process_structure = async (opts: IngressOptions) => {
    if (!opts.rcsb_id) throw new Error("No structure ID provided.")

    let struct_id = opts.rcsb_id.toUpperCase()

    if (opts.acquirePDBRecord) {
        let ribosome = await processPDBRecord(struct_id)
        let filename = await save_struct_profile(ribosome)
        process.stdout.write(`Saved structure profile:\t${filename}`);
    }

    if (opts.downloadCifFile) {
        await download_unpack_place(struct_id)
        process.stdout.write(`Saved structure cif file ${struct_id}.cif`);
    }
    if (opts.renderStructHero) {
        shell.exec(`${process.env["PYTHONBIN"]} ${process.env["RENDER_THUMBNAIL_PY"]} -s ${struct_id}`)
    }

    if (opts.splitRenameChains) {
        shell.exec(`${process.env["PYTHONBIN"]} ${process.env["SPLIT_RENAME_PY"]} -s ${struct_id}`)
    }

    if (opts.extractBindingSites) {
        shell.exec(`${process.env["PYTHONBIN"]}   \
         ${process.env["EXTRACT_BSITES_PY"]}       \
         -s ${struct_id}                           \
         --save`)
    }

    if (opts.commitToNeo4j) {
        var shellstring = `${process.env["COMMIT_STRUCTURE_SH"]} -s ${opts.rcsb_id.toUpperCase()} \
         -d ${process.env.NEO4J_CURRENTDB} \
         -a "${process.env.NEO4J_URI}"`;
        shell.exec(shellstring)
    }

}

main()