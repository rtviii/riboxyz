import { Command, Flags } from '@oclif/core'
import { exec } from 'child_process';
export { run } from '@oclif/core'

/**                               Assumptions for updating the db.
 * ---------------------------------------------------------------------------------------
 * The ingress scripts can be ran from anywhere and need not access the database server.
 * The script has access to
 * -    the internet
 * -    locally installed `cypher-shell`
 * -    the database remote: $NEO4J_URI,$NEO4J_PASSWORD, $NEO4J_USER and $NEO4J_CURRENTDB
 * -    the (local) filesystem containing the parsed structure-profiles: $RIBETL_DATA
 * 
 * ---------------------------------------------------------------------------------------
 * 
 * To access the cypher-shell on a the db host (ex.):
 * 
 * $ echo "match (n:RibosomeStructure) return n.rcsb_id;" | cypher-shell    /
 * -a "neo4j://ribosome.xyz:7687"                                           /
 * -u 'rt'                                                                  /
 * -p 'rrr'                                                                 /
 * --database 'riboauth'                                                    /
 * --format plain 
 * 
 * TODO:1. The script begins by verifying that the structure is not already in the database. 
*/


export abstract class BaseCommand extends Command {


    public neo4j_vars = ["NEO4J_URI", "NEO4J_USER", "NEO4J_PASSWORD", "NEO4J_CURRENTDB", "RIBETL_DATA"]
    constructor(argv: string[], config: any) {
        const       SCRIPTS_DIR            = "/home/rxz/dev/docker_ribxz/cli/scripts"
        process.env["EXTRACT_BSITES_PY"]   = `${SCRIPTS_DIR}/extract_bsites.py`
        process.env["RENDER_THUMBNAIL_PY"] = `${SCRIPTS_DIR}/render_thumbnail.py`
        process.env["COMMIT_STRUCTURE_SH"] = `${SCRIPTS_DIR}/commit_structure.sh`
        process.env["SPLIT_RENAME_PY"]     = `${SCRIPTS_DIR}/split_rename.py`

        super(argv, config);
    }
    static globalFlags = {
        RIBETL_DATA: Flags.string({
            char: 'r',
            multiple: false,           // allow setting this flag multiple times
            env: 'RIBETL_DATA',   // default to value of environment variable
        }),
        //add environment variables need to log into neo4j and use environment by default from system
        NEO4J_URI: Flags.string({
            char: 'a',
            multiple: false,                             // allow setting this flag multiple times
            env: 'NEO4J_URI',                              // default to value of environment variable
        }),

        NEO4J_PASSWORD: Flags.string({
            char: 'p',
            multiple: false,                             // allow setting this flag multiple times
            env: 'NEO4J_PASSWORD',                              // default to value of environment variable
        }),

        NEO4J_USER: Flags.string({
            char: 'u',
            multiple: false,                             // allow setting this flag multiple times
            env: 'NEO4J_USER',                              // default to value of environment variable
        }),
        NEO4J_CURRENTDB: Flags.string({
            char: 'd',
            multiple: false,                             // allow setting this flag multiple times
            env: 'NEO4J_CURRENTDB',                              // default to value of environment variable
        }),
        PYTHONBIN: Flags.string({
            multiple: false,                             // allow setting this flag multiple times
            env: 'PYTHONBIN',                              // default to value of environment variable
        }),
        COMMIT_STRUCTURE_SH: Flags.string({
            multiple: false,                             // allow setting this flag multiple times
            env: 'COMMIT_STRUCTURE_SH',                              // default to value of environment variable
        }),
        EXTRACT_BSITES_PY: Flags.string({
            multiple: false,                             // allow setting this flag multiple times
            env: 'EXTRACT_BSITES_PY',                              // default to value of environment variable
        }),
        RENDER_THUMBNAIL_PY: Flags.string({
            multiple: false,                             // allow setting this flag multiple times
            env: 'RENDER_THUMBNAIL_PY',                              // default to value of environment variable
        }),
        SPLIT_RENAME_PY: Flags.string({
            multiple: false,                             // allow setting this flag multiple times
            env: 'SPLIT_RENAME_PY',                              // default to value of environment variable
        }),
        env: Flags.string({
            char: 'e',
            description: 'Environment variable',
        }),
    };


    async run() {
        let { args, flags } = await this.parse(BaseCommand)

        console.log(`Base args :${args}`, args)
        for (var v of ["NEO4J_URI", "NEO4J_USER", "NEO4J_PASSWORD", "NEO4J_CURRENTDB", "RIBETL_DATA"]) {
            if (Object.keys(flags).includes(v)) {
                process.env[v] = flags[v]
            }
        }
        for (var v of ["PYTHONBIN", "RENDER_THUMBNAIL_PY", "SPLIT_RENAME_PY", "COMMIT_STRUCTURE_SH", "EXTRACT_BSITES_PY"]) {
            if (Object.keys(flags).includes(v)) {
                process.env[v] = flags[v]
            }
        }

        // if (flags.env) {
        //     const varstring = flags.env
        //     process.env[varstring.split("=")[0]] = varstring.split("=")[1]
        // }

        for (var ee of this.neo4j_vars) {
            this.log(`env: ${ee} = ${process.env[ee]}`);
        }
    }

}

