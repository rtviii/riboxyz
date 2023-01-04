// https://search.rcsb.org/#introduction
import axios from "axios";
import { Command, Flags } from '@oclif/core'




export class DbCommand  extends Command {

    public neo4j_vars = ["NEO4J_URI", "NEO4J_USER", "NEO4J_PASSWORD", "NEO4J_CURRENTDB", "RIBETL_DATA"]

    constructor(argv: string[], config: any) {
        super(argv, config);
    }

    static args = []

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
        env: Flags.string({
            char: 'e',
            description: 'Environment variable',
        }),
    };

    async run() {
        let { args, flags } = await this.parse(DbCommand)
        for (var v of ["NEO4J_URI", "NEO4J_USER", "NEO4J_PASSWORD", "NEO4J_CURRENTDB", "RIBETL_DATA"]) {
            if (Object.keys(flags).includes(v)) {
                process.env[v] = flags[v]
            }
        }
        if (flags.env) {
            const varstring = flags.env
            process.env[varstring.split("=")[0]] = varstring.split("=")[1]
        }
        for (var ee of this.neo4j_vars) {
            this.log(`env: ${ee} = ${process.env[ee]}`);
        }
    }

}