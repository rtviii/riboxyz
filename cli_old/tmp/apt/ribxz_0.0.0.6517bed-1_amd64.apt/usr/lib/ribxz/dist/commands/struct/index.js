"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.StructureCommand = void 0;
const core_1 = require("@oclif/core");
class StructureCommand extends core_1.Command {
    constructor(argv, config) {
        super(argv, config);
        this.neo4j_vars = ["NEO4J_URI", "NEO4J_USER", "NEO4J_PASSWORD", "NEO4J_CURRENTDB", "RIBETL_DATA"];
    }
    async run() {
        let { args, flags } = await this.parse(StructureCommand);
        console.log(`struct args :${args}`, args);
        for (var v of ["NEO4J_URI", "NEO4J_USER", "NEO4J_PASSWORD", "NEO4J_CURRENTDB", "RIBETL_DATA"]) {
            if (Object.keys(flags).includes(v)) {
                process.env[v] = flags[v];
            }
        }
        if (flags.env) {
            const varstring = flags.env;
            process.env[varstring.split("=")[0]] = varstring.split("=")[1];
        }
        for (var ee of this.neo4j_vars) {
            this.log(`env: ${ee} = ${process.env[ee]}`);
        }
    }
}
exports.StructureCommand = StructureCommand;
StructureCommand.args = [{
        name: 'rcsb_id',
        required: true
    }];
StructureCommand.globalFlags = {
    RIBETL_DATA: core_1.Flags.string({
        char: 'r',
        multiple: false,
        env: 'RIBETL_DATA', // default to value of environment variable
    }),
    //add environment variables need to log into neo4j and use environment by default from system
    NEO4J_URI: core_1.Flags.string({
        char: 'a',
        multiple: false,
        env: 'NEO4J_URI', // default to value of environment variable
    }),
    NEO4J_PASSWORD: core_1.Flags.string({
        char: 'p',
        multiple: false,
        env: 'NEO4J_PASSWORD', // default to value of environment variable
    }),
    NEO4J_USER: core_1.Flags.string({
        char: 'u',
        multiple: false,
        env: 'NEO4J_USER', // default to value of environment variable
    }),
    NEO4J_CURRENTDB: core_1.Flags.string({
        char: 'd',
        multiple: false,
        env: 'NEO4J_CURRENTDB', // default to value of environment variable
    }),
    env: core_1.Flags.string({
        char: 'e',
        description: 'Environment variable',
    }),
};
