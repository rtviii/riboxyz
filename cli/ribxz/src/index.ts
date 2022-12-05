import { Command, Flags } from '@oclif/core'


export { run } from '@oclif/core'


export abstract class BaseCommand extends Command {
    constructor(argv: string[], config: any) {
        super(argv, config);
        if (!process.env.NEO4J_URI) {

        }
        console.log("Process nev : ", process.env.NEO4J_URI)
        this.log("Process nev : ", process.env.NEO4J_URI)
    }

    static globalFlags = {
        //add environment variables need to log into neo4j and use environment by default from system
        NEO4J_URI: Flags.string({
            char:'a',
            multiple: false,                             // allow setting this flag multiple times
            env: 'NEO4J_URI',                              // default to value of environment variable
            parse: async (input) => (process.env["NEO4J_URI"] = input),              // instead of the user input, return a different value
        }),

        NEO4J_PASSWORD: Flags.string({
            char:'p',
            multiple: false,                             // allow setting this flag multiple times
            env: 'NEO4J_PASSWORD',                              // default to value of environment variable
            parse: async (input) => (process.env["NEO4J_PASSWORD"] = input),              // instead of the user input, return a different value
        }),

        NEO4J_USER: Flags.string({
            char:'u',
            multiple: false,                             // allow setting this flag multiple times
            env: 'NEO4J_USER',                              // default to value of environment variable
            parse: async (input) => (process.env["NEO4J_USER"] = input),              // instead of the user input, return a different value
        }),

        NEO4J_CURRENTDB: Flags.string({
            char:'d',
            multiple: false,                             // allow setting this flag multiple times
            env: 'NEO4J_CURRENTDB',                              // default to value of environment variable
            parse: async (input) => (process.env["NEO4J_CURRENTDB"] = input),              // instead of the user input, return a different value
        }),
        env: Flags.string({
            char: 'e',
            description: 'Environment variable',
            parse: async (input: string) => {
                process.env[input.split("=")[0]] = input.split("=")[1]
                return input
            }
        }),
    };
}
