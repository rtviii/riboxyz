import { Command, Flags } from '@oclif/core'
import { BaseCommand } from '../..'

export default class Db extends BaseCommand {
  static description = 'describe the command here'
  static flags = {
    force    : Flags.boolean({ char: 'f' }),
    // neo4j_uri: Flags.string({
    //   multiple: false,                                                 // allow setting this flag multiple times
    //   env     : 'NEO4J_URI',                                           // default to value of environment variable
    //   parse   : async (input) => (process.env["NEO4J_URI"] = input),   // instead of the user input, return a different value
    // }),
    // neo4j_password: Flags.string({
    //   multiple: false,                                                      // allow setting this flag multiple times
    //   env     : 'NEO4J_PASSWORD',                                           // default to value of environment variable
    //   parse   : async (input) => (process.env["NEO4J_PASSWORD"] = input),   // instead of the user input, return a different value
    // }),
    // neo4j_currentdb: Flags.string({
    //   multiple: false,                                                       // allow setting this flag multiple times
    //   env     : 'NEO4J_CURRENTDB',                                           // default to value of environment variable
    //   parse   : async (input) => (process.env["NEO4J_CURRENTDB"] = input),   // instead of the user input, return a different value
    // }),
    // neo4j_user: Flags.string({
    //   multiple: false,                                                  // allow setting this flag multiple times
    //   env     : 'NEO4J_USER',                                           // default to value of environment variable
    //   parse   : async (input) => (process.env["NEO4J_USER"] = input),   // instead of the user input, return a different value
    // })
  }

  // static args = [{name: 'file'}]

  public async run(): Promise<void> {
    const { args, flags } = await this.parse(Db)

    this.log(`hello from /home/rxz/dev/docker_ribxz/__ingress/ribxz/src/commands/db.ts`)
    if (args.file && flags.force) {
      this.log(`you input --force and --file: ${args.file}`)
    }
  }
}
