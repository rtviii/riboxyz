import 'yargs'
import yargs from 'yargs'
import { CmdStructure } from './cmd_structure'
//create a yargs command "structure" that takes the a positional argument 'rcsb_id' and the following flags (options): '--acquirePDBRecord', '--downloadCifFile', '--splitRenameChains', '--renderStructHero', '--extractBindingSites', '--commitToNeo4j'

export interface CmdStructureLs extends CmdStructure{
    db    : boolean,
    assets: boolean
} 

exports.command  = 'ls <rcsb_id>'
exports.describe = 'list structure'

// introduce a subcommand "ls" that prints the name of the structure
exports.builder = (yargs:CmdStructureLs) => {
    return yargs
        .option('db', {
            alias   : 'db',
            describe: 'list db',
            type    : 'boolean',
            default : false
        })
}

exports.handler = async function handler(argv: CmdStructureLs) {
    console.log("Hellllloo");
    console.log(
        "Got", argv.rcsb_id
    )
}