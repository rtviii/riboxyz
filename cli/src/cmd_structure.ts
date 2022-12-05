import 'yargs'
import yargs, { ArgumentsCamelCase, Argv, BuilderCallback } from 'yargs'
//create a yargs command "structure" that takes the a positional argument 'rcsb_id' and the following flags (options): '--acquirePDBRecord', '--downloadCifFile', '--splitRenameChains', '--renderStructHero', '--extractBindingSites', '--commitToNeo4j'

export interface CmdStructure extends yargs.Argv {
    rcsb_id: string,
    acquirePDBRecord: boolean,
    downloadCifFile: boolean,
    splitRenameChains: boolean,
    renderStructHero: boolean,
    extractBindingSites: boolean,
    commitToNeo4j: boolean
}

exports.command = 'structure <rcsb_id>'

exports.describe = 'Acquire and process PDB structure'


// introduce a subcommand "ls" that prints the name of the structure
exports.builder = (yargs: CmdStructure) => {
    return yargs
        .option('acquirePDBRecord', {
            alias: 'a',
            describe: 'Acquire the PDB record from the RCSB PDB API',
            type: 'boolean',
            default: false
        })
        .option('downloadCifFile', {
            alias: 'd',
            describe: 'Download the CIF file from the RCSB PDB FTP server',
            type: 'boolean',
            default: false
        })
        .option('splitRenameChains', {
            alias: 's',
            describe: 'Split and rename chains',
            type: 'boolean',
            default: false
        })
        .option('renderStructHero', {
            alias: 'r',
            describe: 'Render a StructHero image',
            type: 'boolean',
            default: false
        })
        .option('extractBindingSites', {
            alias: 'e',
            describe: 'Extract binding sites',
            type: 'boolean',
            default: false
        })
        .option('commitToNeo4j', {
            alias: 'c',
            describe: 'Commit to Neo4j',
            type: 'boolean',
            default: false
        })
}



exports.subcommands = 

exports.handler = async function handler(argv: CmdStructure) {
    console.log("Hellllloo");

    console.log(
        "Got", argv.rcsb_id
    )
}