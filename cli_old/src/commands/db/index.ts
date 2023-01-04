import { Command, Flags } from '@oclif/core'
import { BaseCommand } from '../..'
// https://search.rcsb.org/#introduction
import axios from "axios";
import {exec, config} from "shelljs";

export default class Db extends BaseCommand {
  static description = 'describe the command here'
  static flags = {
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


export const missing_structures = async () => {
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
        config.silent = true
        exec(`echo \"match (struct:RibosomeStructure) return struct.rcsb_id\" | cypher-shell -a \"${process.env["NEO4J_URI"]}\" --format plain -u ${process.env.NEO4J_USER} -p ${process.env.NEO4J_PASSWORD} --database ${process.env.NEO4J_CURRENTDB}`,
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