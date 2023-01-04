import { Flags } from '@oclif/core'
import { BaseCommand, queryCypher } from '../..'
import { existsSync, readFileSync } from 'fs'
import { RibosomeStructure } from '../../structure_processing/structure_types'
// https://search.rcsb.org/#introduction
import axios from "axios";
import { ungzip } from 'node-gzip'
import { mkdirSync, writeFileSync } from 'fs'
import { processPDBRecord } from "../../structure_processing/structure_json_profile";
import path = require("path");
import { StructureFolder } from '../struct/show'
import { exec, config } from "shelljs";

export default class Add extends BaseCommand {

    static description = 'Query structure in the database'
    static flags = {
        // repair: Flags.boolean({ char: 'R' }),
    }

    static args = [
        { name: 'rcsb_id', required: true }
    ]

    public async run(): Promise<void> {
        const { args, flags } = await this.parse(Add)

        // console.log(await queryCypher("match (n:RibosomeStructure) return count(n)"))
        let out = await queryCypher("match (n:RibosomeStructure) return count(n)")
    }
}
