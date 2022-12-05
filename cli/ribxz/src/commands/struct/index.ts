import {Command, Flags} from '@oclif/core'
import { BaseCommand } from '../..'

export default class Structure extends BaseCommand {
  static description = 'Query structure in the database'

  static flags = {
    dryrun: Flags.boolean(),
    // flag with no value (-f, --force)
    force: Flags.boolean({char: 'f'}),
  }

  static args = [{name: 'rcsb_id'}]

  public async run(): Promise<void> {
    const {args, flags} = await this.parse(Structure)
    const dryrun        = flags.dryrun
  }

}
