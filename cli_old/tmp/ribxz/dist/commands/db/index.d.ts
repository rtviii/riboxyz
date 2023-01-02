import { BaseCommand } from '../..';
export default class Db extends BaseCommand {
    static description: string;
    static flags: {
        force: import("@oclif/core/lib/interfaces").BooleanFlag<boolean>;
    };
    run(): Promise<void>;
}
export declare const missing_structures: () => Promise<string[] | never[]>;
