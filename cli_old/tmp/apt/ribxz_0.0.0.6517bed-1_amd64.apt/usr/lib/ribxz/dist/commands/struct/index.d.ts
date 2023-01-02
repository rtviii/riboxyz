import { Command } from '@oclif/core';
export declare class StructureCommand extends Command {
    neo4j_vars: string[];
    constructor(argv: string[], config: any);
    static args: {
        name: string;
        required: boolean;
    }[];
    static globalFlags: {
        RIBETL_DATA: import("@oclif/core/lib/interfaces").OptionFlag<string | undefined>;
        NEO4J_URI: import("@oclif/core/lib/interfaces").OptionFlag<string | undefined>;
        NEO4J_PASSWORD: import("@oclif/core/lib/interfaces").OptionFlag<string | undefined>;
        NEO4J_USER: import("@oclif/core/lib/interfaces").OptionFlag<string | undefined>;
        NEO4J_CURRENTDB: import("@oclif/core/lib/interfaces").OptionFlag<string | undefined>;
        env: import("@oclif/core/lib/interfaces").OptionFlag<string | undefined>;
    };
    run(): Promise<void>;
}
