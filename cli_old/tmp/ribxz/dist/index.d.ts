import { Command } from '@oclif/core';
export { run } from '@oclif/core';
/**                               Assumptions for updating the db.
 * ---------------------------------------------------------------------------------------
 * The ingress scripts can be ran from anywhere and need not access the database server.
 * The script has access to
 * -    the internet
 * -    locally installed `cypher-shell`
 * -    the database remote: $NEO4J_URI,$NEO4J_PASSWORD, $NEO4J_USER and $NEO4J_CURRENTDB
 * -    the (local) filesystem containing the parsed structure-profiles: $RIBETL_DATA
 *
 * ---------------------------------------------------------------------------------------
 *
 * To access the cypher-shell on a the db host (ex.):
 *
 * $ echo "match (n:RibosomeStructure) return n.rcsb_id;" | cypher-shell    /
 * -a "neo4j://ribosome.xyz:7687"                                           /
 * -u 'rt'                                                                  /
 * -p 'rrr'                                                                 /
 * --database 'riboauth'                                                    /
 * --format plain
 *
 * TODO:1. The script begins by verifying that the structure is not already in the database.
*/
export declare abstract class BaseCommand extends Command {
    neo4j_vars: string[];
    constructor(argv: string[], config: any);
    static globalFlags: {
        RIBETL_DATA: import("@oclif/core/lib/interfaces").OptionFlag<string | undefined>;
        NEO4J_URI: import("@oclif/core/lib/interfaces").OptionFlag<string | undefined>;
        NEO4J_PASSWORD: import("@oclif/core/lib/interfaces").OptionFlag<string | undefined>;
        NEO4J_USER: import("@oclif/core/lib/interfaces").OptionFlag<string | undefined>;
        NEO4J_CURRENTDB: import("@oclif/core/lib/interfaces").OptionFlag<string | undefined>;
        PYTHONBIN: import("@oclif/core/lib/interfaces").OptionFlag<string | undefined>;
        COMMIT_STRUCTURE_SH: import("@oclif/core/lib/interfaces").OptionFlag<string | undefined>;
        EXTRACT_BSITES_PY: import("@oclif/core/lib/interfaces").OptionFlag<string | undefined>;
        RENDER_THUMBNAIL_PY: import("@oclif/core/lib/interfaces").OptionFlag<string | undefined>;
        SPLIT_RENAME_PY: import("@oclif/core/lib/interfaces").OptionFlag<string | undefined>;
        env: import("@oclif/core/lib/interfaces").OptionFlag<string | undefined>;
    };
    run(): Promise<void>;
}
