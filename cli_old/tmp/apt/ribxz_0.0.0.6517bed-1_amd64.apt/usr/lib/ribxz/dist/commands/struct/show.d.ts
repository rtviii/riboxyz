import { BaseCommand } from '../..';
import { RibosomeStructure } from '../../structure_processing/ribosome_types';
export default class Show extends BaseCommand {
    static description: string;
    static flags: {
        files: import("@oclif/core/lib/interfaces").BooleanFlag<boolean>;
        db: import("@oclif/core/lib/interfaces").BooleanFlag<boolean>;
        repair: import("@oclif/core/lib/interfaces").BooleanFlag<boolean>;
        commit: import("@oclif/core/lib/interfaces").BooleanFlag<boolean>;
        dryrun: import("@oclif/core/lib/interfaces").BooleanFlag<boolean>;
        force: import("@oclif/core/lib/interfaces").BooleanFlag<boolean>;
    };
    static args: {
        name: string;
        required: boolean;
    }[];
    run(): Promise<void>;
}
export declare class StructureFolder {
    folder_path: string;
    private cif_filepath;
    private cif_modified_filepath;
    private json_profile_filepath;
    private chains_folder;
    private png_thumbnail_filepath;
    private rcsb_id;
    __structure: RibosomeStructure;
    private ligands;
    private ligand_like_polymers;
    constructor(rcsb_id: string);
    private __verify_cif;
    private __verify_cif_modified;
    private __verify_json_profile;
    private __verify_png_thumbnail;
    private __verify_chains_folder;
    private __verify_chain_files;
    private __verify_ligands_and_polymers;
    assets_verify(try_fix?: boolean): Promise<void>;
    db_verify(): void;
}
export declare const save_struct_profile: (r: RibosomeStructure) => string;
