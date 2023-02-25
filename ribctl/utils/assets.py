import json
import os
from pydantic import BaseModel, parse_obj_as
import gzip
import os
import requests
from render_thumbnail import render_thumbnail
from ribctl.types_ribosome import RibosomeStructure
from ribctl.utils.process_structure import process_pdb_record
from split_rename import split_rename


def download_unpack_place(struct_id: str) -> None:

    BASE_URL = "http://files.rcsb.org/download/"
    FORMAT = ".cif.gz"

    structid = struct_id.upper()
    url = BASE_URL + structid + FORMAT
    compressed = requests.get(url).content
    decompressed = gzip.decompress(compressed)

    destination_chains = os.path.join(
        os.environ["RIBETL_DATA"],
        structid,
        "CHAINS")

    if not os.path.exists(destination_chains):
        os.mkdir(destination_chains)
        print(f"Created directory {destination_chains}.")

    structfile = os.path.join(
        os.environ["RIBETL_DATA"],
        structid,
        structid + ".cif")

    with open(structfile, "wb") as f:
        f.write(decompressed)


class RibosomeAssets(BaseModel):
    rcsb_id: str

    def __init__(self, rcsb_id: str):
        self.rcsb_id = rcsb_id

    def envcheck(self):
        if not os.environ.get("RIBETL_DATA"):
            raise Exception(
                "RIBETL_DATA environment variable not set. Cannot access assets.")

    def dir_path(self):
        self.envcheck()
        return f"{os.environ.get('RIBETL_DATA')}/{self.rcsb_id}"

    def cif_filepath(self):
        self.envcheck()
        return f"{self.dir_path()}/{self.rcsb_id}.cif"

    def cif_modified_filepath(self):
        self.envcheck()
        return f"{self.dir_path()}/{self.rcsb_id}_modified.cif"

    def json_profile_filepath(self):
        self.envcheck()
        return f"{self.dir_path()}/{self.rcsb_id}.json"

    def chains_folder(self):
        self.envcheck()
        return f"{self.dir_path()}/CHAINS"

    def png_thumbnail_filepath(self):
        self.envcheck()
        return f"{self.dir_path()}/_ray_{self.rcsb_id}.png"
      
    @staticmethod
    def save_json_profile(filepath:str, profile:dict):
        with open(filepath, "w") as f:
            json.dump(profile, f)

    def __verify_cif(self, obtain: bool = False) -> bool:
        if os.path.exists(self.cif_filepath()):
            return True
        else:
            if obtain:
                download_unpack_place(self.rcsb_id)
                return True
            else:
                return False

    def __verify_cif_modified(self, obtain: bool = False) -> bool:
        if os.path.exists(self.cif_modified_filepath()):
            return True
        else:
            if obtain:
                split_rename(self.rcsb_id)
                return True
            else:
                return False

    def __verify_json_profile(self, obtain: bool = False) -> bool:
        if os.path.exists(self.json_profile_filepath()):
            return True
        else:
            if obtain:
                ribosome = process_pdb_record(self.rcsb_id)
                if not parse_obj_as(RibosomeStructure, ribosome):
                    raise Exception("Invalid ribosome structure profile.")
                  
                self.save_json_profile(self.json_profile_filepath(), ribosome)
                print(f"Saved structure profile:\t{self.json_profile_filepath()}")
                return True
            else:
                return False

    def __verify_png_thumbnail(self, obtain: bool = False) -> bool:
        if os.path.exists(self.png_thumbnail_filepath()):
            return True
        else:
            if obtain:
                print("Obtaning thumbnail...")
                render_thumbnail(self.rcsb_id)

                return True
            else:
                return False


# export class RibosomeAssets {
#     rcsb_id: string

#     constructor(
#         rcsb_id: string
#     ) {
#         this.rcsb_id = rcsb_id
#     }
#     envcheck = () => {
#         if (!process.env["RIBETL_DATA"]) {
#             throw Error("RIBETL_DATA environment variable not set. Cannot access assets.")
#         }
#     }

#     folder_path = () => {
#         this.envcheck()
#         return `${process.env["RIBETL_DATA"]}/${this.rcsb_id}`
#     }

#     cif_filepath = () => {
#         this.envcheck()
#         return `${this.folder_path()}/${this.rcsb_id}.cif`
#     }

#     cif_modified_filepath = () => {
#         this.envcheck()
#         return `${this.folder_path()}/${this.rcsb_id}_modified.cif`
#     }

#     json_profile_filepath = () => {
#         this.envcheck()
#         return `${this.folder_path()}/${this.rcsb_id}.json`
#     }

#     chains_folder = () => {
#         this.envcheck()
#         return `${this.folder_path()}/CHAINS`
#     }

#     png_thumbnail_filepath = () => {
#         this.envcheck()
#         return `${this.folder_path()}/_ray_${this.rcsb_id}.png`
#     }


#     async __verify_cif(obtain: boolean = false):Promise<boolean> {
#         if (existsSync(this.cif_filepath())) {
#             return true
#         } else {

#             if (obtain) {
#                 await download_unpack_place(this.rcsb_id); return true
#             } else return false
#         }
#     }
#     async __verify_cif_modified(obtain: boolean = false) {
#         if (existsSync(this.cif_modified_filepath())) { return true } else {
#             if (obtain) {
#                 this.__verify_cif(true)
#                 let y = exec(`${process.env["PYTHONBIN"]} ${process.env["SPLIT_RENAME_PY"]} -s ${this.rcsb_id}`,
#                     (err, stdout, stderr) => {
#                         console.log(err);
#                         console.log(stdout.toString());
#                         console.log(stderr.toString());
#                     })
#                 console.log(y.stdout)
#             } else return false
#         }
#     }
#     async __verify_json_profile(obtain: boolean = false) {
#         if (existsSync(this.json_profile_filepath())) { return true } else {
#             if (obtain) {
#                 let ribosome = await processPDBRecord(this.rcsb_id)
#                 let filename = await save_struct_profile(ribosome)
#                 process.stdout.write(`Saved structure profile:\t${filename}`);
#             } else return false
#         }
#     }
#     async __verify_png_thumbnail(obtain: boolean = false) {
#         if (existsSync(this.png_thumbnail_filepath())) { return true } else {
#             if (obtain) {
#                await new Promise<void>((rs,rj)=>{
#                exec(`${process.env["PYTHONBIN"]} ${process.env["RENDER_THUMBNAIL_PY"]} -s ${this.rcsb_id}`,
#                     (err, stdout, stderr) => {
#                         if (err !==0 ){
#                             rj(stderr)
#                         }
#                         else{rs()}
#                     })
#                })
#             } else return false
#         }
#     }
#     async __verify_chains_folder(obtain: boolean = false) {
#         if (existsSync(this.chains_folder())) { return true } else {
#             if (obtain) {
#                 this.__verify_cif(true)
#                 exec(`${process.env["PYTHONBIN"]} ${process.env["SPLIT_RENAME_PY"]} -s ${this.rcsb_id}`, (err, stdout, stderr) => {
#                     console.log(err);
#                     console.log(stdout.toString());
#                     console.log(stderr.toString());
#                 })
#             } else return false
#         }
#     }
#     //verify that each chain file exists
#     async __verify_chain_files(parent_structure: RibosomeStructure): Promise<boolean> {
#         if (!this.__verify_chains_folder()) {
#             return false
#         }
#         const chain_files_all = [
#             parent_structure.proteins?.map((c) => { return c.auth_asym_id }),
#             (parent_structure.rnas || []).map((c) => { return c.auth_asym_id })
#         ].map(
#             (chain_id) => {
#                 if (!existsSync(`${this.chains_folder()}/${this.rcsb_id}_STRAND_${chain_id}.cif`)) {
#                     console.log(`[${this.rcsb_id}]: NOT FOUND ${this.chains_folder()}/${this.rcsb_id}_STRAND_${chain_id}.cif`)
#                     return false
#                 } else return true
#             }
#         ).reduce((prev, cur) => { return prev && cur }, true)
#         return chain_files_all
#     }

#     async __verify_ligands_and_polymers(struct: RibosomeStructure) {
#         let ligs = struct.ligands && struct.ligands.map((lig_chem_id) => {
#             if (!existsSync(`${this.folder_path()}/LIGAND_${lig_chem_id}.json`)) {
#                 console.log(`[${this.rcsb_id}]: NOT FOUND ${this.folder_path()}/LIGAND_${lig_chem_id}.json (Is it an ION? Expected.)`)
#                 return false
#             } else return true
#         }).reduce((prev, cur) => { return prev && cur }, true)


#         let ligandlike: string[] = []
#         let ligands: string[] = []

#         for (var chain of [...struct.proteins, ...(struct.rnas || [])]) {
#             if (chain.ligand_like) {
#                 ligandlike = [...ligandlike, chain.auth_asym_id]
#             }
#         }
#         if (!struct.ligands) {
#         } else {
#             for (var lig of struct.ligands) {
#                 ligands = [...ligands, lig.chemicalId]
#             }
#         }


#         let polys = ligandlike.map((polymer_id) => {
#             if (!existsSync(`${this.folder_path()}/POLYMER_${polymer_id}.json`)) {
#                 console.log(`[${this.rcsb_id}]: NOT FOUND ${this.folder_path()}/POLYMER_${polymer_id}.json`)
#                 return false
#             } else return true
#         })

#         if ((!polys || !ligs)) {
#             console.log("Some ligands are missing. Calling script:", process.env["PYTHONBIN"], process.env["EXTRACT_BSITES_PY"])
#             exec(`${process.env["PYTHONBIN"]} ${process.env["EXTRACT_BSITES_PY"]} -s ${this.rcsb_id} --save`, (err, stdout, stderr) => {
#                 console.log(err);
#                 console.log(stdout);
#                 console.log(stderr);
#             })
#         }
#     }


#     async init_assets(obtain: boolean = false) {
#     }
# }
