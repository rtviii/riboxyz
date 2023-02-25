



import gzip
import os
import requests


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
 

# export const commit_struct_to_db_sync = (rcsb_id: string): Promise<void> => {
#     console.log(`Commiting ${rcsb_id} to the database`)
#     const commit_script = process.env["COMMIT_STRUCTURE_SH"]
#     let   current_db    = process.env["NEO4J_CURRENTDB"]
#     let   uri           = process.env["NEO4J_URI"]
#     let   invocation    = `${commit_script} -s ${rcsb_id} -d ${current_db} -a "${uri}"`
#     return new Promise<void>((resolve, reject) => {
#         let proc = cp.exec(invocation)
#         if (proc.stderr !== null) {
#             proc.stderr.on("data",
#                 (data) => {
#                     console.log(data)
#                 }
#             )
#         }

#         proc.stdout?.on("data", (data) => { console.log(data) })
#         proc.stderr?.on("data", (data) => { console.log(data); reject(data) })
#         proc.on("exit", () => { process.exit(); resolve() })

#     })

# }