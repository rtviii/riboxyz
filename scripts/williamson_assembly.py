
from concurrent.futures import ThreadPoolExecutor
import os
import urllib.request

def download_file(url, dest):
    try:
        urllib.request.urlretrieve(url, dest)
        print(f"File downloaded and saved to {dest}")
    except Exception as e:
        print(f"An error occurred while downloading {url}: {str(e)}")

assembly_stages = {
"dead-B-a1" : "https://www.ebi.ac.uk/emdb/EMD-40517",
"dead-B-a2" : "https://www.ebi.ac.uk/emdb/EMD-40519",
"dead-B-a3" : "https://www.ebi.ac.uk/emdb/EMD-40520",
"dead-B-a4" : "https://www.ebi.ac.uk/emdb/EMD-40524",
"dead-B-a5" : "https://www.ebi.ac.uk/emdb/EMD-40526",
"dead-B-a6" : "https://www.ebi.ac.uk/emdb/EMD-40528",
"dead-B-b1" : "https://www.ebi.ac.uk/emdb/EMD-40530",
"dead-C-a1" : "https://www.ebi.ac.uk/emdb/EMD-40532",
"dead-C-a2" : "https://www.ebi.ac.uk/emdb/EMD-40534",
"dead-C-a3" : "https://www.ebi.ac.uk/emdb/EMD-40536",
"dead-C-a4" : "https://www.ebi.ac.uk/emdb/EMD-40538",
"dead-C-a5" : "https://www.ebi.ac.uk/emdb/EMD-40540",
"dead-C-b1" : "https://www.ebi.ac.uk/emdb/EMD-40542",
"dead-C-b2" : "https://www.ebi.ac.uk/emdb/EMD-40544",
"dead-E-a1" : "https://www.ebi.ac.uk/emdb/EMD-40546",
"dead-G1"   : "https://www.ebi.ac.uk/emdb/EMD-40548",
"dead-G2"   : "https://www.ebi.ac.uk/emdb/EMD-40550",
"dead-G3"   : "https://www.ebi.ac.uk/emdb/EMD-40552",
"dead-G4"   : "https://www.ebi.ac.uk/emdb/EMD-40555",
"dead-preB1": "https://www.ebi.ac.uk/emdb/EMD-40551",
"dead-preB2": "https://www.ebi.ac.uk/emdb/EMD-40549",
"rl17-B-b1" : "https://www.ebi.ac.uk/emdb/EMD-40309",
"rl17-B-b2" : "https://www.ebi.ac.uk/emdb/EMD-40311",
"rl17-B-b3" : "https://www.ebi.ac.uk/emdb/EMD-40313",
"rl17-B-b4" : "https://www.ebi.ac.uk/emdb/EMD-40511",
"rl17-C-b1" : "https://www.ebi.ac.uk/emdb/EMD-40314",
"rl17-C-b2" : "https://www.ebi.ac.uk/emdb/EMD-40315",
"rl17-C-b3" : "https://www.ebi.ac.uk/emdb/EMD-40317",
"rl17-C-b4" : "https://www.ebi.ac.uk/emdb/EMD-40319",
"rl17-C-b5" : "https://www.ebi.ac.uk/emdb/EMD-40321",
"rl17-C-b6" : "https://www.ebi.ac.uk/emdb/EMD-40323",
"rl17-D-a1" : "https://www.ebi.ac.uk/emdb/EMD-40327",
"rl17-D-a2" : "https://www.ebi.ac.uk/emdb/EMD-40329",
"rl17-D-a3" : "https://www.ebi.ac.uk/emdb/EMD-40331",
"rl17-D-b1" : "https://www.ebi.ac.uk/emdb/EMD-40333",
"rl17-D-b2" : "https://www.ebi.ac.uk/emdb/EMD-40512",
"rl17-D-b3" : "https://www.ebi.ac.uk/emdb/EMD-40514",
"rl17-D-b4" : "https://www.ebi.ac.uk/emdb/EMD-40516",
"rl17-E-a1" : "https://www.ebi.ac.uk/emdb/EMD-40518",
"rl17-E-a2" : "https://www.ebi.ac.uk/emdb/EMD-40521",
"rl17-E-a3" : "https://www.ebi.ac.uk/emdb/EMD-40523",
"rl17-E-a4" : "https://www.ebi.ac.uk/emdb/EMD-40525",
"rl17-E-a5" : "https://www.ebi.ac.uk/emdb/EMD-40527",
"rl17-E-b1" : "https://www.ebi.ac.uk/emdb/EMD-40529",
"rl17-E-b10": "https://www.ebi.ac.uk/emdb/EMD-40531",
"rl17-E-b2" : "https://www.ebi.ac.uk/emdb/EMD-40533",
"rl17-E-b3" : "https://www.ebi.ac.uk/emdb/EMD-40535",
"rl17-E-b4" : "https://www.ebi.ac.uk/emdb/EMD-40537",
"rl17-E-b5" : "https://www.ebi.ac.uk/emdb/EMD-40539",
"rl17-E-b6" : "https://www.ebi.ac.uk/emdb/EMD-40541",
"rl17-E-b7" : "https://www.ebi.ac.uk/emdb/EMD-40543",
"rl17-E-b8" : "https://www.ebi.ac.uk/emdb/EMD-40545",
"rl17-E-b9" : "https://www.ebi.ac.uk/emdb/EMD-40547",
"srmb-B-a1" : "https://www.ebi.ac.uk/emdb/EMD-40316",
"srmb-B-a2" : "https://www.ebi.ac.uk/emdb/EMD-40318",
"srmb-C-a1" : "https://www.ebi.ac.uk/emdb/EMD-40320",
"srmb-C-a2" : "https://www.ebi.ac.uk/emdb/EMD-40322",
"srmb-C-a3" : "https://www.ebi.ac.uk/emdb/EMD-40324",
"srmb-C-a4" : "https://www.ebi.ac.uk/emdb/EMD-40325",
"srmb-C-a5" : "https://www.ebi.ac.uk/emdb/EMD-40328",
"srmb-E-a1" : "https://www.ebi.ac.uk/emdb/EMD-40330",
"srmb-E-a2" : "https://www.ebi.ac.uk/emdb/EMD-40332",
"srmb-E-a3" : "https://www.ebi.ac.uk/emdb/EMD-40513",
"srmb-G1"   : "https://www.ebi.ac.uk/emdb/EMD-4051"}



assembly_intermediates = [
"https://www.ebi.ac.uk/emdb/EMD-40517",
"https://www.ebi.ac.uk/emdb/EMD-40519",
"https://www.ebi.ac.uk/emdb/EMD-40520",
"https://www.ebi.ac.uk/emdb/EMD-40524",
"https://www.ebi.ac.uk/emdb/EMD-40526",
"https://www.ebi.ac.uk/emdb/EMD-40528",
"https://www.ebi.ac.uk/emdb/EMD-40530",
"https://www.ebi.ac.uk/emdb/EMD-40532",
"https://www.ebi.ac.uk/emdb/EMD-40534",
"https://www.ebi.ac.uk/emdb/EMD-40536",
"https://www.ebi.ac.uk/emdb/EMD-40538",
"https://www.ebi.ac.uk/emdb/EMD-40540",
"https://www.ebi.ac.uk/emdb/EMD-40542",
"https://www.ebi.ac.uk/emdb/EMD-40544",
"https://www.ebi.ac.uk/emdb/EMD-40546",
"https://www.ebi.ac.uk/emdb/EMD-40548",
"https://www.ebi.ac.uk/emdb/EMD-40550",
"https://www.ebi.ac.uk/emdb/EMD-40552",
"https://www.ebi.ac.uk/emdb/EMD-40555",
"https://www.ebi.ac.uk/emdb/EMD-40551",
"https://www.ebi.ac.uk/emdb/EMD-40549",
"https://www.ebi.ac.uk/emdb/EMD-40309",
"https://www.ebi.ac.uk/emdb/EMD-40311",
"https://www.ebi.ac.uk/emdb/EMD-40313",
"https://www.ebi.ac.uk/emdb/EMD-40511",
"https://www.ebi.ac.uk/emdb/EMD-40314",
"https://www.ebi.ac.uk/emdb/EMD-40315",
"https://www.ebi.ac.uk/emdb/EMD-40317",
"https://www.ebi.ac.uk/emdb/EMD-40319",
"https://www.ebi.ac.uk/emdb/EMD-40321",
"https://www.ebi.ac.uk/emdb/EMD-40323",
"https://www.ebi.ac.uk/emdb/EMD-40327",
"https://www.ebi.ac.uk/emdb/EMD-40329",
"https://www.ebi.ac.uk/emdb/EMD-40331",
"https://www.ebi.ac.uk/emdb/EMD-40333",
"https://www.ebi.ac.uk/emdb/EMD-40512",
"https://www.ebi.ac.uk/emdb/EMD-40514",
"https://www.ebi.ac.uk/emdb/EMD-40516",
"https://www.ebi.ac.uk/emdb/EMD-40518",
"https://www.ebi.ac.uk/emdb/EMD-40521",
"https://www.ebi.ac.uk/emdb/EMD-40523",
"https://www.ebi.ac.uk/emdb/EMD-40525",
"https://www.ebi.ac.uk/emdb/EMD-40527",
"https://www.ebi.ac.uk/emdb/EMD-40529",
"https://www.ebi.ac.uk/emdb/EMD-40531",
"https://www.ebi.ac.uk/emdb/EMD-40533",
"https://www.ebi.ac.uk/emdb/EMD-40535",
"https://www.ebi.ac.uk/emdb/EMD-40537",
"https://www.ebi.ac.uk/emdb/EMD-40539",
"https://www.ebi.ac.uk/emdb/EMD-40541",
"https://www.ebi.ac.uk/emdb/EMD-40543",
"https://www.ebi.ac.uk/emdb/EMD-40545",
"https://www.ebi.ac.uk/emdb/EMD-40547",
"https://www.ebi.ac.uk/emdb/EMD-40316",
"https://www.ebi.ac.uk/emdb/EMD-40318",
"https://www.ebi.ac.uk/emdb/EMD-40320",
"https://www.ebi.ac.uk/emdb/EMD-40322",
"https://www.ebi.ac.uk/emdb/EMD-40324",
"https://www.ebi.ac.uk/emdb/EMD-40325",
"https://www.ebi.ac.uk/emdb/EMD-40328",
"https://www.ebi.ac.uk/emdb/EMD-40330",
"https://www.ebi.ac.uk/emdb/EMD-40332",
"https://www.ebi.ac.uk/emdb/EMD-40513",
"https://www.ebi.ac.uk/emdb/EMD-40515"]


def extract_emdb_id(emdb_url:str):
    # "https://www.ebi.ac.uk/emdb/EMD-40517",
    emdb_id  = emdb_url.split('/')[-1]
    filename = emdb_id.replace("-","_").lower() + ".map.gz"
    return emdb_id, filename




with ThreadPoolExecutor(max_workers=50) as executor:
    # Submit download tasks for each URL
    for url in assembly_intermediates:
        emdb_id, filename = extract_emdb_id(url)
        try:
            url         = "https://ftp.ebi.ac.uk/pub/databases/emdb/structures/{}/map/{}".format(emdb_id, filename)
            destination = os.path.join('/home/rtviii/dev/riboxyz/scripts/assembly_intermediates', filename)
            if not os.path.exists(destination):
                executor.submit(download_file, url, destination)
            else:
                print("Skipping. EXISTS: {} ".format(destination))
        except:
            print("NO EMD LINK FOR {}".format(emdb_id))