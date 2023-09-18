from ribctl import RIBETL_DATA
from ribctl.etl.ribosome_assets import RibosomeAssets
import os
import urllib.request
from concurrent.futures import ThreadPoolExecutor


# Function to download a file from a URL and save it to the specified directory
def download_file(url, dest):
    try:
        urllib.request.urlretrieve(url, dest)
        print(f"File downloaded and saved to {dest}")
    except Exception as e:
        print(f"An error occurred while downloading {url}: {str(e)}")

# Create a thread pool with 10 workers
with ThreadPoolExecutor(max_workers=50) as executor:
    # Submit download tasks for each URL

    for struct in os.listdir(RIBETL_DATA):
        prof = RibosomeAssets(struct).profile()

        try:
            emd_id      = list(filter(lambda x: "EMD-" in x, prof.rcsb_external_ref_id))[0]
            url         = "https://ftp.ebi.ac.uk/pub/databases/emdb/structures/{}/map/{}.map.gz".format(emd_id,emd_id.replace('-','_').lower() )
            destination = os.path.join(RIBETL_DATA, struct.upper(), "{}_{}.map.gz".format(emd_id.replace('-','_').lower() ,struct.lower()))

            if not os.path.exists(destination):
                executor.submit(download_file, url, destination)
            else:
                print("Skipping. EXISTS: {} ".format(destination))
        except:
            print("NO EMD LINK FOR {}".format(struct))

# Replace "https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0393/map/emd_0393.map.gz" and add more URLs to the urls list as needed. The script will create a thread pool with 10 workers to download the files concurrently and save them to the specified directory.
