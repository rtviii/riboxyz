from ribctl.ribosome_assets import RibosomeAssets


for i in ['3j7z', '5afi','7k00']:
    RCSB_ID = i

    prof = RibosomeAssets(RCSB_ID).profile()
    print(prof.rcsb_external_ref_link)
    print(prof.rcsb_external_ref_id)
    print(prof.rcsb_external_ref_type)



# import os
# import urllib.request
# from concurrent.futures import ThreadPoolExecutor

# # Define the list of URLs
# urls = [
#     "https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0393/map/emd_0393.map.gz",
#     # Add more URLs here...
# ]

# # Define the directory where you want to save the files
# DIR_X = "/path/to/DIR_X"

# # Create the directory if it doesn't exist
# if not os.path.exists(DIR_X):
#     os.makedirs(DIR_X)

# # Function to download a file from a URL and save it to the specified directory
# def download_file(url):
#     try:
#         # Extract the filename from the URL
#         file_name = os.path.join(DIR_X, os.path.basename(url))
        
#         # Download the file and save it to the specified directory
#         urllib.request.urlretrieve(url, file_name)
#         print(f"File downloaded and saved to {file_name}")
#     except Exception as e:
#         print(f"An error occurred while downloading {url}: {str(e)}")

# # Create a thread pool with 10 workers
# with ThreadPoolExecutor(max_workers=10) as executor:
#     # Submit download tasks for each URL
#     for url in urls:
#         executor.submit(download_file, url)

# Replace "https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0393/map/emd_0393.map.gz" and add more URLs to the urls list as needed. The script will create a thread pool with 10 workers to download the files concurrently and save them to the specified directory.
