import requests

# Define the UniProt entry name
entry_name = "A0A176W8K1"

# Define the API URL
api_url = f"https://www.uniprot.org/uniprot/{entry_name}.tab"

# Define the query parameters to include only "90% unirefs" sequences and select specific fields
params = {
    "query": f"reviewed:yes AND identity:0.9 AND id:{entry_name}",
    "format": "tab",
    "columns": "organism-id,sequence"
}

# Make the GET request
response = requests.get(api_url, params=params)
# Check if the request was successful
if response.status_code == 200:
    # Print the TSV data
    print(response.text)
else:
    print(f"Error: Unable to retrieve data. Status code: {response.status_code}")