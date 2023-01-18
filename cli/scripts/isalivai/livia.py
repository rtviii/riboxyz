import csv
import requests
from datetime import datetime

ROR_API_ENDPOINT = "https://api.ror.org/organizations"
INPUT_DIR        = ""
OUTPUT_DIR       = ""

input_file  = "affiliation_names.csv"
now         = datetime.now()
output_file = OUTPUT_DIR + now.strftime("%Y-%m-%d") + "_matched_ror_ids.csv"
fields      = ['input_id', 'ror_id']

with open(input_file) as csv_file:
    reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    with open(output_file, 'w') as outpath:
        writer = csv.writer(outpath, delimiter=',')
        writer.writerow(fields)
        for row in reader:
            input_id = row[0]
            print("Finding ROR for: " + input_id)
            search_term = '"' + input_id + '"'
            params      = {'query': search_term}
            ror_id      = ''
            try:
                _resp = requests.get(ROR_API_ENDPOINT,params=[ ('query', search_term) ] )
                response = _resp.json()
                if response['number_of_results'] == 0:
                    ror_id = ''
                elif response['number_of_results'] == 1:
                    ror_id = response['items'][0]['id']
                else:
                    ror_id = ''
                    for items in response:
                        ror_id = ror_id + ", " + response['items'][0]['id']
            except ValueError:
                ror_id = 'Error'
                pass
            finally:
                writer.writerow([input_id, ror_id])