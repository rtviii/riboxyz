import pprint
import urllib.parse
import argparse
import csv
import requests
import json
import pandas as pd
from os import path
from datetime import datetime

ROR_API_ENDPOINT  = "https://api.dev.ror.org/organizations"
INPUT_DIR         = ""
OUTPUT_DIR        = ""
ES_RESERVED_CHARS = [ "+", "-", "&", "|", "!", "(", ")", "{", "}", "[", "]", "^", '"', "~", "*", "?", ":", "\\", "/" ]

def process_file(input_file):
    now         = datetime.now()
    fields      = ['search_term',
                    'query_results',
                    'has_chosen_match',
                    'match_type',
                    'substring',
                    'chosen_match_names',
                    'has_good_match',
                    'matched_field',
                    'first_match_names',
                    'first_match_aliases',
                    'first_match_acronyms',
                    'first_match_labels']

    output_file = now.strftime("%Y-%m-%d") + "_search_results_affiliation_matching.csv"
    NEW_DF =pd.DataFrame(columns=fields)

    df = pd.read_csv(input_file, delimiter=',', index_col=False, header=0)
    for row in df.iterrows():
        query_term =  row[1][1]
        print("Searching for: " + query_term)
        params    = {'affiliation': query_term}
        resp_     = requests.get(ROR_API_ENDPOINT + '?' + urllib.parse.urlencode(params))
        resp_json = resp_.json()

        query_results = resp_json['number_of_results']
        NEW_ROW       = []

        first_match_names    = []
        first_match_aliases  = []
        first_match_acronyms = []
        first_match_labels   = []
        chosen_match_names   = []
        has_good_match       = False
        has_chosen_match     = False
        match_type           = ''
        substring            = ''
        matched_field        = ''

        if query_results == 0:
           first_match_names    = []
           first_match_aliases  = []
           first_match_acronyms = []
           first_match_labels   = []
        else:
            if len(resp_json['items']) <= 5:
                for item in resp_json['items']:
                    if(item['chosen'] == True):

                        has_chosen_match = True
                        match_type       = item['matching_type']
                        substring        = item['substring']
                        chosen_match_names.append(item['organization']['name'])
                        NEW_ROW.append(item['organization']['name'])

                    first_match_names.append(item['organization']['name'])
                    NEW_ROW.append(item['organization']['name'])
                    if len(item['organization']['aliases']) > 0:
                        for alias in item['organization']['aliases']:
                            first_match_aliases.append(alias)
                            NEW_ROW.append(alias)

                    if len(item['organization']['acronyms']) > 0:
                        for acronym in item['organization']['acronyms']:
                            first_match_acronyms.append(acronym)
                            NEW_ROW.append(acronym)

                    if len(item['organization']['labels']) > 0:
                        for label in item['organization']['labels']:
                            first_match_labels.append(label['label'])
                            NEW_ROW.append(label)
            else:
                for i in range(5):
                    if(resp_json['items'][i]['chosen'] == True):
                        has_chosen_match = True
                        match_type = resp_json['items'][i]['matching_type']
                        substring = resp_json['items'][i]['substring']
                        chosen_match_names.append(resp_json['items'][i]['organization']['name'])
                        NEW_ROW.append(resp_json['items'][i]['organization']['name'])

                    first_match_names.append(resp_json['items'][i]['organization']['name'])
                    NEW_ROW.append(resp_json['items'][i]['organization']['name'])

                    if len(resp_json['items'][i]['organization']['aliases']) > 0:
                        for alias in resp_json['items'][i]['organization']['aliases']:
                            first_match_aliases.append(alias)
                            NEW_ROW.append(alias)
                    if len(resp_json['items'][i]['organization']['acronyms']) > 0:
                        for acronym in resp_json['items'][i]['organization']['acronyms']:
                            first_match_acronyms.append(acronym)
                            NEW_ROW.append(acronym)
                    if len(resp_json['items'][i]['organization']['labels']) > 0:
                        for label in resp_json['items'][i]['organization']['labels']:
                            first_match_labels.append(label['label'])
                            NEW_ROW.append(label['label'])

        if first_match_names:
            if query_term in first_match_names:
                has_good_match = True
                matched_field  = 'name'
            elif query_term in first_match_aliases:
                has_good_match = True
                matched_field  = 'aliases'
            elif query_term in first_match_acronyms:
                has_good_match = True
                matched_field  = 'acronyms'
            elif query_term in first_match_labels:
                has_good_match = True
                matched_field  = 'labels'
            else:
                has_good_match = False


        df_row = pd.DataFrame([NEW_ROW], columns=fields)
        print(df_row)

        # NEW_DF = NEW_DF.append()
        # pd_row.to_csv(output_file)


input_file = "org_names.csv"
if path.exists(input_file):
    process_file(input_file)
else:
    print("File " + input_file + " does not exist. Cannot process file.")