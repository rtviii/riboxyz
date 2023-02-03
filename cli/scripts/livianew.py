import urllib.parse
import argparse
import csv
import requests
import json
import pandas as pd
from os import path
from datetime import datetime

ROR_API_ENDPOINT = "https://api.dev.ror.org/organizations"
INPUT_FILE       = "org_names.csv"
INPUT_ORG_NAMES  = pd.read_csv(INPUT_FILE, delimiter=',')

now              = datetime.now()
output_file      = now.strftime("%Y-%m-%d") + "_search_results_affiliation_matching.csv"
DF_OUT           = pd.DataFrame()
fields           = ['search_term', 'query_results', 'has_chosen_match', 'match_type', 'substring', 'chosen_match_names', 'has_good_match', 'matched_field', 'first_match_names', 'first_match_aliases', 'first_match_acronyms', 'first_match_labels']

with open(output_file, 'w') as csvfile:
    for row in INPUT_ORG_NAMES.itertuples():
        search_term= row.x

        print("Searching for: " + search_term)

        params         = {'affiliation': search_term}
        query_response = requests.get(ROR_API_ENDPOINT + '?' + urllib.parse.urlencode(params))
        query_response = query_response.json()
        query_results  = query_response['number_of_results']

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
            if len(query_response['items']) <= 5:
                for item in query_response['items']:
                    if(item['chosen'] == True):
                        has_chosen_match = True
                        match_type = item['matching_type']
                        substring = item['substring']
                        chosen_match_names.append(item['organization']['name'])
                    first_match_names.append(item['organization']['name'])
                    if len(item['organization']['aliases']) > 0:
                        for alias in item['organization']['aliases']:
                            first_match_aliases.append(alias)
                    if len(item['organization']['acronyms']) > 0:
                        for acronym in item['organization']['acronyms']:
                            first_match_acronyms.append(acronym)
                    if len(item['organization']['labels']) > 0:
                        for label in item['organization']['labels']:
                            first_match_labels.append(label['label'])
            else:
                for i in range(5):
                    if(query_response['items'][i]['chosen'] == True):
                        has_chosen_match = True
                        match_type = query_response['items'][i]['matching_type']
                        substring = query_response['items'][i]['substring']
                        chosen_match_names.append(query_response['items'][i]['organization']['name'])
                    first_match_names.append(query_response['items'][i]['organization']['name'])
                    if len(query_response['items'][i]['organization']['aliases']) > 0:
                        for alias in query_response['items'][i]['organization']['aliases']:
                            first_match_aliases.append(alias)
                    if len(query_response['items'][i]['organization']['acronyms']) > 0:
                        for acronym in query_response['items'][i]['organization']['acronyms']:
                            first_match_acronyms.append(acronym)
                    if len(query_response['items'][i]['organization']['labels']) > 0:
                        for label in query_response['items'][i]['organization']['labels']:
                            first_match_labels.append(label['label'])

        if first_match_names:
            if search_term in first_match_names:
                has_good_match = True
                matched_field = 'name'
            elif search_term in first_match_aliases:
                has_good_match = True
                matched_field = 'aliases'
            elif search_term in first_match_acronyms:
                has_good_match = True
                matched_field = 'acronyms'
            elif search_term in first_match_labels:
                has_good_match = True
                matched_field = 'labels'
            else:
                has_good_match = False

        pd_row= [search_term, query_results, has_chosen_match, match_type, substring, chosen_match_names, has_good_match, matched_field, first_match_names, first_match_aliases, first_match_acronyms, first_match_labels]
        DF_OUT = pd.concat([DF_OUT, pd.DataFrame([pd_row], columns=fields)], ignore_index=True)
        DF_OUT.to_csv(output_file, index=False)

