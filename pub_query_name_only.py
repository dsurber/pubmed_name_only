from Bio import Entrez
from Bio.Entrez import efetch
from Bio.Entrez import read
import regex as re
import pandas as pd
import numpy as np
from datetime import datetime
from redcap import Project
import logging
import time
import chardet
import itertools

import config
import name_only_lib

logger = logging.getLogger(__name__)

#validate data and formats in config.py
val = name_only_lib.validate_config(config.ncbi_api, config.grants)
config.ncbi_api = val[0]
config.grants = val[1]

if len(val[2]) == 0:
    print('Successful validation of config.py with no errors.')
    time.sleep(2)
else:
    print(val[2])
    time.sleep(10)

### Create variables querying pubmed api
Entrez.email = "Your.Name.Here@example.org"
Entrez.api_key = config.ncbi_api

# query pubmed for pmids associated with each grant variation
logger.info("Starting pubmed queries...")
# create set for unique list of all pmids from querying pubmed with name variations
pmids = set()

# read query_table.csv
in_file = "query_table.csv"
data = pd.read_csv(in_file, encoding=chardet.detect(open(in_file, 'rb').read())['encoding'])
# validate data formats in the query table
data = name_only_lib.validate_query_table(data)

#data[pd.isnull(data)] = ''
data.fillna('', inplace=True)

# Create tables that will later be written to csv
# create table for researchers with orcid ids and query pubmed for pmids in start and end date window
orcid_table = data[data.orcid != ''].reset_index()

if len(orcid_table) > 0:
    orcid_table['term'] = orcid_table.apply(lambda x: name_only_lib.orcid_query_term(x['orcid'], x['start'], x['end']), axis = 1)
    orcid_table['pmids'] = orcid_table.apply(lambda x: name_only_lib.get_pmids(x['term']), axis = 1)
else:
    orcid_table = 'none'


# Create tables that will later be written to csv
# create table for researchers with orcid ids and query pubmed for pmids in start and end date window
names_table = pd.DataFrame(columns = ['name_variation'])
names_table.name_variation = data.apply(lambda x: name_only_lib.name_variations(x['lname'], x['fname'], x['mname']), axis = 1)
names_table = name_only_lib.flattenColumn(names_table, 'name_variation')
names_table = data.merge(names_table, left_index = True, right_index = True)

# create column of query terms
names_table['term'] = names_table.apply(lambda x: name_only_lib.name_query_term(x['name_variation'], x['start'], x['end'], x['affiliation']), axis = 1)
# query pubmed for pmids resulting from each term
names_table['pmids'] = names_table.apply(lambda x: name_only_lib.get_pmids(x['term']), axis = 1)

if orcid_table == 'none':
    pmids = list(set(itertools.chain(*names_table.pmids)))
else:
    pmids = list(set(itertools.chain(*names_table.pmids, *orcid_table.pmids)))

### Get table of publication details from pubmed for pmids
# make dataframe of publications
pubs_frame = name_only_lib.summary(pmids, Entrez.api_key, config.grants)

## add a column of name variations for each pmid to the end of pubs_frame
all_variations = []
for pmid in pubs_frame.pmid:
    all_variations.append(names_table.name_variation[names_table['pmids'].astype(str).str.contains(pmid)].tolist())
pubs_frame['name_variations'] = all_variations

## clean up and output the csv tables
pubs_frame = pubs_frame.replace(',', ';', regex=True)
pubs_frame = pubs_frame.apply(lambda x: x.str.slice(0, 30000))
pubs_frame.to_csv('./Reports/pmid_details_table.csv', index=False)
names_table.to_csv('./Reports/names_results_table.csv', index=False)

if orcid_table != 'none':
    orcid_table.to_csv('./Reports/orcid_results_table.csv', index=False)

print('\nQuery complete and reports have been generated in "Reports" folder.')
