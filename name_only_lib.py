from Bio import Entrez
from Bio.Entrez import efetch
from Bio.Entrez import read
import regex as re
from datetime import datetime
import time
import logging
import pandas as pd
import numpy as np

def remove_bad_format(vals, checks, table_name):
    error_messages = []
    #remove any grant with a format that would cause failure later
    for check in checks:
        if check == "delete":
            error_messages.append(str(vals[checks.index(check)])+' has been removed due to format that will cause query failure.  Check '+table_name+' to correct the value and format.')
            vals[checks.index(check)] = ""
    return [vals, error_messages]


#module to check the grant format
def check_grant_format(grant):
	if len(grant) == 11:
		if ((re.search('([A-Za-z][A-Za-z][0-9])', grant[0:3]) is not None)|(re.search('([A-Za-z][0-9][0-9])', grant[0:3]) is not None)) & (re.search('([A-Za-z][A-Za-z][0-9]+)', grant[-8:]) is not None):
			return 'pass'
		else:
			return 'error'
	elif (len(grant) >= 8):
		if re.search('([A-Za-z][A-Za-z][0-9]+$)', grant[-8:]) is not None:
			return 'pass'
	elif len(grant) == 0:
		return 'pass'
	else:
		return 'delete'
	return 'error'


def validate_config(ncbi_api, grants):
	error_messages = []
	grant_check = []

	#check that ncbi_api is a str type and has 36 characters
	if isinstance(ncbi_api, str) == False|len(ncbi_api) != 36:
		error_messages.append('NCBI API token is not correct length or format.  Query will continue without NCBI API token but check "config.py" to correct the value and format.')
		ncbi_api = ""

	if isinstance(grants, list) == False:
		error_messages.append('Grants in the "config.py" file are not in a list.  Query will continue but table showing publications linked to grants may return unexpected results or code failure.  Check "config.py" to correct the value and format.')
		grants = list(grants)

	#run each grant in the list through the module to check the grant format
	checks = list(map(check_grant_format, grants))
	
	#add error message for unexpected format found in grants
	if "error" in checks:
		error_messages.append('1 or more grants are in an unrecognized format.  Query will continue but table showing publications linked to grants may return unexpected results or code failure.  Check "config.py" to correct the value and format.')

	checked = remove_bad_format(grants, checks, 'config.py')
	
	if len(checked[1]) > 0:
		error_messages.append(checked[1])
		grants = checked[0]

	return [ncbi_api, grants, error_messages]


#module to check the Orcid format
def check_orcid_format(orcid):
    orcid = str(orcid)
    if orcid == 'nan':
        return 'pass'
    elif (len(orcid) == 19) & (re.search('([0-9]{4}-[0-9]{4}-[0-9]{4}-[0-9]{4})', orcid) is not None):
        return 'pass'
    else:
        return 'error'


def check_date_format(date_text):
    date_text = str(date_text)
    try:
        datetime.strptime(date_text, '%m/%d/%y')
        return 'pass'
    except ValueError:
        if date_text == 'nan':
            return'delete'
        else:
            raise ValueError("Incorrect data format in query_table.csv, start and end date format should be MM-DD-YY")
            return 'delete'


def validate_query_table(table):
	error_messages = []
	orcid_checks = table.apply(lambda x:check_orcid_format(x['orcid']), axis = 1)
	orcid_checked = remove_bad_format(table['orcid'], orcid_checks, 'config.py')
	
	if len(orcid_checked[1]) > 0:
		error_messages.append(orcid_checked[1])
		table['orcid'] = orcid_checked[0]

	start_date_checks = list(table.apply(lambda x:check_date_format(x['start']), axis = 1))
	start_date_checked = list(remove_bad_format(table['start'], start_date_checks, 'start dates in query_table.csv'))
	if len(start_date_checked[1]) > 0:
		error_messages.append(start_date_checked[1])
		table['start'] = start_date_checked[0]

	end_date_checks = list(table.apply(lambda x:check_date_format(x['end']), axis = 1))
	end_date_checked = list(remove_bad_format(table['end'], end_date_checks, 'end dates in query_table.csv'))
	if len(end_date_checked[1]) > 0:
		error_messages.append(end_date_checked[1])
		table['end'] = end_date_checked[0]

	if len(error_messages) == 0:
		print('Successful validation of query_table.csv with no errors.')
		time.sleep(2)
	else:
		print(error_messages)
		time.sleep(10)
	return table

def name_variations(lname, fname, mname):
	## Set Dev Values...
#selection = 4
#lname = data_frame.lname[selection]
#fname = data_frame.fname[selection]
#mname = data_frame.mname[selection]

#affiliation = data_frame.affiliation[selection]
	if mname == np.nan:
		mname = ''
	else:
		mname = str(mname)

	# break names into lists at white space or punctuation like a hyphen
	lnames = re.findall(r'\w+', lname)
	# add unbroken last name to list of last name variations
	lnames.append(lname)
	fnames = re.findall(r'\w+', fname)
	mnames = re.findall(r'\w+', mname)

	# all combinations of last name variations with first initial variations
	author_str = [x+' '+y[0] for x in lnames for y in fnames]

	# if middle name was listed, add all combinations for lastname with first and middle 
	# initial variations as well as last names with just middle initial variations
	if len(mnames) > 0:
	    author_str.extend([x+' '+y[0]+z[0] for x in lnames for y in fnames for z in mnames])
	    author_str.extend([x+' '+z[0] for x in lnames for z in mnames])

	# remove duplicates due to initial duplication e.g. Mary-Margaret
	author_str = list(set(author_str))
	return author_str


def flattenColumn(input, column):
    '''
    column is a string of the column's name.
    for each value of the column's element (which might be a list),
    duplicate the rest of columns at the corresponding row with the (each) value.
    '''
    column_flat = pd.DataFrame(
        [
            [i, c_flattened]
            for i, y in input[column].apply(list).iteritems()
            for c_flattened in y
        ],
        columns=['I', column]
    )
    column_flat = column_flat.set_index('I')
    return (
        input.drop(column, 1)
             .merge(column_flat, left_index=True, right_index=True)
    )


def name_query_term(auth_name, start, end, affiliation):
	## Set Dev Values...
#selection = 4
#start = data_frame.start[selection]
#end = data_frame.end[selection]

	start = str(start)
	end = str(end)

	# format start, end, and affiliation for pubmed query
	start = str(datetime.strptime(start, "%m/%d/%y").strftime("%Y/%m/%d"))
	if end == '': 
	    end = '3000'
	else:
	    end = str(datetime.strptime(end, "%m/%d/%y").strftime("%Y/%m/%d"))

	## create list of pubmed query terms using new author_str list
	if affiliation == '':
		term = '("'+auth_name+'"[Author]) AND ("'+start+'"[Date - Publication] : '+end+'[Date - Publication])'
	else:
		term = '("'+auth_name+'"[Author]) AND ("'+start+'"[Date - Publication] : '+end+'[Date - Publication]) AND ("'+affiliation+'"[Affiliation])'

	#query_frame = pd.DataFrame({'author': author_str, 'start': start, 'end': end, 'affiliation': affiliation, 'term':term},
	#					columns = ['author', 'index', 'start', 'end', 'affiliation', 'term'])
	return term


def orcid_query_term(orcid, start, end):
	## Set Dev Values...
#selection = 4
#start = data_frame.start[selection]
#end = data_frame.end[selection]
	start = str(start)
	end = str(end)

	# format start, end, and affiliation for pubmed query
	start = str(datetime.strptime(start, "%m/%d/%y").strftime("%Y/%m/%d"))
	if end == '': 
	    end = '3000'
	else:
	    end = str(datetime.strptime(end, "%m/%d/%y").strftime("%Y/%m/%d"))

	## create list of pubmed query terms using new author_str list
	term = '("'+orcid+'"[Identifier]) AND ("'+start+'"[Date - Publication] : '+end+'[Date - Publication])'

	#query_frame = pd.DataFrame({'author': author_str, 'start': start, 'end': end, 'affiliation': affiliation, 'term':term},
	#					columns = ['author', 'index', 'start', 'end', 'affiliation', 'term'])
	return term	


def get_pmids(term):
	pmids=''
	attempt = 0
	while attempt <= 3:
	    try:
	        handle = Entrez.esearch(db='pubmed', 
	                                #term='"'+name+'"',
	                                term=term,
	                                #field='author', #or 'orcid', #or'identifier' 
	                                retmax=5000,
	                                usehistory='y', 
	                                retmode='xml')
	        record = Entrez.read(handle)
	        handle.close()
	        if int(record['Count']) > 0:
	            pmids = record['IdList']
#	            logger.info('Entrez ESearch returns %i Ids for %s' % (int(record['Count']), str(term)))
	        else:
	            pmids = ('')
	        attempt = 4
	    except Exception as e:
#	        logger.warning('Received error from server: %s' % str(e))
#	        logger.warning('Attempt %i of 3 for %s.' % (attempt, str(term)))
	        attempt += 1
	        time.sleep(2)
#	logger.debug('Name %s queried.' % str(term))

	## Add code to write out a .csv table of terms ?even pass in author value? with resulting pmids 
	return pmids


def details(pub, variations):
    # remove all white space and \n to help regex function
    pub = ''.join(pub.split('\n'))

    pmid = re.search('<PMID.*?>(.*?)</PMID>', pub).group(1)
    # if pmc exists
    if re.search('pmc\">PMC(.*?)</ArticleId>', pub) is not None:
        pmcid = re.search('pmc\">PMC(.*?)</ArticleId>', pub).group(1)
    else:
        pmcid = ''
    # if nihms exists
    if re.search('mid\">NIHMS(.*?)</ArticleId>', pub) is not None:
        nihmsid = re.search('mid\">NIHMS(.*?)</ArticleId>', pub).group(1)
    else:
        nihmsid = ''
    # if nctid exists
    nctid = []
    if re.search('NCT(.*?)</AccessionNumber>', pub) is not None:
        nctid = re.findall('<AccessionNum.*?(NCT[0-9].*?)</AccessionNum', pub)

    nctid = ', '.join(nctid)

    if re.search('<ArticleTitle>(.*?)</ArticleTitle>', pub) is not None:
        pub_title = re.search('<ArticleTitle>(.*?)</ArticleTitle>', pub).group(1)
    else:
        pub_title = ''

    ## loop to get all author info
    # initialize lists
    authors_lnames = []
    authors_initials = []
    authors_fnames = []
    authors_affil = []

    # split into xml batches of author info
    author_list = re.split('<Author Valid', pub)

    # loop through to get author info
    for x in range(1,len(author_list)):
        if re.search('<LastName>(.*?)</LastName>', author_list[x]) is not None:
            authors_lnames.append(re.search('<LastName>(.*?)</LastName>',
                                            author_list[x]).group(1))
        else:
            authors_lnames.append('Unknown')
        if re.search('<Initials>(.*?)</Initials>', author_list[x]) is not None:
            authors_initials.append(re.search('<Initials>(.*?)</Initials>',
                                              author_list[x]).group(1))
        else:
            authors_initials.append('Unknown')
        if re.search('<ForeName>(.*?)</ForeName>', author_list[x]) is not None:
            authors_fnames.append(re.search('<ForeName>(.*?)</ForeName>',
                                            author_list[x]).group(1))
        else:
            authors_fnames.append('Unknown')
        if re.search('Affiliation>', author_list[x]) is not None:
            authors_affil.append(re.search('Affiliation>(.*?)</Affiliation',
                                           author_list[x]).group(1))
        else:
            authors_affil.append('')
    # combine fname and lname to get full list of author names
    authors = [i+' '+j for i, j in zip(authors_fnames, authors_lnames)]

    authors_lnames = ', '.join(authors_lnames)
    authors_fnames = ', '.join(authors_fnames)
    authors_initials = ', '.join(authors_initials)
    authors_affil = ', '.join(authors_affil)
    authors = ', '.join(authors)


    ## get pub_date from when journal was published
    if re.search('<JournalIssue.*?<PubDate>(.*?)</PubDate>', pub) is not None:
        publish_date = re.search('<JournalIssue.*?<PubDate>(.*?)</PubDate>',
                                 pub).group(1)
    else:
        publish_date = ''

    # clean up pub_date
    if re.search('<Medline', publish_date) is None:
        if re.search('<Year>[0-9]{4}', publish_date) is None:
            year = '2099'
        else:
            year = re.search('<Year>([0-9]{4})</Year>', publish_date).group(1)
        if re.search('<Month>', publish_date) is None:
            month = '01'
        elif re.search('<Month>([A-Za-z].*?)-.*</Month>',
                       publish_date) is not None:
            month = re.search('<Month>([A-Za-z].*?)-.*</Month>',
                              publish_date).group(1)
        else:
            month = re.search('<Month>(.*)</Month>', publish_date).group(1)
        if re.search('<Day>', publish_date) is None:
            day = '01'
        else:
            day = re.search('<Day>(.*)</Day>', publish_date).group(1)
    else:
        medline = re.search('<Medline.*?>(.*?)</Medline.*?>',
                            publish_date).group(1)
        if re.search('[A-Za-z]', medline) is not None:
            year = re.search('^.*?([0-9]{4}).*?$', medline).group(1)
            
            if re.search('Summer', medline) is not None:
                month = '06'
            else:
                try:
                    month = re.search('^.*?([A-Za-z]{3}).*?[-|/].*$', medline).group(1)
                except Exception as err:
                    print(medline)
                    month = '01'
            day = '01'
        else:
            year = re.search('^([0-9]{4}).*?$', medline).group(1)
            month = '01'
            day = '01'

    # remove all whitespace from 'month'
    month = ''.join(month.split())

    # combine month day and year
    pub_date = year+ '-' + month + '-' + day
    # check if month is letters and change into numerical date type
    if re.search('[a-zA-Z]', pub_date) is not None:
        pub_date = datetime.strptime(pub_date,
                                              "%Y-%b-%d").strftime("%Y-%m-%d")

    if re.search('</JournalIssue>.*?<ISOAbbreviation>(.*?)</ISOAbbre',
                              pub) is not None:
        journal_short = re.search('</JournalIssue>.*?<ISOAbbreviation>(.*?)</ISOAbbre',
                                  pub).group(1)
    else:
        journal_short = 'Unknown'

    if re.search('/JournalIssue>.*?<Title>(.*?)</Title>',
                             pub) is not None:
        journal_full = re.search('/JournalIssue>.*?<Title>(.*?)</Title>',
                                 pub).group(1)
    else:
        journal_full = 'Unknown'

    ## get grant list to clean up and compare with variations to get pubmed tags
    pubmed_tags = []
    grant_list = re.split('<Grant>', pub)
    for x in range(1, len(grant_list)):
        if re.search('<GrantID>', grant_list[x]) is not None:
            if re.search('<GrantID>(.*?)</GrantID>',
                         grant_list[x]).group(1) in variations:
                pubmed_tags.append(re.search('<GrantID>(.*?)</GrantID>',
                                             grant_list[x]).group(1))

    pubmed_tags = ', '.join(pubmed_tags)

    row = [pmid, pmcid, nihmsid,  nctid, pub_title, authors,
            authors_lnames, authors_initials, authors_affil,
            pub_date, journal_short, journal_full, pubmed_tags]

    return row


def summary(pmids, ncbi_key, grants):
	#***!!! developing !!!***
	Entrez.email = "Your.Name.Here@example.org"
	Entrez.api_key = ncbi_key
	logger = logging.getLogger(__name__)
	#***!!! dev variables... ***!!!

	#pmids = '22339280,24860250,27801552,28570838,23544412,27239297,28760970,30368762,21451738'
	#grants = variations
	#Entrez.api_key = keys.ncbi_key
	#logger = logging.getLogger(__name__)
	#***!!!

	try:
			from urllib.error import HTTPError # for Python 3
	except ImportError:
			from urllib2 import HTTPError # for Python 2

	count = len(pmids)
	records = []
	attempt = 0

	while attempt < 3:
		attempt += 1
		logger.info('Going to Epost pmid list results')
		try:
			# query pubmed with pmids and post results with ePost
			post_xml = Entrez.epost('pubmed', id=','.join(pmids))
			# read results
			search_results = Entrez.read(post_xml)
			# close the link
			post_xml.close()
			attempt = 4
		except HTTPError as err:
			if 500 <= err.code <= 599:
				logger.warning('Received error from server: %s' % err)
				logger.warning('Attempt %i of 3' % attempt)
				time.sleep(10)
			else:
				raise

	# set paramater values from ePost location to get xml with eFetch
	webenv = search_results['WebEnv']
	query_key = search_results['QueryKey']

	batch_size = 500

	for start in range(0, count, batch_size):
		end = min(count, start+batch_size)
		logger.info('Going to fetch record %i to %i' % (start+1, end))
		attempt = 0
		while attempt < 3:
			attempt += 1
			try:
				# use eFetch to get xml information out of ePost results
				fetch_handle = Entrez.efetch(db='pubmed',
				                             retstart=start, retmax=batch_size,
				                             webenv=webenv, query_key=query_key,
				                             retmode='xml')
				records.extend(fetch_handle.read())
				fetch_handle.close
				attempt = 4
			except HTTPError as err:
				if 500 <= err.code <= 599:
					logger.warning('Received error from server: %s' % err)
					logger.warning('Attempt %i of 3' % attempt)
					time.sleep(10)
				else:
					raise

	pub_list = re.split('<PubmedArticle>', ''.join(records))
	rows = []
	for x in range(1, len(pub_list)):
		# assemble list of publication details
		rows.append(details(pub_list[x], grants))

	pubs_frame = pd.DataFrame(rows, columns=[
	                          'pmid', 'pmcid', 'nihmsid',  'nctid', 'pub_title',
	                          'authors', 'authors_lnames', 'authors_initials',
	                          'authors_affil', 'pub_date', 'journal_short',
	                          'journal_full', 'pubmed_tags'])
	return pubs_frame
