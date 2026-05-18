📚 PubMed Researcher Query Tool
A Python-based tool for querying PubMed to retrieve publications associated with researchers using name variations and/or ORCID IDs, with optional grant filtering.
This project automates:
	Generating name variations for authors
	Querying PubMed via the NCBI Entrez API
	Collecting publication metadata
	Producing clean CSV reports for analysis


🚀 Features
	🔎 Query PubMed using:
		Author name variations
		ORCID identifiers
	📅 Filter publications by date range
	🏫 Optional affiliation filtering
	🎯 Grant matching support (NIH-style grant numbers)
	📊 Output detailed publication reports as CSV files
	✅ Built-in validation for input data and configuration


📦 Installation
1. Clone repository
	git clone https://github.com/your-username/pubmed-query-tool.gitcd pubmed-query-toolShow more lines
2. Install dependencies
	pip install biopython pandas numpy regex python-redcap chardetShow more lines

⚙️ Configuration
Edit the config.py file:
	ncbi_api = "YOUR_NCBI_API_KEY"
	grants = [""] ##comma delimited if listing

ncbi_api: Your NCBI API key (optional but recommended for higher rate limits)
grants: List of grant numbers to match against publications


📁 Input Data
Provide a CSV file named:
query_table.csv

Required columns:
column 		description
lname 		Last name
fname 		First name
mname 		Middle name (optional)
orcid 		ORCID ID (optional)
start 		Start date (MM/DD/YYYY)
end 		End date (MM/DD/YYYY or blank)
affiliation Institution (optional)

Example:
lname,fname,mname,orcid,start,end,affiliation
brasier,allan,R,,1/1/2023,6/1/2025,University of Wisconsin

▶️ Usage
Run the main script:
	python pub_query_name_only.py

📊 Output
Results are written to a Reports/ directory:

File 					Description
pmid_details_table.csv	Full publication metadata
names_results_table.csv	Name-based query results
orcid_results_table.csv	ORCID-based query results (if applicable)

🧠 How It Works
1. Validates configuration (config.py)
2. Loads and validates input data
3. Generates:
	Name variations
	PubMed search terms
4. Queries PubMed via Entrez API
5. Collects PMIDs and retrieves publication details
6. Matches publications to grants (optional)
7. Outputs cleaned CSV reports

🛠 Project Structure
.
├── pub_query_name_only.py   # Main execution script
├── name_only_lib.py         # Core logic and utilities
├── config.py                # Configuration settings
├── query_table.csv          # Input data
├── Reports/                 # Output files
└── README.md


⚠️ Notes
Ensure date format is MM/DD/YYYY
ORCID must follow format: 0000-0000-0000-0000
Invalid inputs are automatically cleaned or flagged
Large queries may take time due to API limits


🤝 Contributing
Contributions are welcome!

Fork the repository
Create a feature branch
Submit a pull request


📄 License (MIT)
MIT License

Copyright (c) 2019 David Surber

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to do so, subject to the
following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


📬 Contact
David Surber
dsurber@wisc.edu
https://github.com/dsurber/