# About
* In an effort to scrape GenBank Genomes for hosts, where host data is very lacking in HEV, we've made a tool that scrapes host data from GenBank entires and tries to figure out which entries have which hosts. 
* We have a finite set of possible hosts that we have predetermined (but not limited to): 
    * Human 
    * Avian
    * Boar 
    * Deer
    * Primate
    * Mongoose
    * Leporidae
    * Camel

* The major challenge in doing so is manually assigning each category (the vast language people use in GenBank) into our finite set.

* Currently, there are a few places that we look for where a host mapping could be in a GenBank entry if it exists, namely the fields: 'host', 'lab_host', 'isolate', and 'isolation_source' in this order. 


# Setup
For setting up python dependencies:
``` 
pip install -R ./requirements.txt 
```
Before using a python file, you may need to set     
```
Entrez.email = 'JohnDoe@gmail.com'
Entrez.apikey = "1234"
```
See https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/ for more info.


# Usage
* Run the pipeline on a GenBank accession list and output the log to logs/VIPBRC-NCBI.txt
```
./query_hepA_nuc.py -p ./accs/ViprBRC-07-04-2022-accession-genes.acc | tee logs/VIPBRC-NCBI.txt
```


# Improving the amount of scraped hosts
* The log files will output entries that it was unable to map. This debug output can help you make manual adjustments to our manual adjustments. There are three places where we store data for our maps: manualMap and species2group in Mapping.py, and hosts.csv.

* manualMap (Mapping.py)
    * The most straightforward way to assign a mapping is by accession number in manualMap. If you see explicitly in the study linked by the GenBank entry that they were used on X host, simply assign it to one of our hosts. However, this is inefficent and should only be used if there is nothing useful in the GenBank entry.

* species2group (Mapping.py)
    * In our approach to parsing the GenBank entries, we look in 'host', 'lab_host', 'isolation_source', and 'isolation_source' for useful strings. From this, there are a lot of unique strings, for which we can assign a host per string.

* hosts.csv
    * For HEV in specific, there is a specific study that has many GenBank entries that are not assignable via our heuristic; however, they tagged each of their entries with a specific "isolate" tag in their GenBank entries, in which they provided a csv file (hosts.csv) which maps isolate tags to host strings that can be mapped using species2group. (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7875884/#SM2 supplemental data)
    * If there is an isolate field in a GenBank entry that you know the host to, add it as a new row in hosts.csv

* Debug output:
    * There are a few warnings that tell you which GenBank entries need assignment and some general advice for each.
	* For each entry marked as "Unable to map to host because we are missing the following entry in mapping.py:", we may simply add a mapping in species2group in Mappings.py to resolve the issue.
	* If a entry has "Unable to map to host because no valid fields existed", examine the GenBank entry if it contains any useful information. You may need to add more functionality in the "Extraction of host features" loop in query_hepA_nuc.py if it's a field we do not check. If there is no useful information, you will need to look through the publication associated with the entry or contact the authors to make a mapping in manualMap.




















