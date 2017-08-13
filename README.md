# Reference_Proteome_Manager project
There are two main programs to allow researchers to download and maintain protein databases.
 
"UniProt_reference_proteome_manager.py" allows users to select from the UniProt
[reference proteomes](http://www.uniprot.org/help/reference_proteome) available for [FTP download](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/). There are many more completed proteomes (over 150 thousand) than
reference proteomes (less than 10 thousand). Only reference proteomes are supported.

Reference proteomes consist of "canonical" sequences (one protein from each gene) and
(in many cases) additional protein isoforms. Canonical downloads are named accordingly
and isoforms (if present) are added to canonical sequences and combined databases are
created (they contain "all" in their names).

"Ensembl_proteome_manager.py" lets users download the main Ensembl vertebrate protein
databases. There are over 70 species available (mid 2017). Ensembl protein databases
are high quality, complete proteomes. Ensembl is gene-centric and the protein sequences are easy to map to their respective genes. Ensembl genomes are often used in transcriptomics and comparing trasncriptomic and proteomic results may be easier if Ensembl protein databases are used. Protein description strings in Ensembl databases are 
a little verbose and a reformatting step is performed to make the descriptions more 
readable.

Both programs support file/folder creation and naming to help manage your protein
databases. The programs automatically create folders with hold the downloaded database
files. Release dates and versions are added to folder names for each species of interest.
There are options to add common contaminants and to add reversed decoy sequences. Default
lists of species to download are supported so that updates to commonly used databases are
easy to do.    

## User Guide
A PowerPoint User Guide is included to help users. Please check that for more detailed information on how to use the program. Some screen shot graphics may not exactly match
current program versions but most of the information should be up-to-date.

## Getting Started
### Prerequisites
The programs were written using Python v.3.6.2.
They have been tested with Windows 7, Windows 10, and macOS 10.12.6.
The dependent files are "fasta_lib_Py3.py", "reverse_fasta.py", and "Thermo_contams_fixed.fasta" that are required for either program to run with all features. Many program features will work without a contaminants FASTA file. A different contaninants FASTA file can be used provided the name is kept as
"Thermo_contams_fixed.fasta" or if the "reverse_fasta.py" script source file is modified.

Note: if the FASTA file accessions will be processed by a regular expression, the
accession format for contaminants should be consistent with the main FASTA file.
The contaminant database may need modification. The format of decoy accessions can also
be an issue. Parsing protein accessions with regular expressions is not a good idea
unless thorough testing is performed. Running a search engine is not proper testing!

### Installing
Click [here](https://www.python.org/downloads/release/python-362/) to download Python 3.6.2

Download this repo as a .zip file or clone it to access the files.

## Authors
- Delan Huang (Oregon Health and Sciences Universtity)
- Phil Wilmarth (Oregon Health and Sciences Universtity)

## License
This project is licensed under the MIT License.
