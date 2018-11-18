# Proteome Manager project
There are two main programs to allow proteomics researchers to download and maintain protein databases. Both programs support file/folder creation and naming to help manage your protein databases. The programs automatically create folders with hold the downloaded database files. Release dates and versions are added to folder names for each species of interest. There are options to add common contaminants and to add reversed decoy sequences. Default lists of species to download are supported so that updates to commonly used databases are easy to do.

### UniProt_reference_proteome_manager.py

Allows users to select from the UniProt [reference proteomes](http://www.uniprot.org/help/reference_proteome) available for [FTP download](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/). There are many more completed proteomes (over 150 thousand) than
reference proteomes (less than 10 thousand). Only reference proteomes are supported via this GUI application.

Reference proteomes consist of "canonical" sequences (one protein from each gene) and
(in many cases) additional protein isoforms. Canonical downloads are named accordingly
and isoforms (if present) are added to canonical sequences and combined databases are
created (they contain "all" in their names).

### Ensembl_proteome_manager.py

Lets users download the main Ensembl vertebrate protein
databases. There are over 70 species available (mid 2017). Ensembl protein databases
are high quality, complete proteomes. Ensembl is gene-centric and the protein sequences are easy to map to their respective genes. Ensembl genomes are often used in transcriptomics and comparing trasncriptomic and proteomic results may be easier if Ensembl protein databases are used. Protein description strings in Ensembl databases are a little verbose and a reformatting step is performed to make the descriptions more
readable.

---
## Getting Started

### Dependancies

The programs were written using Python v.3.6 and later.
They have been tested with Windows 7, Windows 10, and macOS 10.14.
The dependent files are "fasta_lib_Py3.py", "reverse_fasta.py", and "Thermo_contams_fixed.fasta" that are required for either program to run with all features. Many program features will work without a contaminants FASTA file. A different contaninants FASTA file can be used provided the name is kept as
"Thermo_contams_fixed.fasta" or if the "reverse_fasta.py" script source file is modified.

Note: if the FASTA file accessions will be processed by a regular expression, the
accession format for contaminants should be consistent with the main FASTA file.
The contaminant database may need modification. The format of decoy accessions can also
be an issue. Parsing protein accessions with regular expressions is not a good idea
unless thorough testing is performed. Running a search engine is not proper testing!

### Installation

Click [here](https://www.python.org/downloads/) to download the latest Python 3, or you can install a scientific Python distribution like [Anaconda](https://www.anaconda.com/distribution/). Download this repository as a .zip file or clone it to access the files.

### Authors

- Delan Huang (Oregon Health and Sciences University)
- Phil Wilmarth (Oregon Health and Sciences University)

### Additional FASTA tools
There are additional FASTA utilities at [this repository](https://github.com/pwilmart/fasta_utilities). The repository also has some more general background reading on protein FASTA files.

### License

This project is licensed under the MIT License.

---
## Ensembl Proteome Manager User Guide

- Ensembl_proteome_manager.py
- fasta_lib.py
- Ensembl_current_release.pickle
- Ensembl_fixer.py
- default_Ensembl_species.txt

This script uses a GUI window (in addition to some console output) to show you the list of vertebrate species in the current Ensembl release (149 proteomes as of 11/12/2018). There are options to filter the list of proteomes to find those of interest. And options to add contaminants and/or decoys to the downloaded databases. Different contaminant databases can be used. The list of downloaded proteomes can be saved so that those species can be updated more easily. File and folder naming is done automatically to append release information and keep the downloaded databases organized. The FASTA header lines in Ensembl databases are not very friendly for typical researches (my opinion) and they are reformatted and shortened to be more useful.

![Ensembl Main GUI window](/images/Ensembl_1_main_edited.jpeg)

**Ensembl main window.** The GUI has a frame at the top to facilitate searching for proteomes and for specifying how to process the downloaded FASTA files. The lower frame has a left side with the available proteomes and a right side with the desired proteomes to download. The center set of buttons manage the left and right lists and the downloading. There is a status bar at the bottom.

---

![Ensembl filtering controls](/images/Ensembl_2_top_edited.jpeg)

**Filtering the proteome list and processing options.** The left list of proteomes can be filtered based on species names or taxonomy numbers. The searching is not case sensitive and does a simple "in" test. Substrings will return results and the general idea is to make the left list a little shorter so the species of interest can be found more easily. The downloaded databases will be compressed. During decompression, common contaminants can be added from a specified contaminants FASTA file. Sequence reversed decoys can also be added. The two check box options are independent and both can be checked.

---

![Ensembl filter for mouse](/images/Ensembl_3_mouse_edited.jpeg)

**Example of how to get mouse proteomes.** The taxonomy number for mouse is 10090. If we enter that in the TAxonomy ID field, and click the Show Filtered List button, we will get 13 mouse proteomes. Ensembl has specific proteomes for 13 common mouse strains. The top one in the left list is the typical mouse genome of the most commonly used strain (C57BL/6J, I think).

---

![Ensembl selecting downloads](/images/Ensembl_4_add_edited.jpeg)

**Adding mouse to the download list with human.** If we select the first mouse line on the left, then click the Add Proteome(s) button, that proteome is added to the right window. We can click the Download button to download and process these databases.

---

![Ensembl download dialog](/images/Ensembl_5_download_edited.jpeg)

**A dialog box lets you select the location for Ensembl databases on your computer.** The script will take care of creating release version named subfolders. You want to pick a "container" folder for your collection of Ensembl databases. Examination of the subfolders and their contents will give you an idea of the general organization and naming scheme. Some information in the filenames is redundant with information in the folder names on purpose. When adding FASTA files to data repositories or as Supplemental files, all of the release information should be contained in the filename because the file is usually taken out of it folder path context.

---

![Ensembl file organization](/images/Ensembl_6_files.png)

**Subfolder organization.** The compressed downloaded files from the Ensembl FTP site are located in the folders with species information. The decompressed databases have the ".fasta" file extensions. We did not select any processing options, so we just have the target databases without any common contaminants. There is also a log file with the information that was shown in the console window when the script ran. This window also has one of the mouse strains (left over from an earlier testing).

---

## UniProt Reference Proteome Manager User Guide

- UniProt_reference_proteome_manager.py
- fasta_lib.py
- UniProt_current_release.pickle
- default_UniProt_species.txt
- reverse_fasta.py

UniProt has several ways to find and download databases. The main web site options are the easiest to find and use. They have limitations, however. There is a [UniProt FTP](https://www.uniprot.org/downloads) site that is often overlooked. There is a reduced list of higher quality reference proteomes, for example. There are (as of 11/15/2018) 439 archaea, 8895 bacteria, 1184 eukaryota, and 6178 virus reference proteomes available via FTP. The sequence collections for each species are split into a canonical set (sort of a one gene one protein idea) and (optionally) any additional isoforms of canonical proteins.

Protein databases available through the main web site are split into reviewed sequences (Swiss-Prot entries) and unreviewed entries (TrEMBL entries). Swiss-Prot (and only Swiss-Prot) entries can have optional annotated isoforms. The canonical sequence collections contain both Swiss-Prot and TrEMBL entries to make up "complete" proteomes. The higher eukaryotic canonical proteomes all have around 21000 sequences, for example. These canonical databases are, therefore, reasonably complete with minimal peptide redundancy. These databases are particularly good choices for shotgun quantitative proteomics data.

There are README files that provide the mappings from the species names and taxonomy numbers to the UniProt proteome numbers. The actual FTP file listings only have the proteome numbers. This can make finding the right databases to download a little challenging. This script gets information from the FTP site and presents it in a more human friendly format. There is automatic folder creation and file naming logic to help keep databases organized and make sure that the UniProt release information is captured. There are also some convenience options to add common contaminants and decoy sequences. _**Note:** The folder naming options are not yet implemented._      

![UniProt main window](/images/UniProt_1_main_edited.jpeg)

**UniProt reference proteome manager main window.** There is a top pane for controlling what proteomes are presented in the left list of the lower pane. Different kingdoms can be selected, species names can be restricted to specific species names, and taxonomy numbers can also be restricted to those of interest. After downloading databases, they can be processed to add contaminants or decoys (and contaminants). The user can select different contaminant databases if desired.

The bottom pane has available proteomes listed on the left, and the desired databases to download on the right. The right list can be saved as a default list that loads when the GUI launches. There are controls to move proteomes from the left to the right, to drop proteomes from the right list, and download the databases. Databases can be downloaded as canonical sequences only, or canonical sequences and isoform sequences. If isoforms are downloaded, they will be automatically added to the canonical sequences to make a single combined protein FASTA database.

---

![UniProt bovine species search](/images/UniProt_2A_bovine_edited.jpeg)

**Filtering for bovine (cow) sequences.** We can restrict the kingdom to Eukaryota by unchecking the other kingdom boxes. We can require that the species name contain "bovine" (a case insensitive "in" test), and click the Show Filtered List button.

---

![UniProt left list filtered](/images/UniProt_2B_bovine_edited.jpeg)

**Left list will update.** We now have just two possible proteomes on the left. We can select the bovine proteome (taxonomy 9913) and click the Add Proteome(s) button.

---

![UniProt moved to right](/images/UniProt_2C_bovine_edited.jpeg)

**Bovine proteome added to right list.** The bovine proteome has been added to our download list. The right list has some species loaded from our default list that we do not need.

---

![UniProt select to drop](/images/UniProt_3A_drop_edited.jpeg)

**Drop any unneeded proteomes.** We can select the mouse, rat, yeast, and E. coli rows and then click the Drop Proteome(s) button to remove them.

---

![UniProt right updated and download](/images/UniProt_3B_drop_edited.jpeg)

**Ready to download databases.** We are ready to download some protein databases to test if there are human proteins that make us behave more like a herd of cattle, a flock of sheep, or if we really are just a bunch of dirty little pigs. We will download just the canonical sequences and we will add decoys and contaminants so the databases are ready to use with a search engine program, such as [Comet](http://comet-ms.sourceforge.net/).

---

![UniProt download dialog](/images/UniProt_4A_download_edited.jpeg)

**Specify the download location.** We want to select a folder where we will keep all of our UniProt protein database. Subfolder creation, naming, and file naming will be taken care of by the script.

---

![UniProt console](/images/UniProt_4B_download_edited.jpeg)

**Console window also has information.** Download progress is logged to the console window with details on filenames, locations, and sequence counts.

---

![UniProt after download and quit](/images/UniProt_4C_download_edited.jpeg)

**Quit after downloads finish.** The status bar and a dialog alert box will let you know when downloads have finished. If you do not have any additional databases to download, it is time to quit. Click the quit button or close the GUI window. You may also want to quit your Python 3 shell.

---

![UniProt update defaults](/images/UniProt_5_defaults_edited.jpeg)

**Right list can be saved.** The right list might be a collection of species that you want to download on a regular basis. If the current right list differs from the previously saved list, you will be asked if you want to save the changes.

---

![UniProt files](/images/UniProt_6_files_edited.jpeg)

**Example of what downloaded files/folder look like.** The compressed download files are saved in nicely named folders. The downloaded files are decompressed, descriptively named, and any selected processing performed. A FASTA file of the downloaded database is always created in addition to any desired processed versions (with decoys or contaminants). A log file is also present.

---

### Updates and bug fixes

- Feb. 28, 2018 (PW): UniProt FTP README file layout changed and entries were not getting read.

- Nov. 8, 2018 (PW): Some minor GUI improvements and support for different contaminant databases. Has some retry loops for better FTP connectivity.

- Nov. 15, 2018 (PW): Some debugging in Ensembl GUI and updated README.
