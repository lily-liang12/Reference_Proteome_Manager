"""'reference_proteome_manager.py' written by Delan Huang, OHSU, July 2017.

The MIT License (MIT)

Copyright (c) 2017 OHSU

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

Direct questions to:
Technology & Research Collaborations, Oregon Health & Science University,
Ph: 503-494-8200, FAX: 503-494-4729, Email: techmgmt@ohsu.edu.
"""

"""TODO: Inspect importing behavior"""
# Built-in module imports
from tkinter import *
from tkinter.ttk import *
from tkinter import messagebox, filedialog
import os
import sys
import ftplib
import datetime
import re
import copy
import pickle
from pathlib import PureWindowsPath

# Imports dependent on other files
# This python file only uses built-in modules, no external downloads required
try:
    import fasta_lib_Py3 as fasta_lib
    import reverse_fasta as add_rev
except ImportError:
    print("Could not import all files.")
    sys.exit("Imports failed!")

# Helper Classes
class Checkboxes(Frame):
    """Creates and packs a set of checkboxes."""
    def __init__(self, parent=None, checkboxes=[], side=LEFT):
        """Constructor creates the checkboxes in the checkboxes list."""
        Frame.__init__(self, parent)
        self.vars = []
        for checkbox in checkboxes:
            var = IntVar()
            check = Checkbutton(self, text=checkbox, variable=var)
            check.pack(side=side, fill=X, expand=YES, padx=10)
            self.vars.append(var)

    def get_state(self):
        """Returns the state of the check boxes"""
        return map((lambda var: var.get()), self.vars)

    def check_all(self):
        """Sets all check boxes to checked."""
        for var in self.vars:
            var.set(1)

    def uncheck_all(self):
        """Unchecks all checkboxes."""
        for var in self.vars:
            var.set(0)

class ReadMeEntry:
    """Container for data parsed from README table rows."""
    def __init__(self, line_entry):
        """Create placeholders for variables and then parse the line"""
        self.kingdom = ""                           # Major phylogenic categories
        self.proteome_ID = ""                       # UniProt refence proteome designation
        self.tax_ID = ""                            # NCBI taxonomy number
        self.oscode = ""                            # UniProt OSCODE string
        self.main_fasta = ""                        # Number of entries in the main fasta file
        self.additional_fasta = ""                  # Number of entries in the additional fasta file
        self.gene2acc = ""                          # Number of entries in the gene2acc file
        self.species_name = ""                      # Latin species name
        self.short_name = ""                        # Shortened species name with underscores
        self.ftp_download_list = []                 # FTP downloadable files for each species
        self.ftp_file_path = ""                     # Kingdom branch path at FTP site
        self.download_folder_name = ""              # More descriptive folder name to hold download files

        # Regular expression for parsing README table rows       
        self.parser = re.compile('^(\S+)\s([0-9]+)\s(.+?)\s+([0-9]+)\s+([0-9]+)\s+([0-9]+)\s+(.*)$')
        
        # List of characters that cannot be in folder names
        self.illegal_pattern = r"[\\#%&{}/<>*?:]"

        self.setAttributes(line_entry)  # Populate object attributes
        self.makeShortName()            # Makes some shorter species names

    # Parse README table line
    def setAttributes(self, line):
        """Parse attributes from table line."""
        m = self.parser.match(line)
        if not m:   # This can be used to skip over rows before or after the main table
            raise ValueError('Invalid line')

        # Get the matching groups and load attributes
        groups = m.groups()
        self.proteome_ID = groups[0]
        self.tax_ID = groups[1]
        self.oscode = groups[2]
        self.main_fasta = groups[3]
        self.additional_fasta = groups[4]
        self.gene2acc = groups[5]
        self.species_name = groups[6]

    def makeShortName(self):
        """Clunky attempt to get a shorter species name to add to download filenames"""
        m = re.match(r"([A-Z][a-z]+\s)+[a-z]+", self.species_name)
        if m:
            self.short_name = re.sub(r"\s", "_", m.group())
            if self.short_name.endswith('_sp'):
                self.short_name = self.short_name[:-3]
                
    def makeFolderName(self, date, dash=True):
        """ This function will remove any characters found in the remove chars list from the
            target string.
        """
        # Remove invalid folder name characters
        fixed_name = re.sub(self.illegal_pattern, "", self.species_name).strip()
        if dash:
            fixed_name = fixed_name.replace(" ", "-")

        # Make the local download folder name
        self.download_folder_name = '_'.join([date, self.proteome_ID, fixed_name])

    # Define all getter/setter methods
    def getProtID(self):
        return self.proteome_ID
    def setProtID(self, protID):
        self.proteome_ID = protID
    def getTaxID(self):
        return self.tax_ID
    def setTaxID(self, taxID):
        self.tax_ID = taxID
    def getOscode(self):
        return self.oscode
    def setOscode(self, oscode):
        self.oscode = oscode
    def getMainFasta(self):
        return self.main_fasta
    def setMainFasta(self, main_fasta):
        self.main_fasta = main_fasta
    def getAdditionalFasta(self):
        return self.additional_fasta
    def setAdditionalFasta(self, additional_fasta):
        self.additional_fasta = additional_fasta
    def getGene2Acc(self):
        return self.gene2acc
    def setGene2Acc(self, gene2acc):
        self.gene2acc = gene2acc
    def getSpeciesName(self):
        return self.species_name
    def setSpeciesName(self, species_name):
        self.species_name = species_name
    def getKingdom(self):
        return self.kingdom
    def setKingdom(self, kingdom):
        self.kingdom = kingdom

    def snoop(self):
        """Diagnostic dump."""
        print('kingdom:', self.kingdom)
        print('proteome ID:', self.proteome_ID)
        print('tax ID:', self.tax_ID)
        print('Oscode:', self.oscode)
        print('Fasta entries:', self.main_fasta)
        print('Additional entries:', self.additional_fasta)
        print('Gene To Acc entries:', self.gene2acc)
        print('species name:', self.species_name)
        print('short name:', self.short_name)
        print('download list:', self.ftp_download_list)
        print('ftp file path:', self.ftp_file_path)
        print('download folder name:', self.download_folder_name)         
    
# Build GUI
class GUI:
    """Main GUI class for application."""
    def __init__(self, URL, ref_prot_path, kingdom_paths, headers, banned_list):
        """Create object and set some state attributes."""
        self.url = URL                          # Url of UniProt FTP site
        self.ftp = None                         # FTP object (set in login method)
        self.ref_prot_path = ref_prot_path      # Specifies top level directory of the Uniprot ftp database
        self.kingdom_paths = kingdom_paths      # List of directory names where files are located (kingdoms)
        self.kingdom_selections = []            # List of subpaths user specified
        self.all_entries = []                   # List of selected entry object attributes
        self.banned_list = banned_list          # List of file identifiers to be omitted when downloading
        self.date = ""                          # This should be a UniProt version (i.e. 2017.07 for July release)        
        self.headers = headers                  # Needed for columns in tables
        self.proteome_IDs = []                  # List of unique proteome IDs
        self.import_root = os.getcwd()          # Path location of defaults.pickle
        self.data = None                        # Data from pickle file
        self.quit_save_state = "not triggered"  # Specifically refers to event when user wants to save database after quitting program
                
        # List of characters that cannot be in folder names
        self.illegal_characters = r"[\\#%&{}/<>*?:]"

    # Helper Class Functions
    def login(self):
        """Open an FTP connection and login."""
        self.ftp = ftplib.FTP()
        self.ftp.connect(str(self.url))
        self.ftp.login()

    def logout(self):
        """Close the FTP connection."""
        try:
            self.ftp.quit()
        except:
            pass # Catch error if there is no FTP connection to close
        
    def loadAllEntries(self):
        self.data = self.unpickle_entries()
        date = self.data["Date"]
        entries = self.data["Entries"]

        exitParse = False

        # If defaults file date matches current database version, then load entries from file
        self.login()
        self.ftp.cwd(self.ref_prot_path)  # move into README file
        listing = []
        self.ftp.retrlines('RETR README', listing.append)

        # Get the release version information
        for line in listing:
            if "release" in line.lower():
                version = line.replace(',', '')
                version = version.replace('_', '.')
                self.date = version.split()[1]

        if self.date == date:
            for entry in entries:
                self.all_entries.append(entry)
            exitParse = True

        return exitParse
        
    def parseReadMe(self):
        """Fetches the README file and parses the table in "ReadMeEntry" objects."""
        self.ftp.cwd(self.ref_prot_path)  # move into README file
        listing = []
        self.ftp.retrlines('RETR README', listing.append)

        # Try to load entry objects from pickle file unless first time running, then user needs to save defaults
        try:
            exitParse = self.loadAllEntries()
            if exitParse:
                return None
        except FileNotFoundError:
            # Get the release version information
            for line in listing:
                if "release" in line.lower():
                    version = line.replace(',', '')
                    version = version.replace('_', '.')
                    self.date = version.split()[1]
        
        # Find and parse the table
        header_index = listing.index('Proteome_ID Tax_ID  OSCODE     #(1)    #(2)    #(3)  Species Name')
        last_index = listing.index('Gene mapping files (*.gene2acc)')
        for line in listing[header_index:last_index]:
            try:
                entry = ReadMeEntry(line)
            except ValueError:
                continue
            self.all_entries.append(entry)

        # Add the kingdom categories and download file lists
        self.get_kingdoms()

    def get_kingdoms(self):
        """Walks the kingdom FTP pages and sets additional entry attributes."""
        for kingdom in self.kingdom_paths:
            kingdom_proteome = {}
            kingdom_path = self.ref_prot_path + kingdom
            self.ftp.cwd(kingdom_path)  # Move into category location

            listing = []    # To hold file listing
            self.ftp.retrlines('LIST', listing.append)   # Get the listing and save

            # Count the number of proteomes (each has several files)
            for line in listing:
                line = line.strip() # Want last item, so strip EOL
                fname = line.split()[-1] # Get the file name
                if fname.split('_')[0].startswith('UP'):
                    key = fname.split('_')[0]   # Parse the reference proteome string

                    # Save all filenames for each species
                    if key in kingdom_proteome:
                        kingdom_proteome[key].append(fname)
                    else:
                        kingdom_proteome[key] = [fname]

            kingdom_keys = list(kingdom_proteome.keys())
            for entry in self.all_entries:
                if entry.proteome_ID in kingdom_keys:
                    entry.ftp_download_list = list(kingdom_proteome[entry.proteome_ID]) # makes copy of list
                    entry.kingdom = kingdom
                    entry.ftp_file_path = kingdom_path  # save file path in new variable

            print(kingdom, 'count is', len(kingdom_keys))

        return
    
    def filterEntries(self):
        """ Checks values from checkboxes and search fields, filters all proteome IDs associated with
        selected kingdoms, taxon numbers, and/or species names, then returns a list with all matching entries.
        """
        # Grab values from checkboxes and assign them to their associated kingdom
        self.checkbox_values = list(self.checkboxes.get_state())
        kingdoms = dict(zip(self.kingdom_paths, self.checkbox_values))

        # Get the species and taxonomy substring filters
        species_entry = self.searchSpecies.get().lower()
        tax_entry = self.searchTax.get()

        # Filter for Kingdoms that were selected
        self.kingdom_selections = [key for key in kingdoms if kingdoms[key] == 1]        
        self.selected_entries = [entry for entry in self.all_entries if entry.kingdom in self.kingdom_selections]

        # Filter on taxonomy number substring
        self.selected_entries = [entry for entry in self.selected_entries if tax_entry in entry.getTaxID()]

        # Filter on species name substring
        self.selected_entries = [entry for entry in self.selected_entries if species_entry in entry.getSpeciesName().lower()]

    def get_filtered_proteome_list(self):
        """ Calls relevant methods to create filtered lists, then finds intersection of the lists, 
        and outputs relevant info to user
        """
        self.filterEntries()        

        if len(self.selected_entries) == 0:
            # Ask if user wants all entries shown if no filters are selected
            answer = messagebox.askyesno("Are you sure?",
                                         "No filters were selected and/or found. Would you like to show all databases?")
            if answer:
                self.selected_entries = self.all_entries
            else:
                return None
                    
        # Only show relevant info to user in entries
        entries = [[entry.getTaxID(), int(entry.getMainFasta()), int(entry.getAdditionalFasta()),
                    entry.getKingdom(), entry.getSpeciesName()]
                    for entry in self.selected_entries]

        # Clear entries before importing
        for row in self.tree_left.get_children():
            self.tree_left.delete(row)
        for entry in sorted(entries):
            self.tree_left.insert('', 'end', values=entry)

        self.update_status_bar("List updated with %s entries" % len(self.selected_entries))        
        
    def reset_filters(self):
        """Resets filters to defaults."""
        self.checkboxes.check_all()
        self.searchSpecies.delete(0, END)
        self.searchTax.delete(0, END)
        self.rb_var.set(2)
        
    def sort_text_column(self, tv, col, reverse=False):
        l = [(tv.set(k, col), k) for k in tv.get_children('')]
        l.sort(key=lambda x: x[0].lower(), reverse=reverse)

        # Rearrange items in sorted positions
        for index, (val, k) in enumerate(l):
            tv.move(k, '', index)

        # Reverse sort next time
        tv.heading(col, command=lambda col_=col: self.sort_text_column(tv, col_, not reverse))
    
    def sort_num_column(self, tv, col, reverse=False):
        l = [(tv.set(k, col), k) for k in tv.get_children('')]
        l.sort(key=lambda x: int(x[0]), reverse=reverse)

        # Rearrange items in sorted positions
        for index, (val, k) in enumerate(l):
            tv.move(k, '', index)

        # Reverse sort next time
        tv.heading(col, command=lambda col_=col: self.sort_num_column(tv, col_, not reverse))
    
    def move_to_left(self):
        selection = self.tree_right.selection()  # Creates sets with elements "I001", etc.
        
        for selected in selection:
            selected_copy = self.tree_right.item(selected)  # Creates a set of dicts
            self.tree_right.delete(selected)
            self.tree_left.insert('', 'end', values=selected_copy['values'])
        try:
            self.update_status_bar("{} dropped".format(selected_copy['values'][-1]))
        except UnboundLocalError:
            print("User tried to remove a proteome even though none was selected!")

    def move_to_right(self):
        selection = self.tree_left.selection()  
        
        for selected in selection:
            selected_copy = self.tree_left.item(selected)
            self.tree_left.delete(selected)
            self.tree_right.insert('', 'end', values=selected_copy['values'])
        try:
            self.update_status_bar("{} added".format(selected_copy['values'][-1]))  # Species name should be last
        except UnboundLocalError:
            print("User tried to add a proteome even though none was selected!")
        
    def pickle_entries(self, databases):
        """Saves list of all entry objects into defaults.pickle"""
        text = {"Databases":databases, "Date":self.date, "Entries":self.all_entries}
        # Make sure we are in the correct folder
        try:
            os.chdir(self.import_root)
        except OSError:
            print("OSError occurred during pickling. Cwd: {}".format(os.getcwd()))
        
        with open('defaults_UniProt.pickle', 'wb') as file:
            pickle.dump(text, file)

    def unpickle_entries(self):
        """Loads list of all entry objects into defaults.pickle"""
        with open('defaults_UniProt.pickle', 'rb') as file:
            return pickle.load(file)

    def save_to_defaults(self):
        answer = True
        # Throw a warning to overwrite if user is trying to save to defaults
        if self.quit_save_state == "not triggered":  # Check to make sure quit state wasn't triggered
            if os.path.isfile("defaults_UniProt.pickle"):
                answer = messagebox.askyesno("File Detected!",
                                             "A defaults.pickle file was already found. Would you like to overwrite?")
        if answer:
            if self.import_root:
                save_path = self.import_root
            else:
                save_path = filedialog.askdirectory()
                self.import_root = save_path
            try:
                os.chdir(save_path)

                # Attempt to write data struct in byte form to defaults.pickle
                items = self.tree_right.get_children()
                databases = [self.tree_right.item(item)['values'] for item in items]
                # Remove duplicates
                db_set = set(tuple(x) for x in databases)
                databases = [list(x) for x in db_set]
                self.pickle_entries(databases)
                # Re-import right tree
                for row in self.tree_right.get_children():
                    self.tree_right.delete(row)
                self.data = self.unpickle_entries()
                for database in self.data["Databases"]:
                    self.tree_right.insert('', 'end', values=database)

                self.update_status_bar("Databases saved to defaults_UniProt.pickle")
            except OSError:
                messagebox.showwarning("Invalid Filename!",
                                       "Cannot save defaults_UniProt.pickle to selected folder!")
                return None
        else:
            return None
        
    def import_defaults(self, initial=False):
        try:
            if initial:
                # Load in entries from databases
                databases = self.data["Databases"]
                # Save cwd path for save_to_defaults()
                self.update_status_bar("defaults_UniProt.pickle imported.")
            else:
                self.import_root = filedialog.askopenfilename().rsplit("/", 1)[0] + "/"  
                os.chdir(self.import_root)

                self.data = self.unpickle_entries()
                databases = self.data["Databases"]
                self.update_status_bar("defaults_UniProt.pickle imported.")
        except FileNotFoundError:
            self.update_status_bar("No defaults imported/defaults could not be found")
            return None
        except OSError:
            messagebox.showwarning("Invalid File!", "Invalid file selection!")
            return None
        except TypeError:
            self.update_status_bar("No defaults imported/defaults could not be found")
            # print("If self.data is None, self.data hasn't been initialized yet: ", type(self.data))
            return None

        
        # Clear selected databases before importing
        for row in self.tree_right.get_children():
                self.tree_right.delete(row)
    
        for database in databases:
            # Make a new zombie list to parse kingdom from species name
            tax = database[0]
            canonical = database[1]
            additional = database[2]
            kingdom = database[3]
            species_name = database[4]

            clean_database = [tax, canonical, additional, kingdom, species_name]
            self.tree_right.insert('', 'end', values=clean_database)

    def addRevSequences(self, fasta_file, contam_location):
        """Gets selection value from radiobuttons and then passes those values to imported fasta_reverse function.
        More documentation on how fasta_reverse works can be found in the reverse_fasta.py file.
        """
        rb_values = self.rb_var.get()
        # Initially set everything to false
        forward = False
        reverse = False
        both = False
        if rb_values == 0:
            forward = True
        elif rb_values == 1:
            reverse = True
        elif rb_values == 2:
            both = True
        else:
            print("Error occurred in determining radiobutton values!")
        add_rev.fasta_reverse(fasta_file, forward, reverse, both, contam_location)
            
    def download_databases(self):
        """Fetches the database files for the selected species."""
        self.login()    # Refresh the FTP connection
        
        # Throw warning if no databases selected
        if len(self.tree_right.get_children()) == 0:
            messagebox.showwarning("Empty Selection", "No databases were selected for download!")
            return None  # Exit function
            
        # Get parent folder location for database download
        db_default = os.getcwd()
        abs_path = filedialog.askdirectory(parent=self.root, initialdir=db_default,
                                           title='Select container for UniProt downloads')
        if not abs_path:
            return None

        # Make a separate folder to contain all files
        uniprot_dir_name = r"UniProt_{}".format(self.date)
        uniprot_dir_path = os.path.join(abs_path, uniprot_dir_name)
        try:
            os.mkdir(uniprot_dir_path)
        except FileExistsError:
            pass
        os.chdir(uniprot_dir_path)

        # Get taxonomy ID numbers for right (download) list
        tax_id_list = [self.tree_right.item(entry)['values'][0] for entry in self.tree_right.get_children()]
        set_tax_id_list = list(set(tax_id_list))  # remove duplicates (if any)
        if len(tax_id_list) != len(set_tax_id_list):
            messagebox.showwarning("Duplicates found!", "Duplicate databases were selected and will be ignored!")

        # Get the entry objects for the right taxonomy numbers
        download_entries = [entry for entry in self.all_entries if int(entry.tax_ID) in set_tax_id_list]

        # Add normalized folder name attribute
        [entry.makeFolderName(self.date) for entry in download_entries]

        for entry in download_entries:
            # Move to the FTP site branch where files are located
            self.ftp.cwd(entry.ftp_file_path)
                
            # Set local location for the download
            try:
                os.mkdir(os.path.join(uniprot_dir_path, entry.download_folder_name))
                os.chdir(os.path.join(uniprot_dir_path, entry.download_folder_name))
            except FileExistsError:
                os.chdir(os.path.join(uniprot_dir_path, entry.download_folder_name))
            except OSError:
                print("OSError")
                print('Download for this entry failed:')
                entry.snoop()
                continue

            # Download reference proteome database(s)
            for file in entry.ftp_download_list:
                # Skip any files that we do not want to download                    
                if self.banned_file(file):
                    continue
                
                # Download the file
##                #(skip if already dowloaded)
##                if os.path.exists(os.path.join(uniprot_dir_path, entry.download_folder_name, file)):
##                    continue
                self.update_status_bar("Downloading {} file".format(file))
                self.ftp.retrbinary('RETR {}'.format(file), open('{}'.format(file), 'wb').write)
                print("{} is done downloading".format(file))

            self.process_fasta_files(uniprot_dir_path, entry)

        messagebox.showinfo("All Downloads Completed!", "Downloads Finished!")

    def banned_file(self, fname):
        """False if fname in banned list."""
        skip = False
        for ban in self.banned_list:
            if ban.lower() in fname.lower():
                skip = True
        return skip

    def process_fasta_files(self, uniprot_dir_path, entry):
        """Uncompresses canonical FASTA file and does some analysis. Also
        combines fasta and additional fasta files with decompression.
        """
        # Get the list of protein fasta files
        temp_files = [x for x in entry.ftp_download_list if 'fasta' in x.lower()]
        fasta_files = []
        combined_files = []
        for f in temp_files:
            if not self.banned_file(f):
                fasta_files.append(f)
        fasta_files.sort()

        fasta_file = fasta_files[0].replace('.fasta.gz', '')
        fasta_file = fasta_file + '_' + entry.short_name + '_canonical.fasta'
        combined_files.append(fasta_file)
        fasta_obj_list = [open(os.path.join(uniprot_dir_path, fasta_file), 'w')]
        if len(fasta_files) == 2:
            fasta_file = fasta_files[1].replace('_additional.fasta.gz', '')
            fasta_file = fasta_file + '_' + entry.short_name + '_all.fasta'
            fasta_obj_list.append(open(os.path.join(uniprot_dir_path, fasta_file), 'w'))
            combined_files.append(fasta_file)

        # Set up to read the fasta file entries and init counters
        print('proteome:', entry.proteome_ID, 'species:', entry.species_name)
        p = fasta_lib.Protein()

        # Read entries and write to new file
        for i, fasta in enumerate(fasta_files):
            sp_count = 0
            iso_count = 0
            tr_count = 0
            p_count = 0
            f = fasta_lib.FastaReader(os.path.join(uniprot_dir_path, entry.download_folder_name, fasta))
            while f.readNextProtein(p, False):
                p_count += 1
                if p.accession.startswith('sp|'):
                    sp_count += 1
                if p.accession.startswith('tr|'):
                    tr_count += 1
                if ('-' in p.accession) or ('Isoform of' in p.description):
                    iso_count += 1
                if i == 0:
                    for obj in fasta_obj_list:
                        p.printProtein(obj)
                else:
                    p.printProtein(fasta_obj_list[i])

            # Print stats
            print('..database:', fasta)
            print('....tot_count:', p_count, 'sp count:', sp_count, 'tr count:', tr_count, 'isoform count:', iso_count)

        # Close output file(s)
        for obj in fasta_obj_list:
            obj.close()

        # chdir into correct folder and make sure all file paths are set up correctly
        cwd = PureWindowsPath(os.getcwd())
        contam_location = cwd.parents[1]  # Contams file is 2 directories up
        os.chdir(cwd.parents[0])
        for file in combined_files:
            # Add forward/reverse/contams
            self.addRevSequences(file, contam_location)
            

##        # print stats
##        print('proteome:', entry.proteome_ID, 'species:', entry.species_name)
##        print('tot_count:', p_count, 'sp count:', sp_count, 'tr count:', tr_count, 'isoform count:', iso_count)

    def update_defaults(self):
        """If the entries in right tree do not match original defaults file, ask user to save updated list"""
        tree_items = [self.tree_right.item(entry)['values'] for entry in self.tree_right.get_children()]

        # Get data from current defaults file
        try:
            self.data = self.unpickle_entries()
        except FileNotFoundError:
            # print("First time defaults file has been created.")
            self.pickle_entries([])

        # Match current selected database to defaults
        try:
            if tree_items != self.data["Databases"]:
                answer = messagebox.askyesno("Unsaved Progress",
                                             "Database selections are different than defaults! Would you like to save?")
                if answer:
                    self.quit_save_state = "triggered"
                    self.save_to_defaults()
        except TypeError:  # Triggers when self.data hasn't been initialized yet
            answer = messagebox.askyesno("Unsaved Progress",
                                         "Would you like to save currently selected databases to defaults?")
            if answer:
                self.quit_save_state = "triggered"
                self.save_to_defaults()
        
    def update_status_bar(self, _text):
        """Updates status bar with new text"""
        self.status_bar.config(text=_text)
        self.status_bar.update_idletasks()
        self.root.after(100)
           
    def quit_gui(self):
        """Quits the GUI application."""
        self.logout()   # Close the FTP connection
        self.update_defaults()
        self.root.withdraw()
        self.root.update_idletasks()
        self.root.destroy()
        sys.exit()


    # Main Create GUI Function
    def create_gui(self):
        """Creates the main GUI window and starts the event loop."""
        self.root = Tk()
        self.root.title("UniProt Reference Proteome Downloader")
        self.root.geometry("1250x750+150+50")
        self.root.minsize(1250, 650)

        # Check boxes and Import button Frame
        ## Main Frame
        optionFrame = LabelFrame(self.root, text="Options")
        optionFrame.pack(side=TOP, padx=5, pady=5)

        ## Kingdom Frame
        kingdomFrame = LabelFrame(optionFrame, text="Kingdoms")
        kingdomFrame.pack(side=TOP, fill=BOTH, expand=YES, padx=5, pady=5)

        ## Generate checkboxes
        self.checkboxes = Checkboxes(kingdomFrame, self.kingdom_paths)
        self.checkboxes.pack(side=LEFT, fill=X)
        self.checkboxes.check_all()

        # Search Window
        ## Main Frame
        searchWindowFrame = LabelFrame(optionFrame, text="Additional Filters")
        searchWindowFrame.pack(side=BOTTOM, fill=BOTH, expand=YES, padx=5, pady=5)

        # Create search bars/buttons
        # Species Search Bar
        species_frame = Frame(searchWindowFrame)
        species_frame.pack(fill=X, padx=5, pady=5)
        species_label = Label(species_frame, text="Species Name:")
        species_label.pack(side=LEFT, padx=5, pady=5)
        self.searchSpecies = Entry(species_frame)
        self.searchSpecies.pack(side=RIGHT, fill=X, expand=YES, padx=5, pady=5)        

        # Taxonomy ID Search Bar
        tax_frame = Frame(searchWindowFrame)
        tax_frame.pack(fill=X, padx=5, pady=5)
        tax_label = Label(tax_frame, text="Taxonomy ID:")
        tax_label.pack(side=LEFT, padx=5, pady=5)
        self.searchTax = Entry(tax_frame)
        self.searchTax.pack(side=RIGHT, fill=X, expand=YES, padx=5, pady=5)

        # Radiobuttons for contams and/or decoy databases
        self.rb_var = IntVar()
        rbutton_frame = LabelFrame(optionFrame, text="Additional Database Types")
        rbutton_frame.pack(fill=X, padx=5, pady=5)
        forward_rb = Radiobutton(rbutton_frame, text="Forward Sequences Only", variable=self.rb_var, value=0)
        rev_rb = Radiobutton(rbutton_frame, text="Reverse Sequences Only", variable=self.rb_var, value=1)
        for_rev_rb = Radiobutton(rbutton_frame, text="Forward and Reverse Sequences", variable=self.rb_var, value=2)
        forward_rb.pack(padx=5, pady=5, anchor=W)
        rev_rb.pack(padx=5, pady=5, anchor=W)
        for_rev_rb.pack(padx=5, pady=5, anchor=W)
        self.rb_var.set(2)  # Set BOTH to default

        ## Show filtered list button and reset filters button
        filter_button = Button(searchWindowFrame, text="Show Filtered List", command=self.get_filtered_proteome_list)
        filter_button.pack(side=LEFT, padx=10, pady=10)
        clear_button = Button(searchWindowFrame, text="Reset Filters", command=self.reset_filters)
        clear_button.pack(side=RIGHT, padx=10, pady=10)

        # Entry mover-thingy Frame
        ## Main Frame
        entryFrame = LabelFrame(self.root, text="Entries")
        entryFrame.pack(side=TOP, fill=BOTH, expand=YES, padx=5, pady=5)

        ## Left Window
        leftWindowFrame = LabelFrame(entryFrame, text="Reference Proteomes")
        leftWindowFrame.pack(fill=BOTH, expand=YES, side=LEFT, padx=5, pady=10)

        # Create TreeView
        self.tree_left = Treeview(leftWindowFrame, columns=self.headers, show="headings")
        self.tree_left.pack(fill=BOTH, expand=YES, side=LEFT, padx=5, pady=5)
        for col in self.headers:
            if col in ["TAX ID", "CANONICAL ENTRIES", "ADDITIONAL ENTRIES"]:
                self.tree_left.heading(col, text=col.title(),
                                       command=lambda col_=col: self.sort_num_column(self.tree_left, col_))
                self.tree_left.column(col, minwidth=25, width=105, stretch=NO, anchor=E)
            else:
                self.tree_left.heading(col, text=col.title(),
                                       command=lambda col_=col: self.sort_text_column(self.tree_left, col_))
                self.tree_left.column(col, minwidth=25, width=105, stretch=NO)
        self.tree_left.heading(self.headers[-1], anchor=W)
        self.tree_left.column(self.headers[-1], minwidth=25, width=650, stretch=YES)  # assumes species name is always last

    
        # Add scrollbars to the TreeView 
        leftScrollY = Scrollbar(leftWindowFrame, orient=VERTICAL)
        leftScrollY.pack(side=RIGHT, fill=Y)
        
        leftScrollX = Scrollbar(self.tree_left, orient=HORIZONTAL)
        leftScrollX.pack(side=BOTTOM, fill=X)    

        self.tree_left.config(yscrollcommand=leftScrollY.set, xscrollcommand=leftScrollX.set)
        leftScrollY.config(command = self.tree_left.yview)
        leftScrollX.config(command = self.tree_left.xview)
        
        
        ## Menu Buttons
        buttonFrame = LabelFrame(entryFrame, text="Menu Buttons")
        buttonFrame.pack(side=LEFT)

        addButton = Button(buttonFrame, text="Add Proteome(s)", command=self.move_to_right)
        addButton.pack()
        addButton.config(width=15)

        removeButton = Button(buttonFrame, text="Drop Proteome(s)", command=self.move_to_left)
        removeButton.pack()
        removeButton.config(width=15)
        
        saveButton = Button(buttonFrame, text="Save Defaults", command=self.save_to_defaults)  
        saveButton.pack()
        saveButton.config(width=15)

        importButton = Button(buttonFrame, text="Import Defaults", command=self.import_defaults)
        importButton.pack()
        importButton.config(width=15)

        downloadButton = Button(buttonFrame, text="Download", command=self.download_databases)
        downloadButton.pack()
        downloadButton.config(width=15)

        quitButton = Button(buttonFrame, text="Quit", command=self.quit_gui)
        quitButton.pack()
        quitButton.config(width=15)

        ## Right Window
        rightWindowFrame = LabelFrame(entryFrame, text="Selected Proteomes")
        rightWindowFrame.pack(fill=BOTH, expand=YES, side=RIGHT, padx=5, pady=10)
        
        self.tree_right = Treeview(rightWindowFrame, columns=self.headers, show="headings")
        self.tree_right.pack(fill=BOTH, expand=YES, side=LEFT, padx=5, pady=5)
        for col in self.headers:
            if col in ["TAX ID", "CANONICAL ENTRIES", "ADDITIONAL ENTRIES"]:
                self.tree_right.heading(col, text=col.title(),
                                       command=lambda col_=col: self.sort_num_column(self.tree_right, col_))
                self.tree_right.column(col, minwidth=25, width=105, stretch=NO, anchor=E)
            else:
                self.tree_right.heading(col, text=col.title(), 
                                       command=lambda col_=col: self.sort_text_column(self.tree_right, col_))
                self.tree_right.column(col, minwidth=25, width=105, stretch=NO)
        self.tree_right.heading(self.headers[-1], anchor=W)
        self.tree_right.column(self.headers[-1], width=650, stretch=YES) # Assumes species names are last
        
        rightScrollX = Scrollbar(self.tree_right, orient=HORIZONTAL)
        rightScrollX.pack(side=BOTTOM, fill=X)

        rightScrollY = Scrollbar(rightWindowFrame, orient=VERTICAL)
        rightScrollY.pack(side=RIGHT, fill=Y)

        self.tree_right.config(yscrollcommand=rightScrollY.set, xscrollcommand=rightScrollX.set)
        rightScrollY.config(command = self.tree_right.yview)
        rightScrollX.config(command = self.tree_right.xview)
        
        # Miscellaneous Frame
        miscFrame = Frame(self.root)
        miscFrame.pack(side=BOTTOM, fill=X, padx=5, pady=5)

        # Status Bar
        status_frame = LabelFrame(miscFrame, text="Status")
        status_frame.pack(side=TOP, fill=X, padx=5, pady=5)
        self.status_bar = Label(status_frame, text="", relief=SUNKEN)
        self.status_bar.pack(fill=X, padx=5, pady=5)
        
        # open the FTP connection
        self.login()
        self.parseReadMe()      # Create Entry objects if there are no entry objects 
        self.import_defaults(True)  # Initial import of defaults
        self.root.protocol("WM_DELETE_WINDOW", self.quit_gui)  # Override window close event
        self.root.mainloop()


# Main Function
if __name__ == '__main__':
    # Global Variables
    URL = 'ftp.uniprot.org'
    REF_PROT_PATH = '/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/'
    KINGDOM_PATHS = ('Archaea', 'Bacteria', 'Eukaryota', 'Viruses')
    HEADERS = ["TAX ID", "CANONICAL ENTRIES", "ADDITIONAL ENTRIES", "KINGDOM", "SPECIES NAME"]
    BANNED = ["DNA", "gene2acc", "idmapping"]
    
    gui = GUI(URL, REF_PROT_PATH, KINGDOM_PATHS, HEADERS, BANNED)
    gui.create_gui()

# End
