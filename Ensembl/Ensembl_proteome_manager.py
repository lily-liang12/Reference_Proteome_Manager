"""'Ensembl_proteome_manager.py' written by Delan Huang, OHSU, July 2017.

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

Delan Huang, 2017-07-12
TODO:
 - More Error Checking (adding/dropping when nothing selected, )
 - Aesthetic/housekeeping changes to both UI and code
 - Overall, program is very rough but functional
"""
# debugging and edits -PW 8/10/2017

# Built-in module imports
from tkinter import *
from tkinter.ttk import *
from tkinter import messagebox
from tkinter import filedialog
import os
import sys
import ftplib
##import datetime
import re
import urllib.request
import pickle
from datetime import datetime

# Imports dependent on other files
# This python file only uses built-in modules, no external downloads required
try:
    import fasta_lib_Py3 as fasta_lib
    import Ensembl_fixer
    import reverse_fasta as add_rev
except ImportError:
    print("Could not import all files.")
    sys.exit()

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
            
class AnimalEntry:
    def __init__(self, c_n, l_n, taxid, e_a, acc, g_m, v_d, r_d, p_a):
        self.common_name = c_n          # Species Common Name (string)
        self.latin_name = l_n           # Species Latin Name (string)
        self.taxID = taxid              # Taxonomy ID Number (int)
        self.ensembl_assembly = e_a     # Ensembl assembly (string?)
        self.accession = acc            # 
        self.genebuild_method = g_m     # 
        self.variation_database = v_d   # 
        self.reg_database = r_d         # 
        self.pre_assembly = p_a         #
        self.folder_name = ""           # Folder Name for each species
        self.ftp_file_path = ""         # Species ftp download path

    # Define getter/setter methods (Can include more as necessary)
    def getCommonName(self):
        return self.common_name
    
    def setCommonName(self, c_n):
        self.common_name = c_n
        
    def getLatinName(self):
        return self.latin_name
    
    def setLatinName(self,l_n):
        self.latin_name = l_n
        
    def getTaxID(self,):
        return self.taxID
    
    def setTaxID(self, taxid):
        self.taxID = taxid
        
    def getEnsemblAssembly(self):
        return self.ensembl_assembly
    
    def setEnsemblAssembly(self, e_a):
        self.ensembl_assembly = e_a
        
    def getFolderName(self):
        return self.folder_name
    
    def setFolderName(self, folder_name):
        self.folder_name = folder_name
        
    def getFTPFile(self):
        return self.ftp_file_path
    
    def setFTPFile(self, file_path):
        self.ftp_file_path = file_path

# Build GUI
class GUI:
    """Main GUI class for application."""
    def __init__(self, url, prot_path, text, headers, banned_list, script_location):
        """Create object and set some state attributes."""
        self.url = url                          # Url of Ensembl FTP site
        self.ensembl_prot_path = prot_path      # Location of Ensembl databases
        self.ensembl_ftp = os.path.dirname(prot_path)   # top FTP level
        self.ftp = None                         # FTP object (set in login method)
        self.text = text                        # HTML text of webpage
        self.raw_table = []                     # HTML text of just animals table
        self.selected_entries = []              # List of selected AnimalEntry objects
        self.animal_list = []                   # List of all AnimalEntry objects
        self.banned_list = banned_list          # List of file identifiers to be omitted when downloading
        self.date = ""                          # This should be the date uploaded(i.e. 2017.07 for July release)
        self.release = ''                       # Ensembl release number
        self.headers = headers                  # Needed for columns in tables
        self.proteome_IDs = []                  # List of unique proteome IDs
        self.script_location = script_location  # Script path location
        self.selected_default = os.path.join(script_location, 'default_Ensembl_species.txt')     # typical default species file path
        self.data = None                        # Holds unpickled information
        self.quit_save_state = "not triggered"  # Trigger for updating defaults file on quit status
        
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
            pass # we will get error if there is no FTP connection to close
        
    def cleanCommonName(self, name):
        p = re.compile(r"alt=\"(.*?)\"")
        m = p.search(name)
        return m.groups()[0]

    def cleanLatinName(self, name):
        p = re.compile(r"<i\b[^>]*>(.*?)</i>")
        m = p.search(name)
        return m.groups()[0]

    def createRawTable(self):
        # Setup html file to find required information 
        # Find start and end of h3 header block
        TEXT = self.text
        if "<td" in TEXT:
            start_ind = TEXT.index("<td")
        if "</table>" in TEXT:
           end_ind = TEXT.index("</table>")

        # Text Block that needs to be parsed
        self.raw_table = TEXT[start_ind:end_ind]

    def loadAllEntries(self):
        """Loads Ensembl proteome entries from pickle file.
        If file does not exist or file is out-of-date, return False.
        """
        # get the ccontents of current_README file
        self.login()
        self.ftp.cwd(self.ensembl_ftp)  # move into current_README file location
        listing = []
        self.ftp.retrlines('RETR current_README', listing.append)

        # Get the current release version from current_README
        for line in listing:
            if "Ensembl Release" in line:
                items = line.split()
                release = int(items[items.index('Release') + 1])

        # see if pickle file exists
        if not os.path.exists(os.path.join(self.script_location, 'Ensembl_current_release.pickle')):
            print('pickle file not present')
            self.release = release
            return False

        # get data from pickle file
        self.data = self.unpickle_entries()
        self.date = self.data["Date"]
        self.release = self.data["Release"]
        self.animal_list = self.data["Entries"]

        # if pickled version matches current database version, then load entries from pickle file
        if self.release == release:
            return True
        else:
            print('saved release is out-of-date')
            self.release = release  # set this to the current release version
            return False

    def parseRawTable(self):
        """Gets Ensembl proteome entries. Looks for pickle file first and checks if current, if not fetches from web."""
        if self.loadAllEntries():
            return  # pickled entries were read in and were current            
        else:
            print('fetching data from web')
            # Parse header into animal list
            # Need an alternative path for missing entries where gene build method is "import"
            parser = re.compile(r'<td\b[^>]*>(.*?)</td>|</span\b[^>]*>(.*?)</span>')
            matched_groups = parser.findall(self.raw_table)
            parsed = []
            for i in range(0, len(matched_groups), 9):  # Split 1D list into 2D so that each animal has 9 attributes
                animal = matched_groups[i:i+9]
                parsed.append(animal)
                
            # We want to remove the empty space produced by alternative path in regex
            for animal in parsed:
                for i in range(len(animal)):
                    for path in animal[i]:
                        if path:
                            animal[i] = path
                common_name = self.cleanCommonName(animal[0])
                latin_name = self.cleanLatinName(animal[1])
                tax_id = animal[2]
                if not str(tax_id).isdigit():  # In case tax_id is something other than a number
                    tax_id = "000"

                # Create main animal entry
                animal_obj = AnimalEntry(common_name, latin_name, tax_id, animal[3], animal[4],
                                         animal[5], animal[6], animal[7], animal[8])

                # Set animal object's folder name and ftp download path
                folder_name = "{}_{}_{}".format(animal_obj.getCommonName(), animal_obj.getLatinName(), animal_obj.getTaxID())
                folder_name = folder_name.replace(" ", "-")
                animal_obj.setFolderName(folder_name)
                
                download_latin_name = latin_name.replace(" ", "_").lower()
                # download_path = os.path.join(self.ensembl_prot_path, download_latin_name, "pep", "")
                download_path = r"{}/{}/pep/".format(self.ensembl_prot_path, download_latin_name)
                animal_obj.setFTPFile(download_path)
                self.animal_list.append(animal_obj)
                
            self.removeInvalidAnimals()
            self.getDate()

            # save the fetched species information
            self.pickle_entries()

    def removeInvalidAnimals(self):
        self.login()
        # If we cant find the animal directory, remove it from animal list
        for animal in self.animal_list:
            try:
                self.ftp.cwd(animal.getFTPFile())
            except ftplib.error_perm:
                self.animal_list.remove(animal)

    def getDate(self):
        self.login()
        # Just hope that this animal actually exists in ftp database, because aardvark doesn't
        self.ftp.cwd(self.animal_list[1].getFTPFile())  
        
        # Create a list of all files in each species folder
        listing = []
        self.ftp.retrlines('LIST', listing.append)
        
        # Download each selected entry's fasta file
        line = listing[0]
        line = line.strip()  # Want last item, so strip EOL
        fname = line.split()[-1]
        modifiedTime = self.ftp.sendcmd('MDTM {}'.format(fname))
        modifiedTime = datetime.strptime(modifiedTime[4:], "%Y%m%d%H%M%S").strftime("%d %B %Y %H:%M:%S")
        # Prints something like "07 May 2017 22:22:07"
        month = modifiedTime.split()[1]
        year = modifiedTime.split()[2]
        self.date = "{}.{}".format(month, year)
        
    def filterEntries(self):
        """Checks values search fields, filters all animals associated with
        taxon numbers, and/or species names, then returns a list with all matching entries.
        """
        # get the species and taxonomy substring filters
        species_entry = self.searchSpecies.get().lower()
        tax_entry = self.searchTax.get()

        # filter on taxonomy number substring
        self.selected_entries = [entry for entry in self.animal_list if tax_entry in entry.getTaxID()]

        # filter on species name substring
        self.selected_entries = [entry for entry in self.selected_entries
                                 if species_entry in entry.getCommonName().lower()
                                 or species_entry in entry.getLatinName().lower()]
        
    def get_filtered_proteome_list(self):
        """Calls relevant methods to create filtered lists, then finds intersection of the lists, 
        and outputs relevant info to user
        """
        self.filterEntries()

        if len(self.selected_entries) == 0:
            # Ask if user wants all entries shown if no filters are selected
            answer = messagebox.askyesno("Are you sure?",
                                         "No databases found. Would you like to show all databases?")
            if answer:
                self.selected_entries = self.animal_list
            else:
                return None
                    
        # Only show relevant info to user in entries
        entries = [[entry.getCommonName(), entry.getLatinName(),
                    entry.getTaxID(), entry.getEnsemblAssembly()]
                    for entry in self.selected_entries]

        # clear entries before importing
        for row in self.tree_left.get_children():
            self.tree_left.delete(row)
        for entry in sorted(entries):
            self.tree_left.insert('', 'end', values=entry)

        self.update_status_bar("List updated with %s entries" % len(self.selected_entries))
        
    def reset_filters(self):
        """Resets filters to defaults."""
        self.searchSpecies.delete(0, END)
        self.searchTax.delete(0, END)
        self.reverse_contams.uncheck_all()
        
    def sort_text_column(self, tv, col, reverse=False):
        """Sorts entries in treeview tables alphabetically."""
        l = [(tv.set(k, col), k) for k in tv.get_children('')]
        l.sort(key=lambda x: x[0].lower(), reverse=reverse)

        # rearrange items in sorted positions
        for index, (val, k) in enumerate(l):
            tv.move(k, '', index)

        # reverse sort next time
        tv.heading(col, command=lambda col_=col: self.sort_text_column(tv, col_, not reverse))
    
    def sort_num_column(self, tv, col, reverse=False):
        """Sorts entries in treeview tables numerically."""
        l = [(tv.set(k, col), k) for k in tv.get_children('')]
        l.sort(key=lambda x: int(x[0]), reverse=reverse)

        # rearrange items in sorted positions
        for index, (val, k) in enumerate(l):
            tv.move(k, '', index)

        # reverse sort next time
        tv.heading(col, command=lambda col_=col: self.sort_num_column(tv, col_, not reverse))
    
    def move_to_left(self):
        """Movies entry(ies) from right treeview to left."""
        selection = self.tree_right.selection()  # creates sets with elements "I001", etc.
        
        for selected in selection:
            selected_copy = self.tree_right.item(selected)  # creates a set of dicts
            self.tree_right.delete(selected)
            self.tree_left.insert('', 'end', values=selected_copy['values'])
        self.update_status_bar("{} dropped".format(selected_copy['values'][0]))

    def move_to_right(self):
        """Movies entry(ies) from left treeview to right."""
        selection = self.tree_left.selection()  
        
        for selected in selection:
            selected_copy = self.tree_left.item(selected)
            self.tree_left.delete(selected)
            self.tree_right.insert('', 'end', values=selected_copy['values'])
        self.update_status_bar("{} added".format(selected_copy['values'][0]))  # Species name should be first

    def pickle_entries(self):
        """Saves full left display list to make subsequent launches faster."""
        text = {"Date": self.date, "Release": self.release, "Entries": self.animal_list}

        # make sure we are in the location with the script
        try:
            os.chdir(self.script_location)
        except OSError:
            print("OSError occurred during pickling. Cwd: {}".format(os.getcwd()))

        with open('Ensembl_current_release.pickle', 'wb') as file:
            pickle.dump(text, file)

    def unpickle_entries(self):
        """Loads saved full left display list of species."""
        with open('Ensembl_current_release.pickle', 'rb') as file:
            return pickle.load(file)

    def save_defaults(self, overwrite=False):
        """Saves species in the right display box to a default species text file"""
        desired_file = self.selected_default
        if not overwrite:
            print('should be asking for save file name')
            desired_file = fasta_lib.save_file(self.script_location, [('Text files', '*.txt')],
                                               default_file=self.selected_default,
                                               title_string='Specify a default species file name')
        if desired_file:
            try:
                # write default species list to file
                items = self.tree_right.get_children()
                databases = [self.tree_right.item(item)['values'] for item in items]
                for database in databases:
                    database[-1] = database[-1].rstrip(r"""\'"*""") # seem to accumulate EOL characters
                
                # Remove duplicates
                db_set = set(tuple(x) for x in databases)
                databases = sorted([list(x) for x in db_set], key=lambda y: y[0]) # sort DBs by common name

                with open(desired_file, "w") as defaults_txt:
                    self.selected_default = desired_file
                    for database in databases:
                        defaults_txt.write("{}\n".format(database))

                self.status_bar.config(text="Databases saved to species text file")
            except OSError:
                messagebox.showwarning("Invalid Filename!", "Cannot save species list to selected folder!")
        
    def select_defaults_and_load(self):
        """Let user browse to a defaults file and load the species."""
        self.selected_default = fasta_lib.get_file(self.script_location,
                                                   [('Text files', '*.txt')],
                                                   'Select a default Ensembl species list file')
        self.load_defaults()
                        
    def load_defaults(self, display=True):
        try:
            with open(self.selected_default, "r") as defaults_txt:
                databases = defaults_txt.readlines()
            self.status_bar.config(text="default species list imported.")
                
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
        if display:
            for row in self.tree_right.get_children():
                    self.tree_right.delete(row)

        loaded_databases = []    
        for database in databases:
            # load the right list from the defaults
            database = database[1:-1] # trim square brackets
            common_name = database.split(', ')[0][1:-1] # trim quotes
            latin_name = database.split(', ')[1][1:-1] # trim quotes
            tax_id = int(database.split(', ')[2])
            e_a = database.split(', ')[3][1:-1]
            e_a = e_a.rstrip(r"""'\"*""") # seem to get extra chracters to remove
            loaded_databases.append([common_name, latin_name, tax_id, e_a])

        loaded_databases = sorted(loaded_databases, key=lambda x: x[0]) # sort DBs by common name

        if display:
            for database in loaded_databases:
                self.tree_right.insert('', 'end', values=database)

        return loaded_databases

    def update_defaults(self):
        """If the entries in right tree do not match original defaults file, ask user to save updated list"""
        right_tree_items = [self.tree_right.item(entry)['values'] for entry in self.tree_right.get_children()]

        # Remove duplicates
        db_set = set(tuple(x) for x in right_tree_items)
        right_tree_items = sorted([list(x) for x in db_set], key=lambda y: y[0]) # sort DBs by common name
                
        # compare current right-side databases to stored defaults
        if right_tree_items != self.load_defaults(display=False):
            if os.path.exists(self.selected_default):
                answer = messagebox.askyesno("Unsaved Progress",
                                             "Right species list differs from defaults! Would you like to save?")
                if answer:
                    self.quit_save_state = True
                    self.save_defaults(overwrite=True)
            else:              
                answer = messagebox.askyesno("Unsaved Progress",
                                             "Save right species list for next time?")
                if answer:
                    self.quit_save_state = True
                    self.save_defaults(overwrite=True)
            
    def download_databases(self):
        """Fetches the database files for the selected species."""
        self.login()    # refresh the FTP connection
        
        # throw warning if no databases selected
        if len(self.tree_right.get_children()) == 0:
               messagebox.showwarning("Empty Selection", "No databases were selected for download!")
               return None  # exit function
            
        # get parent folder location for database download
        db_default = os.getcwd()
        self.abs_dl_path = filedialog.askdirectory(parent=self.root, initialdir=db_default,
                                           title='Select container for Ensembl downloads')
        if not self.abs_dl_path:
            return None

        # Make a separate folder to contain all files
        ensembl_dir_name = r"Ensembl_v{}".format(self.release)
        ensembl_dir_path = os.path.join(self.abs_dl_path, ensembl_dir_name)
        try:
            os.mkdir(ensembl_dir_path)
        except FileExistsError:
            pass
        os.chdir(ensembl_dir_path)

        # Grab entries from right tree view
        download_taxid = [self.tree_right.item(entry)['values'][2] for entry in self.tree_right.get_children()]
        set_download_taxid = list(set(download_taxid))
        if len(download_taxid) != len(set_download_taxid):
            messagebox.showwarning("Duplicates found!", "Duplicate databases were selected and will be ignored!")
            
        # Create a list of selected animal objects from list of tax id's selected
        download_entries = [entry for taxid in download_taxid for entry in self.animal_list
                            if int(taxid) == int(entry.getTaxID())]

        # Change ftp directory for each species
        for entry in download_entries:
            self.ftp.cwd(entry.getFTPFile())

            # Create a folder for each species
            try:
                os.mkdir(os.path.join(ensembl_dir_path, entry.getFolderName()))
                os.chdir(os.path.join(ensembl_dir_path, entry.getFolderName()))
            except FileExistsError:
                os.chdir(os.path.join(ensembl_dir_path, entry.getFolderName()))
                
            # Create a list of all files in each species folder
            listing = []
            self.ftp.retrlines('LIST', listing.append)
            
            # Download each selected entry's fasta file
            for line in listing:
                line = line.strip() # Want last item, so strip EOL
                # self.date = line.split()[-4]  # Gives the month uploaded
                fname = line.split()[-1] # Get the file name
                
                # Skip any files that we do not want to download
                if self.banned_file(fname):
                    continue
                self.update_status_bar("Downloading {} file".format(fname))
                self.ftp.retrbinary('RETR {}'.format(fname), open('{}'.format(fname), 'wb').write)
                print("{} is done downloading".format(fname))
                self.process_fasta_files(os.path.join(ensembl_dir_path, entry.getFolderName(), fname), entry)

        messagebox.showinfo("All Downloads Completed!", "Downloads Finished!")

    def process_fasta_files(self, file_location, entry):
        """Uncompresses FASTA file, reformats descriptions, and does some analysis.
        """
        # analyze and fix descriptions (also uncompresses the file)
        new_fasta_file = Ensembl_fixer.main(file_location)

        # chdir into correct folder and make sure all file paths are set up correctly
        contam_location = self.script_location
        ensembl_dir_name = r"Ensembl_v{}".format(self.release)
        os.chdir(os.path.join(self.abs_dl_path, ensembl_dir_name))
        
        # Add forward/reverse/contams, as specified by checkboxes
        self.addRevSequences(new_fasta_file, contam_location)

    def addRevSequences(self, fasta_file, contam_location):
        """Gets selection value from radiobuttons and then passes those values to imported fasta_reverse function.
        More documentation on how fasta_reverse works can be found in the reverse_fasta.py file.
        """
        reverse_values = list(self.reverse_contams.get_state())
        # Initially set everything to false
        forward = False
        reverse = False
        both = False
        decoy = reverse_values[0]
        contams = reverse_values[1]
        
        if decoy:
            both = True
        if contams and not decoy:
            forward = True
        if not contams:
            contam_location = os.path.join(contam_location, "block")  # Prevent script from finding contams file

        if contams or decoy:        
            add_rev.fasta_reverse(fasta_file, forward, reverse, both, contam_path=contam_location)
        
    def banned_file(self, fname):
        """False if fname in banned list."""
        skip = False
        for ban in BANNED:
            if ban.lower() in fname.lower():
                skip = True
        return skip
    
    def update_status_bar(self, _text):
        """Updates status bar with new text"""
        self.status_bar.config(text=_text)
        self.root.update_idletasks()
        
    def quit_gui(self):
        """Quits the GUI application."""
        self.logout()   # close the FTP connection
        self.update_defaults()
        self.root.withdraw()
        self.root.update_idletasks()
        self.root.destroy()
        sys.exit()


    # Main Create GUI Function
    def create_gui(self):
        """Creates the main GUI window and starts the event loop."""
        self.root = Tk()
        self.root.title("Ensembl Reference Proteome Downloader")
        self.root.geometry("1250x700+250+150")
        self.root.minsize(1250, 650)

        # Check boxes and Import button Frame
        ## Main Frame
        optionFrame = LabelFrame(self.root, text="Options")
        optionFrame.pack(side=TOP, padx=5, pady=5)

        # Additiona database types Frame
        ## Main Frame
        revFrame = LabelFrame(optionFrame, text="Additional Database Types")
        revFrame.pack(fill=BOTH, expand=YES, padx=5, pady=5)
        self.reverse_contams = Checkboxes(revFrame, ["Decoy Database(s)", "Contaminants"])
        self.reverse_contams.pack(side = LEFT, fill=X, padx=5, pady=5)
        
        # Search Window
        ## Main Frame
        searchWindowFrame = LabelFrame(optionFrame, text="Filters")
        searchWindowFrame.pack(side=BOTTOM, fill=BOTH, expand=YES, padx=5, pady=5)
        
        # Create search bars/buttons
        species_frame = Frame(searchWindowFrame)
        species_frame.pack(fill=X, padx=5, pady=5)
        species_label = Label(species_frame, text="Species Name:")
        species_label.pack(side=LEFT, padx=5, pady=5)
        self.searchSpecies = Entry(species_frame)
        self.searchSpecies.pack(side=RIGHT, fill=X, expand=YES, padx=5, pady=5)        

        tax_frame = Frame(searchWindowFrame)
        tax_frame.pack(fill=X, padx=5, pady=5)
        tax_label = Label(tax_frame, text="Taxonomy ID:")
        tax_label.pack(side=LEFT, padx=5, pady=5)
        self.searchTax = Entry(tax_frame)
        self.searchTax.pack(side=RIGHT, fill=X, expand=YES, padx=5, pady=5)       

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
            if col in ["TAX ID"]:
                self.tree_left.heading(col, text=col.title(),
                                       command=lambda col_=col: self.sort_num_column(self.tree_left, col_))
                self.tree_left.column(col, minwidth=25, width=100, stretch=NO, anchor=E)
            else:
                self.tree_left.heading(col, text=col.title(), anchor=W,
                                       command=lambda col_=col: self.sort_text_column(self.tree_left, col_))
                self.tree_left.column(col, minwidth=25, width=150, stretch=NO)
    
        # Add scrollbars to the TreeView 
        leftScrollY = Scrollbar(leftWindowFrame, orient=VERTICAL)
        leftScrollY.pack(side=RIGHT, fill=Y)
        
        leftScrollX = Scrollbar(self.tree_left, orient=HORIZONTAL)
        leftScrollX.pack(side=BOTTOM, fill=X)    

        self.tree_left.config(yscrollcommand=leftScrollY.set, xscrollcommand=leftScrollX.set)
        leftScrollY.config(command = self.tree_left.yview)
        leftScrollX.config(command = self.tree_left.xview)
        
        
        ## Menu Buttons
        button_width = 19
        buttonFrame = LabelFrame(entryFrame, text="Menu Buttons")
        buttonFrame.pack(side=LEFT)

        addButton = Button(buttonFrame, text="Add Proteome(s)", command=self.move_to_right)
        addButton.pack()
        addButton.config(width=button_width)

        removeButton = Button(buttonFrame, text="Drop Proteome(s)", command=self.move_to_left)
        removeButton.pack()
        removeButton.config(width=button_width)
        
        saveButton = Button(buttonFrame, text="Save Defaults", command=self.save_defaults)  
        saveButton.pack()
        saveButton.config(width=button_width)

        importButton = Button(buttonFrame, text="Import Defaults", command=self.load_defaults)
        importButton.pack()
        importButton.config(width=button_width)
        
        downloadButton = Button(buttonFrame, text="Download Databases", command=self.download_databases)
        downloadButton.pack()
        downloadButton.config(width=button_width)

        quitButton = Button(buttonFrame, text="Quit", command=self.quit_gui)
        quitButton.pack()
        quitButton.config(width=button_width)
        

        ## Right Window
        rightWindowFrame = LabelFrame(entryFrame, text="Selected Proteomes")
        rightWindowFrame.pack(fill=BOTH, expand=YES, side=RIGHT, padx=5, pady=10)
        
        self.tree_right = Treeview(rightWindowFrame, columns=self.headers, show="headings")
        self.tree_right.pack(fill=BOTH, expand=YES, side=LEFT, padx=5, pady=5)
        for col in self.headers:
            if col in ["TAX ID"]:
                self.tree_right.heading(col, text=col.title(),
                                       command=lambda col_=col: self.sort_num_column(self.tree_right, col_))
                self.tree_right.column(col, minwidth=25, width=100, stretch=NO, anchor=E)
            else:
                self.tree_right.heading(col, text=col.title(), anchor=W,
                                       command=lambda col_=col: self.sort_text_column(self.tree_right, col_))
                self.tree_right.column(col, minwidth=25, width=150, stretch=NO)
        
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
        self.load_defaults(True)  # initial import of defaults
        self.createRawTable()
        self.parseRawTable()  # Create Entry objects
        self.root.protocol("WM_DELETE_WINDOW", self.quit_gui)  # Override window close event
        self.get_filtered_proteome_list()   # show the full left list to start
        self.root.mainloop()

# Main Function
if __name__ == '__main__':
    # Global Variables
    FTP_URL = 'ftp.ensembl.org'
    PROT_PATH = '/pub/current_fasta'
    HEADERS = ["COMMON NAME", "LATIN NAME", "TAX ID", "ENSEMBL ASSEMBLY"]
    BANNED = ["README", "CHECKSUMS", "abinitio.fa.gz"]
    SCRIPT_LOCATION = os.path.dirname(os.path.realpath(__file__))

    # Get HTML page from Ensembl for parsing
    PARSE_URL = r'http://www.ensembl.org/info/about/species.html'
    RESPONSE = urllib.request.urlopen(PARSE_URL)
    DATA = RESPONSE.read()
    TEXT = DATA.decode('utf-8')
    
    gui = GUI(FTP_URL, PROT_PATH, TEXT, HEADERS, BANNED, SCRIPT_LOCATION)
    gui.create_gui()
