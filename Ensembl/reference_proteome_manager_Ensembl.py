"""'reference_proteome_manager_Ensembl.py' written by Delan Huang, OHSU, July 2017.

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
import urllib.request
import pickle
import gzip
from datetime import datetime

# Imports dependent on other files
# This python file only uses built-in modules, no external downloads required
try:
    import fasta_lib_Py3 as fasta_lib
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
        self.common_name = c_n          # Species Common Name
        self.latin_name = l_n           # Species Latin Name
        self.taxID = taxid              # Taxonomy ID Number
        self.ensembl_assembly = e_a     # 
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
    def __init__(self, url, ref_prot_path, text, headers, banned_list, script_location):
        """Create object and set some state attributes."""
        self.url = url                          # Url of UniProt FTP site
        self.ref_prot_path = ref_prot_path      # Location of databases
        self.ftp = None                         # FTP object (set in login method)
        self.text = text                        # HTML text of webpage
        self.raw_table = []                     # HTML text of just animals table
        self.selected_entries = []              # List of selected AnimalEntry objects
        self.animal_list = []                   # List of all AnimalEntry objects
        self.banned_list = banned_list          # List of file identifiers to be omitted when downloading
        self.date = ""                          # This should be the date uploaded(i.e. 2017.07 for July release)        
        self.headers = headers                  # Needed for columns in tables
        self.proteome_IDs = []                  # List of unique proteome IDs
        self.script_location = script_location  # Script path location
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

    def parseRawTable(self):
        try:
            # Load Entries if already saved from defaults and make sure its up to date
            self.animal_list = self.data['Entries']
            self.getDate()
            if self.date == self.data['Date']:
                return None
            
        except (IndexError, TypeError) as err:  # These errors generally occur if a defaults file hasn't been created yet
            print("Error: ", err)
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
                animal_obj.setFolderName("{}_{}_{}".format(animal_obj.getCommonName(),
                                                           animal_obj.getLatinName(),
                                                           animal_obj.getTaxID()))
                
                download_latin_name = latin_name.replace(" ", "_").lower()
                # download_path = os.path.join(self.ref_prot_path, download_latin_name, "pep", "")
                download_path = r"{}/{}/pep/".format(self.ref_prot_path, download_latin_name)
                animal_obj.setFTPFile(download_path)
                self.animal_list.append(animal_obj)
            self.removeInvalidAnimals()
            self.getDate()

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
        
    # TODO: Species name search is working, but taxID is not
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
        l = [(tv.set(k, col), k) for k in tv.get_children('')]
        l.sort(key=lambda x: x[0].lower(), reverse=reverse)

        # rearrange items in sorted positions
        for index, (val, k) in enumerate(l):
            tv.move(k, '', index)

        # reverse sort next time
        tv.heading(col, command=lambda col_=col: self.sort_text_column(tv, col_, not reverse))
    
    def sort_num_column(self, tv, col, reverse=False):
        l = [(tv.set(k, col), k) for k in tv.get_children('')]
        l.sort(key=lambda x: int(x[0]), reverse=reverse)

        # rearrange items in sorted positions
        for index, (val, k) in enumerate(l):
            tv.move(k, '', index)

        # reverse sort next time
        tv.heading(col, command=lambda col_=col: self.sort_num_column(tv, col_, not reverse))
    
    def move_to_left(self):
        selection = self.tree_right.selection()  # creates sets with elements "I001", etc.
        
        for selected in selection:
            selected_copy = self.tree_right.item(selected)  # creates a set of dicts
            self.tree_right.delete(selected)
            self.tree_left.insert('', 'end', values=selected_copy['values'])
        self.update_status_bar("{} dropped".format(selected_copy['values'][0]))

    def move_to_right(self):
        selection = self.tree_left.selection()  
        
        for selected in selection:
            selected_copy = self.tree_left.item(selected)
            self.tree_left.delete(selected)
            self.tree_right.insert('', 'end', values=selected_copy['values'])
        self.update_status_bar("{} added".format(selected_copy['values'][0]))  # Species name should be first

    def pickle_entries(self, databases):
        text = {"Databases":databases, "Date":self.date, "Entries":self.animal_list}
        try:
            os.chdir(self.script_location)
        except OSError:
            print("OSError occurred during pickling. Cwd: {}".format(os.getcwd()))

        with open('defaults_Ensembl.pickle', 'wb') as file:
            pickle.dump(text, file)

    def unpickle_entries(self):
        with open('defaults_Ensembl.pickle', 'rb') as file:
            return pickle.load(file)

    def save_to_defaults(self):
        answer = True
        # Throw a warning to overwrite if user is trying to save to defaults
        if self.quit_save_state == "not triggered":  # Check to make sure quit state wasn't triggered
            if os.path.isfile("defaults_Ensembl.pickle"):
                answer = messagebox.askyesno("File Detected!",
                                             "A defaults.pickle file was already found. Would you like to overwrite?")
        if answer:
            if self.script_location:
                save_path = self.script_location
            else:
                save_path = filedialog.askdirectory()
                self.script_location = save_path
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

                self.update_status_bar("Databases saved to defaults_Ensembl.pickle")
            except OSError:
                messagebox.showwarning("Invalid Filename!",
                                       "Cannot save defaults_Ensembl.pickle to selected folder!")
                return None
        else:
            return None
        
    def import_defaults(self, initial=False):
        try:
            if initial:
                self.data = self.unpickle_entries()
                # Load in entries from databases
                databases = self.data["Databases"]
                # Save cwd path for save_to_defaults()
                self.update_status_bar("defaults_Ensembl.pickle imported.")
            else:
                self.script_location = filedialog.askopenfilename().rsplit("/", 1)[0] + "/"  
                os.chdir(self.script_location)

                self.data = self.unpickle_entries()
                databases = self.data["Databases"]
                self.update_status_bar("defaults_Ensembl.pickle imported.")
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
            common_name = database[0]
            latin_name = database[1]
            tax_id = database[2]
            e_a = database[3]

            clean_database = [common_name, latin_name, tax_id, e_a]
            self.tree_right.insert('', 'end', values=clean_database)

    def update_defaults(self):
        """If the entries in right tree do not match original defaults file, ask user to save updated list"""
        tree_items = [self.tree_right.item(entry)['values'] for entry in self.tree_right.get_children()]

        # Get data from current defaults file
        try:
            os.chdir(self.script_location)
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
            
    """TODO: This function """        
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
        ensembl_dir_name = r"Ensembl_{}".format(self.date)
        ensembl_dir_path = os.path.join(self.abs_dl_path, ensembl_dir_name)
        try:
            os.mkdir(ensembl_dir_path)
        except FileExistsError:
            pass
        os.chdir(ensembl_dir_path)

        ## TODO: Create new program flow for Ensembl
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
        """Uncompresses canonical FASTA file and does some analysis. Also
        combines fasta and additional fasta files with decompression.
        """
        # Get the list of protein fasta files
        with gzip.open(file_location, 'rb') as in_file:
            file_content = in_file.read()
            
        fasta_file = file_location.replace(".gz", "")
        with open(fasta_file, 'wb') as out_file:
            out_file.write(file_content)
        
        # chdir into correct folder and make sure all file paths are set up correctly
        contam_location = self.script_location
        ensembl_dir_name = r"Ensembl_{}".format(self.date)
        os.chdir(os.path.join(self.abs_dl_path, ensembl_dir_name))
        # Add forward/reverse/contams
        self.addRevSequences(fasta_file, contam_location)

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
        
        if decoy == 1 and contams == 1:
            both = True
        elif decoy == 1 and contams == 0:
            reverse = True
            contam_location = os.path.join(contam_location, "block")  # Prevent script from finding contams file
        elif decoy ==0 and contams == 1:
            forward = True
        else:
            print("Error occurred in determining checkbox values or no selection made!")
        add_rev.fasta_reverse(fasta_file, forward, reverse, both, contam_location)
        
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
        
        saveButton = Button(buttonFrame, text="Save Defaults", command=self.save_to_defaults)  
        saveButton.pack()
        saveButton.config(width=button_width)

        importButton = Button(buttonFrame, text="Import Defaults", command=self.import_defaults)
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
        self.import_defaults(True)  # initial import of defaults
        self.createRawTable()
        self.parseRawTable()  # Create Entry objects
        self.root.protocol("WM_DELETE_WINDOW", self.quit_gui)  # Override window close event
        self.root.mainloop()

# Main Function
if __name__ == '__main__':
    # Global Variables
    FTP_URL = 'ftp.ensembl.org'
    REF_PROT_PATH = '/pub/current_fasta'
    HEADERS = ["COMMON NAME", "LATIN NAME", "TAX ID", "ENSEMBL ASSEMBLY"]
    BANNED = ["README", "CHECKSUMS", "abinitio.fa.gz"]
    SCRIPT_LOCATION = os.path.dirname(os.path.realpath(__file__))

    # Get HTML page from Ensembl for parsing
    PARSE_URL = r'http://www.ensembl.org/info/about/species.html'
    RESPONSE = urllib.request.urlopen(PARSE_URL)
    DATA = RESPONSE.read()
    TEXT = DATA.decode('utf-8')
    
    gui = GUI(FTP_URL, REF_PROT_PATH, TEXT, HEADERS, BANNED, SCRIPT_LOCATION)
    gui.create_gui()
