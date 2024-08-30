#package for GUI setup
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog, messagebox
from tkinter import scrolledtext
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import sys
import subprocess
import threading
from glob import glob
import pkg_resources


#package for astro analyis
from glob import glob
from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clip
import os 
import pyregion
import re
from scipy.interpolate import interp1d, splrep, splev
from astropy.coordinates import SkyCoord
import astropy.units as u
from typing import List

from .utils import *



#DRP required package
from scipy.ndimage import shift
import math
import ref_index



"""
Author: Zhuyun Zhuang, with the help from chatGPT! (free version!)
Date: 05/28/2024
Usage: python kcwi_zap_viewer.py

Description: Load data from a directory and allow user to interactively
perform sky subtraction of KCWI/KCRM data with ZAP, flux calibration and telluric correction.

This version is still under development. If you run into any issue, please drop a message to zzhuang at astro.caltech.edu 
"""
####Some user-defined setup. People should update it based on their own needs######
#only set it for test purporse. When it browse the directory, it will start from your favorite directory storing the data :)
# initial_dir = '/scr/zzhuang/keck_obs/kcwi'
initial_dir = os.getcwd()

#Set it to the place where you put the TelFit file from pypeit. Can download it via "pypeit_install_telluric TelFit_MaunaKea_3100_26100_R20000.fits"
#Please do not download the TelPCA file (the default of Pypeit). The updated TelPCA file would cause weird shape in the telluric model so please stick to TelFit_MaunaKea!
#telgridfile = '/Users/yuguangchen/.pypeit/cache/download/url/5f17ecc1fcc921d6ec01e18d931ec2f8/content'

#pick this region to generate the white-lighted image because sky lines are much stronger elsewhere. 
# TODO: Can also make it as an input or variable parameter in the GUI
wlimg_wave_range_red = [6380, 7200] 
wlimg_wave_range_blue = [3600, 5500] 

#the reference wavelength for DAR correction, set it to be the one used for M_BL_4500 to avoid second alignment. 
# If you don't want to use the function, please set it to None
DAR_ref_wave = 4577.858334891979
DAR_shift_order = 1


        
class KCWIViewerApp:
    def __init__(self, root):
        """
        Initialize the GUI window
        """
        self.root = root
        self.root.title("KCWI Data Viewer")
        self.last_focused_entry = None

        # Create a custom style
        style = ttk.Style()

        # Create a new style for the Notebook tabs
        style.configure('TNotebook.Tab', font=('Helvetica', '10', 'bold'), padding=[10, 4])

        # Ensure the selected tab has a different foreground color
        style.map('TNotebook.Tab',
                foreground=[('selected', 'black')],
                background=[('selected', 'lightgrey')])

        ############ Input directory #################
        self.input_frame = tk.Frame(root)
        self.input_frame.grid(row=0, column=0, columnspan=3, sticky='nsew', padx=5, pady=5)

        self.input_dir_label = tk.Label(self.input_frame, text="Input Directory:")
        self.input_dir_label.grid(row=0, column=0, sticky='ew')
        self.input_dir_entry = tk.Entry(self.input_frame)
        self.input_dir_entry.insert(0, initial_dir)
        self.input_dir_entry.grid(row=0, column=1, sticky='ew')
        self.browse_button = tk.Button(self.input_frame, text="Browse", command=self.browse_input_directory)
        self.browse_button.grid(row=0, column=2, sticky='ew')

        ############# Output directory############
        self.output_frame = tk.Frame(root)
        self.output_frame.grid(row=0, column=3, columnspan=3, sticky='nsew', padx=5, pady=5)

        self.output_dir_label = tk.Label(self.output_frame, text="Output Directory:")
        self.output_dir_label.grid(row=0, column=0, sticky='ew')
        self.output_dir_entry = tk.Entry(self.output_frame)
        self.output_dir_entry.insert(0, os.path.join(initial_dir, 'kskywizard'))
        self.output = self.output_dir_entry.get()
        self.output_dir_entry.grid(row=0, column=1, sticky='ew')
        self.browse_button = tk.Button(self.output_frame, text="Browse", command=self.browse_output_directory)
        self.browse_button.grid(row=0, column=2, sticky='ew')

        self.input_frame.columnconfigure(0, weight=1)
        self.input_frame.columnconfigure(1, weight=1)
        self.input_frame.columnconfigure(2, weight=1)
        self.output_frame.columnconfigure(0, weight=1)
        self.output_frame.columnconfigure(1, weight=1)
        self.output_frame.columnconfigure(2, weight=1)

        # Create a Notebook widget (for tabs)
        self.notebook = ttk.Notebook(root, style='TNotebook')
        self.notebook.grid(row=1, column=0, columnspan=6, sticky="nsew")

        # Create two frames to hold the content of each tab
        self.tab1 = ttk.Frame(self.notebook)
        self.tab2 = ttk.Frame(self.notebook)

        # Add the tabs to the notebook
        self.notebook.add(self.tab1, text="Invsens")
        self.notebook.add(self.tab2, text="Science")
        
        
        ###### tab 1 #####
        ############# Setup the flux calibration part for stds ###############
        #self.stddir = re.sub('py/kcwi_tools.py', 'data/stds', kcwi_tools.__file__)  #the base directory to read in the flux-calibrated spectrum of a given std
        self.stddir = pkg_resources.resource_filename(__name__, 'data/stds')
        self.std_index_label = tk.Label(self.tab1, text="Standard star invsens (DRP): ")
        self.std_index_label.grid(row = 1, column = 0, sticky='ew')
        self.std_entry = tk.Entry(self.tab1)
        self.std_entry.grid(row = 1, column = 1, sticky='ew')
        self.load_std_button = tk.Button(self.tab1, text = 'Browse DRP invsens', command = lambda: self.load_invsens('DRP'))
        self.load_std_button.grid(row = 1, column = 2, sticky='ew')
        self.load_std_update_button = tk.Button(self.tab1, text = 'Browse updated invsens', command = lambda: self.load_invsens('updated'))
        self.load_std_update_button.grid(row = 1, column = 3, sticky='ew')
        self.save_std_button = tk.Button(self.tab1, text = 'Save updated invsens', command = self.save_updated_invsens)
        self.save_std_button.grid(row = 1, column = 4, sticky='ew')
        # self.std_bspline = tk.Label(root, text = 'B-Spline pars [breakpoints, polyorder]:')
        # self.std_bspline.grid(row = 4, column = 3)
        # self.std_bspline_entry = tk.Entry(root)
        # self.std_bspline_entry.grid(row = 4, column = 5)
        # self.std_bspline_entry.bind("<Return>", self.update_bspline_pars)

        self.tab1.columnconfigure(0, weight=1)
        self.tab1.columnconfigure(1, weight=1)
        self.tab1.columnconfigure(2, weight=1)
        self.tab1.columnconfigure(3, weight=1)
        self.tab1.columnconfigure(4, weight=1)

        ######## tab 2 #############
        ############# Science index input############
        self.index_label = tk.Label(self.tab2, text="Science Frame No.:")
        self.index_label.grid(row=2, column=0, sticky='ew')
        self.index_entry = tk.Entry(self.tab2)
        self.index_entry.grid(row=2, column=1, sticky='ew')
        self.index_entry.bind("<Return>", self.update_index)

        ############# Buttons for increasing and decreasing index for science frame############
        self.increase_button = tk.Button(self.tab2, text="Previous", command=self.decrease_index)
        self.increase_button.grid(row=2, column=2, sticky='ew')
        self.decrease_button = tk.Button(self.tab2, text="Next", command=self.increase_index)
        self.decrease_button.grid(row=2, column=3, sticky='ew')

        #############Check box to select whether an off-field sky is used############
        self.use_index2_var = tk.BooleanVar()
        self.use_index2_checkbox = tk.Checkbutton(self.tab2, text="Use Off-field Sky Frame No.", variable=self.use_index2_var, command=self.toggle_index2_entry)
        self.use_index2_checkbox.grid(row=2, column=4, sticky='ew')

        ############# sky index input############
        self.index2_entry = tk.Entry(self.tab2)
        self.index2_entry.grid(row=2, column=5, sticky='ew')
        self.index2_entry.bind("<Return>", self.update_index2)

        ############# Naming structure selection - used for loading the data cube for science data
        # self.structure_label = tk.Label(root, text="Science Data Product type:")
        # self.structure_label.grid(row=2, column=0)
        self.structure_var = tk.StringVar()
        self.structure_options = ["icubed", "icube", "icubes"]
        self.structure_var.set('icubed')
        self.structure_menu = tk.OptionMenu(self.tab2, self.structure_var, *self.structure_options)
        self.structure_menu.grid(row=3, column=0, sticky='ew')


        #############  Load raw data (DRP-reduced cube) button############
        self.load_button = tk.Button(self.tab2, text="Load Raw Cube", command=lambda: self.load_data('raw'))
        self.load_button.grid(row=3, column=1, sticky='ew')

        #############  Save cropped data (good wavelength region) button############
        self.load_button = tk.Button(self.tab2, text="Save Cropped Cube", command=self.save_cropped_data)
        self.load_button.grid(row=3, column=2, sticky='ew')

        ############# Load cropped data button############
        self.load_crop_button = tk.Button(self.tab2, text="Load Cropped Cube", command=lambda: self.load_data('cropped'))
        self.load_crop_button.grid(row=3, column=3, sticky='ew')

        ############# input the redshift of a given source #########
        self.redshift_label = tk.Label(self.tab2, text = 'Redshift:')
        self.redshift_label.grid(row = 3, column =4, sticky='ew')
        self.redshift_entry = tk.Entry(self.tab2)
        self.redshift_entry.grid(row = 3, column = 5, sticky='ew')
        self.redshift_entry.bind("<Return>", self.update_redshift)


        ############# Input the mask frame No. (can be different from the science one since we usally take the at least three frames) ############
        # This function is only used for convenience, so that for each pair we only need to create one 
        self.mask_index_label = tk.Label(self.tab2, text="ZAP Mask Frame No.:")
        self.mask_index_label.grid(row=4, column=0, sticky='ew')
        self.mask_entry = tk.Entry(self.tab2)
        self.mask_entry.grid(row = 4, column=1, sticky='ew')
        self.mask_entry.bind("<Return>", self.update_mindex)
        self.update_mask_button = tk.Button(self.tab2, text = 'Update ZAP mask', command=self.update_zap_mask)
        self.update_mask_button.grid(row = 4, column = 2, sticky='ew')

        #############Check box to select whether multiple skysegment used, and whether have an additional sky seg near Halpha############
        self.use_multi_skyseg = tk.BooleanVar()
        self.use_multi_skyseg_checkbox = tk.Checkbutton(self.tab2, text="Use multiple skyseg in ZAP", variable=self.use_multi_skyseg, 
                                                        onvalue = True, offvalue = False)
        self.use_multi_skyseg_checkbox.grid(row=4, column=3, sticky='ew')
        self.use_multi_skyseg.set(True)
        self.use_Ha_seg = tk.BooleanVar()
        self.use_Ha_seg_checkbox = tk.Checkbutton(self.tab2, text='Additional Sky Seg near Halpha', variable=self.use_Ha_seg,
                                                  onvalue = True, offvalue = False)
        self.use_Ha_seg_checkbox.grid(row = 4, column = 4, sticky='ew')
        self.use_Ha_seg.set(True)


        #############Run button for the ZAP ############
        self.run_zap_button = tk.Button(self.tab2, text = 'Run ZAP', command = self.run_zap_precondition)
        self.run_zap_button.grid(row =4, column =5, sticky='ew')

        self.tab2.columnconfigure(0, weight=1)
        self.tab2.columnconfigure(1, weight=1)
        self.tab2.columnconfigure(2, weight=1)
        self.tab2.columnconfigure(3, weight=1)
        self.tab2.columnconfigure(4, weight=1)
        self.tab2.columnconfigure(5, weight=1)

        self.tab2.rowconfigure(2, weight=1)
        self.tab2.rowconfigure(3, weight=1)
        self.tab2.rowconfigure(4, weight=1)
        

        
        ######################## Plotting area
        self.figure = Figure(figsize = (12, 4))
        self.ax = self.figure.add_subplot(111)
        self.figure.tight_layout(rect=[0.05, 0.03, 1, 0.92])
        # Global font settings
        plt.rc('font', family='serif', size=12)

        self.canvas = FigureCanvasTkAgg(self.figure, master=root)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=5, column=0, columnspan=6, sticky='nsew')
        self.canvas.get_tk_widget().focus_force()


        # Add navigation toolbar
        self.toolbarFrame = tk.Frame(master=root)
        self.toolbarFrame.grid(row=6,column=2, columnspan = 2, sticky='ew')
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)
        # self.toolbar.home_toggle(True)

        #add the interactive 
        # self.canvas.get_tk_widget().bind('<Button-1>', lambda event: self.canvas.focus_set()) 
        self.canvas.get_tk_widget().focus_set()#only activate the interactive mode when left click the left canvas first
        self.canvas.mpl_connect('button_press_event', self.on_click) #add all the interactive functions with mouse button 
        self.canvas.mpl_connect('key_press_event', self.on_type) #add all the interactive functions with keyboard press

        #Add the entry to update the indices of x and y coordinates in the DS9 region used for spectrum extraction
        self.region_box_label = tk.Label(root, text = 'Region for spectrum extraction from cube (x1, y1, x2, y2):')
        self.region_box_label.grid(row = 6, column = 4, sticky='ew')
        self.region_box_entry = tk.Entry(root)
        self.region_box_entry.grid(row = 6, column = 5, sticky='ew')
        self.region_box_entry.bind("<Return>", self.update_region_idx)



        # Output text widget
        # self.output_text = tk.Text(root, height=5, width=100, font=("Arial", 12))
        self.output_text = scrolledtext.ScrolledText(root, height=10, width=140, font=("Arial", 12), wrap="word")
        self.output_text.grid(row=7, column=0, columnspan=6, sticky='nsew')
        sys.stderr = RedirectText(self.output_text)
        sys.stdout = RedirectText(self.output_text)
        self.redirect_text = RedirectText(self.output_text)
        self.output_text.bind("<Return>", self.update_zap_skyseg)

        #text widget for adjusting sky segment
        # self.output_text = scrolledtext.ScrolledText(root, height = 5, width = 40, font = ('Arial', 12))
        # self.output_text.grid(row =7, column = 4, columnspan = 2)
        # self.output_text.bind('<Return>', self.ontype_skyseg)
        # self.skyseg_input_active = False

        # Initialize some parameters
        self.index = -1 #science frame number
        self.index2 = -1 #sky frame number; if specified
        self.redshift = 0.0
        self.redshift_entry.insert(tk.END, '0.0')
        self.index2_entry.config(state=tk.DISABLED)
        self.update_index_entries()
        self.zap = {}
        self.ctype = self.structure_var.get()


        ######### final configuration
        # self.menubar.add_cascade(label="File", menu=self.filemenu)
        # root.config(menu = self.menubar)
        root.columnconfigure(0, weight=1)
        root.columnconfigure(1, weight=1)
        root.columnconfigure(2, weight=1)
        root.columnconfigure(3, weight=1)
        root.columnconfigure(4, weight=1)
        root.columnconfigure(5, weight=1)
        
        root.rowconfigure(1, weight=1)
        root.rowconfigure(5, weight=1)
        root.rowconfigure(7, weight=1)




        #read in the line table
        line_list_fn = pkg_resources.resource_filename(__name__, 'data/lines.dat')
        self.line_tbl = Table.read(line_list_fn, format = 'csv')

    def insert_text(self, text, newline = True, color="black"):

        # Configure the tag with the specified text color
        self.output_text.tag_configure("colored_text", foreground=color)

        #finish this line
        if newline:
            self.output_text.insert(tk.END, text + '\n', "colored_text")

        #allow other texts to be input into the current line
        else:
            self.output_text.insert(tk.END, text, "colored_text")


        self.output_text.see(tk.END)

    def focus_in(self, event):
        """
        A variable to keep track of which entry widget was last focused, so that the code would only update one variable at one time.
        """
        self.last_focused_entry = event.widget


    def toggle_index2_entry(self):
        """Decide whether index2 (off-field sky frame should be loaded)"""

        self.index2_entry.config(state=tk.NORMAL if self.use_index2_var.get() else tk.DISABLED)
        #when the off-field sky option is turned off
        if not self.use_index2_var.get():
            self.index2_entry.delete(0, tk.END)
            self.index2 = -1
        else:
            self.index2_entry.delete(0, tk.END)

    def increase_index(self):
        """"Increase the science frame number"""
        self.index += 1
        # self.update_index_entries()
        self.index_entry.delete(0, tk.END)
        self.index_entry.insert(tk.END, str(self.index))
        self.insert_text(f"[INFO] Set the science frame: {self.prefix}_{self.index:05d}")

        #drop the attribute "zap" to reinitialize the ZAP configuration
        self.zap = {}

    def decrease_index(self):
        """Decrease the science frame number"""
        self.index = max(0, self.index - 1)
        # self.update_index_entries()
        self.index_entry.delete(0, tk.END)
        self.index_entry.insert(tk.END, str(self.index))
        self.insert_text(f"[INFO] Set the science frame: {self.prefix}_{self.index:05d}")

        #drop the attribute "zap" to reinitialize the ZAP configuration
        self.zap = {}

    def update_index_entries(self):
        #if update the science frame number
        if self.last_focused_entry == self.index_entry and self.index > 0:
            self.index_entry.delete(0, tk.END)
            self.index_entry.insert(tk.END, str(self.index))
            self.insert_text(f"[INFO] Set the science frame: {self.prefix}_{self.index:05d}")
            # self.index_entry.selection_clear()

        #if update the sky frame number
        if self.last_focused_entry == self.index2_entry and self.index2 > 0:
            self.index2_entry.delete(0, tk.END)
            self.index2_entry.insert(tk.END, str(self.index2))
            self.insert_text(f"[INFO] Set the sky frame for science frame {self.prefix}_{self.index:05d}: {self.prefix}_{self.index2:05d} ")
            # self.index2_entry.selection_clear()

        #if update the mask frame number
        if self.last_focused_entry == self.mask_entry:
            self.mask_entry.delete(0, tk.END)
            self.mask_entry.insert(tk.END, str(self.mindex))
            self.insert_text(f"[INFO] Set the ZAP mask frame for science frame {self.prefix}_{self.index:05d}: {self.prefix}_{self.mindex:05d} ")

        #update the redshift
        if self.last_focused_entry == self.redshift_entry:
            self.redshift_entry.delete(0, tk.END)
            self.redshift_entry.insert(tk.END, str(self.redshift))
            self.insert_text(f"[INFO] Set the redshift to z = {self.redshift:0.3f} ")

        #update the region box
        if self.last_focused_entry == self.region_box_entry:
            self.region_box_entry.delete(0, tk.END)
            self.region_box_entry.insert(tk.END, ','.join([str(idx + 1) for idx in self.region_idx]))
            self.insert_text(f"[INFO] Extract the spectrum from (x1, y1, x2, y2) = " + ','.join([str(idx + 1) for idx in self.region_idx]))

            # ', '.join(map(str, self.zap['skyseg']))

         # Set focus back to the master window to allow canvas to capture keys
        self.root.focus_set()


    def update_index(self, event):
        try:
            self.index = int(self.index_entry.get())

            #drop the attribute "zap" to reinitialize the ZAP configuration
            self.zap = {}

        except ValueError:
            pass
        self.update_index_entries()

    def update_index2(self, event):
        """
        Update the sky frame number if the "Use off-field sky" box is selected
        """
        try:
            self.index2 = int(self.index2_entry.get()) #get the sky frame number
        except ValueError:
            pass #no sky frame
        self.update_index_entries()

    def update_mindex(self, event):
        """
        Update the ZAP mask frame number
        """
        try:
            #get the mask frame number        
            self.mindex = int(self.mask_entry.get())
        except ValueError:
            self.insert_text(f"[ERROR] Need to specify the frame used for generating the object mask!")
        self.update_index_entries()

    def update_redshift(self, event):
        """
        Update the input redshift
        """
        try:
            #get the mask frame number        
            self.redshift = float(self.redshift_entry.get())
        except ValueError:
            self.redshift = 0.0
            # self.insert_text(f"[INFO] Redshift not set. Adopt z = 0")
        self.update_index_entries()

    def update_region_idx(self, event):
        """
        Read in the input and update the regions in ds9 to extract spectrum from datacube
        """
        try:
            region_arr = self.region_box_entry.get().split(',')
            self.region_idx = [int(i)-1 for i in region_arr] #minus 1 as ds9 used 1-based index while python is 0-based index
        except:
            self.insert_text('[ERROR] Invalid input of the region. Need to be "x1, y1, x2, y2" in ds9 coordinates. Use the default coordinates')
            self.region_idx = [12, 43, 22, 53]

        self.plot_spec_dict['xrange'] = [self.region_idx[0], self.region_idx[2]]
        self.plot_spec_dict['yrange'] = [self.region_idx[1], self.region_idx[3]]
        self.plot_spec_dict['restore_limit'] = False

        self.plot_spectrum(**self.plot_spec_dict)

            
        # self.plot_spectrum(hdu = self.scihdu[0].data, xrange = [self.region_idx[0], self.region_idx[2]], yrange = [self.region_idx[1], self.region_idx[3]])
        self.update_index_entries()

            
            

    def browse_input_directory(self):
        """Select the input directory"""

        #specify the input directory
        directory = filedialog.askdirectory(initialdir=initial_dir)

        # Check if the directory is not empty (i.e., the user didn't click "Cancel")
        if directory:
            self.input_dir_entry.delete(0, tk.END)
            self.input_dir_entry.insert(tk.END, directory)
            self.base = self.input_dir_entry.get()
            self.insert_text(f"[INFO] Loading the data from {self.base}")

            # Find the common string name for a given date
            try:
                self.prefix = os.path.basename(glob(f'{directory}/k*.fits')[0])[:8]
            except IndexError:
                self.prefix = None

    def browse_output_directory(self):
        """Select the output directory"""
        directory = filedialog.askdirectory(initialdir=initial_dir)

        if directory:    
            self.output_dir_entry.delete(0, tk.END)
            self.output_dir_entry.insert(tk.END, directory)
            self.output = self.output_dir_entry.get()
            # print(f'[INFO] Set the output directory to {self.output}')
            self.insert_text(f"[INFO] Set the output directory to {self.output}")

    def load_data(self, datatype):
        """
        Load the DRP-reduced data cube
        """
        #get the cube type
        self.ctype = self.structure_var.get()

        #determine where to load the data
        if datatype == 'raw':
            base = self.base
        elif datatype == 'cropped':
            base = self.output
        else:
            self.insert_text(f'[ERROR] datatype not recognized! Need to be "raw" or "cropped" ')

        # load the science frame
        if self.index > 0:
            # self.scipath = 
            self.scihdu = fits.open(f'{base}/{self.prefix}_{self.index:05d}_{self.ctype}.fits')
            self.objname = self.scihdu[0].header['OBJECT']
            if datatype == 'raw':
                self.scihdu = self.crop_cube(self.scihdu) #crop the data cube to good wavelength region
                self.insert_text(f"[INFO] Loading the DRP-reduced science frame {self.prefix}_{self.index:05d}")
            elif datatype == 'cropped':
                self.insert_text(f"[INFO] Loading the cropped DRP-reduced science frame {self.prefix}_{self.index:05d}")

            self.scihdr = self.scihdu[0].header
            self.obswave = (np.arange(self.scihdr['NAXIS3']) + 1 - self.scihdr['CRPIX3']) * self.scihdr['CD3_3'] + self.scihdr['CRVAL3']
            if self.obswave[-1] > wlimg_wave_range_red[0]:
                self.wlimg_wave_range = wlimg_wave_range_red
            else:
                self.wlimg_wave_range = wlimg_wave_range_blue



            #replace the bad pixels (flags >0) with NaNs
            # self.scihdu[0].data[self.scihdu['FLAGS'].data > 0] = np.nan 
            # self.scihdu['UNCERT'].data[self.scihdu['FLAGS'].data > 0] = np.nan 

        else:
            self.insert_text(f"[ERROR] Wrong science frame! Need to set it to a positive integer. Check the KCWI log for the frame number!")

        # self.z = 0
        if self.index2 > 0:
            self.skypath = f'{base}/{self.prefix}_{self.index2:05d}_{self.ctype}.fits'
            self.skyhdu = fits.open(self.skypath)
            if datatype == 'raw':
                self.skyhdu = self.crop_cube(self.skyhdu) #crop the data cube to good wavelength region
                self.insert_text(f"[INFO] Loading the DRP-reduced sky frame {self.prefix}_{self.index2:05d}")
            else:
                self.insert_text(f"[INFO] Loading the cropped DRP-reduced sky frame {self.prefix}_{self.index2:05d} for {self.prefix}_{self.index:05d}")
            # self.plot_spectrum(self.scihdu[0].data, hdu_sky=self.skyhdu[0].data, 
            #                     yunit = self.scihdr['BUNIT']) #plot the spectrum of the central region (x = [12, 22], y = [43, 53]; 10x10 box) for a quick look. 
            self.plot_spec_dict = {'datacube': self.scihdu[0].data, 'errcube': self.scihdu['UNCERT'].data,
                                    'flagcube': self.scihdu['FLAGS'].data, 
                                    'z': 0.0,
                                   'yunit': self.scihdr['BUNIT'], 'skycube': self.skyhdu[0].data}
            

            #replace the bad pixels (flags >0) with NaNs
            # self.skyhdu[0].data[self.skyhdu['FLAGS'].data > 0] = np.nan  
            # self.skyhdu['UNCERT'].data[self.skyhdu['FLAGS'].data > 0] = np.nan 

            
        else:
            self.skyhdu = None
            # self.plot_spectrum(self.scihdu[0].data, yunit = self.scihdr['BUNIT']) #plot the spectrum of the central region (x = [12, 22], y = [43, 53]; 10x10 box) for a quick look. 
            self.plot_spec_dict = {'datacube': self.scihdu[0].data, 'errcube': self.scihdu['UNCERT'].data, 
                                    'flagcube': self.scihdu['FLAGS'].data, 'z': 0.0,
                                   'yunit': self.scihdr['BUNIT'], 'skycube': None} 

        self.plot_spectrum(**self.plot_spec_dict)

        self.region_box_entry.delete(0, tk.END)
        self.region_box_entry.insert(tk.END, '13, 44, 23, 54')

    def save_cropped_data(self):
        """Save the cropped datacube within the good wavelength region, white-lighted image, and preliminary mask for ZAP"""

        #for science frame
        hdu_list = []
        indices = []
        if self.index > 0:
            hdu_list.append(self.scihdu)
            indices.append(self.index)
        if self.index2 > 0:
            
            #scale the sky spectrum if the exposure time between the science and sky doesn't match.
            scaling_factor = self.scihdr['XPOSURE'] / self.skyhdu[0].header['XPOSURE']
            self.skyhdu[0].data *= scaling_factor
            self.skyhdu[0].header['XPOSURE'] = self.scihdr['XPOSURE']
            self.insert_text(f"[INFO] Scale the sky cube by {scaling_factor:0.2f}")

            hdu_list.append(self.skyhdu)
            indices.append(self.index2)
        
        if len(hdu_list) > 0:
            for i, hdu in  enumerate(hdu_list):
                index = indices[i]

                self.insert_text(f"[INFO] Saving the cropped datacube for {self.prefix}_{index:05d}")

                #cropped datacube
                hdu.writeto(f'{self.output}/{self.prefix}_{index:05d}_{self.ctype}.fits', overwrite = True)

                #white-lighted image
                wlimg_index = np.where((self.obswave >= self.wlimg_wave_range[0]) & (self.obswave <= self.wlimg_wave_range[1]))[0]
                wlimg = np.sum(hdu[0].data[wlimg_index], axis = 0)
                hdr2d = collapse_header(self.scihdu[0].header)
                wlhdu = fits.PrimaryHDU(wlimg, header = hdr2d)
                mask = np.mean(self.scihdu['FLAGS'].data, axis = 0)
                mhdu = fits.ImageHDU(mask, header = hdr2d)
                hdulist = fits.HDUList([wlhdu, mhdu])
                check_dir(self.output)
                hdulist.writeto(f'{self.output}/{self.prefix}_{index:05d}_wlimg.fits', overwrite = True)

                #preliminary ZAP mask - only for the sky cube
                if i > 0:
                    edgemask = mask > 1
                    allmask = np.zeros_like(mask)
                    allmask[edgemask] = 1
                    maskhdu = fits.PrimaryHDU(allmask, header = hdr2d)
                    maskhdu.writeto(f'{self.output}/{self.prefix}_{index:05d}_zap_mask.fits', overwrite = True)

        if len(hdu_list) == 1:
            self.insert_text(f'[INFO] No off-field sky frame set for {self.prefix}_{self.index:05d}.' + \
                             f'Remember to check the white-lighteded image and create a ds9 region file masking out the source and save it as {self.prefix}_{self.index:05d}.reg in the output directory!')

    def update_zap_mask(self):
        """Update the ZAP mask for a given science frame"""

        mindex = self.mindex

        mhdu = fits.open(f'{self.output}/{self.prefix}_{self.index:05d}_wlimg.fits')
        edgemask = mhdu[1].data > 1

        region_path = f'{self.output}/{self.prefix}_{mindex:05d}.reg'
        if not os.path.exists(region_path):
            self.insert_text(f"[ERROR] Region file {self.prefix}_{mindex:05d}.reg not exists!")
        else:
            self.insert_text(f"[INFO] Reading region file {self.prefix}_{mindex:05d}.reg for {self.prefix}_{self.index:05d}")
            r = pyregion.open(region_path)
            region_mask  = r.get_mask(hdu=mhdu[0])
            allmask = np.zeros_like(mhdu[0].data)
            # allmask[edgemask | region_mask] = 1
            allmask[region_mask] = 2
            allmask[edgemask] = 1
            maskhdu = fits.PrimaryHDU(allmask, header = mhdu[0].header)

            #save ZAP mask
            filename = f'{self.output}/{self.prefix}_{self.mindex:05d}_zap_mask.fits'
            if os.path.exists(filename):
                confirm = messagebox.askyesno("File Exists",
                                              f'ZAP mask {self.prefix}_{self.mindex:05d}_zap_mask.fits already exists. Do you want to overwrite it?')
                if not confirm:
                    return
            self.insert_text(f"[INFO] Saving ZAP mask to {self.prefix}_{self.mindex:05d}_zap_mask.fits")
            maskhdu.writeto(filename, overwrite = True)

    
    def crop_cube(self, hdu):
        """Crop the HDUList to only consist of the region"""
        
        hdr = hdu[0].header
        wave = (np.arange(hdr['NAXIS3']) + 1 - hdr['CRPIX3']) * hdr['CD3_3'] + hdr['CRVAL3']
        wavegood_idx = np.where((wave >= hdr['WAVGOOD0']) & (wave <= hdr['WAVGOOD1']))[0]
        
        for h in hdu:
            h.data = h.data[wavegood_idx]
        hdu[0].header['CRVAL3'] = wave[wavegood_idx][0]
        return hdu


    def plot_spectrum(self, datacube, errcube, flagcube, z = 0, xrange = [12, 22], yrange = [43, 53], skycube = None, restore_limit = False,
                      yunit = 'electron', unzapped_skycube = False, show_lines = False):
        """
        plot the spectrum of a given region for
        
        Args:
            datacube (3D array)
            errcube (3D arr): the uncertainty cube of datacube
            z (float): the redshift of the object
            xrange (list): the start and ending of the x-axis from which the spec is extracted
            yrange (list): the start and ending of the x-axis from which the spec is extracted
            skycube (3D arr): the off-field sky cube, sky model cube, or the unzapped cube
            unzapped_skycube (bool): if True, the skycube should be the datacube pre ZAP subtraction
        """
        #restore the limit of the original plot, so that the plot would stay in the zoomed-in version
        if restore_limit:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()

        self.ax.clear()

        datacube = datacube.copy()
        errcube = errcube.copy()
        # datacube[flagcube > 0] = np.nan
        # errcube[flagcube > 0] = np.nan

        # spec = np.nanmean(datacube[:, yrange[0]:yrange[1], xrange[0]: xrange[1]], axis = (1,2))
        spec = np.mean(datacube[:, yrange[0]:yrange[1], xrange[0]: xrange[1]], axis = (1,2))
        err = np.sqrt(np.nansum(errcube[:, yrange[0]:yrange[1], xrange[0]: xrange[1]]**2, axis = (1,2))) / np.sum(np.isfinite(errcube[:, yrange[0]:yrange[1], xrange[0]: xrange[1]]), axis = (1,2))
        self.ax.step(self.obswave / (1+z), spec, color ='k', lw = 1, label = 'sci spec')
        self.ax.step(self.obswave / (1+z), err, color ='grey', lw = 1, label = 'sci err')
        self.ax.fill_between(self.obswave / (1+z), spec - err, spec + err, color = 'lightgrey', step = 'pre')
        
        if skycube is not None:
            skyspec = np.nanmean(skycube[:, yrange[0]:yrange[1], xrange[0]: xrange[1]], axis = (1,2))

            # rescale it to the similar level
            if unzapped_skycube:
                skyspec = skyspec / (np.median(skyspec)/ np.median(spec))
                label = 'pre-ZAP spec (rescaled)'
            else:
                label = 'sky spec'

            self.ax.step(self.obswave / (1+z), skyspec, color ='lightskyblue', lw = 1, label = label, alpha = 0.5)
        self.ax.legend()
        self.ax.set_title(self.scihdr['OFNAME'] + f'  - {self.objname}')
        # self.figure.tight_layout()
        if np.abs(z) < 1e-10:
            self.ax.set_xlabel('Obs. Wavelength [A]')
        else:
            self.ax.set_xlabel('Rest. Wavelength [A]')
        self.ax.set_ylabel('Flux [%s]'%yunit)

        if show_lines:
            line_tbl = self.line_tbl[(self.line_tbl['vacuum'] >= self.obswave[0] / (1+z)) & (self.line_tbl['vacuum'] <= self.obswave[-1]/ (1+z))]
            ymin, ymax = self.ax.get_ylim()
            ytext = 0.1*(ymax - ymin) + ymin
            for row in line_tbl:
                if np.abs(z) > 1e-10:                    
                    self.ax.axvline(row['vacuum'], color = 'lightpink', ymin = 0.15, ls = 'dotted', alpha = 0.5, zorder = 1)
                    self.ax.text(row['vacuum'], ytext, row['line'], ha ='center', va = 'top', rotation = 90, fontsize = 8, color = 'darkorchid')

        if restore_limit:
            self.ax.set_xlim(xlim)
            self.ax.set_ylim(ylim)
        else:
            self.toolbar.update()
            
        self.canvas.draw()


    def run_zap_precondition(self):
        """
        Initialize the ZAP setup, used for the connection with ZAP button
        """

        ########Reinitialize the zap dictionary if the multi_sky_seg or use_Ha_seg has been updated
        if len(self.zap) > 0:
            if (self.use_multi_skyseg.get() != self.zap['use_multi_skyseg']) or (self.use_Ha_seg.get() != self.zap['use_Ha_seg']):
                self.zap = {}

        ####initialize the zap dictionary if the self.zap is an empty dictionary
        if not self.zap:
            ### Only use one sky segment
            if not self.use_multi_skyseg.get():
                self._zap_skyseg = [] #use one segment, default setting of ZAP
                self._zap_cfwidth = 300 #cfwidth = 300, default setting of ZAP

            else:
                #Use the sky segment from the old version of MUSE. See Table 1 in Soto+16 for details 
                skyseg0 = [0, 5400, 5850, 6440, 6750, 7200, 7700, 8265, 8602, 8731, 9275, 10000] 
                #remove the extra sky segement falls outside of the spectral region in case ZAP runs into problem
                skyseg0 = np.array(skyseg0)
                idx_remove_low = np.where(skyseg0 < self.obswave[0])[0]
                idx_remove_up = np.where(skyseg0 > self.obswave[-1])[0]
                idx_remove = np.concatenate((idx_remove_low[1:], idx_remove_up[1:]))
                self._zap_skyseg = np.delete(skyseg0, idx_remove) #the first index is zero; need to keep
                self._zap_skyseg = self._zap_skyseg.tolist()
                self._zap_cfwidth = [300] * (len(self._zap_skyseg) - 1) #cfwidth = 300, default setting of ZAP

                if self.use_Ha_seg.get():
                    self.insert_text(f'[INFO] Add additional sky component near Halpha+NII complex, assuming z = {self.redshift:0.4f}.')
                    Ha_comp = [6540, 6595] #spectral range of the Halpha+NII complex. 
                    self._zap_skyseg, self._zap_cfwidth = reassign_skyseg(self._zap_skyseg, self._zap_cfwidth, 
                                                            (Ha_comp[0]*(1+self.redshift), Ha_comp[1]*(1+self.redshift)), 5)
            
            self.zap['skyseg'] = self._zap_skyseg.copy()
            if np.isscalar(self._zap_cfwidth):
                self.zap['cfwidth'] = self._zap_cfwidth
            else:
                self.zap['cfwidth'] = self._zap_cfwidth.copy()

            self.zap['use_multi_skyseg'] = self.use_multi_skyseg.get()
            self.zap['use_Ha_seg'] = self.use_Ha_seg.get()

        #TODO: add the line to indicate different sky segments

        #print the default sky segment                                                
        self.print_sky_seg(self.zap['skyseg'], self.zap['cfwidth'])
        self.skyseg_input_active = True

    def update_zap_skyseg(self, event):
        """
        Main program to run ZAP
        """
        if self.skyseg_input_active:
            # get the current content of the last line
            current_content = self.output_text.get(1.0, tk.END)
            lines = current_content.splitlines()[-2:] #the last three lines are usable

            #the reset the sky segment
            if 'r' in lines[-1]:
                self.zap['skyseg'] = self._zap_skyseg
                self.zap['cfwidth'] = self._zap_cfwidth
                self.print_sky_seg(self.zap['skyseg'], self.zap['cfwidth'])


            #running ZAP
            elif 'q' in lines[-1]:
                self.skyseg_input_active = False
                self.run_zap()

            elif 'u' in lines[-1]:
                skyseg = []
                cfwidth = []

                #multiple sky segment
                if ';' in lines[0]:
                    tuple_list = lines[0].split(';')[:-1]
                    for t in tuple_list:
                        t = re.sub('[()]','',t)
                        start, end, width = t.split(',')
                        #the first element
                        if len(skyseg) == 0:
                            skyseg.extend([int(start), int(end)])
                        else:
                            if skyseg[-1] != int(start):
                                self.insert_text('\n[ERROR] Invalid input!')
                                self.print_sky_seg(self.zap['skyseg'], self.zap['cfwidth'])
                                return "break"
                            skyseg.append(int(end))
                        cfwidth.append(int(width))

                #one sky segment
                else:
                    start, end, width = lines[0][1:-1].split(',')
                    skyseg = []
                    cfwidth = int(width)                    

                self.zap['skyseg'] = skyseg
                self.zap['cfwidth'] = cfwidth
                self.print_sky_seg(self.zap['skyseg'], self.zap['cfwidth'])
            
            else:
                self.insert_text('\n[ERROR] Invalid input!')
                # self.insert_text(f"[INSTRUCTION] Type 'q' if satisfied with the segment, type 'r' to reset the default segmengt, or edit the following skyseg and cfwidth directy and type 'u'!")
                # self.insert_text('skyseg = [' + ', '.join(map(str, self.zap['skyseg'])) + ']')
                # self.insert_text('cfwidth = [' + ', '.join(map(str, self.zap['cfwidth'])) + ']')
                # self.insert_text('-->: ', newline = False)
                self.print_sky_seg(self.zap['skyseg'], self.zap['cfwidth'])
                return "break"


            # Prevent the default behavior of the Return key
            return "break"
            
        
    
    def print_sky_seg(self, skyseg, cfwidth = 300, skip_title = False):
        """
        print the sky segment into the text widget
        Args:
            skyseg: If not using multiple sky segment for the PCA analysis, it is set as [] by default (same as the official MUSE ZAP);
                    If using multiple sky segment, it takes a list of the endpoints in each interval (including the start and the end wavelength)  - N elements          
            cfwidth: A scalar (for one sky segment) or a list with (N-1 elements representing the cfwidth in each segment)
        """
        nseg = len(skyseg)

        self.insert_text('\n ##########Sky segment for ZAP##########')
        self.insert_text(f"[INSTRUCTION] Type 'q' if satisfied with the segment, type 'r' to reset the default segmengt, or edit the following skyseg and cfwidth directy and type 'u'!\n")

        #one cfwidth all segments
        if np.isscalar(cfwidth):

            if nseg == 0:
                self.insert_text(f'Using 1 sky segment (start, end, cfwidth):\n')
                self.insert_text('(0, 10000, %i)'%cfwidth)

            else:
                self.insert_text(f'Using {nseg:d} sky segment (start, end, cfwidth):\n')
                for i in range(nseg - 1):
                    if i == nseg - 2:
                        self.insert_text(f'(%i, %i, %i); '%(skyseg[i], skyseg[i+1], cfwidth))
                    else:
                        self.insert_text(f'(%i, %i, %i); '%(skyseg[i], skyseg[i+1], cfwidth), newline = False)
                # self.insert_text(f'Using {nseg:d} sky segment!')


        #multiple sky segments
        else:
            if nseg - 1 != len(cfwidth):
                self.insert_text(f'Length mismatch! The length of cfwidth needs to be N-1 for N sky segments!')
                return

            self.insert_text(f'Using {nseg:d} sky segment (start, end, cfwidth):\n')
            for i in range(nseg - 1):
                if i == nseg - 2:
                    self.insert_text(f'(%i, %i, %i); '%(skyseg[i], skyseg[i+1], cfwidth[i]))
                else:
                    self.insert_text(f'(%i, %i, %i); '%(skyseg[i], skyseg[i+1], cfwidth[i]), newline = False)
            # self.insert_text('\n')
            # self.insert_text(f'Using {nseg:d} sky segment!')
        
        #print out the current sky seg
        # self.insert_text('skyseg = [' + ', '.join(map(str, self.zap['skyseg'])) + ']')
        # self.insert_text('cfwidth =    [' + ', '.join(map(str, self.zap['cfwidth'])) + ']')
        self.insert_text('-->: ', newline = False)


    def run_zap(self):
        """
        The main function running the ZAP
        """
        if self.index < 0:
            self.insert_text('\n [ERROR] Set the science frame number first!')
            return

        self.insert_text('\n' +f'[INFO] Running ZAP on {self.prefix}_{self.index:05d}_{self.ctype}.fits...')

        try:
            self.scihdu = fits.open(f'{self.output}/{self.prefix}_{self.index:05d}_{self.ctype}.fits')
            self.scihdr = self.scihdu[0].header
        except FileNotFoundError:
            self.insert_text(f'[ERROR] Cropped data cubes not exist! Load the raw DRP-reduced cubes and save the cropped cubes first!')
            return
        self.obswave = (np.arange(self.scihdr['NAXIS3']) + 1 - self.scihdr['CRPIX3']) * self.scihdr['CD3_3'] + self.scihdr['CRVAL3']
        

        #In-field sky:
        if self.index2 < 0:
            self.insert_text(f'[INFO] Use in-field sky based on {self.prefix}_{self.mindex:05d}_zap_mask.fits')
            self.skyhdu = None
            maskpath = f'{self.output}/{self.prefix}_{self.mindex:05d}_zap_mask.fits'
            zobj = zap.process(f'{self.output}/{self.prefix}_{self.index:05d}_{self.ctype}.fits',
                  mask = maskpath, interactive = True, 
                  cfwidthSP = self.zap['cfwidth'], cfwidthSVD = self.zap['cfwidth'], skyseg = self.zap['skyseg'], zlevel = 'median')

        #Off-field sky
        if self.index2 > 0:
            self.insert_text(f'[INFO] Use off-field sky {self.prefix}_{self.index2:05d}. Mask file -  {self.prefix}__{self.index2:05d}_zap_mask.fits')
            self.skyhdu = fits.open(f'{self.output}/{self.prefix}_{self.index2:05d}_{self.ctype}.fits')
            maskpath = f'{self.output}/{self.prefix}_{self.index2:05d}_zap_mask.fits'
            extSVD = zap.SVDoutput(f'{self.output}/{self.prefix}_{self.index2:05d}_{self.ctype}.fits',
                       mask = maskpath,
                        skyseg = self.zap['skyseg'], zlevel = 'median')
            zobj = zap.process(f'{self.output}/{self.prefix}_{self.index:05d}_{self.ctype}.fits', extSVD=extSVD, interactive = True,
                  cfwidthSP = self.zap['cfwidth'], 
                   skyseg = self.zap['skyseg'])

        nsig = 3
        # skycube = zobj.cube - zobj.cleancube
        # ### if the object pixels are either over-subtracted or under-subtracted near Halpha (>3sigma), replace the sky pixel value with the median
        # mask = fits.getdata(maskpath) #get the mask to avoid the edge pixel. Should be similar if using the off-field sky
        # use = np.abs(mask - 1) > 1e-6 #mask = 1 for edge mask
        # skycube_clipped = sigma_clip(skycube[:,use], sigma = nsig, axis = 1)
        # # skycube_clipped = sigma_clip(skycube[:,use], axis = 1)

        # median_cube = np.zeros_like(skycube) #3D median cube
        # std_cube = np.zeros_like(skycube) #3D STD cube
        # npix = np.sum(use)
        # median_cube[:, use] = np.repeat(np.nanmedian(skycube_clipped, axis = 1).data, npix).reshape(len(skycube), npix)
        # std_cube[:, use] = np.repeat(np.nanstd(skycube_clipped, axis = 1).data, npix).reshape(len(skycube), npix)
        # #find and replace the 3sigma outlier
        # bpm = np.where((np.abs(skycube - median_cube) > nsig*std_cube) & (median_cube > 0))
        # # bpm = np.where(np.abs(skycube - median_cube) > nsig*std_cube)
        # skycube[bpm] = median_cube[bpm]
        # cleancube = zobj.cube - skycube

        skycube0 = zobj.cube - zobj.cleancube
        skycube = skycube0.copy()

        mask = fits.getdata(maskpath) #get the mask to avoid the edge pixel. Should be similar if using the off-field sky
        use = np.abs(mask - 1) > 1e-6 #mask = 1 for edge mask
        skycube[:,~use] = np.nan
        skycube_clipped = sigma_clip(skycube, sigma = nsig, axis = (1,2))
        median_sky = np.ma.median(skycube_clipped, axis = (1,2)).data
        median_cube = median_sky[:, np.newaxis, np.newaxis] * np.ones((1, np.shape(skycube)[1], np.shape(skycube)[2]))
        skycube[skycube_clipped.mask] = median_cube[skycube_clipped.mask]
        skycube[:, ~use] = skycube0[:, ~use]
        cleancube = zobj.cube - skycube

        # save the output cube
        self.cleanhdu = self.scihdu.copy()
        self.cleanhdu.append(self.scihdu[0])
        self.cleanhdu[-1].name = 'UNZAPPED'
        self.cleanhdu[0].data = cleancube
        # bad = np.where(self.cleanhdu['FLAGS'].data >=8)
        # self.cleanhdu[0].data[bad] = 0
        skyhdu = fits.ImageHDU(data=skycube, header=self.scihdu[0].header)
        skyhdu.name = 'SKYMODEL_ZAP'
        self.cleanhdu.append(skyhdu)
        self.cleanhdu.writeto(f'{self.output}/{self.prefix}_{self.index:05d}_zap_{self.ctype}.fits', overwrite = True)

        #Need to run the DAR correction for the unrectified cube
        if self.ctype == 'icube':

            self.insert_text(f'[INFO] Correcting for DAR!')
            image_size = self.cleanhdu[0].data.shape

            # get wavelengths
            w0 = self.scihdr['CRVAL3']
            dw = self.scihdr['CD3_3']
            waves = w0 + np.arange(image_size[0]) * dw
            wgoo0 = self.scihdr['WAVGOOD0']
            wgoo1 = self.scihdr['WAVGOOD1']
            if DAR_ref_wave is None: 
                wref = self.scihdr['WAVMID']
            else:
                wref = DAR_ref_wave
            self.insert_text("Ref WL = %.1f, good WL range = (%.1f - %.1f)" %
                            (wref, wgoo0, wgoo1))
            
            # spatial scales in arcsec/item
            y_scale = self.scihdr['PXSCL'] * 3600.
            x_scale = self.scihdr['SLSCL'] * 3600.

            # padding depends on grating
            if 'H' in self.scihdr['RGRATNAM']:
                padding_as = 2.0
            elif 'M' in self.scihdr['RGRATNAM']:
                padding_as = 3.0
            else:
                padding_as = 4.0

            padding_x = int(padding_as / x_scale)
            padding_y = int(padding_as / y_scale)

            # update WCS
            crpix1 = self.scihdr['CRPIX1']
            crpix2 = self.scihdr['CRPIX2']
            self.scihdr['CRPIX1'] = crpix1 + float(padding_x)
            self.scihdr['CRPIX2'] = crpix2 + float(padding_y)

            # airmass
            airmass = self.scihdr['AIRMASS']
            self.insert_text("Airmass: %.3f" % airmass)

            # IFU orientation
            ifu_pa = self.scihdr['IFUPA']

            # Parallactic angle
            parallactic_angle = self.scihdr['PARANG']

            # Projection angle in radians
            projection_angle_deg = ifu_pa - parallactic_angle
            projection_angle = math.radians(projection_angle_deg)

            self.insert_text("DAR Angles: ifu_pa, parang, projang (deg): "
                         "%.2f, %.2f, %.2f" % (ifu_pa, parallactic_angle,
                                               projection_angle_deg))
            
             # dispersion over goo wl range in arcsec
            dispersion_max_as = atm_disper(wgoo1, wgoo0, airmass)

            # projected onto IFU
            xdmax_as = dispersion_max_as * math.sin(projection_angle)
            ydmax_as = dispersion_max_as * math.cos(projection_angle)
            self.insert_text("DAR over GOOD WL range: total, x, y (asec): "
                            "%.2f, %.2f, %.2f" % (dispersion_max_as, xdmax_as,
                                                ydmax_as))
            # now in pixels
            xdmax_px = xdmax_as / x_scale
            ydmax_px = ydmax_as / y_scale
            dmax_px = math.sqrt(xdmax_px**2 + ydmax_px**2)
            self.insert_text("DAR over GOOD WL range: total, x, y (pix): "
                            "%.2f, %.2f, %.2f" % (dmax_px, xdmax_px, ydmax_px))
            
            # prepare output cubes
            output_image = np.zeros((image_size[0], image_size[1]+2*padding_y,
                                    image_size[2]+2*padding_x), dtype=np.float64)
            output_stddev = output_image.copy()
            output_mask = np.zeros((image_size[0], image_size[1]+2*padding_y,
                                    image_size[2]+2*padding_x), dtype=np.uint8)
            output_flags = np.zeros((image_size[0], image_size[1] + 2*padding_y,
                                    image_size[2] + 2 * padding_x), dtype=np.uint8)
            output_skymodel_zap = np.zeros((image_size[0], image_size[1]+2*padding_y,
                                    image_size[2]+2*padding_x), dtype=np.float64)
            output_unzapped = np.zeros((image_size[0], image_size[1]+2*padding_y,
                                    image_size[2]+2*padding_x), dtype=np.float64)
            
            if 'NOSKYSUB' in [hdu.name for hdu in self.cleanhdu]:
                output_noskysub = np.zeros((image_size[0],
                                            image_size[1] + 2*padding_y,
                                            image_size[2] + 2 * padding_x),
                                        dtype=np.float64)
            else:
                output_noskysub = None

            # DAR padded pixel flag
            output_flags += 128

            output_image[:, padding_y:(padding_y+image_size[1]),
                        padding_x:(padding_x+image_size[2])] = self.cleanhdu[0].data

            output_stddev[:, padding_y:(padding_y+image_size[1]),
                        padding_x:(padding_x+image_size[2])] = self.cleanhdu['UNCERT'].data

            output_mask[:, padding_y:(padding_y+image_size[1]),
                        padding_x:(padding_x+image_size[2])] = self.cleanhdu['MASK'].data

            output_flags[:, padding_y:(padding_y+image_size[1]),
                        padding_x:(padding_x+image_size[2])] = self.cleanhdu['FLAGS'].data
            
            output_unzapped[:, padding_y:(padding_y+image_size[1]),
                        padding_x:(padding_x+image_size[2])] = self.cleanhdu['UNZAPPED'].data
            
            output_skymodel_zap[:, padding_y:(padding_y+image_size[1]),
                        padding_x:(padding_x+image_size[2])] = self.cleanhdu['SKYMODEL_ZAP'].data
            
            if output_noskysub is not None:
                output_noskysub[:, padding_y:(padding_y + image_size[1]),
                                padding_x:(padding_x + image_size[2])] = self.cleanhdu['NOSKYSUB'].data
                
            self.insert_text(f"Image cube DAR order = {DAR_shift_order}")
            self.insert_text(f"Std. Dev. cube DAR order = {DAR_shift_order}")
            self.insert_text(f"Mask cube DAR order = 1 (constant)")
            self.insert_text(f"Flag cube DAR order = 1 (constant)")
            # Perform correction
            for j, wl in enumerate(waves):
                dispersion_correction = atm_disper(wref, wl, airmass)
                x_shift = dispersion_correction * \
                    math.sin(projection_angle) / x_scale
                y_shift = dispersion_correction * \
                    math.cos(projection_angle) / y_scale
                output_image[j, :, :] = shift(output_image[j, :, :],
                                            (y_shift, x_shift), order = DAR_shift_order)
                output_stddev[j, :, :] = shift(output_stddev[j, :, :],
                                            (y_shift, x_shift))
                output_mask[j, :, :] = np.ceil(shift(output_mask[j, :, :], (y_shift,
                                                                x_shift), order=1, mode = 'constant', cval=128))
                output_flags[j, :, :] = np.ceil(shift(output_flags[j, :, :], (y_shift,
                                                                  x_shift), order=1, mode = 'constant', cval=128))
                output_unzapped[j, :, :] = shift(output_unzapped[j, :, :],
                                            (y_shift, x_shift), order = DAR_shift_order)
                output_skymodel_zap[j, :, :] = shift(output_skymodel_zap[j, :, :],
                                            (y_shift, x_shift), order = DAR_shift_order)
                if output_noskysub is not None:
                    output_noskysub[j, :, :] = shift(output_noskysub[j, :, :],
                                                    (y_shift, x_shift), order = DAR_shift_order)
                    
            #update the data
            self.cleanhdu[0].data = output_image
            self.cleanhdu['UNCERT'].data = output_stddev
            self.cleanhdu['MASK'].data = output_mask
            self.cleanhdu['FLAGS'].data = output_flags
            self.cleanhdu['UNZAPPED'].data = output_unzapped
            self.cleanhdu['SKYMODEL_ZAP'].data = output_skymodel_zap

            if output_noskysub is not None:
                self.cleanhdu['NOSKYSUB'].data = output_noskysub

            # update header
            self.scihdr['DARCOR'] = (True, 'DAR corrected?')
            self.scihdr['DARANG'] = (projection_angle_deg,
                                                        'DAR projection angle '
                                                        '(deg)')
            self.scihdr['DARPADX'] = (padding_x,
                                                        'DAR X padding (pix)')
            self.scihdr['DARPADY'] = (padding_y,
                                                        'DAR Y padding (pix)')
            self.scihdr['DAREFWL'] = (wref, 'DAR reference wl (Ang)')

            self.cleanhdu[0].header = self.scihdr

            self.cleanhdu.writeto(f'{self.output}/{self.prefix}_{self.index:05d}_zap_icubed.fits', overwrite = True)

                    

        #no available std
        if not hasattr(self, 'std'):
            if self.ctype == 'icubes':
                self.insert_text(f'[INFO] The input datacube has been flux calibrated! Skip the flux calibration.')
                self.insert_text(f'[WARNING] No telluric model from standard star available. Skip the telluric correction.')
                self.cleanhdu_flux = self.cleanhdu.copy()
            else:
                self.insert_text(f'[ERROR] No sensitive function and telluric model available! Cannot apply the flux calibration!')
                return

        else:
            use = (self.std['wave']>= self.obswave[0]) & (self.std['wave'] <= self.obswave[-1])
            tellmodel = self.std['tellmodel'][use]**(self.cleanhdu[0].header['AIRMASS']) #convert the telluric model at AM=1.0 to the real AM
            
            #flux calibration and telluric correction
            if self.ctype == 'icubes':
                self.insert_text(f'[INFO] The input datacube has been flux calibrated! Skip the flux calibration. Running telluric correction...')
                mscal = 1. / tellmodel
            else:
                self.insert_text(f'[INFO] Running the flux calibration and telluric correction...') 
                mscal = self.std['invsens_model'] * 1e16 / self.cleanhdu[0].header['XPOSURE'] #normalize by the exposure time
                mscal, self.cleanhdu[0].header = kcwi_correct_extin(mscal, self.cleanhdu[0].header)
                mscal = mscal[use] / tellmodel #include the telluric correction

            #reshpae the 1D mscal to 3D 
            mscal = mscal[:, np.newaxis, np.newaxis]
            #also need to crop the mscal as the wavelength of the std goes a little beyond WAVGOOD
            use = ()
            self.cleanhdu_flux = self.cleanhdu.copy() #flux-calibrated cube
            self.cleanhdu_flux[0].data *= mscal
            self.cleanhdu_flux['UNCERT'].data *= mscal
            self.cleanhdu_flux['UNZAPPED'].data *= mscal
            self.cleanhdu_flux['SKYMODEL_ZAP'].data *= mscal
            self.cleanhdu_flux[0].header['BUNIT'] = '1e-16 erg / (Angstrom cm2 s)'
            self.cleanhdu_flux[0].header['STDCOR'] = (True, 'std corrected?')
            self.cleanhdu_flux[0].header['MSFILE'] = ('%s_invsens_updated.fits'%self.std['frame'], 'Master std filename')
            self.cleanhdu_flux[0].header['MSIMNO'] = (int(self.std['frame'][-5:]), 'master std image number')

            try:
                self.cleanhdu_flux['NOSKYSUB'].data *= mscal
            except KeyError:
                pass

            #write the new cubes
            self.cleanhdu_flux.writeto(f'{self.output}/{self.prefix}_{self.index:05d}_zap_icubes.fits', overwrite = True)
                

        #save the white-lighted image of the clean cube
        wlimg_index = np.where((self.obswave >= self.wlimg_wave_range[0]) & (self.obswave <= self.wlimg_wave_range[1]))[0]
        wlimg = np.sum(self.cleanhdu_flux[0].data, axis = 0)
        hdr2d = collapse_header(self.cleanhdu_flux[0].header)
        wlhdu = fits.PrimaryHDU(wlimg, header = hdr2d)
        wlhdu.header['WAVWLIMG0'] = self.wlimg_wave_range[0]
        wlhdu.header['WAVWLIMG1'] = self.wlimg_wave_range[1]
        mask = np.mean(self.cleanhdu_flux['FLAGS'].data, axis = 0)
        mhdu = fits.ImageHDU(mask, header = hdr2d)
        hdulist = fits.HDUList([wlhdu, mhdu])
        hdulist.writeto(f'{self.output}/{self.prefix}_{self.index:05d}_zapclean_wlimg.fits', overwrite = True)

        #TODO make the errcube_combined as the official cleanhdu_flux['UNCERT'].data
        # mask = fits.getdata(maskpath) #get the mask to avoid the edge pixel. Should be similar if using the off-field sky
        # use = np.abs(mask - 1) > 1e-6 #mask = 1 for edge mask
        # skystd = np.std(self.cleanhdu_flux['SKYMODEL_ZAP'].data[:, use], axis = 1) #std in sky model of each wavelength slice
        # errcube = self.cleanhdu_flux['UNCERT'].data.copy()
        # errcube_combined = np.sqrt(errcube**2 + skystd[:, np.newaxis, np.newaxis]**2)


        #plot the spectrum
        self.plot_spec_dict = {'datacube': self.cleanhdu_flux[0].data, 'errcube': self.cleanhdu_flux['UNCERT'].data,
                                'flagcube': self.cleanhdu_flux['FLAGS'].data, 
                                'z': self.redshift, 'yunit': self.cleanhdu_flux[0].header['BUNIT'], 
                                'skycube': self.cleanhdu_flux['UNZAPPED'].data, 'unzapped_skycube': True,
                                'restore_limit': False, 'show_lines': True}

        self.insert_text(f'[INFO] ZAP Done for {self.prefix}_{self.index:05d}!')
        self.plot_spectrum(**self.plot_spec_dict)
        inst = ("\n[INSTRUCTION] Press 'n' to turn off the pre-ZAP spec;"
                "press 's' to display the pre-ZAP spec;"
                "press 'r' to plot the spec in the rest frame;"
                "press 'o' to plot the spec in the observed frame."
                "If you want to plot the spectrum of the other box region, please update the region index in the box above as (x1, y1, x2, y2) [lower left + upper right].")

        self.insert_text(inst)

        


    def load_invsens(self, type):
        """
        Load the  _invsens.fits of a standard star.
        """
        if type == 'DRP':
            filename = filedialog.askopenfilename(initialdir=self.base) #REAL code

        else:
            filename = filedialog.askopenfilename(initialdir=self.output) #REAL code
        # filename = '/scr/zzhuang/keck_obs/kcwi/2023sep23/red/redux/kr230923_00174_invsens.fits' #for test purpose only
        self.std_entry.delete(0, tk.END)
        self.std_entry.insert(tk.END, filename)

        #load the invsens file 
        hdu = fits.open(self.std_entry.get())
        self.std = {}

        hdr = hdu[0].header
        wvl = (np.arange(hdr['NAXIS1']) + 1 - hdr['CRPIX1']) * hdr['CDELT1'] + hdr['CRVAL1']
        #Add 3A on each side to give the spline-fit some buffer 
        good = np.where((wvl >= hdr['WAVGOOD0'] - 3) & (wvl<= hdr['WAVGOOD1'] + 3))[0]
        self.std['invsens_hdr'] = hdr

        #raw invsens from the DRP
        if type == 'DRP':
            self.std['wave'] = wvl[good]
            #the first row is the raw ratio of real flux and count to be fitted for; second row the best-fit invsens; third row is the raw electron/s 
            self.std['invsens_data'] = hdu[0].data[0, good]
            self.std['invsens_model_drp'] = hdu[0].data[1, good]
            self.std['counts'] = hdu[0].data[2, good]
            self.std['name'] = hdr['OBJECT'].lower()
            self.std['name'] = re.sub('[\ _]', '', self.std['name'])
            self.std['frame'] = re.sub('_invsens.fits', '', os.path.basename(self.std_entry.get()))
            flag = np.full(len(self.std['wave']), 1, dtype = int)
            self.std['flag'] = self.mask_skyline_region(self.std['wave'], flag)

            #sens func and telluric model to be updated
            self.std['invsens_model'] = None
            self.std['tellmodel'] = None

        #the updated version, so no need to crop the data
        elif type == 'updated':
            self.std['wave'] = wvl
            #the first row is the raw ratio of real flux and count to be fitted for; second row the best-fit invsens; third row is the raw electron/s 
            self.std['invsens_data'] = hdu[0].data[0]
            self.std['invsens_model_drp'] = hdu[0].data[1]
            self.std['counts'] = hdu[0].data[2]
            self.std['name'] = hdr['OBJECT'].lower()
            self.std['invsens_model'] = hdu[0].data[3]
            self.std['flag'] = hdu[0].data[4]
            if len(hdu[0].data) > 5:
                self.std['tellmodel'] = hdu[0].data[5]
            else:
                self.std['tellmodel'] = None
            self.std['frame'] = re.sub('_invsens_updated.fits', '', os.path.basename(self.std_entry.get()))



        #setup the B-Spline fit parameters
        self.std['bspline_bkpt'] = 150 #breakpoints
        self.std['bspline_polyorder'] = 3 #polynomial order between interval
        
        self.insert_text(f"[INFO] Loading the {self.std_entry.get()}")


        #load the flux-calibrated spec for comparison
        try:
            std = fits.getdata('{0}/{1}.fits'.format(self.stddir, self.std['name']))
        except ValueError:
            self.insert_text(f"Cannot find the std spec for {self.std_name} in {self.stddir}")

        use = np.where((std['WAVELENGTH'] >= hdr['WAVGOOD0'] - 5) & (std['WAVELENGTH'] <= hdr['WAVGOOD1'] + 5))[0] #only load the spectrum in the good wavelength region
        self.std['spec_calib'] = np.column_stack((std['WAVELENGTH'][use], std['FLUX'][use]))

        # self.std['flag'] = self.mask_skyline_region(self.std['wave']) #flag indicated if a region is masked out
        # self.std['use_ind'] = np.array([], dtype = int) #a list of indices for the single points used for fitting
        #Initialize the region for later interactive selection

        self.std['region_start'] = None
        self.std['region_end'] = None

        #set the focuse to the canvas page
        self.canvas.get_tk_widget().focus_set()

        self.plot_std()
        inst = ("[INSTRUCTION] Press 'i' twice on each side of the region you want to include in the fit;"
                "press 'e' twice on each side of the region you want to exclude from the fit;"
                "press 'a' on the single data point you want to include in the fit;"
                "and press 'd' on the data point you want to exclude from the fit. "
                "First, to fit the sensitive function, press 'f'. Second, 'press' t to fit the telluric model [this may take a while]")
        self.insert_text(inst)

    def save_updated_invsens(self):

        """
        Save the updated invsens file to the output directory
        """
        
        ##########save the updated invsens to the output directory
        newhdr = self.std['invsens_hdr'].copy()
        newhdr['NAXIS1'] = len(self.std['wave'])
        newhdr['CRVAL1'] = self.std['wave'][0]
        #still similar to the format of the original invsens.fits; add the updated invsens model and flag to the last two rows
        newdata = np.vstack((self.std['invsens_data'], self.std['invsens_model_drp'], self.std['counts'],
                             self.std['invsens_model'], self.std['flag']))

        #add the tellmodel @ AM =1 to the file
        if self.std['tellmodel'] is not None:
            newdata = np.vstack((newdata, self.std['tellmodel']))

        newhdu = fits.PrimaryHDU(newdata, header = newhdr)

        frame = self.std['frame']
        filename = f'{self.output}/{frame}_invsens_updated.fits'
        if os.path.exists(filename):
            confirm = messagebox.askyesno("File Exists",
                                            f'Updated invsens file {frame}_invsens_updated.fits already exists. Do you want to overwrite it?')
            if not confirm:
                return
        self.insert_text(f"[INFO] Saving the updated invsens file to {frame}_invsens_updated.fits")
        check_dir(os.path.dirname(filename))
        newhdu.writeto(filename, overwrite = True)
        

    def plot_std(self, restore_limit = False):
        """
        Plot the best-fit STD 
        """
        #restore the limit of the original plot, so that the plot would stay in the zoomed-in version
        if restore_limit:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
        
         #plot the flux-calibrated spec
        self.ax.clear()
        self.ax.step(self.std['wave'], self.std['counts'] * self.std['invsens_model_drp'], 
                        color = 'r', lw =1, label = 'best-fit model from DRP', where = 'mid', alpha = 0.5) #raw count x invsens = flux-calibrated spec
        self.ax.plot(self.std['spec_calib'][:,0], self.std['spec_calib'][:,1], color = 'k', label = 'flux-calibrated spec of std')


        use_region = np.where(self.std['flag'] == 1)[0]
        use_point = np.where(self.std['flag'] == 2)[0]

        #plot the updated flux-calibrated model
        if self.std['invsens_model'] is not None:
            self.ax.step(self.std['wave'], self.std['counts'] * self.std['invsens_model'], color = 'cyan', 
                         lw =1, label = 'best-fit new model', where = 'mid') #raw count x invsens = flux-calibrated spec
            self.ax.plot(self.std['wave'][use_region], (self.std['counts'] * self.std['invsens_model'])[use_region], 'x', 
                        color = 'lightgreen', ms = 5, label = 'selected for fitting') #selected regions 
            if len(use_point) > 0:
                self.ax.plot(self.std['wave'][use_point], (self.std['counts'] * self.std['invsens_model'])[use_point], 'o', 
                            color = 'darkgreen', ms = 10, label = 'selected for fitting') #selected pixels 

             #plot the telluric-corrected model
            if self.std['tellmodel'] is not None:
                telluric = self.std['tellmodel']**(self.std['invsens_hdr']['AIRMASS']) #convert the model at AM=1.0 to the real AM
                self.ax.step(self.std['wave'], self.std['counts'] * self.std['invsens_model'] / telluric, color = 'royalblue',
                            where = 'mid', lw = 1, label = 'telluric-corrected std spec')
        
        #plot the DRP-reduced flux-calibrated model
        else:
            self.ax.plot(self.std['wave'][use_region], (self.std['counts'] * self.std['invsens_model_drp'])[use_region], 'x', 
                    color = 'lightgreen', ms = 5, label = 'selected for fitting') #selected regions 
            if len(use_point) > 0:
                self.ax.plot(self.std['wave'][use_point], (self.std['counts'] * self.std['invsens_model_drp'])[use_point], 'o', 
                            color = 'darkgreen', ms = 10, label = 'selected for fitting') #selected pixels 


        self.ax.set_title('%s - %s'%(self.std['frame'], self.std['name']))
        self.ax.set_xlabel('Obs. Wavelength [A]')
        self.ax.set_ylabel('flam [erg/s/cm2/A]')
        self.ax.legend()

        if restore_limit:
            self.ax.set_xlim(xlim)
            self.ax.set_ylim(ylim)
        else:
            # Capture the current view as the new home view
            self.toolbar.update()

        self.canvas.draw()
        

    def mask_skyline_region(self, wave, flag):
        """
        Mask out the region with dense sky line and telluric absorption from flux calibration fitting
        """
        regions = [[6274, 6302],[6864.00, 6950.00], [7160.00, 7385.00], [7590, 7691],[8102, 8375],
                   [8943, 9225], [9300, 9400],
                   [6514.00, 6600.00] #Halpha region
                   ]

        mask = False
        for r in regions:
            mask |= (wave >= r[0]) & (wave <= r[1])

        flag[mask] = 0

        return flag
        # [6864, 6935], [7164, 7345], [7591, 7694], [8131]]

    def run_telluric(self): 

        """
        Run the telluric correction on the standard star
        """

        if not hasattr(self, 'std'):
            self.insert_text(f"[ERROR] The updated invsens file not exists! Run the sensfunc fit first with the key'f' first!")
            return
        
        ################# create the "fake" Pypeit output for telluric correction ############
        frame = self.std['frame']
        self.insert_text(f"[INFO] Running telluric correction on {frame}...")

        newhdr = self.std['invsens_hdr'].copy()
        newhdr['NAXIS1'] = len(self.std['wave'])
        newhdr['CRVAL1'] = self.std['wave'][0]
        std_flux = self.std['counts'] * self.std['invsens_model']  / 1e-17 #pypeit takes the unit as 1e-17 erg/s/cm2/A
        SNR = 1000 #fake SNR of STD (required by Pypeit input)
        std_ivar = 1. / (std_flux / SNR)**2
        col1 = fits.Column(name = 'wave', format = '1D', array = self.std['wave'])
        col2 = fits.Column(name = 'flux', format = '1D', array = std_flux )
        col3 = fits.Column(name = 'ivar', format = '1D', array = std_ivar)

        #mask used for telluric correction, bad pixels with mask = 0
        mask = np.full(len(std_flux), 1, dtype = int)
        #the std spec in both DRP and pypeit seems to have some problems here; mask it out for g19b2b, should check for other stds
        # mask[(self.std['wave'] >=6310) & self.std['wave'] <= 6380] = 0 
        col4 = fits.Column(name = 'mask', format = '1K', array = mask)
        hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4])
        # hdu.header['PYP_SPEC'] = 'keck_lris_red'
        hdulist = fits.HDUList([fits.PrimaryHDU(header = newhdr), hdu])
        hdulist[1].header['PYP_SPEC'] = 'keck_kcrm'

        #need to update the header for pypeit input
        keys_1 = { 'DMODCLS': 'OneSpec ', 'DMODVER': '1.0.2   ', 'FLUXED': True, 
        'CHECKSUM': 'CG5BCD59CD5ACD59', 'DATASUM': '60558086'
       }
        keys_2 = {'PYP_SPEC': 'keck_kcrm', 'PYPELINE': 'SlicerIFU',  'TARGET':newhdr['OBJECT'],
        'DISPNAME': self.std['invsens_hdr']['RGRATNAM'], 'decker': self.std['invsens_hdr']['IFUNAM'],  
        'binning': self.std['invsens_hdr']['BINNING'], 'FILENAME': '%s.fits'%frame,
        'NCOMP': 1,
        # 'AIRMASS':1.0 #modify the airmass to one because the std has been extinction corrected by DRP
        }

        for k in keys_1.keys():
            hdulist[1].header[k] = keys_1[k]
        for k in keys_2.keys():
            hdulist[0].header[k] = keys_2[k]
        std_spec1d_path = f'{self.output}/{frame}_standard_spec1d.fits'
        hdulist.writeto(f'{std_spec1d_path}', overwrite = True)

        #Generate the pypeit input file
        coord = SkyCoord(ra = self.std['invsens_hdr']['RA'], dec = self.std['invsens_hdr']['DEC'], unit = (u.hourangle, u.deg))
        # lines = ['[telluric]',
        #          f'\ttelgridfile = {telgridfile}',
        #          '\tobjmodel = star',
        #          f'\tstar_ra = {coord.ra.value:0.5f}',
        #          f'\tstar_dec = {coord.dec.value:0.5f}',
        #          '\tteltype = grid']

        # tellfile_out = (f'{self.output}/{frame}_standard.tell')
        # with open(tellfile_out, 'w') as file:
        #     file.writelines(lines)

        # Run the pypeit_tellfit command for telluric correction
        # command = f'pypeit_tellfit -t {tellfile_out} {std_spec1d_path}'
        # process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
        # # self.read_tellfit_output()
        # while process.poll() is None: # check whether process is still running
        #     msg = process.stdout.readline().strip() # read a line from the process output
        #     if msg:
        #         self.redirect_text.write(msg)
        # print('finished')
        # command = ['pypeit_tellfit', '-t', f'{tellfile_out}', f'{std_spec1d_path}']
        # # run_command(command, self.output_text)
        # process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,)
        # process.wait()
        telluric_correct(infile_path = std_spec1d_path, star_ra = coord.ra.value, star_dec = coord.dec.value)

        # #the output files would be in the current directory. Need to move them to the output directory
        # files = ['telluric.par',  f'{frame}_standard_spec1d_tellcorr.fits', f'{frame}_standard_spec1d_tellmodel.fits']
        # for f in files:
        #     os.rename(f, f'{self.output}/{f}')

        #load the telluric model
        tellfile = fits.getdata(re.sub('.fits', '_tellcorr.fits', std_spec1d_path))
        tellmodel = tellfile['telluric'] #best-fit telluric model at the airmass of the standard star
        tellmodel = tellmodel**(1. / newhdr['AIRMASS']) #correct the telluric spectrum to be the one at airmass of 1
        self.std['tellmodel'] = tellmodel
        self.insert_text(f"[INFO] Telluric correction on {frame} DONE")


    def read_tellfit_output(self):
        if self.tellfit_process.poll() is None:  # Process is still running
            output_line = self.tellfit_process.stdout.readline()
            if output_line:
                self.redirect_text.write(output_line)
        self.root.after(100, self.read_tellfit_output)  # Check for new output every 100 milliseconds


    def on_click(self, event):
        """
        Mouse button used for general plotting analysis
        """
        # Right click to resume the home window
        if event.button == 3:  
            self.toolbar.home()

    def on_type(self, event):

        """
        Key interaction used for std analysis
        """
        # exclude/include the regions used for standard star analysis
        if event.key == 'e' or event.key == 'i':

            #select the starting point
            if self.std['region_start'] is None:
                self.std['region_start'] = event.xdata

            #select the ending point
            else:
                self.std['region_end'] = event.xdata   

                #flag the region to True to be masked out from fitting
                if event.key == 'e':
                    self.std['flag'][(self.std['wave']>= self.std['region_start']) & (self.std['wave']<= self.std['region_end'])] = 0
                
                #flag the region to False to be included in the fitting
                else:
                    self.std['flag'][(self.std['wave']>= self.std['region_start']) & (self.std['wave']<= self.std['region_end'])] = 1
                    
                self.std['region_start'] = None #reset the starting point for the next input
                self.plot_std(restore_limit = True) #update the std plot

        # add single one continuum data point closest to the mouse location; should be used in the regions with dense sky features
        if event.key == 'a':
            #find the index of the point closest to the mouse location
            # idx = np.argmin((self.std['wave'] - event.xdata)**2 + (self.std['counts'] * self.std['invsens_model_drp'] - event.ydata)**2) 
            idx = np.argmin(np.abs(self.std['wave'] - event.xdata))
            self.std['flag'][idx] = 2
            # print(self.std['use_ind'])
            self.plot_std(restore_limit = True)
        
        # delete single one continuum data point closest to the mouse location; should be used in the regions with dense sky features
        if event.key == 'd':
            if np.sum(self.std['flag'] ==2) > 0:
                # idx = np.argmin((self.std['wave'][self.std['use_ind']] - event.xdata)**2 + ((self.std['counts'] * self.std['invsens_model_drp'])[self.std['use_ind']] - event.ydata)**2) 
                ind_point = np.where(self.std['flag'] ==2 )[0]
                idx = np.argmin(np.abs(self.std['wave'][ind_point] - event.xdata))
                self.std['flag'][ind_point[idx]] = 0
                # self.std['use_ind'] = np.delete(self.std['use_ind'], idx)
                self.plot_std(restore_limit = True)

        #running the invsens fitting
        if event.key == 'f':
            use = self.std['flag'] > 0
            # use[self.std['use_ind']] = True
            # print(len(use), self.std['wave'])
            self.std['invsens_model'] = self.fit_bspline(self.std['wave'], self.std['invsens_data'], self.std['bspline_bkpt'], self.std['bspline_polyorder'], use)
            
            self.plot_std(restore_limit = True)

        #running the telluric correction
        if event.key == 't':
            self.run_telluric()
            self.plot_std()

        #plot the spectrum into the observed wavelength
        if event.key == 'o':
            self.plot_spec_dict['z'] = 0.0
            self.plot_spec_dict['restore_limit'] = False
            self.plot_spectrum(**self.plot_spec_dict)

        #plot the spectrum into the rest wavelength
        if event.key == 'r':
            self.plot_spec_dict['z'] = self.redshift
            self.plot_spec_dict['restore_limit'] = False
            self.plot_spectrum(**self.plot_spec_dict)

        #plot the sky model cube
        if event.key == 's':
            self.plot_spec_dict['skycube'] = self.cleanhdu_flux['UNZAPPED'].data
            self.plot_spec_dict['unzapped_skycube'] = True
            self.plot_spec_dict['restore_limit'] = False
            self.plot_spectrum(**self.plot_spec_dict)

        #unplot the sky model cube
        if event.key == 'n':
            self.plot_spec_dict['skycube'] = None
            self.plot_spec_dict['unzapped_skycube'] = False
            self.plot_spec_dict['restore_limit'] = False
            self.plot_spectrum(**self.plot_spec_dict)


            # self.std['invsens_model'] = bspline()

    def fit_bspline(self, x_full, y_full, bkpt, k, use):
        """
        A simple BSpline-fit wrapper for flux calibration

        Args:
            x (1D arr): x data
            y (1D arr): y data
            bkpt (int): breakpoints (end point of each interval)
            use (boolen arry, same shape as x and y): True for pixels used for fitting
        """
        x = x_full[use]
        y = y_full[use]
        # x_idx = np.arange(x.size)
        # n_knots = int(x.size // bkpt)        
        # knots_idx = np.linspace(np.min(x[use])+ 0.1, np.max(x[use])-0.1, n_knots)
        # print(n_knots, knots)
        # print(x[use][0], x[use][-1])
        knots_idx = np.arange(bkpt, x.size, bkpt)
        # knots = x[knots_idx]

        # #insert the first and the last index to the knots
        # if x[knots_idx][-1] -1 < x[-1]:
        #     knots = np.append(knots, x[-1] - 0.5)
        # # print
        # knots = np.insert(knots, 0, x[0]+0.5)
        bspline = splrep(x, y, k = k, task=-1, t=x[knots_idx])

        return splev(x_full, bspline)


    # def set_zap_skyseg(self):
    #     if not self.use_multi_skyseg.get():



    
class RedirectText:
    """
    A class redirecting the system output to the text widget in the GUI
    """
    def __init__(self, text_widget):
        self.text_widget = text_widget

    def write(self, string):
        self.text_widget.insert(tk.END, string)
        self.text_widget.see(tk.END)

    def flush(self):
        pass


def update_text_box(text_box: tk.Text, proc: subprocess.Popen):
    """Create a closure to capture proc and text_box which will be called to update the text box"""

    # closure to capture proc and text_box
    def _update_text_box():
        for line in proc.stdout:
            text_box.insert(tk.END, line.decode())
            text_box.yview(tk.END)

    return _update_text_box


def run_command(command: List[str], text_box: tk.Text):
    """Run a command and redirect output to a text box in real-time"""
    proc = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    # Start a separate thread to update the text box with output
    thread = threading.Thread(target=update_text_box(text_box, proc))
    thread.start()
    proc.wait()

def main():
    root = tk.Tk()
    root.geometry("1200x800")
    root.resizable(True, True)
    app = KCWIViewerApp(root)
    root.bind("<FocusIn>", app.focus_in)
    root.mainloop()

if __name__ == "__main__":
    main()