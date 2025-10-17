#!/usr/bin/env python3
"""
GECKO-A Python bindings

This module provides Python interfaces to the GECKO-A Fortran library
using ctypes to call C-compatible functions.
"""

import ctypes
import os
import sys

# Get the path to the shared library
LIB_PATH = os.path.join(os.path.dirname(__file__), '..', 'OBJ', 'libgecko.so')

# Load the shared library
try:
    libgecko = ctypes.CDLL(LIB_PATH)
except OSError as e:
    print(f"Error loading GECKO-A shared library from {LIB_PATH}: {e}")
    sys.exit(1)

# Define C function signatures

# Initialize default parameters
gecko_define_defaults = libgecko.gecko_define_defaults
gecko_define_defaults.argtypes = []
gecko_define_defaults.restype = None

# Read namelist
gecko_rd_nml = libgecko.gecko_rd_nml
gecko_rd_nml.argtypes = []
gecko_rd_nml.restype = None

# Initialize dictionaries and stacks
gecko_initdictstack = libgecko.gecko_initdictstack
gecko_initdictstack.argtypes = []
gecko_initdictstack.restype = None

# Write log
gecko_wrtlog = libgecko.gecko_wrtlog
gecko_wrtlog.argtypes = []
gecko_wrtlog.restype = None

# Read reaction info
gecko_rdreac_info = libgecko.gecko_rdreac_info
gecko_rdreac_info.argtypes = []
gecko_rdreac_info.restype = None

# Load C1 mechanism
gecko_loadc1mch = libgecko.gecko_loadc1mch
gecko_loadc1mch.argtypes = []
gecko_loadc1mch.restype = None

# Write kc1
gecko_wrt_kc1 = libgecko.gecko_wrt_kc1
gecko_wrt_kc1.argtypes = []
gecko_wrt_kc1.restype = None

# Load database
gecko_loaddb = libgecko.gecko_loaddb
gecko_loaddb.argtypes = []
gecko_loaddb.restype = None

# Read chemical input
gecko_rdchemin = libgecko.gecko_rdchemin
gecko_rdchemin.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.POINTER(ctypes.c_int)]
gecko_rdchemin.restype = None

# Check parenthesis
gecko_check_parenthesis = libgecko.gecko_check_parenthesis
gecko_check_parenthesis.argtypes = [ctypes.c_char_p, ctypes.c_int]
gecko_check_parenthesis.restype = None

# Process input chemical
gecko_in1chm = libgecko.gecko_in1chm
gecko_in1chm.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_int]
gecko_in1chm.restype = None

# Sort name list
gecko_sort_namlst = libgecko.gecko_sort_namlst
gecko_sort_namlst.argtypes = []
gecko_sort_namlst.restype = None

# Compute dictionary elements
gecko_dictelement = libgecko.gecko_dictelement
gecko_dictelement.argtypes = []
gecko_dictelement.restype = None

# Write dictionary
gecko_wrt_dict = libgecko.gecko_wrt_dict
gecko_wrt_dict.argtypes = []
gecko_wrt_dict.restype = None

# Write max yields
gecko_wrt_mxyield = libgecko.gecko_wrt_mxyield
gecko_wrt_mxyield.argtypes = []
gecko_wrt_mxyield.restype = None

# Write RO2 species
gecko_wrt_ro2 = libgecko.gecko_wrt_ro2
gecko_wrt_ro2.argtypes = []
gecko_wrt_ro2.restype = None

# Change phase
gecko_changephase = libgecko.gecko_changephase
gecko_changephase.argtypes = []
gecko_changephase.restype = None

# Write size
gecko_wrt_size = libgecko.gecko_wrt_size
gecko_wrt_size.argtypes = []
gecko_wrt_size.restype = None

# Write vapor pressure
gecko_wrt_psat = libgecko.gecko_wrt_psat
gecko_wrt_psat.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int]
gecko_wrt_psat.restype = None

# Write Henry's law coefficient
gecko_wrt_henry = libgecko.gecko_wrt_henry
gecko_wrt_henry.argtypes = []
gecko_wrt_henry.restype = None

# Write deposition
gecko_wrt_depo = libgecko.gecko_wrt_depo
gecko_wrt_depo.argtypes = []
gecko_wrt_depo.restype = None

# Write heat of formation
gecko_wrt_heatf = libgecko.gecko_wrt_heatf
gecko_wrt_heatf.argtypes = []
gecko_wrt_heatf.restype = None

# Write Tg
gecko_wrt_Tg = libgecko.gecko_wrt_Tg
gecko_wrt_Tg.argtypes = []
gecko_wrt_Tg.restype = None

# Write diffusion volume
gecko_diffusion_vol = libgecko.gecko_diffusion_vol
gecko_diffusion_vol.argtypes = []
gecko_diffusion_vol.restype = None

# Get number of records
gecko_get_nrec = libgecko.gecko_get_nrec
gecko_get_nrec.argtypes = []
gecko_get_nrec.restype = ctypes.c_int

# Get number of VOCs in stack
gecko_get_nhldvoc = libgecko.gecko_get_nhldvoc
gecko_get_nhldvoc.argtypes = []
gecko_get_nhldvoc.restype = ctypes.c_int

# Get number of radicals in stack
gecko_get_nhldrad = libgecko.gecko_get_nhldrad
gecko_get_nhldrad.argtypes = []
gecko_get_nhldrad.restype = ctypes.c_int

# Get output directory
gecko_get_dirout = libgecko.gecko_get_dirout
gecko_get_dirout.argtypes = [ctypes.c_char_p, ctypes.c_int]
gecko_get_dirout.restype = None

# Run simple processing
gecko_run_simple = libgecko.gecko_run_simple
gecko_run_simple.argtypes = []
gecko_run_simple.restype = None


# High-level Python interface
class GeckoA:
    """High-level interface to GECKO-A library"""
    
    def __init__(self):
        """Initialize GECKO-A"""
        pass
    
    def define_defaults(self):
        """Initialize default parameters and flags"""
        gecko_define_defaults()
    
    def read_namelist(self):
        """Read namelist configuration file"""
        gecko_rd_nml()
    
    def initialize(self):
        """Initialize dictionaries and stacks"""
        gecko_initdictstack()
    
    def write_log(self):
        """Write log information"""
        gecko_wrtlog()
    
    def load_reaction_info(self):
        """Read reaction info/references"""
        gecko_rdreac_info()
    
    def load_c1_mechanism(self):
        """Load C1 mechanism"""
        gecko_loadc1mch()
    
    def write_kc1(self):
        """Write koh for C1 species"""
        gecko_wrt_kc1()
    
    def load_database(self):
        """Load database"""
        gecko_loaddb()
    
    def read_chemical_input(self, filename):
        """Read chemical input from file
        
        Args:
            filename (str): Path to cheminput.dat file
            
        Returns:
            int: Number of primary species
        """
        ninp = ctypes.c_int()
        filename_bytes = filename.encode('utf-8')
        gecko_rdchemin(filename_bytes, len(filename_bytes), ctypes.byref(ninp))
        return ninp.value
    
    def check_parenthesis(self, chem):
        """Check parenthesis in chemical formula
        
        Args:
            chem (str): Chemical formula
            
        Returns:
            str: Checked chemical formula
        """
        chem_bytes = bytearray(chem.encode('utf-8'))
        gecko_check_parenthesis(ctypes.c_char_p(bytes(chem_bytes)), len(chem_bytes))
        return chem_bytes.decode('utf-8').rstrip('\x00')
    
    def process_input_chemical(self, chem):
        """Process input chemical species
        
        Args:
            chem (str): Chemical formula
            
        Returns:
            str: ID name of the species
        """
        chem_bytes = chem.encode('utf-8')
        idnam_bytes = bytearray(10)  # 6 chars + null terminator + padding
        gecko_in1chm(chem_bytes, len(chem_bytes), 
                     ctypes.c_char_p(bytes(idnam_bytes)), len(idnam_bytes))
        return idnam_bytes.decode('utf-8').rstrip('\x00 ')
    
    def sort_namelst(self):
        """Sort name list"""
        gecko_sort_namlst()
    
    def compute_dictionary_elements(self):
        """Compute molar mass and atoms for species"""
        gecko_dictelement()
    
    def write_dictionary(self):
        """Write dictionary"""
        gecko_wrt_dict()
    
    def write_max_yields(self):
        """Write maximum yields"""
        gecko_wrt_mxyield()
    
    def write_ro2(self):
        """Write peroxy species"""
        gecko_wrt_ro2()
    
    def change_phase(self):
        """Change phase (gas to particle/wall)"""
        gecko_changephase()
    
    def write_size(self):
        """Write mechanism size information"""
        gecko_wrt_size()
    
    def write_vapor_pressure(self, pvap_sar):
        """Write vapor pressure
        
        Args:
            pvap_sar (int): SAR selector (1=JR-MY, 2=Nannoolal, 3=SIMPOL-1)
        """
        if pvap_sar == 1:
            gecko_wrt_psat(1, 0, 0)
        elif pvap_sar == 2:
            gecko_wrt_psat(0, 1, 0)
        elif pvap_sar == 3:
            gecko_wrt_psat(0, 0, 1)
    
    def write_henry(self):
        """Write Henry's law coefficients"""
        gecko_wrt_henry()
    
    def write_deposition(self):
        """Write deposition parameters"""
        gecko_wrt_depo()
    
    def write_heat_formation(self):
        """Write heat of formation"""
        gecko_wrt_heatf()
    
    def write_Tg(self):
        """Write glass transition temperature"""
        gecko_wrt_Tg()
    
    def write_diffusion_volume(self):
        """Write diffusion volume"""
        gecko_diffusion_vol()
    
    def get_nrec(self):
        """Get number of records in dictionary
        
        Returns:
            int: Number of records
        """
        return gecko_get_nrec()
    
    def get_nhldvoc(self):
        """Get number of VOCs in stack
        
        Returns:
            int: Number of VOCs
        """
        return gecko_get_nhldvoc()
    
    def get_nhldrad(self):
        """Get number of radicals in stack
        
        Returns:
            int: Number of radicals
        """
        return gecko_get_nhldrad()
    
    def get_dirout(self):
        """Get output directory
        
        Returns:
            str: Output directory path
        """
        dirout_bytes = bytearray(256)
        gecko_get_dirout(ctypes.c_char_p(bytes(dirout_bytes)), len(dirout_bytes))
        return dirout_bytes.decode('utf-8').rstrip('\x00 ')
    
    def run_simple(self):
        """Run simplified processing (write outputs without full chemistry loop)"""
        gecko_run_simple()
