#!/usr/bin/env python3
"""
GECKO-A Main Program (Python Version)

This is a Python reimplementation of the main GECKO-A program logic.
It uses the compiled Fortran shared library through Python bindings.

Note: This is a simplified version that handles the initialization and 
basic workflow. The complex chemistry loop logic remains in Fortran for now,
as it involves intricate interactions with many chemistry modules.
"""

import os
import sys

# Add the LIB directory to Python path
lib_path = os.path.join(os.path.dirname(__file__), 'LIB')
sys.path.insert(0, lib_path)

from gecko_python import GeckoA

def main():
    """Main GECKO-A program logic"""
    
    print("=" * 70)
    print("GECKO-A - Python Version")
    print("Generator for Explicit Chemistry and Kinetics of Organics in the Atmosphere")
    print("=" * 70)
    print()
    
    # Initialize GECKO-A interface
    gecko = GeckoA()
    
    # ------------------------------------------
    # INITIALIZATION
    # ------------------------------------------
    print("Initializing...")
    
    # Get default values for user-settable parameters & flags
    gecko.define_defaults()
    
    # Read namelist
    gecko.read_namelist()
    
    # Initialize data in the dictionaries and stacks
    gecko.initialize()
    
    # Get output directory
    dirout = gecko.get_dirout()
    print(f"Output directory: {dirout}")
    
    # ------------------------------------------
    # OPEN OUTPUT FILES
    # ------------------------------------------
    print("Opening output files...")
    
    # Note: File opening is handled by Fortran code when functions are called
    # The Fortran library manages file units and writes to them
    
    # Open log file and write log
    gecko.write_log()
    
    # ------------------------------------------
    # READ ALL DATA
    # ------------------------------------------
    print("Reading data...")
    
    # Read comment/references dictionary
    gecko.load_reaction_info()
    
    # Read the dictionaries and mechanisms for C1, inorg
    gecko.load_c1_mechanism()
    gecko.write_kc1()
    
    # Read all data to be stored in the database module
    gecko.load_database()
    
    # Load species from cheminput.dat
    print("Reading the list of primary species...")
    filename = './cheminput.dat'
    ninp = gecko.read_chemical_input(filename)
    print(f"Number of primary species: {ninp}")
    
    # ------------------------------------------
    # NOTE ABOUT CHEMISTRY LOOPS
    # ------------------------------------------
    # The main chemistry loop is complex and involves:
    # 1. Loop over primary species
    # 2. For each species, loop over the chemical stacks (radicals and VOCs)
    # 3. For each species in stack:
    #    - Parse chemical structure
    #    - Apply various chemistry transformations (OH, O3, NO3, photolysis, etc.)
    #    - Generate products and add to stacks
    #    - Continue until stacks are empty or generation limit reached
    #
    # This requires access to many internal Fortran data structures and
    # chemistry modules that are not easily exposed through simple C bindings.
    # 
    # For a complete Python implementation, you would need to either:
    # 1. Expose many more Fortran data structures and functions
    # 2. Rewrite the entire chemistry logic in Python
    # 3. Create a hybrid approach with key functions exposed
    #
    # For now, this demonstrates the basic workflow and shows how to
    # call the Fortran library from Python.
    
    print()
    print("=" * 70)
    print("Note: The chemistry loop processing is handled by Fortran modules.")
    print("A full Python implementation would require extensive additional")
    print("bindings or a complete rewrite of the chemistry modules.")
    print("=" * 70)
    print()
    
    # ------------------------------------------
    # END OF REACTIONS - WRITE DATA OUT
    # ------------------------------------------
    print("Writing output files...")
    
    # Get status
    nrec = gecko.get_nrec()
    print(f"Number of records in dictionary: {nrec}")
    
    # Compute molar mass and atoms for species
    print("Computing molar masses...")
    gecko.compute_dictionary_elements()
    
    # Write the dictionaries
    print("Writing dictionary...")
    gecko.write_dictionary()
    
    # Run simplified output generation
    # This writes dictionaries, peroxy species, mass transport, and size info
    gecko.run_simple()
    
    print()
    print("=" * 70)
    print("GECKO-A processing completed!")
    print("=" * 70)
    print()
    print("For a complete simulation, use the original Fortran executable")
    print("or extend the Python bindings to include chemistry processing.")
    

if __name__ == '__main__':
    main()
