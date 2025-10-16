"""
GECKO-A Python Module
=====================

Python interface to GECKO-A C library.

This module provides Python bindings to the GECKO-A atmospheric chemistry
mechanism generator, converted from Fortran to C.

Modules
-------
- keyparameter: Core constants and parameters
- keyflag: Configuration flags and settings
- minidict: Mini dictionary for species lookup
- toolbox: Common utility functions

Example
-------
>>> import geckoa
>>> 
>>> # Initialize default parameters
>>> geckoa.define_defaults()
>>> 
>>> # Get configuration
>>> config = geckoa.get_config()
>>> print(f"Temperature: {config['TK']} K")
>>> 
>>> # Use mini dictionary
>>> geckoa.clean_minid()
>>> geckoa.add_fo("CH4", "CH4")
>>> formula = geckoa.get_fo("CH4")
>>> print(f"Formula for CH4: {formula}")
>>> 
>>> # Calculate Arrhenius rate constant
>>> k = geckoa.kval([1.5e-12, 0.0, 0.0], 298.0)
>>> print(f"Rate constant: {k}")

Constants
---------
MXLCO : int
    Maximum length of species names (code)
MXLFO : int
    Maximum length of a formula
MXNODE : int
    Maximum number of nodes allowed
MXPS : int
    Maximum number of "primary" species
MXPD : int
    Maximum number of products per reaction
"""

# This will be the compiled C extension module
# The actual implementation is in geckoa_module.c
