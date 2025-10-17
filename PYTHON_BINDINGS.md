# GECKO-A Python Bindings

This document describes the Python bindings for GECKO-A and how to use them.

## Overview

The GECKO-A Fortran code has been compiled into a shared library (`libgecko.so`) that can be called from Python using ctypes. This allows for:

1. **Modular Usage**: Use specific GECKO-A functions in Python workflows
2. **Integration**: Integrate GECKO-A with other Python tools
3. **Flexibility**: Build custom analysis pipelines

## Files Added

### Fortran Side

- **LIB/gecko_wrapper.f90**: Fortran module providing C-compatible interfaces (using iso_c_binding) to GECKO-A functions
- **OBJ/makefile**: Updated with `shared` target to compile shared library with `-fPIC` flag

### Python Side

- **LIB/gecko_python.py**: Python module with ctypes bindings and high-level `GeckoA` class
- **main_python.py**: Python reimplementation of main GECKO-A workflow (demonstration)

## Building the Shared Library

### Prerequisites

- gfortran compiler
- GNU Make

### Compilation

```bash
cd OBJ
make clean
make shared
```

This creates `libgecko.so` in the OBJ directory.

## Using Python Bindings

### Basic Example

```python
from gecko_python import GeckoA

# Initialize GECKO-A
gecko = GeckoA()

# Set up defaults and read configuration
gecko.define_defaults()
gecko.read_namelist()
gecko.initialize()

# Load data
gecko.load_reaction_info()
gecko.load_c1_mechanism()
gecko.load_database()

# Read input species
ninp = gecko.read_chemical_input('./cheminput.dat')
print(f"Loaded {ninp} primary species")

# Get status
nrec = gecko.get_nrec()
print(f"Dictionary has {nrec} records")
```

### Running the Python Main Script

```bash
cd RUN
cp ../INPUT/gecko.nml ./
cp ../INPUT/cheminput.dat ./
python3 ../main_python.py
```

**Note**: The namelist file `gecko.nml` must be edited to set correct paths for `dirgecko` and `dirout`.

## Available Functions

The `GeckoA` class provides the following methods:

### Initialization
- `define_defaults()`: Initialize default parameters and flags
- `read_namelist()`: Read namelist configuration file
- `initialize()`: Initialize dictionaries and stacks

### Data Loading
- `load_reaction_info()`: Read reaction info/references
- `load_c1_mechanism()`: Load C1 mechanism
- `load_database()`: Load database
- `read_chemical_input(filename)`: Read primary species from file

### Processing
- `check_parenthesis(chem)`: Check parenthesis in chemical formula
- `process_input_chemical(chem)`: Process input chemical species
- `sort_namelst()`: Sort name list

### Output
- `compute_dictionary_elements()`: Compute molar mass and atoms
- `write_dictionary()`: Write dictionary
- `write_max_yields()`: Write maximum yields
- `write_ro2()`: Write peroxy species
- `change_phase()`: Change phase (gas to particle/wall)
- `write_size()`: Write mechanism size information
- `write_vapor_pressure(pvap_sar)`: Write vapor pressure
- `write_henry()`: Write Henry's law coefficients
- `write_deposition()`: Write deposition parameters
- `write_heat_formation()`: Write heat of formation
- `write_Tg()`: Write glass transition temperature
- `write_diffusion_volume()`: Write diffusion volume
- `run_simple()`: Run simplified output generation

### Status Queries
- `get_nrec()`: Get number of records in dictionary
- `get_nhldvoc()`: Get number of VOCs in stack
- `get_nhldrad()`: Get number of radicals in stack
- `get_dirout()`: Get output directory

## Limitations

The current implementation provides a demonstration of how to interface Python with GECKO-A. The main chemistry loop (species processing with OH, O3, NO3 reactions, etc.) is not yet fully exposed through Python bindings because it requires:

1. Access to complex internal Fortran data structures (stacks, dictionaries)
2. Many interdependent chemistry modules
3. Extensive conditional logic

### Future Improvements

To create a fully functional Python version, consider:

1. **Expose More Data Structures**: Add C bindings for stack and dictionary access
2. **Module-by-Module Wrapping**: Wrap individual chemistry modules (hochem, no3chem, etc.)
3. **Complete Rewrite**: Rewrite chemistry logic in Python (significant effort)
4. **Hybrid Approach**: Keep chemistry in Fortran, expose higher-level workflow functions

## Technical Details

### C Bindings

The Fortran wrapper uses `iso_c_binding` to create C-compatible interfaces:

```fortran
SUBROUTINE c_initdictstack() BIND(C, name="gecko_initdictstack")
  USE maintool, ONLY: initdictstack
  CALL initdictstack()
END SUBROUTINE c_initdictstack
```

### Python ctypes

The Python module uses ctypes to call these functions:

```python
gecko_initdictstack = libgecko.gecko_initdictstack
gecko_initdictstack.argtypes = []
gecko_initdictstack.restype = None
```

### Compilation Flags

The shared library is compiled with:
- `-fPIC`: Position-independent code (required for shared libraries)
- `-shared`: Create shared library
- `-fdefault-real-8`: Use 8-byte reals (consistent with original)
- `-mcmodel=large`: Support large memory model

## Troubleshooting

### Library Not Found

If you get "library not found" errors, ensure:
1. The shared library exists: `ls OBJ/libgecko.so`
2. You're running from the correct directory
3. The path in `gecko_python.py` is correct

### Segmentation Faults

If you encounter segfaults:
1. Ensure you call functions in the correct order (initialize before processing)
2. Check that input files exist and paths are correct
3. Verify the Fortran code was compiled with same flags as when testing

### File Not Found Errors

GECKO-A expects certain files and directories:
- Input files: `cheminput.dat`, `gecko.nml`
- Data directory: `DATA/` (should exist in GECKO-A root)
- Output directory: Must exist before running (specified in `gecko.nml`)

## Example Use Cases

### 1. Batch Processing

```python
from gecko_python import GeckoA

species_list = ["CH3CH3", "CH3CH2CH3", "CH3CH2CH2CH3"]

for species in species_list:
    # Write to cheminput.dat
    with open('cheminput.dat', 'w') as f:
        f.write(f"{species}\nEND\n")
    
    # Run GECKO-A
    gecko = GeckoA()
    gecko.define_defaults()
    gecko.read_namelist()
    gecko.initialize()
    # ... process species
```

### 2. Custom Analysis

```python
from gecko_python import GeckoA

gecko = GeckoA()
# ... initialize and process

# Query status during processing
nrec = gecko.get_nrec()
nvoc = gecko.get_nhldvoc()
nrad = gecko.get_nhldrad()

print(f"Status: {nrec} species, {nvoc} VOCs, {nrad} radicals in stack")
```

## Contributing

To add more Python bindings:

1. Add C-compatible wrapper in `LIB/gecko_wrapper.f90`
2. Add function signature in `LIB/gecko_python.py`
3. Add high-level method to `GeckoA` class
4. Rebuild shared library: `cd OBJ && make shared`
5. Test the new bindings

## References

- GECKO-A Wiki: https://gitlab.in2p3.fr/ipsl/lisa/geckoa/public/gecko-a/-/wikis/home
- Fortran-Python Interface: https://docs.python.org/3/library/ctypes.html
- iso_c_binding: https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fC_005fBINDING.html
