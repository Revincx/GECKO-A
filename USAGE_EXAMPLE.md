# GECKO-A Usage Examples

This document provides examples of using both the original Fortran executable and the new Python bindings.

## Table of Contents

1. [Original Fortran Usage](#original-fortran-usage)
2. [Python Bindings Usage](#python-bindings-usage)
3. [Comparison](#comparison)

## Original Fortran Usage

### Building the Executable

```bash
cd OBJ
make clean
make
```

This creates the `cm` executable.

### Running a Simulation

```bash
cd RUN
cp ../INPUT/gecko.nml ./
cp ../INPUT/cheminput.dat ./
cp ../OBJ/cm ./

# Edit gecko.nml and cheminput.dat as needed

./cm

# Or use the provided script:
./gecko.sh
```

### Example: Processing n-heptane

1. Edit `cheminput.dat`:
```
CH3CH2CH2CH2CH2CH2CH3
END
```

2. Edit `gecko.nml` (set correct paths):
```fortran
&dir
dirgecko='/path/to/GECKO-A/'
dirout='/path/to/GECKO-A/RUN/OUT/'
/
&thresholds
critvp = -13
brcut  = 0.05
yldcut = 1E-3
rxloss = 1E-10
/
&env_cond
TK     = 298.
/
&reductions
maxgen = 1
isomerfg = .TRUE.
highnoxfg = .FALSE.
rx_ro2_multiclass = .TRUE.
/
&sar
pvap_sar = 2
/
&process
g2pfg = .TRUE.
g2wfg = .FALSE.
/
```

3. Run:
```bash
./cm
```

4. Check output in `OUT/` directory:
   - `scheme.log`: Log file with settings
   - `reactions.dum`: Generated reactions
   - `gasspe.dum`: Gas phase species
   - `dictionnary.dat`: Species dictionary
   - And more...

## Python Bindings Usage

### Building the Shared Library

```bash
cd OBJ
make clean
make shared
```

This creates `libgecko.so`.

### Test the Bindings

```bash
python3 test_python_bindings.py
```

All tests should pass:
```
======================================================================
GECKO-A Python Bindings Test Suite
======================================================================
Test 1: Loading Python module...
✓ Module loaded successfully

Test 2: Loading shared library...
✓ Shared library loaded successfully

Test 4: Creating GeckoA object...
✓ GeckoA object created successfully

Test 3: Calling basic functions...
  Testing define_defaults()...
  ✓ define_defaults() executed
  Testing initialize()...
  ✓ initialize() executed
  Testing get_nrec()...
  ✓ get_nrec() returned: 0

======================================================================
Test Results: 4/4 passed
======================================================================

✓ All tests passed!
```

### Using the Python Interface

#### Example 1: Basic Initialization

```python
#!/usr/bin/env python3
import sys
sys.path.insert(0, 'LIB')
from gecko_python import GeckoA

# Initialize GECKO-A
gecko = GeckoA()
gecko.define_defaults()
gecko.read_namelist()
gecko.initialize()

print(f"Initialized with {gecko.get_nrec()} records")
```

#### Example 2: Loading Data

```python
#!/usr/bin/env python3
import sys
sys.path.insert(0, 'LIB')
from gecko_python import GeckoA

# Setup
gecko = GeckoA()
gecko.define_defaults()
gecko.read_namelist()
gecko.initialize()

# Load all databases
print("Loading reaction info...")
gecko.load_reaction_info()

print("Loading C1 mechanism...")
gecko.load_c1_mechanism()
gecko.write_kc1()

print("Loading database...")
gecko.load_database()

# Read input species
print("Reading input species...")
ninp = gecko.read_chemical_input('./cheminput.dat')
print(f"Loaded {ninp} primary species")

# Get status
nrec = gecko.get_nrec()
print(f"Dictionary has {nrec} records")
```

#### Example 3: Batch Processing (Concept)

```python
#!/usr/bin/env python3
"""
Process multiple species in a loop
Note: Full chemistry processing not yet implemented
"""
import sys
sys.path.insert(0, 'LIB')
from gecko_python import GeckoA

species_list = [
    "CH3CH3",           # ethane
    "CH3CH2CH3",        # propane
    "CH3CH2CH2CH3",     # butane
]

for species in species_list:
    print(f"\nProcessing {species}...")
    
    # Write to cheminput.dat
    with open('cheminput.dat', 'w') as f:
        f.write(f"{species}\nEND\n")
    
    # Initialize and process
    gecko = GeckoA()
    gecko.define_defaults()
    gecko.read_namelist()
    gecko.initialize()
    
    # Load databases (only need to do once if we restructure)
    gecko.load_reaction_info()
    gecko.load_c1_mechanism()
    gecko.load_database()
    
    # Read input
    ninp = gecko.read_chemical_input('./cheminput.dat')
    
    # Process (would need full implementation)
    # ... chemistry loop here ...
    
    # Write outputs
    gecko.compute_dictionary_elements()
    gecko.write_dictionary()
    gecko.write_ro2()
    gecko.write_size()
    
    print(f"Completed processing {species}")
```

### Running the Demo Python Script

```bash
cd RUN
cp ../INPUT/gecko.nml ./
cp ../INPUT/cheminput.dat ./

# Edit gecko.nml to set correct paths
# Edit cheminput.dat to set species

python3 ../main_python.py
```

Example output:
```
======================================================================
GECKO-A - Python Version
Generator for Explicit Chemistry and Kinetics of Organics in the Atmosphere
======================================================================

Initializing...
Output directory: ./OUT/
Opening output files...
Reading data...
Reading the list of primary species...
Number of primary species: 1

======================================================================
Note: The chemistry loop processing is handled by Fortran modules.
A full Python implementation would require extensive additional
bindings or a complete rewrite of the chemistry modules.
======================================================================

Writing output files...
Number of records in dictionary: 0
Computing molar masses...
Writing dictionary...
Writing mass transfer equations (if any)...

======================================================================
GECKO-A processing completed!
======================================================================

For a complete simulation, use the original Fortran executable
or extend the Python bindings to include chemistry processing.
```

## Comparison

### Original Fortran Executable

**Advantages:**
- Complete implementation with all chemistry
- Tested and validated
- Fast execution
- Production-ready

**Disadvantages:**
- Not easily integrated with Python workflows
- Limited flexibility for custom analysis
- Standalone tool

### Python Bindings

**Advantages:**
- Can be integrated into Python workflows
- Allows custom analysis and processing
- Modular - use only what you need
- Easy to extend for new use cases

**Disadvantages:**
- Not all features exposed yet
- Main chemistry loop not fully wrapped
- Requires both Fortran library and Python
- Still under development

### When to Use Which

**Use the Fortran Executable when:**
- Running complete GECKO-A simulations
- Production work
- Standard mechanism generation
- Batch processing with shell scripts

**Use the Python Bindings when:**
- Integrating GECKO-A into Python pipelines
- Building custom tools on top of GECKO-A
- Automating analysis workflows
- Prototyping new features
- Needing programmatic access to GECKO-A functions

## Current Limitations

The Python bindings currently support:

✓ Initialization and configuration
✓ Database loading
✓ Input species reading
✓ Basic output generation
✓ Status queries

Not yet supported:

✗ Full chemistry loop processing
✗ Species stack manipulation from Python
✗ Direct access to internal data structures
✗ Individual chemistry module calls (OH, O3, NO3, etc.)

## Future Development

To make Python bindings fully functional:

1. **Expose Stack Operations**: Add bindings to manipulate species stacks
2. **Wrap Chemistry Modules**: Create bindings for hochem, no3chem, etc.
3. **Data Structure Access**: Expose dictionaries and internal data
4. **Complete Workflow**: Implement full chemistry loop in Python or Fortran wrapper

## Getting Help

- Documentation: See `PYTHON_BINDINGS.md` for detailed API reference
- Test suite: Run `python3 test_python_bindings.py` to verify setup
- Original docs: See `README.md` and GECKO-A wiki

## Examples Directory Structure

```
GECKO-A/
├── OBJ/
│   ├── makefile          (updated with 'shared' target)
│   ├── cm                (Fortran executable)
│   └── libgecko.so       (shared library)
├── LIB/
│   ├── gecko_wrapper.f90 (C bindings)
│   └── gecko_python.py   (Python interface)
├── RUN/
│   ├── gecko.sh          (original run script)
│   └── OUT/              (output directory)
├── INPUT/
│   ├── cheminput.dat     (input species)
│   └── gecko.nml         (configuration)
├── main_python.py        (Python demo script)
├── test_python_bindings.py (test suite)
├── PYTHON_BINDINGS.md    (API documentation)
└── USAGE_EXAMPLE.md      (this file)
```
