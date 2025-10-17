# GECKO-A Python Bindings

Python interface to the GECKO-A atmospheric chemistry mechanism generator.

## What's New

GECKO-A can now be compiled as a shared library and called from Python! This enables:

- 🐍 **Python Integration**: Use GECKO-A functions in Python scripts
- 🔧 **Modular Access**: Call individual functions without running full simulation  
- 📊 **Custom Workflows**: Build analysis pipelines with Python
- 🔌 **Interoperability**: Integrate with other Python atmospheric chemistry tools

## Quick Links

- [Quick Start Guide](QUICK_START.md) - Get started in 5 minutes
- [Python API Reference](PYTHON_BINDINGS.md) - Complete API documentation
- [Usage Examples](USAGE_EXAMPLE.md) - Fortran vs Python examples
- [Implementation Details](IMPLEMENTATION_SUMMARY.md) - Technical details
- [Original README](README.md) - Original GECKO-A documentation

## Installation

### 1. Build the Shared Library

```bash
cd OBJ
make clean
make shared
```

This creates `libgecko.so` in the OBJ directory.

### 2. Test the Installation

```bash
cd ..
python3 test_python_bindings.py
```

You should see:
```
✓ All tests passed!
```

### 3. Try the Demo

```bash
cd RUN
cp ../INPUT/gecko.nml ./
cp ../INPUT/cheminput.dat ./
# Edit gecko.nml to set correct paths
python3 ../main_python.py
```

## Usage Example

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

# Load databases
gecko.load_reaction_info()
gecko.load_c1_mechanism()
gecko.load_database()

# Read input species
ninp = gecko.read_chemical_input('./cheminput.dat')
print(f"Loaded {ninp} primary species")

# Query status
nrec = gecko.get_nrec()
print(f"Dictionary has {nrec} records")

# Write outputs
gecko.compute_dictionary_elements()
gecko.write_dictionary()
gecko.write_ro2()
```

## Available Functions

The `GeckoA` class provides 30+ methods including:

### Initialization
- `define_defaults()` - Set default parameters
- `read_namelist()` - Read configuration
- `initialize()` - Initialize data structures

### Data Loading  
- `load_database()` - Load chemistry database
- `load_c1_mechanism()` - Load C1 chemistry
- `read_chemical_input(filename)` - Read input species

### Processing
- `check_parenthesis(chem)` - Validate formula
- `process_input_chemical(chem)` - Process species
- `sort_namelst()` - Sort name list

### Output
- `write_dictionary()` - Write species dictionary
- `write_ro2()` - Write peroxy species
- `write_vapor_pressure(sar)` - Write Pvap estimates
- Many more...

### Status Queries
- `get_nrec()` - Number of species
- `get_nhldvoc()` - VOCs in stack
- `get_nhldrad()` - Radicals in stack

See [PYTHON_BINDINGS.md](PYTHON_BINDINGS.md) for complete API reference.

## Project Structure

```
GECKO-A/
├── OBJ/
│   ├── makefile          # Updated with 'shared' target
│   ├── cm                # Original Fortran executable
│   └── libgecko.so       # New shared library
│
├── LIB/
│   ├── *.f90             # Original Fortran modules
│   ├── gecko_wrapper.f90 # NEW: C bindings
│   └── gecko_python.py   # NEW: Python interface
│
├── main_python.py        # NEW: Python demo
├── test_python_bindings.py  # NEW: Test suite
│
└── Documentation (NEW):
    ├── PYTHON_BINDINGS.md      # API reference
    ├── USAGE_EXAMPLE.md        # Usage examples
    ├── IMPLEMENTATION_SUMMARY.md  # Implementation details
    ├── QUICK_START.md          # Quick start
    └── README_PYTHON.md        # This file
```

## Features

### ✅ What's Implemented

- **Shared Library Compilation**: Full GECKO-A as shared library
- **C Bindings**: 30+ functions with C-compatible interfaces
- **Python Wrapper**: High-level `GeckoA` class
- **Type Safety**: Proper type definitions for all functions
- **Documentation**: Comprehensive docs (1400+ lines)
- **Test Suite**: Automated testing
- **Examples**: Working examples and demos

### ⚠️ Current Limitations

The main chemistry loop (OH, O3, NO3 reactions, etc.) is not yet fully exposed because it requires:
- Access to internal Fortran data structures
- Many interdependent chemistry modules
- Complex control flow

**Workaround**: Use the original Fortran executable (`cm`) for full chemistry processing, or extend the Python bindings as needed.

### 🚀 Future Enhancements

- Expose stack manipulation functions
- Add chemistry module wrappers
- Implement full chemistry loop in wrapper
- Add data export to JSON/CSV
- Real-time progress monitoring

## Compatibility

### Tested With
- ✓ Python 3.12.3
- ✓ gfortran 13.3.0
- ✓ Ubuntu 24.04 LTS
- ✓ Linux x86_64

### Requirements
- Python 3.x (standard library only, no external dependencies)
- gfortran compiler
- GNU Make
- GECKO-A data files

## Comparison: Fortran vs Python

### Use Fortran Executable When:
- Running complete GECKO-A simulations
- Production mechanism generation
- Standard workflows
- Maximum performance needed

### Use Python Bindings When:
- Integrating with Python pipelines
- Building custom analysis tools
- Prototyping new features
- Needing programmatic access

**Best of Both Worlds**: Use Python for workflow orchestration and the Fortran library for computation!

## Examples

### Example 1: Batch Processing

```python
species_list = ["CH3CH3", "CH3CH2CH3", "CH3CH2CH2CH3"]

for species in species_list:
    with open('cheminput.dat', 'w') as f:
        f.write(f"{species}\nEND\n")
    
    gecko = GeckoA()
    gecko.define_defaults()
    gecko.read_namelist()
    # ... process species ...
```

### Example 2: Custom Analysis

```python
gecko = GeckoA()
# ... initialize and load data ...

# Monitor processing
nrec = gecko.get_nrec()
nvoc = gecko.get_nhldvoc()
nrad = gecko.get_nhldrad()
print(f"Status: {nrec} species, {nvoc} VOCs, {nrad} radicals")
```

### Example 3: Integration with Other Tools

```python
import numpy as np
import pandas as pd
from gecko_python import GeckoA

# Use GECKO-A with other Python tools
gecko = GeckoA()
# ... generate mechanism ...

# Analyze with pandas/numpy
# Export to other formats
# Visualize with matplotlib
# Integrate with other models
```

## Troubleshooting

See [QUICK_START.md](QUICK_START.md#common-issues) for common issues and solutions.

## Contributing

To add new Python bindings:

1. Add C-compatible wrapper in `LIB/gecko_wrapper.f90`:
```fortran
SUBROUTINE c_my_function() BIND(C, name="gecko_my_function")
  USE my_module, ONLY: my_function
  CALL my_function()
END SUBROUTINE
```

2. Add Python binding in `LIB/gecko_python.py`:
```python
gecko_my_function = libgecko.gecko_my_function
gecko_my_function.argtypes = []
gecko_my_function.restype = None
```

3. Add method to `GeckoA` class:
```python
def my_function(self):
    """My function description"""
    gecko_my_function()
```

4. Rebuild: `cd OBJ && make shared`

5. Test: Add test to `test_python_bindings.py`

## Performance

The Python bindings have minimal overhead:
- Function call overhead: ~100-1000 ns
- String conversion: ~microseconds  
- Core computation: Native Fortran speed

For GECKO-A's use case, this overhead is negligible.

## License

Same as original GECKO-A. See [LICENSE](LICENSE).

## Citation

If you use GECKO-A in your research, please cite the original papers (see GECKO-A wiki).

If you use the Python bindings specifically, please also mention:
```
GECKO-A Python Bindings (2024)
GitHub: https://github.com/Revincx/GECKO-A
```

## Support

- **Documentation**: See docs in this repository
- **GECKO-A Help**: See GECKO-A wiki (link in README.md)
- **Issues**: Open an issue on GitHub
- **Questions**: Check documentation first, then ask

## Version History

- **2024**: Python bindings added
  - Shared library compilation
  - C bindings with iso_c_binding
  - Python ctypes wrapper
  - Comprehensive documentation
  - Test suite

- **Original**: Fortran implementation
  - See original README.md and AUTHORS.md

## Acknowledgments

- Original GECKO-A developers (see AUTHORS.md)
- Python bindings implementation: GitHub Copilot (2024)

## Related Projects

- GECKO-A Web Interface: https://geckoa.ineris.fr/
- GECKO-A Wiki: https://gitlab.in2p3.fr/ipsl/lisa/geckoa/public/gecko-a/-/wikis/home

## Learn More

- [Quick Start](QUICK_START.md) - Start here!
- [API Reference](PYTHON_BINDINGS.md) - Complete function list
- [Examples](USAGE_EXAMPLE.md) - Detailed examples
- [Implementation](IMPLEMENTATION_SUMMARY.md) - How it works

---

**Ready to get started?** See [QUICK_START.md](QUICK_START.md)!
