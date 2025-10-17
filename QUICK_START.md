# GECKO-A Quick Start Guide

## For Original Fortran Users

### Build
```bash
cd OBJ
make clean
make
```

### Run
```bash
cd RUN
cp ../INPUT/gecko.nml ./
cp ../INPUT/cheminput.dat ./
cp ../OBJ/cm ./
# Edit gecko.nml and cheminput.dat as needed
./cm
```

## For Python Users

### Build Shared Library
```bash
cd OBJ
make clean
make shared
```

### Test Installation
```bash
cd ..
python3 test_python_bindings.py
```

Expected output:
```
======================================================================
GECKO-A Python Bindings Test Suite
======================================================================
...
Test Results: 4/4 passed
======================================================================
✓ All tests passed!
```

### Run Python Demo
```bash
cd RUN
cp ../INPUT/* ./
# Edit gecko.nml to set correct paths:
#   dirgecko='/path/to/GECKO-A/'
#   dirout='/path/to/GECKO-A/RUN/OUT/'
python3 ../main_python.py
```

## Simple Python Example

Create a file `my_gecko.py`:

```python
#!/usr/bin/env python3
import sys
sys.path.insert(0, 'LIB')

from gecko_python import GeckoA

# Initialize
gecko = GeckoA()
gecko.define_defaults()
gecko.read_namelist()
gecko.initialize()

# Load data
gecko.load_database()
gecko.load_c1_mechanism()

# Get status
print(f"Number of records: {gecko.get_nrec()}")
print(f"VOCs in stack: {gecko.get_nhldvoc()}")
print(f"Radicals in stack: {gecko.get_nhldrad()}")
```

Run it:
```bash
cd RUN
cp ../INPUT/gecko.nml ./
python3 ../my_gecko.py
```

## File Structure

```
GECKO-A/
├── OBJ/
│   ├── makefile        # Build system
│   ├── cm              # Fortran executable (after 'make')
│   └── libgecko.so     # Shared library (after 'make shared')
│
├── LIB/
│   ├── *.f90           # Fortran source files
│   ├── gecko_wrapper.f90   # C bindings
│   └── gecko_python.py     # Python interface
│
├── RUN/
│   ├── gecko.sh        # Run script
│   └── OUT/            # Output directory
│
├── INPUT/
│   ├── gecko.nml       # Configuration
│   └── cheminput.dat   # Input species
│
├── DATA/               # Database files
│
└── Documentation:
    ├── README.md               # Original README
    ├── PYTHON_BINDINGS.md      # API reference
    ├── USAGE_EXAMPLE.md        # Detailed examples
    ├── IMPLEMENTATION_SUMMARY.md   # Implementation details
    └── QUICK_START.md          # This file
```

## Common Issues

### Issue: Library not found
```
OSError: libgecko.so: cannot open shared object file
```

**Solution**: Make sure you're in the correct directory and the library exists:
```bash
ls OBJ/libgecko.so  # Should show the file
```

### Issue: Module not found
```
ModuleNotFoundError: No module named 'gecko_python'
```

**Solution**: Add LIB directory to Python path:
```python
import sys
sys.path.insert(0, 'LIB')
from gecko_python import GeckoA
```

### Issue: Segmentation fault
```
Segmentation fault (core dumped)
```

**Solution**: Call functions in the correct order:
1. `define_defaults()`
2. `read_namelist()` or `initialize()`
3. Other functions

### Issue: File not found
```
Error: File not found: cheminput.dat
```

**Solution**: Make sure you're running from the RUN directory:
```bash
cd RUN
cp ../INPUT/cheminput.dat ./
```

## Build Targets

### Fortran Executable
```bash
make          # Build cm executable
make clean    # Remove build artifacts
```

### Shared Library
```bash
make shared   # Build libgecko.so
```

### Both
```bash
make clean && make && make shared
```

## Getting Help

- **API Reference**: See `PYTHON_BINDINGS.md`
- **Examples**: See `USAGE_EXAMPLE.md`
- **Implementation**: See `IMPLEMENTATION_SUMMARY.md`
- **Original docs**: See `README.md`
- **GECKO-A wiki**: https://gitlab.in2p3.fr/ipsl/lisa/geckoa/public/gecko-a/-/wikis/home

## Next Steps

1. **Try the test suite**: `python3 test_python_bindings.py`
2. **Run the demo**: `python3 main_python.py` (from RUN/)
3. **Read the docs**: Start with `PYTHON_BINDINGS.md`
4. **Write custom scripts**: Use `GeckoA` class
5. **Extend bindings**: Add more functions as needed

## Tips

- Always call `define_defaults()` first
- Use `initialize()` before processing
- Check return values from functions
- Run test suite after rebuilding library
- Keep backup of original `cm` executable

## Version Info

- GECKO-A: Original Fortran version
- Python Bindings: 2024 implementation
- Tested with:
  - Python 3.12.3
  - gfortran 13.3.0
  - Ubuntu 24.04 LTS
