# GECKO-A C Version

This directory contains the C conversion of GECKO-A Fortran code.

## Conversion Status

### Completed Modules (8 of 81)

1. **keyparameter** - Core constants and parameters
2. **keyflag** - Configuration flags and settings
3. **minidict** - Mini dictionary for species lookup
4. **tempoci** - Temporary OCI data structures
5. **references** - Reference management
6. **tempflag** - Temporary computation flags
7. **sortstring** - String sorting utilities (quicksort + insertion sort)
8. **toolbox** - Common utility functions

### Build System

- **Makefile.c** - Makefile for C compilation
- Successfully compiles and links all converted modules
- Produces executable `cm` in OBJ directory

## Building

```bash
cd OBJ
make -f Makefile.c clean
make -f Makefile.c all
```

## Running

```bash
cd OBJ
./cm
```

The current version demonstrates:
- Module initialization
- Configuration parameter management
- Mini dictionary operations
- Utility functions

## C Code Structure

Each Fortran module is converted to:
- `.h` header file with declarations
- `.c` implementation file

Key design decisions:
- Fortran `MODULE` → C header + implementation
- Fortran arrays → C arrays (0-indexed)
- Fortran strings → C char arrays
- Fortran `LOGICAL` → C `bool` (stdbool.h)
- Fortran `REAL` → C `double` (default)
- Fortran `INTEGER` → C `int`

## Modules Remaining

Approximately 73 modules still need conversion, including:
- **High priority**: database, dictstackdb, searching, mapping
- **Chemistry modules**: Various chemistry calculation modules
- **I/O modules**: File reading/writing utilities
- **Main program**: Full mechanism generation logic

## Python Bindings (EXTRA - In Progress)

Python bindings using the Python C API will be added to enable:
- Direct calls to C functions from Python
- Type-safe interfaces with proper error handling
- Shared library compilation target

See `python/` directory for Python wrapper development.

## Testing

Currently includes basic functional demonstration in main.c.
Comprehensive test suite to be added as more modules are converted.
