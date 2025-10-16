# GECKO-A Python Bindings

Python interface to the GECKO-A C library.

## Installation

### Build from source

```bash
cd python
python3 setup.py build_ext --inplace
```

Or install:

```bash
cd python
pip install .
```

## Usage

```python
import geckoa

# Initialize default parameters
geckoa.define_defaults()

# Get configuration
config = geckoa.get_config()
print(f"Temperature: {config['TK']} K")
print(f"Max generations: {config['maxgen']}")

# Use mini dictionary
geckoa.clean_minid()
geckoa.add_fo("CH4", "CH4")
geckoa.add_fo("C2H6", "CC")

formula = geckoa.get_fo("CH4")
print(f"Formula for CH4: {formula}")

# Count substring occurrences
line = "CH3-CH2-CH3"
count = geckoa.countstring(line, "CH")
print(f"'CH' appears {count} times in '{line}'")

# Calculate Arrhenius rate constant
# k = A * T^n * exp(-Ea/R/T)
# arrh = [A, n, Ea/R]
k = geckoa.kval([1.5e-12, 0.0, 0.0], 298.0)
print(f"Rate constant at 298K: {k}")
```

## Available Functions

### Configuration
- `define_defaults()` - Initialize default parameters
- `get_config()` - Get configuration dictionary

### Mini Dictionary
- `clean_minid()` - Clear the mini dictionary
- `add_fo(name, formula)` - Add a species formula
- `get_fo(name)` - Get formula for a species

### Utilities
- `countstring(line, substring)` - Count substring occurrences
- `kval(arrh, T)` - Calculate Arrhenius rate constant

### Constants
- `MXLCO` - Maximum species name length (6)
- `MXLFO` - Maximum formula length (120)
- `MXNODE` - Maximum nodes (35)
- `MXPS` - Maximum primary species (600)
- `MXPD` - Maximum products per reaction (38)

## Type Signatures

```python
def define_defaults() -> None: ...
def get_config() -> dict[str, Union[str, int, float]]: ...
def clean_minid() -> None: ...
def add_fo(name: str, formula: str) -> None: ...
def get_fo(name: str) -> str: ...
def countstring(line: str, substring: str) -> int: ...
def kval(arrh: list[float], T: float) -> float: ...
```

## Testing

```bash
cd python
python3 -c "import geckoa; geckoa.define_defaults(); print(geckoa.get_config())"
```
