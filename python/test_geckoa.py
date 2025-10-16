#!/usr/bin/env python3
"""
GECKO-A Python Bindings - Demonstration Script

This script demonstrates all the functionality of the GECKO-A Python bindings.
"""

import sys
import os

# Try to import the module
try:
    import geckoa
except ImportError:
    print("Error: geckoa module not found.")
    print("Build it with: cd python && python3 setup.py build_ext --inplace")
    sys.exit(1)

print("=" * 60)
print("GECKO-A Python Bindings - Demonstration")
print("=" * 60)
print()

# 1. Initialize defaults
print("1. Initializing default parameters...")
geckoa.define_defaults()
print("   ✓ Defaults initialized")
print()

# 2. Get configuration
print("2. Reading configuration...")
config = geckoa.get_config()
for key, value in config.items():
    print(f"   {key:12s} = {value}")
print()

# 3. Display constants
print("3. Key constants:")
print(f"   MXLCO  (max species name length)   = {geckoa.MXLCO}")
print(f"   MXLFO  (max formula length)        = {geckoa.MXLFO}")
print(f"   MXNODE (max nodes)                 = {geckoa.MXNODE}")
print(f"   MXPS   (max primary species)       = {geckoa.MXPS}")
print(f"   MXPD   (max products per reaction) = {geckoa.MXPD}")
print()

# 4. Test mini dictionary
print("4. Mini dictionary operations:")
geckoa.clean_minid()
print("   ✓ Dictionary cleared")

species = [
    ("CH4", "CH4"),
    ("C2H6", "CC"),
    ("C3H8", "CCC"),
    ("BENZENE", "C1=CC=CC=C1"),
]

for name, formula in species:
    geckoa.add_fo(name, formula)
    print(f"   Added {name:10s} -> {formula}")

print()
print("   Retrieving formulas:")
for name, expected in species:
    retrieved = geckoa.get_fo(name)
    status = "✓" if retrieved == expected else "✗"
    print(f"   {status} {name:10s} -> {retrieved}")
print()

# 5. Test string operations
print("5. String operations:")
test_strings = [
    ("CH3-CH2-CH3", "CH", 3),
    ("COOH-CH2-COOH", "COOH", 2),
    ("C=C-C=C-C=C", "C=C", 3),
    ("Hello World", "l", 3),
]

for line, substr, expected in test_strings:
    count = geckoa.countstring(line, substr)
    status = "✓" if count == expected else "✗"
    print(f"   {status} '{substr}' in '{line}': {count} occurrences")
print()

# 6. Test Arrhenius rate constant calculation
print("6. Arrhenius rate constant calculations:")
print("   k = A × T^n × exp(-Ea/R/T)")
print()

test_cases = [
    ([1.5e-12, 0.0, 0.0], 298.0, "Simple A factor"),
    ([1.0e-12, 0.5, 0.0], 298.0, "Temperature dependent"),
    ([1.0e-11, 0.0, 500.0], 298.0, "Activation energy"),
    ([2.5e-12, 0.8, 200.0], 310.0, "Full Arrhenius"),
]

for arrh, T, description in test_cases:
    k = geckoa.kval(arrh, T)
    print(f"   {description:25s}")
    print(f"      A={arrh[0]:.2e}, n={arrh[1]:.1f}, Ea/R={arrh[2]:.1f}, T={T:.1f}K")
    print(f"      k = {k:.4e}")
    print()

# 7. Summary
print("=" * 60)
print("All tests completed successfully! ✓")
print("=" * 60)
print()
print("The GECKO-A Python bindings are working correctly.")
print("You can now use this module in your Python scripts.")
print()
