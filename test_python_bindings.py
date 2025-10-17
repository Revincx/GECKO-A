#!/usr/bin/env python3
"""
Test script for GECKO-A Python bindings

This script tests basic functionality of the Python bindings
without requiring full data files.
"""

import os
import sys

# Add the LIB directory to Python path
lib_path = os.path.join(os.path.dirname(__file__), 'LIB')
sys.path.insert(0, lib_path)

def test_module_loading():
    """Test that the module can be loaded"""
    print("Test 1: Loading Python module...")
    try:
        from gecko_python import GeckoA
        print("✓ Module loaded successfully")
        return True
    except Exception as e:
        print(f"✗ Failed to load module: {e}")
        return False

def test_library_loading():
    """Test that the shared library can be loaded"""
    print("\nTest 2: Loading shared library...")
    try:
        from gecko_python import libgecko
        print("✓ Shared library loaded successfully")
        return True
    except Exception as e:
        print(f"✗ Failed to load shared library: {e}")
        return False

def test_basic_functions():
    """Test that basic functions can be called"""
    print("\nTest 3: Calling basic functions...")
    try:
        from gecko_python import GeckoA
        gecko = GeckoA()
        
        # Test a simple function that doesn't require initialization
        print("  Testing define_defaults()...")
        gecko.define_defaults()
        print("  ✓ define_defaults() executed")
        
        # Test initialization
        print("  Testing initialize()...")
        gecko.initialize()
        print("  ✓ initialize() executed")
        
        # Test status query
        print("  Testing get_nrec()...")
        nrec = gecko.get_nrec()
        print(f"  ✓ get_nrec() returned: {nrec}")
        
        return True
    except Exception as e:
        print(f"✗ Failed during function calls: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_object_creation():
    """Test that GeckoA object can be created"""
    print("\nTest 4: Creating GeckoA object...")
    try:
        from gecko_python import GeckoA
        gecko = GeckoA()
        print("✓ GeckoA object created successfully")
        return True
    except Exception as e:
        print(f"✗ Failed to create GeckoA object: {e}")
        return False

def main():
    """Run all tests"""
    print("=" * 70)
    print("GECKO-A Python Bindings Test Suite")
    print("=" * 70)
    
    tests = [
        test_module_loading,
        test_library_loading,
        test_object_creation,
        test_basic_functions,
    ]
    
    results = []
    for test in tests:
        results.append(test())
    
    print("\n" + "=" * 70)
    print(f"Test Results: {sum(results)}/{len(results)} passed")
    print("=" * 70)
    
    if all(results):
        print("\n✓ All tests passed!")
        return 0
    else:
        print("\n✗ Some tests failed.")
        return 1

if __name__ == '__main__':
    sys.exit(main())
