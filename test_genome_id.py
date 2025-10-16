#!/usr/bin/env python
"""Test the genome_id validation in uSeq2Tracks"""

import pandas as pd
import yaml

# Test 1: Valid genome_id
print("=== Test 1: Valid genome_id ===")
try:
    config = {"genome_id": "galGal6", "outdir": "./test_results"}
    
    if not config.get("genome_id") or config["genome_id"].strip() == "":
        raise ValueError("genome_id validation failed")
    
    GENOME_ID = config["genome_id"].strip().replace(" ", "_")
    BASE_OUTDIR = config["outdir"]
    GENOME_OUTDIR = f"{BASE_OUTDIR}/{GENOME_ID}"
    
    print(f"✓ GENOME_ID: {GENOME_ID}")
    print(f"✓ Output directory: {GENOME_OUTDIR}")
    print("✓ Test 1 PASSED")
except Exception as e:
    print(f"✗ Test 1 FAILED: {e}")

print()

# Test 2: Empty genome_id (should fail)
print("=== Test 2: Empty genome_id (should fail) ===")
try:
    config = {"genome_id": "", "outdir": "./test_results"}
    
    if not config.get("genome_id") or config["genome_id"].strip() == "":
        raise ValueError(
            "ERROR: 'genome_id' is required in config.yaml but not set!\n"
            "Please add a genome identifier, e.g.:\n"
            "  genome_id: 'galGal6'\n"
            "This will be used to tag all output files and directories."
        )
    
    print("✗ Test 2 FAILED: Should have raised an error")
except ValueError as e:
    print("✓ Test 2 PASSED: Correctly caught missing genome_id")
    print(f"✓ Error message: {str(e).split()[0]}...")

print()

# Test 3: Missing genome_id (should fail)
print("=== Test 3: Missing genome_id key (should fail) ===")
try:
    config = {"outdir": "./test_results"}  # No genome_id key
    
    if not config.get("genome_id") or config["genome_id"].strip() == "":
        raise ValueError("genome_id validation failed")
    
    print("✗ Test 3 FAILED: Should have raised an error")
except ValueError as e:
    print("✓ Test 3 PASSED: Correctly caught missing genome_id key")

print()

# Test 4: Genome ID with spaces (should be cleaned)
print("=== Test 4: Genome ID with spaces ===")
try:
    config = {"genome_id": " gal Gal6 ", "outdir": "./test_results"}
    
    if not config.get("genome_id") or config["genome_id"].strip() == "":
        raise ValueError("genome_id validation failed")
    
    GENOME_ID = config["genome_id"].strip().replace(" ", "_")
    BASE_OUTDIR = config["outdir"]
    GENOME_OUTDIR = f"{BASE_OUTDIR}/{GENOME_ID}"
    
    print(f"✓ Original: '{config['genome_id']}'")
    print(f"✓ Cleaned GENOME_ID: '{GENOME_ID}'")
    print(f"✓ Output directory: {GENOME_OUTDIR}")
    print("✓ Test 4 PASSED")
except Exception as e:
    print(f"✗ Test 4 FAILED: {e}")

print("\n=== All tests completed ===")
