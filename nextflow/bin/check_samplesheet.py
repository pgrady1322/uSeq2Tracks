#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
uSeq2Tracks v1.0.0

Samplesheet validation for the Nextflow pipeline.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

import sys
import csv
import argparse
from pathlib import Path


def parse_args(args=None):
    parser = argparse.ArgumentParser(
        description="Validate uSeq2Tracks samplesheet and check for required fields"
    )
    parser.add_argument("samplesheet", type=Path, help="Input samplesheet file")
    parser.add_argument("output", type=Path, help="Output validated samplesheet")
    return parser.parse_args(args)


def is_valid_value(value):
    """Check if value is not empty or NA"""
    if value is None:
        return False
    value_str = str(value).strip().upper()
    return value_str not in ['', 'NA', 'NAN', 'NONE', 'NULL', '<NA>', 'N/A']


def check_samplesheet(samplesheet_path, output_path):
    """
    Validate samplesheet format and contents
    """
    required_cols = ['sample_id', 'type']
    valid_types = ['atacseq', 'chipseq', 'cutrun', 'rnaseq', 'wgs', 
                   'nanopore', 'pacbio', 'ancientdna']
    
    seen_ids = set()
    samples = []
    
    with open(samplesheet_path, 'r') as f:
        reader = csv.DictReader(f)
        
        # Check required columns
        if not all(col in reader.fieldnames for col in required_cols):
            missing = [col for col in required_cols if col not in reader.fieldnames]
            print(f"ERROR: Missing required columns: {', '.join(missing)}", file=sys.stderr)
            sys.exit(1)
        
        for i, row in enumerate(reader, start=2):
            # Check for empty sample_id
            if not is_valid_value(row.get('sample_id')):
                print(f"ERROR: Line {i}: sample_id is empty or invalid", file=sys.stderr)
                sys.exit(1)
            
            sample_id = row['sample_id'].strip()
            
            # Check for duplicate sample_ids
            if sample_id in seen_ids:
                print(f"ERROR: Line {i}: Duplicate sample_id '{sample_id}'", file=sys.stderr)
                sys.exit(1)
            seen_ids.add(sample_id)
            
            # Check assay type
            assay_type = row.get('type', '').strip().lower()
            if assay_type not in valid_types:
                print(f"ERROR: Line {i}: Invalid type '{assay_type}'. " 
                      f"Must be one of: {', '.join(valid_types)}", file=sys.stderr)
                sys.exit(1)
            
            # Check data source
            has_sra = is_valid_value(row.get('sra_id'))
            has_local = is_valid_value(row.get('read1'))
            
            if not has_sra and not has_local:
                print(f"ERROR: Line {i}: Sample '{sample_id}' has neither SRA ID "
                      "nor local files (read1)", file=sys.stderr)
                sys.exit(1)
            
            # Check read2 exists if read1 exists (for paired-end)
            if has_local:
                read1 = Path(row['read1'])
                if not read1.exists():
                    print(f"WARNING: Line {i}: read1 file does not exist: {read1}", 
                          file=sys.stderr)
                
                if is_valid_value(row.get('read2')):
                    read2 = Path(row['read2'])
                    if not read2.exists():
                        print(f"WARNING: Line {i}: read2 file does not exist: {read2}", 
                              file=sys.stderr)
            
            samples.append(row)
    
    # Write validated samplesheet
    if samples:
        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=reader.fieldnames)
            writer.writeheader()
            writer.writerows(samples)
        
        print(f"✓ Validated {len(samples)} samples")
        print(f"  Assay types: {', '.join(set(s['type'] for s in samples))}")
    else:
        print("ERROR: No valid samples found in samplesheet", file=sys.stderr)
        sys.exit(1)


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.samplesheet, args.output)


if __name__ == '__main__':
    main()

# uSeq2Tracks v1.0.0
# Any usage is subject to this software's license.
