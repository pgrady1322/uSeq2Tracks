#!/usr/bin/env python3
"""
Diagnostic script to check uSeq2Tracks path generation
Run this on your server to see what paths the pipeline is actually generating
"""

import pandas as pd
import yaml
import sys
import os

print("🔍 uSeq2Tracks Path Diagnostic")
print("=" * 50)

# Check if files exist
files_to_check = ['config.yaml', 'Snakefile', 'sample_sheet.csv']
for f in files_to_check:
    exists = os.path.exists(f)
    print(f"File {f}: {'✅ EXISTS' if exists else '❌ MISSING'}")
    if exists:
        print(f"  Size: {os.path.getsize(f)} bytes")

print()

try:
    # Load config
    with open('config.yaml', 'r') as f:
        config = yaml.safe_load(f)
    
    print("📋 Configuration:")
    print(f"  genome_id: {repr(config.get('genome_id'))}")
    print(f"  outdir: {repr(config.get('outdir'))}")
    print(f"  samplesheet: {repr(config.get('samplesheet'))}")
    print(f"  rapid_mode: {config.get('rapid_mode', False)}")
    
    # Calculate paths like Snakefile does
    GENOME_ID = config.get('genome_id', 'Unknown')
    BASE_OUTDIR = config["outdir"]
    GENOME_OUTDIR = f"{BASE_OUTDIR}/{GENOME_ID}"
    RAPID_MODE = config.get('rapid_mode', False)
    
    print()
    print("🎯 Expected Paths:")
    print(f"  BASE_OUTDIR: {BASE_OUTDIR}")
    print(f"  GENOME_OUTDIR: {GENOME_OUTDIR}")
    print(f"  RAPID_MODE: {RAPID_MODE}")
    
    # Load samples
    samples_df = pd.read_csv(config['samplesheet'])
    print()
    print("📊 Sample Data:")
    print(f"  Total samples: {len(samples_df)}")
    print(f"  Sample types: {samples_df['type'].value_counts().to_dict()}")
    
    # Check sample types like pipeline does
    SAMPLES = samples_df.set_index('sample_id').to_dict('index')
    SAMPLE_TYPES = {}
    for sample_id, sample_info in SAMPLES.items():
        sample_type = sample_info.get('type', '').lower()
        if sample_type not in SAMPLE_TYPES:
            SAMPLE_TYPES[sample_type] = []
        SAMPLE_TYPES[sample_type].append(sample_id)

    RNASEQ_SAMPLES = SAMPLE_TYPES.get('rnaseq', [])
    print(f"  RNA-seq samples: {len(RNASEQ_SAMPLES)}")
    
    # Show expected target paths
    print()
    print("🎯 Expected Target Files:")
    if RNASEQ_SAMPLES:
        target = f"{GENOME_OUTDIR}/rnaseq/tracks_complete.txt"
        print(f"  RNA-seq completion: {target}")
    
    target = f"{GENOME_OUTDIR}/ucsc/hub.txt"
    print(f"  UCSC hub: {target}")
    
    if not RAPID_MODE:
        target = f"{GENOME_OUTDIR}/qc/multiqc_report.html"
        print(f"  MultiQC report: {target}")
    
    # Check a few sample files to see what paths they'd generate
    print()
    print("📁 Sample File Path Examples:")
    for sample in list(RNASEQ_SAMPLES)[:3]:
        print(f"  Sample: {sample}")
        if not RAPID_MODE:
            # These paths should NOT be in rule all anymore, but show what they would be
            fastqc_path = f"{GENOME_OUTDIR}/qc/fastqc/{sample}_fastqc.html"
            print(f"    FastQC: {fastqc_path}")
        
        bam_path = f"{GENOME_OUTDIR}/rnaseq/bam/{sample}.Aligned.sortedByCoord.out.bam"
        print(f"    BAM: {bam_path}")

except Exception as e:
    print(f"❌ Error: {e}")
    sys.exit(1)

print()
print("✅ Diagnostic complete!")
print()
print("🔧 If paths look wrong:")
print("1. Check genome_id in config.yaml")
print("2. Check outdir path in config.yaml") 
print("3. Make sure you're running from the correct directory")
