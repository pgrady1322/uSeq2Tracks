# ============================================================
# Common utility functions and rules
# ============================================================

import os
import pandas as pd
from pathlib import Path

def is_valid_sra_id(sra_id):
    """
    Robust check for valid SRA ID that handles all edge cases:
    - None, NaN, empty strings, whitespace
    - Various pandas representations of missing data
    - String representations of NaN
    """
    if sra_id is None:
        return False
    
    # Handle pandas NaN values
    if pd.isna(sra_id):
        return False
        
    # Convert to string and check various empty/invalid representations
    sra_str = str(sra_id).strip()
    invalid_values = {'', 'nan', 'NaN', 'none', 'None', 'null', 'NULL', '<NA>', 'N/A', 'n/a'}
    
    return sra_str not in invalid_values and len(sra_str) > 0

def has_local_files(sample_info):
    """Check if sample has valid local file paths"""
    read1 = sample_info.get('read1', '')
    read2 = sample_info.get('read2', '')
    
    # At minimum, read1 should exist and be a valid path
    if not read1 or pd.isna(read1) or str(read1).strip() in {'', 'nan', 'NaN', 'none'}:
        return False
    
    return True

def get_sample_reads(wildcards):
    """Get input reads for a sample - auto-detect layout for SRA samples"""
    sample_info = SAMPLES[wildcards.sample]
    
    # Check data source priority: local files first, then SRA
    if has_local_files(sample_info):
        # Use local files
        reads = {"r1": sample_info["read1"]}
        if 'read2' in sample_info and sample_info['read2']:
            read2 = sample_info['read2']
            # Check if read2 is valid (not empty/nan)
            if not pd.isna(read2) and str(read2).strip() not in {'', 'nan', 'NaN', 'none'}:
                reads["r2"] = read2
        return reads
    
    elif is_valid_sra_id(sample_info.get('sra_id', '')):
        # Use SRA data - assume paired-end and auto-detect actual layout
        r2_file = f"data/raw/{wildcards.sample}_R2.fastq.gz"
        
        # Check if R2 file is a single-end placeholder (created when SRA data is actually single-end)
        if os.path.exists(r2_file):
            try:
                import gzip
                with gzip.open(r2_file, 'rt') as f:
                    first_line = f.readline().strip()
                    if first_line == "# SINGLE_END_PLACEHOLDER":
                        # Actually single-end, return only R1
                        return {"r1": f"data/raw/{wildcards.sample}_R1.fastq.gz"}
            except:
                # If R2 file is small (< 100 bytes), likely a placeholder
                if os.path.getsize(r2_file) < 100:
                    return {"r1": f"data/raw/{wildcards.sample}_R1.fastq.gz"}
                pass  # If we can't read it, assume it's a real paired file
        
        # Default: assume paired-end (will be auto-detected during download)
        return {
            "r1": f"data/raw/{wildcards.sample}_R1.fastq.gz",
            "r2": r2_file
        }
    
    else:
        raise ValueError(f"Sample {wildcards.sample} has neither valid local files nor valid SRA ID. "
                        f"read1: {sample_info.get('read1', 'missing')}, "
                        f"sra_id: {sample_info.get('sra_id', 'missing')}")

def get_genome_indices():
    """Get all genome index files"""
    genome_base = Path(config["genome"]).stem
    indices = {
        "fai": f"{GENOME_OUTDIR}/genome/{genome_base}.fa.fai",
        "dict": f"{GENOME_OUTDIR}/genome/{genome_base}.dict"
    }
    
    # Add mapper-specific indices based on what samples we have
    if WGS_SAMPLES or CHIPSEQ_SAMPLES or CUTRUN_SAMPLES or ATACSEQ_SAMPLES:
        indices["bwa_mem2"] = f"{GENOME_OUTDIR}/genome/bwa_mem2/{genome_base}.fa.0123"
        indices["bowtie2"] = f"{GENOME_OUTDIR}/genome/bowtie2/{genome_base}.1.bt2"
    
    if RNASEQ_SAMPLES:
        indices["star"] = f"{GENOME_OUTDIR}/genome/star/SA"
        indices["hisat2"] = f"{GENOME_OUTDIR}/genome/hisat2/{genome_base}.1.ht2"
    
    if NANOPORE_SAMPLES or PACBIO_SAMPLES:
        indices["minimap2"] = f"{GENOME_OUTDIR}/genome/minimap2/{genome_base}.mmi"
    
    if ADNA_SAMPLES:
        indices["bwa_aln"] = f"{GENOME_OUTDIR}/genome/bwa_aln/{genome_base}.fa.bwt"
    
    return indices

def get_samples_by_type(assay_type):
    """Get all samples of a specific assay type"""
    return [s for s, info in SAMPLES.items() if info.get('type') == assay_type]

def get_control_for_sample(wildcards):
    """Get control sample for ChIP-seq/CUT&RUN analysis"""
    sample_info = SAMPLES[wildcards.sample]
    condition = sample_info.get('condition', '')
    control_tag = config.get('chipseq', {}).get('control_tag', 'input')
    
    # Look for control sample with same condition + control_tag
    control_condition = f"{condition}_{control_tag}"
    
    for sample_id, info in SAMPLES.items():
        if info.get('condition') == control_condition:
            if config.get("chipseq", {}).get("markdup", False):
                return f"{config['outdir']}/chipseq/bam/{sample_id}.dedup.bam"
            else:
                return f"{config['outdir']}/chipseq/bam/{sample_id}.sorted.bam"
    
    # No control found
    return []

def get_fastq_for_fastqc(wildcards):
    """Get FASTQ files for FastQC analysis"""
    reads = get_sample_reads(wildcards)
    if "r2" in reads:
        return [reads["r1"], reads["r2"]]
    else:
        return [reads["r1"]]

def get_trimmed_reads(wildcards):
    """Get trimmed reads for a sample"""
    sample_info = SAMPLES[wildcards.sample]
    
    # Return trimmed reads
    reads = {"r1": f"{config['outdir']}/trimmed/{wildcards.sample}_R1_trimmed.fastq.gz"}
    
    # Check if we should expect R2 (paired-end data)
    if has_local_files(sample_info):
        # Local data - check if read2 exists
        if 'read2' in sample_info and sample_info['read2']:
            read2 = sample_info['read2']
            if not pd.isna(read2) and str(read2).strip() not in {'', 'nan', 'NaN', 'none'}:
                reads["r2"] = f"{config['outdir']}/trimmed/{wildcards.sample}_R2_trimmed.fastq.gz"
    elif is_valid_sra_id(sample_info.get('sra_id', '')):
        # SRA data - will be auto-detected during trimming
        # Always include R2 path, but actual presence depends on SRA layout
        reads["r2"] = f"{config['outdir']}/trimmed/{wildcards.sample}_R2_trimmed.fastq.gz"
    
    return reads

def get_input_reads(wildcards):
    """Smart function that returns trimmed reads if trimming is enabled, otherwise raw reads"""
    if config.get('adapter_trimming', {}).get('enabled', False):
        return get_trimmed_reads(wildcards)
    else:
        return get_sample_reads(wildcards)
