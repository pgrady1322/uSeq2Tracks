# ============================================================
# uSeq2Tracks – Snakefile (Snakemake)  •  v0.5.0
# ============================================================
# Universal Sequencing-to-Genome Pipeline
# Converts heterogeneous sequencing data to UCSC-ready tracks

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd

# ───────── CONFIGURATION ─────────
configfile: "config.yaml"

# Validate config against JSON-schema
validate(config, schema="schemas/config.schema.yaml")

# ───────── CONFIGURATION VALIDATION ─────────
# Ensure required genome_id is set
if not config.get("genome_id") or config["genome_id"].strip() == "":
    raise ValueError(
        "ERROR: 'genome_id' is required in config.yaml but not set!\n"
        "Please add a genome identifier, e.g.:\n"
        "  genome_id: 'galGal6'\n"
        "This will be used to tag all output files and directories."
    )

# Clean genome_id (remove problematic / shell-unsafe characters)
import re as _re
GENOME_ID = _re.sub(r"[^\w.-]", "_", config["genome_id"].strip())
print(f"Pipeline running with genome ID: {GENOME_ID}")

# Check for rapid mode
RAPID_MODE = config.get("rapid_mode", False)
if RAPID_MODE:
    print("🚀 RAPID MODE ENABLED: Skipping QC, generating essential tracks only")

# Update output directory to include genome ID
BASE_OUTDIR = config["outdir"]
GENOME_OUTDIR = f"{BASE_OUTDIR}/{GENOME_ID}"
print(f"Output directory: {GENOME_OUTDIR}")

# Set rapid mode output directory (always define for consistency)
RAPID_OUTDIR = f"{GENOME_OUTDIR}/rapid"
if RAPID_MODE:
    print(f"Rapid mode output: {RAPID_OUTDIR}")

# ───────── SAMPLE SHEET PARSING ─────────
samples_df = pd.read_csv(config["samplesheet"], comment='#')  # Skip comment lines starting with #
samples_df = samples_df.dropna(subset=['sample_id'])  # Remove rows with missing sample_id
samples_df["outdir"] = GENOME_OUTDIR

# Backwards compatibility: fill missing optional columns
if 'experiment_group' not in samples_df.columns:
    samples_df['experiment_group'] = 'exp1'  # Default experiment group

if 'replicate_group' not in samples_df.columns:
    # Use condition if available, otherwise use sample type
    if 'condition' in samples_df.columns:
        samples_df['replicate_group'] = samples_df['condition']
    else:
        samples_df['replicate_group'] = samples_df['type']

# No layout column needed - all SRA samples processed as paired-end with auto-detection

# Create dictionaries for easy access
SAMPLES = samples_df.set_index("sample_id").to_dict('index')
SAMPLE_TYPES = samples_df.groupby("type")["sample_id"].apply(list).to_dict()

# Create replicate groups for composite tracks
REPLICATE_GROUPS = {}
EXPERIMENT_GROUPS = {}
if 'replicate_group' in samples_df.columns:
    for sample_type in samples_df['type'].unique():
        type_samples = samples_df[samples_df['type'] == sample_type]
        
        # Group by experiment first, then by replicate group
        if 'experiment_group' in samples_df.columns:
            EXPERIMENT_GROUPS[sample_type] = {}
            for exp_group in type_samples['experiment_group'].unique():
                exp_samples = type_samples[type_samples['experiment_group'] == exp_group]
                EXPERIMENT_GROUPS[sample_type][exp_group] = exp_samples.groupby('replicate_group')['sample_id'].apply(list).to_dict()
        
        # Legacy: simple replicate grouping (all experiments together)
        REPLICATE_GROUPS[sample_type] = type_samples.groupby('replicate_group')['sample_id'].apply(list).to_dict()

# ───────── CONFIG DEFAULTS ─────────
# Data-driven defaults per assay — avoids 50+ repetitive setdefault lines.
_MAPPING_THREADS = config.get('resources', {}).get('mapping_threads', 8)

_ASSAY_DEFAULTS: dict[str, dict] = {
    "atacseq":    {"threads": _MAPPING_THREADS, "bowtie2_opts": "--very-sensitive",
                   "bwa_mem2_opts": "-M", "shift": -75, "extsize": 150,
                   "macs3_opts": "", "bw_norm": "CPM"},
    "wgs":        {"threads": _MAPPING_THREADS, "mapper": "bwa_mem2",
                   "bowtie2_opts": "--very-sensitive", "bwa_mem2_opts": "-M",
                   "bw_norm": "CPM"},
    "cutrun":     {"threads": _MAPPING_THREADS, "bowtie2_opts": "--very-sensitive",
                   "bwa_mem2_opts": "-M", "shift": 0, "extsize": 160,
                   "macs3_opts": "", "bw_norm": "CPM"},
    "rnaseq":     {"threads": _MAPPING_THREADS,
                   "star_opts": "--outFilterMultimapNmax 20", "bw_norm": "CPM"},
    "nanopore":   {"threads": _MAPPING_THREADS,
                   "minimap2_opts": "--secondary=no", "bw_norm": "CPM"},
    "pacbio":     {"threads": _MAPPING_THREADS,
                   "minimap2_opts": "--secondary=no", "bw_norm": "CPM"},
    "ancientdna": {"threads": _MAPPING_THREADS, "mapper": "bwa_aln",
                   "bwa_aln_opts": "-l 1024 -n 0.01 -o 2",
                   "bwa_mem_opts": "-M", "min_mapq": 30, "min_baseq": 20,
                   "bw_norm": "CPM"},
}

for _assay, _defaults in _ASSAY_DEFAULTS.items():
    config.setdefault(_assay, {})
    for _key, _val in _defaults.items():
        config[_assay].setdefault(_key, _val)

# Get samples by type (with empty defaults)
WGS_SAMPLES = SAMPLE_TYPES.get("wgs", [])
CHIPSEQ_SAMPLES = SAMPLE_TYPES.get("chipseq", [])
CUTRUN_SAMPLES = SAMPLE_TYPES.get("cutrun", [])
ATACSEQ_SAMPLES = SAMPLE_TYPES.get("atacseq", [])
RNASEQ_SAMPLES = SAMPLE_TYPES.get("rnaseq", [])
NANOPORE_SAMPLES = SAMPLE_TYPES.get("nanopore", [])
PACBIO_SAMPLES = SAMPLE_TYPES.get("pacbio", [])
ADNA_SAMPLES = SAMPLE_TYPES.get("ancientdna", [])

# ───────── INCLUDE RULES ─────────
include: "rules/common.smk"
include: "rules/genome_prep.smk"
include: "rules/fastqc.smk"
include: "rules/adapter_trimming.smk"
include: "rules/sra_download.smk"

# Include workflow-specific rules only if samples exist
if WGS_SAMPLES:
    include: "rules/wgs.smk"
if CHIPSEQ_SAMPLES:
    include: "rules/chipseq.smk"
if CUTRUN_SAMPLES:
    include: "rules/cutrun.smk"
if ATACSEQ_SAMPLES:
    include: "rules/atacseq.smk"
if RNASEQ_SAMPLES:
    include: "rules/rnaseq.smk"
if NANOPORE_SAMPLES:
    include: "rules/nanopore.smk"
if PACBIO_SAMPLES:
    include: "rules/pacbio.smk"
if ADNA_SAMPLES:
    include: "rules/ancientdna.smk"

include: "rules/ucsc_tracks.smk"
include: "rules/multiqc.smk"
include: "rules/replicate_merging.smk"
include: "rules/parameter_sweep.smk"
include: "rules/genrich.smk"

# NOTE: is_valid_sra_id() and has_local_files() are defined in rules/common.smk

# ───────── PARAMETER SWEEP FUNCTIONS ─────────
def get_parameter_sweep_outputs() -> list[str]:
    """Return peak file paths for every (sample, q-value) combination.

    Returns:
        List of expected narrowPeak file paths (empty when sweeps are
        disabled in the config).
    """
    if not config.get('parameter_sweep', {}).get('enabled', False):
        return []
    
    qvalues = config.get('parameter_sweep', {}).get('qvalues', [0.05, 0.01, 0.005, 0.001])
    outputs = []
    
    # Convert q-values to strings for file naming (replace dots with underscores)  
    qval_strings = [str(q).replace('.', '') for q in qvalues]
    
    # ATAC-seq parameter sweep outputs
    if ATACSEQ_SAMPLES:
        for sample in ATACSEQ_SAMPLES:
            for qval in qval_strings:
                outputs.append(f"{GENOME_OUTDIR}/atacseq/peaks_sweep/{sample}_q{qval}_peaks.narrowPeak")
    
    # CUT&RUN parameter sweep outputs (excluding control samples)
    if CUTRUN_SAMPLES:
        cutrun_peak_samples = [s for s in CUTRUN_SAMPLES if not any(keyword in s.lower() for keyword in ['igg', 'input', 'control'])]
        for sample in cutrun_peak_samples:
            for qval in qval_strings:
                outputs.append(f"{GENOME_OUTDIR}/cutrun/peaks_sweep/{sample}_q{qval}_peaks.narrowPeak")
    
    # ChIP-seq parameter sweep outputs
    if CHIPSEQ_SAMPLES:
        for sample in CHIPSEQ_SAMPLES:
            for qval in qval_strings:
                outputs.append(f"{GENOME_OUTDIR}/chipseq/peaks_sweep/{sample}_q{qval}_peaks.narrowPeak")
    
    return outputs


def get_genrich_outputs() -> list[str]:
    """Return Genrich peak-calling output paths.

    Returns:
        List of expected Genrich output paths (empty when Genrich is
        disabled in the config).
    """
    outputs = []
    
    # Standard Genrich outputs
    if config.get('genrich', {}).get('enabled', False):
        outputs.append(f"{GENOME_OUTDIR}/genrich/genrich_complete.txt")
    
    # Genrich parameter sweep outputs
    if config.get('parameter_sweep', {}).get('genrich_enabled', False):
        qvalues = config.get('parameter_sweep', {}).get('genrich_qvalues', [0.05, 0.01, 0.005, 0.001])
        
        for assay in ['atacseq', 'cutrun', 'chipseq']:
            if assay not in REPLICATE_GROUPS:
                continue
                
            _type_map = {"atacseq": ATACSEQ_SAMPLES, "cutrun": CUTRUN_SAMPLES, "chipseq": CHIPSEQ_SAMPLES}
            assay_samples = _type_map.get(assay, [])
            if not assay_samples:
                continue
                
            for group in REPLICATE_GROUPS[assay].keys():
                for qval in qvalues:
                    # Convert qvalue to filename format (0.05 -> q005)
                    qval_str = f"q{str(qval).replace('0.', '').zfill(3)}"
                    outputs.append(f"{GENOME_OUTDIR}/{assay}/genrich_sweep/{group}_genrich_{qval_str}.narrowPeak")
        
        outputs.append(f"{GENOME_OUTDIR}/genrich_sweep/genrich_sweep_summary.txt")
    
    return outputs

# ───────── SRA DOWNLOAD FUNCTIONS AND RULES ─────────
def get_sra_download_inputs() -> list[str]:
    """Return expected FASTQ paths for all SRA-sourced samples.

    Returns:
        List of ``data/raw/{sample}_R{1,2}.fastq.gz`` paths for every
        sample that uses SRA and lacks local files.
    """
    inputs = []
    for sample in SAMPLES.keys():
        sample_info = SAMPLES[sample]
        sra_id = sample_info.get('sra_id', '')
        
        # Only include if we have a valid SRA ID AND no local files
        if is_valid_sra_id(sra_id) and not has_local_files(sample_info):
            # All SRA samples use paired-end naming (auto-detection handles actual layout)
            inputs.extend([
                f"data/raw/{sample}_R1.fastq.gz",
                f"data/raw/{sample}_R2.fastq.gz"
            ])
    return inputs

# ───────── TARGET RULE ─────────
rule all:
    input:
        # SRA downloads (if any samples use SRA)
        (get_sra_download_inputs() if any(
            is_valid_sra_id(SAMPLES[s].get('sra_id', '')) and not has_local_files(SAMPLES[s])
            for s in SAMPLES.keys()
        ) else []) +
        
        # Essential pipeline completion files (these will trigger all dependencies)
        ([f"{GENOME_OUTDIR}/wgs/tracks_complete.txt"] if WGS_SAMPLES else []) +
        ([f"{GENOME_OUTDIR}/chipseq/tracks_complete.txt"] if CHIPSEQ_SAMPLES else []) +
        ([f"{GENOME_OUTDIR}/cutrun/tracks_complete.txt"] if CUTRUN_SAMPLES else []) +
        ([f"{GENOME_OUTDIR}/atacseq/tracks_complete.txt"] if ATACSEQ_SAMPLES else []) +
        ([f"{GENOME_OUTDIR}/rnaseq/tracks_complete.txt"] if RNASEQ_SAMPLES else []) +
        ([f"{GENOME_OUTDIR}/nanopore/tracks_complete.txt"] if NANOPORE_SAMPLES else []) +
        ([f"{GENOME_OUTDIR}/pacbio/tracks_complete.txt"] if PACBIO_SAMPLES else []) +
        ([f"{GENOME_OUTDIR}/ancientdna/tracks_complete.txt"] if ADNA_SAMPLES else []) +
        
        # UCSC tracks (both modes)
        [f"{GENOME_OUTDIR}/ucsc/hub.txt"] +
        
        # MultiQC report (skip in rapid mode) - this will trigger FastQC and other QC steps
        ([] if RAPID_MODE else [f"{GENOME_OUTDIR}/qc/multiqc_report.html"]) +
        
        # Rapid mode completion marker
        ([f"{GENOME_OUTDIR}/rapid/rapid_tracks_complete.txt"] if RAPID_MODE else [])

# ───────── RAPID MODE COMPLETION ─────────
rule rapid_tracks_complete:
    """Creates a summary of essential tracks for rapid mode"""
    input:
        lambda wildcards: [
            f"{GENOME_OUTDIR}/ucsc/hub.txt"
        ] + ([f"{GENOME_OUTDIR}/wgs/tracks_complete.txt"] if WGS_SAMPLES else []) + \
          ([f"{GENOME_OUTDIR}/chipseq/tracks_complete.txt"] if CHIPSEQ_SAMPLES else []) + \
          ([f"{GENOME_OUTDIR}/cutrun/tracks_complete.txt"] if CUTRUN_SAMPLES else []) + \
          ([f"{GENOME_OUTDIR}/atacseq/tracks_complete.txt"] if ATACSEQ_SAMPLES else []) + \
          ([f"{GENOME_OUTDIR}/rnaseq/tracks_complete.txt"] if RNASEQ_SAMPLES else []) + \
          ([f"{GENOME_OUTDIR}/nanopore/tracks_complete.txt"] if NANOPORE_SAMPLES else []) + \
          ([f"{GENOME_OUTDIR}/pacbio/tracks_complete.txt"] if PACBIO_SAMPLES else []) + \
          ([f"{GENOME_OUTDIR}/ancientdna/tracks_complete.txt"] if ADNA_SAMPLES else [])
    output:
        f"{GENOME_OUTDIR}/rapid/rapid_tracks_complete.txt"
    run:
        import os
        os.makedirs(f"{GENOME_OUTDIR}/rapid", exist_ok=True)
        with open(output[0], 'w') as f:
            f.write("🚀 Rapid mode track generation complete!\n")
            f.write(f"📊 Essential tracks available in: {GENOME_OUTDIR}\n")
            f.write(f"🌐 UCSC Hub: {GENOME_OUTDIR}/ucsc/hub.txt\n")
            f.write(f"📁 Rapid mode summary directory: {GENOME_OUTDIR}/rapid/\n")

# ───────── WORKFLOW VALIDATION ─────────
def validate_configuration() -> None:
    """Validate pipeline configuration and print a data-source summary.

    Raises:
        ValueError: When required columns are missing, sample IDs are
            duplicated, or samples lack both SRA and local data sources.
    """
    issues = []
    warnings = []
    
    # Check required columns
    required_cols = ['sample_id', 'type']
    for col in required_cols:
        if col not in samples_df.columns:
            issues.append(f"Missing required column: {col}")
    
    # Check for empty sample IDs
    if samples_df['sample_id'].isna().any():
        issues.append("Found empty sample_id values")
    
    # Check for duplicate sample IDs
    if samples_df['sample_id'].duplicated().any():
        duplicates = samples_df[samples_df['sample_id'].duplicated()]['sample_id'].tolist()
        issues.append(f"Duplicate sample_id values: {duplicates}")
    
    # Robust validation of data sources for each sample
    problematic = []
    sra_samples = []
    local_samples = []
    
    for idx, row in samples_df.iterrows():
        sample_id = row['sample_id']
        sample_info = row.to_dict()
        
        has_valid_sra = is_valid_sra_id(sample_info.get('sra_id', ''))
        has_valid_local = has_local_files(sample_info)
        
        if has_valid_sra and has_valid_local:
            warnings.append(f"Sample {sample_id} has both SRA ID and local files - will prioritize local files")
            local_samples.append(sample_id)
        elif has_valid_local:
            local_samples.append(sample_id)
        elif has_valid_sra:
            sra_samples.append(sample_id)
        else:
            problematic.append(sample_id)
    
    if problematic:
        issues.append(f"Samples missing both valid SRA ID and local files: {problematic}")
    
    # Summary information
    print(f"📊 Data source summary:")
    print(f"   • Local files: {len(local_samples)} samples")
    print(f"   • SRA data: {len(sra_samples)} samples") 
    if warnings:
        print(f"   • Mixed sources: {len([w for w in warnings if 'both SRA ID and local files' in w])} samples")
    
    if warnings:
        print("⚠️  Configuration warnings:")
        for warning in warnings:
            print(f"  • {warning}")
    
    if issues:
        error_msg = "❌ Configuration validation failed:\n" + "\n".join([f"  • {issue}" for issue in issues])
        raise ValueError(error_msg)
    else:
        print("✅ Configuration validation passed")

# Validate configuration at pipeline start
validate_configuration()

onstart:
    print("🧬 Starting uSeq2Tracks pipeline...")
    print(f"📊 Found {len(SAMPLES)} samples across {len(SAMPLE_TYPES)} assay types")
    for assay_type, sample_list in SAMPLE_TYPES.items():
        print(f"   • {assay_type}: {len(sample_list)} samples")
    print(f"📁 Output directory: {GENOME_OUTDIR}")
    print(f"🧬 Reference genome: {config['genome']}")

onsuccess:
    print("✅ uSeq2Tracks pipeline completed successfully!")
    print(f"📊 UCSC hub file: {GENOME_OUTDIR}/ucsc/hub.txt")
    print(f"📊 MultiQC report: {GENOME_OUTDIR}/qc/multiqc_report.html")

onerror:
    print("❌ uSeq2Tracks pipeline failed!")
    print("📋 Check the log files in the output directory for details.")
