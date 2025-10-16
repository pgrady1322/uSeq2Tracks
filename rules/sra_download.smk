# ============================================================
# SRA download rules - assume paired-end with auto-detection
# ============================================================

# All SRA samples go through paired-end download with fallback to single-end
rule download_sra_paired:
    """Download SRA data assuming paired-end with automatic single-end detection."""
    output:
        r1="data/raw/{sample}_R1.fastq.gz",
        r2="data/raw/{sample}_R2.fastq.gz"
    params:
        sra_id=lambda wildcards: pd.read_csv(config["samplesheet"], comment='#').set_index('sample_id').loc[wildcards.sample, 'sra_id']
    threads: 4
    resources:
        mem_mb=8000,
        runtime=240
    run:
        import pandas as pd
        
        # Get the SRA ID and validate it robustly
        sra_id = params.sra_id
        
        # Use the same robust validation as in other functions
        def is_valid_sra_id(sra_id):
            if sra_id is None:
                return False
            if pd.isna(sra_id):
                return False
            sra_str = str(sra_id).strip()
            invalid_values = {'', 'nan', 'NaN', 'none', 'None', 'null', 'NULL', '<NA>', 'N/A', 'n/a'}
            return sra_str not in invalid_values and len(sra_str) > 0
        
        if not is_valid_sra_id(sra_id):
            raise ValueError(f"Sample {wildcards.sample} has no valid SRA ID but SRA download rule was triggered. "
                           f"SRA ID value: {repr(sra_id)}. "
                           f"This indicates a configuration error - samples with local files should not trigger SRA download.")
        
        shell("""
        # Create output directory
        mkdir -p $(dirname {output.r1})
        
        # Clean up any existing compressed files from previous runs
        rm -f $(dirname {output.r1})/{params.sra_id}*.fastq.gz
        
        # Download and convert SRA data
        fasterq-dump --threads {threads} {params.sra_id} -O $(dirname {output.r1})
        
        # Check what files were actually created
        if [[ -f $(dirname {output.r1})/{params.sra_id}_1.fastq && -f $(dirname {output.r1})/{params.sra_id}_2.fastq ]]; then
            # True paired-end data
            echo "✅ Found paired-end files for {params.sra_id}"
            gzip $(dirname {output.r1})/{params.sra_id}_1.fastq
            gzip $(dirname {output.r1})/{params.sra_id}_2.fastq
            mv $(dirname {output.r1})/{params.sra_id}_1.fastq.gz {output.r1}
            mv $(dirname {output.r1})/{params.sra_id}_2.fastq.gz {output.r2}
        elif [[ -f $(dirname {output.r1})/{params.sra_id}.fastq ]]; then
            # Actually single-end data - create a special marker for downstream processing
            echo "⚠️  WARNING: {params.sra_id} appears to be single-end - auto-detected and handling appropriately"
            echo "⚠️  Creating placeholder R2 file for downstream auto-detection"
            gzip $(dirname {output.r1})/{params.sra_id}.fastq
            mv $(dirname {output.r1})/{params.sra_id}.fastq.gz {output.r1}
            # Create R2 file with special marker indicating it's actually single-end
            echo "# SINGLE_END_PLACEHOLDER" > {output.r2}
        else
            echo "❌ No FASTQ files found for {params.sra_id}"
            exit 1
        fi
        """)
