# ============================================================
# Adapter trimming rules using fastp
# Trims adapters from Illumina reads for all data types
# ============================================================

# Trim adapters from raw reads using fastp
rule trim_adapters:
    input:
        unpack(get_sample_reads)
    output:
        r1=f"{config['outdir']}/trimmed/{{sample}}_R1_trimmed.fastq.gz",
        r2=f"{config['outdir']}/trimmed/{{sample}}_R2_trimmed.fastq.gz",
        json=f"{config['outdir']}/qc/trimming/{{sample}}_fastp_trim.json",
        html=f"{config['outdir']}/qc/trimming/{{sample}}_fastp_trim.html"
    params:
        # Adapter trimming options from config
        min_length = config.get('adapter_trimming', {}).get('min_length', 20),
        quality_cutoff = config.get('adapter_trimming', {}).get('quality_cutoff', 20),
        adapter1 = config.get('adapter_trimming', {}).get('adapter1', ''),
        adapter2 = config.get('adapter_trimming', {}).get('adapter2', ''),
        extra_opts = config.get('adapter_trimming', {}).get('extra_opts', ''),
        poly_g_min_len = config.get('adapter_trimming', {}).get('poly_g_min_len', 10)
    threads: config.get('adapter_trimming', {}).get('threads', 4)
    run:
        # Check if we have paired-end or single-end data
        has_r2 = hasattr(input, 'r2') and input.r2 and os.path.exists(input.r2)
        
        # For SRA data, check if R2 is a placeholder (single-end data)
        if has_r2 and input.r2:
            try:
                # Check if R2 file is tiny (placeholder file)
                if os.path.getsize(input.r2) < 100:
                    has_r2 = False
                else:
                    # Try to read first line to check for placeholder
                    import gzip
                    with gzip.open(input.r2, 'rt') as f:
                        first_line = f.readline().strip()
                        if first_line == "# SINGLE_END_PLACEHOLDER":
                            has_r2 = False
            except:
                pass  # If we can't read it, assume it's real paired data
        
        shell("mkdir -p $(dirname {output.r1}) $(dirname {output.json})")
        
        # Build adapter arguments
        adapter_args = ""
        if params.adapter1:
            adapter_args += f" --adapter_sequence {params.adapter1}"
        if params.adapter2 and has_r2:
            adapter_args += f" --adapter_sequence_r2 {params.adapter2}"
        
        if has_r2:
            # Paired-end trimming
            shell("""
            fastp -i {input.r1} -I {input.r2} \
                  -o {output.r1} -O {output.r2} \
                  --thread {threads} \
                  --length_required {params.min_length} \
                  --cut_mean_quality {params.quality_cutoff} \
                  --trim_poly_g \
                  --poly_g_min_len {params.poly_g_min_len} \
                  --trim_poly_x \
                  --correction \
                  --detect_adapter_for_pe \
                  {adapter_args} \
                  {params.extra_opts} \
                  --json {output.json} \
                  --html {output.html} \
                  2>/dev/null
            """)
        else:
            # Single-end trimming (create empty R2 file for consistency)
            shell("""
            fastp -i {input.r1} \
                  -o {output.r1} \
                  --thread {threads} \
                  --length_required {params.min_length} \
                  --cut_mean_quality {params.quality_cutoff} \
                  --trim_poly_g \
                  --poly_g_min_len {params.poly_g_min_len} \
                  --trim_poly_x \
                  {adapter_args} \
                  {params.extra_opts} \
                  --json {output.json} \
                  --html {output.html} \
                  2>/dev/null
            """)
            
            # Create empty R2 placeholder for single-end data
            shell("touch {output.r2}")
