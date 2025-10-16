# ============================================================
# CUT&RUN analysis rules
# CUT&RUN: Cleavage Under Targets and Release Using Nuclease
# ============================================================

def get_cutrun_peak_samples():
    """Get CUT&RUN samples that should have peaks called (exclude controls)"""
    control_keywords = ['igg', 'input', 'control']
    peak_samples = []
    
    for sample in CUTRUN_SAMPLES:
        sample_info = SAMPLES[sample]
        condition = str(sample_info.get('condition', '')).lower()
        replicate_group = str(sample_info.get('replicate_group', '')).lower()
        sample_id = str(sample).lower()
        
        # Skip if any field contains control keywords
        is_control = any(keyword in field for keyword in control_keywords 
                        for field in [condition, replicate_group, sample_id])
        
        if not is_control:
            peak_samples.append(sample)
    
    return peak_samples

# Map CUT&RUN reads with bowtie2
rule cutrun_map_bowtie2:
    input:
        unpack(get_input_reads_cutrun),
        idx = f"{GENOME_OUTDIR}/genome/bowtie2/genome.1.bt2"
    output:
        f"{config['outdir']}/cutrun/bam/{{sample}}.raw.bam"
    params:
        bowtie2_opts = config.get('cutrun', {}).get('bowtie2_opts', '--very-sensitive'),
        genome_prefix = f"{GENOME_OUTDIR}/genome/bowtie2/genome"
    threads: config.get('cutrun', {}).get('threads', config['resources']['mapping_threads'])
    run:
        # Check if we have paired-end or single-end data
        has_r2 = hasattr(input, 'r2') and input.r2
        
        shell("mkdir -p $(dirname {output})")
        
        if has_r2:
            # Paired-end
            shell("""
            bowtie2 {params.bowtie2_opts} -p {threads} -x {params.genome_prefix} \
                    -1 {input.r1} -2 {input.r2} \
              | samtools view -bS -o {output} -
            """)
        else:
            # Single-end
            shell("""
            bowtie2 {params.bowtie2_opts} -p {threads} -x {params.genome_prefix} \
                    -U {input.r1} \
              | samtools view -bS -o {output} -
            """)

# Map CUT&RUN reads with bwa-mem2
rule cutrun_map_bwa_mem2:
    input:
        unpack(get_input_reads_cutrun),
        idx = f"{GENOME_OUTDIR}/genome/bwa_mem2/genome.bwt.2bit.64"
    output:
        f"{config['outdir']}/cutrun/bam/{{sample}}.raw.bam"
    params:
        bwa_mem2_opts = config.get('cutrun', {}).get('bwa_mem2_opts', '-M'),
        genome_prefix = f"{GENOME_OUTDIR}/genome/bwa_mem2/genome"
    threads: config.get('cutrun', {}).get('threads', config['resources']['mapping_threads'])
    run:
        # Check if we have paired-end or single-end data
        has_r2 = hasattr(input, 'r2') and input.r2
        
        shell("mkdir -p $(dirname {output})")
        
        if has_r2:
            # Paired-end
            shell("""
            bwa-mem2 mem {params.bwa_mem2_opts} -t {threads} {params.genome_prefix} \
                     {input.r1} {input.r2} \
              | samtools view -bS -o {output} -
            """)
        else:
            # Single-end
            shell("""
            bwa-mem2 mem {params.bwa_mem2_opts} -t {threads} {params.genome_prefix} \
                     {input.r1} \
              | samtools view -bS -o {output} -
            """)

# Sort and index CUT&RUN BAM files
rule cutrun_sort_and_index:
    input:
        f"{config['outdir']}/cutrun/bam/{{sample}}.raw.bam"
    output:
        bam = f"{config['outdir']}/cutrun/bam/{{sample}}.sorted.bam",
        bai = f"{config['outdir']}/cutrun/bam/{{sample}}.sorted.bam.bai"
    params:
    threads: config.get('cutrun', {}).get('threads', config['resources']['mapping_threads'])
    resources:
        mem_mb=config.get('resources', {}).get('sort_memory_mb', 8000),
        tmpdir=config.get('resources', {}).get('tmpdir', '/tmp')
    shell:
        """
        # Set memory limit per thread for samtools sort (in MB)
        mem_per_thread=$((({resources.mem_mb} * 80 / 100) / {threads}))
        
        # Use temporary directory to avoid disk space issues
        mkdir -p {resources.tmpdir}/samtools_sort_{wildcards.sample}
        
        samtools sort -@ {threads} -m ${{mem_per_thread}}M \
            -T {resources.tmpdir}/samtools_sort_{wildcards.sample}/tmp \
            -o {output.bam} {input}
        
        samtools index -@ {threads} {output.bam}
        
        # Clean up temporary directory
        rm -rf {resources.tmpdir}/samtools_sort_{wildcards.sample}
        """

# Mark duplicates in CUT&RUN BAM files
rule cutrun_markdup:
    input:
        f"{config['outdir']}/cutrun/bam/{{sample}}.sorted.bam"
    output:
        bam = f"{config['outdir']}/cutrun/bam/{{sample}}.dedup.bam",
        metrics = f"{config['outdir']}/cutrun/qc/{{sample}}.markdup.metrics"
    threads: config.get('cutrun', {}).get('threads', config['resources']['mapping_threads'])
    shell:
        """
        mkdir -p $(dirname {output.metrics})
        # Mark duplicates but don't remove them (add -s flag for CUT&RUN)
        samtools markdup -@ {threads} -f {output.metrics} -s {input} {output.bam}
        samtools index {output.bam}
        """

# Call peaks with MACS3 for CUT&RUN
rule cutrun_call_peaks:
    input:
        f"{config['outdir']}/cutrun/bam/{{sample}}.dedup.bam" if config.get('cutrun', {}).get('markdup', True) else f"{config['outdir']}/cutrun/bam/{{sample}}.sorted.bam"
    output:
        f"{config['outdir']}/cutrun/peaks/{{sample}}_peaks.narrowPeak"
    wildcard_constraints:
        sample="(?!.*(?:igg|input|control)).*"
    params:
        shift = config.get('cutrun', {}).get('shift', 0),
        extsize = config.get('cutrun', {}).get('extsize', 160),
        macs3_opts = config.get('cutrun', {}).get('macs3_opts', '--qval 0.05'),
        outdir = f"{config['outdir']}/cutrun/peaks",
        prefix = "{sample}"
    shell:
        """
        mkdir -p {params.outdir}
        
        # Ensure MACS3 is available
        if ! command -v macs3 &> /dev/null; then
            echo "ERROR: macs3 command not found"
            echo "Available Python packages:"
            python3 -c "import sys; print('\\n'.join(sys.path))" 2>/dev/null || echo "Python3 not available"
            exit 1
        fi
        
        echo "Using MACS3 version: $(macs3 --version)"
        
        cd {params.outdir}
        macs3 callpeak -t {input} -n {params.prefix} \
                       --nomodel --shift {params.shift} --extsize {params.extsize} \
                       -f BAMPE --keep-dup all \
                       {params.macs3_opts}
        """

# Parameter sweep: Call peaks with multiple q-values for CUT&RUN
rule cutrun_call_peaks_sweep:
    input:
        f"{config['outdir']}/cutrun/bam/{{sample}}.dedup.bam" if config.get('cutrun', {}).get('markdup', True) else f"{config['outdir']}/cutrun/bam/{{sample}}.sorted.bam"
    output:
        f"{config['outdir']}/cutrun/peaks_sweep/{{sample}}_q{{qval}}_peaks.narrowPeak"
    wildcard_constraints:
        sample="(?!.*(?:igg|input|control)).*"
    params:
        shift = config.get('cutrun', {}).get('shift', 0),
        extsize = config.get('cutrun', {}).get('extsize', 160),
        outdir = f"{config['outdir']}/cutrun/peaks_sweep",
        prefix = "{sample}_q{qval}",
        qval = "{qval}"
    shell:
        """
        mkdir -p {params.outdir}
        
        # Ensure MACS3 is available
        if ! command -v macs3 &> /dev/null; then
            echo "ERROR: macs3 command not found"
            echo "Available Python packages:"
            python3 -c "import sys; print('\\n'.join(sys.path))" 2>/dev/null || echo "Python3 not available"
            exit 1
        fi
        
        echo "Using MACS3 version: $(macs3 --version)"
        echo "Parameter sweep: Running with q-value = {params.qval}"
        
        cd {params.outdir}
        macs3 callpeak -t {input} -n {params.prefix} \
                       --nomodel --shift {params.shift} --extsize {params.extsize} \
                       -f BAMPE --keep-dup all \
                       --qval {params.qval}
        """

# Convert CUT&RUN BAM to BigWig for visualization
rule cutrun_bam_to_bigwig:
    input:
        bam = f"{config['outdir']}/cutrun/bam/{{sample}}.dedup.bam" if config.get('cutrun', {}).get('markdup', True) else f"{config['outdir']}/cutrun/bam/{{sample}}.sorted.bam",
        fai = f"{GENOME_OUTDIR}/genome/genome.fa.fai"
    output:
        f"{config['outdir']}/cutrun/bigwig/{{sample}}.bw"
    params:
        normalize = config.get('cutrun', {}).get('bw_norm', 'CPM')
    threads: config.get('cutrun', {}).get('threads', config['resources']['mapping_threads'])
    shell:
        """
        mkdir -p $(dirname {output})
        bamCoverage -b {input.bam} -o {output} \
                    --normalizeUsing {params.normalize} \
                    --numberOfProcessors {threads}
        """

# CUT&RUN quality control profile around genes (optional)
rule cutrun_profile_qc:
    input:
        bam = f"{config['outdir']}/cutrun/bam/{{sample}}.dedup.bam" if config.get('cutrun', {}).get('markdup', True) else f"{config['outdir']}/cutrun/bam/{{sample}}.sorted.bam",
        gff = config.get('cutrun', {}).get('qc', {}).get('gff', '')
    output:
        profile = f"{config['outdir']}/cutrun/qc/{{sample}}_gene_profile.{config.get('cutrun', {}).get('qc', {}).get('plot_format', 'png')}",
        matrix = f"{config['outdir']}/cutrun/qc/{{sample}}_gene_profile.gz"
    params:
        upstream = config.get('cutrun', {}).get('qc', {}).get('upstream', 2000),
        downstream = config.get('cutrun', {}).get('qc', {}).get('downstream', 2000),
        bin_size = config.get('cutrun', {}).get('qc', {}).get('bin_size', 50)
    threads: config.get('cutrun', {}).get('threads', config['resources']['mapping_threads'])
    shell:
        """
        mkdir -p $(dirname {output.profile})
        computeMatrix reference-point -S {input.bam} \
                                     -R {input.gff} \
                                     --referencePoint TSS \
                                     -b {params.upstream} -a {params.downstream} \
                                     --binSize {params.bin_size} \
                                     -p {threads} \
                                     -o {output.matrix}
        
        plotProfile -m {output.matrix} \
                   -out {output.profile} \
                   --plotTitle "CUT&RUN Signal around TSS - {wildcards.sample}"
        """

# Complete CUT&RUN pipeline marker
rule cutrun_complete:
    input:
        expand(f"{config['outdir']}/cutrun/bam/{{sample}}.dedup.bam", sample=CUTRUN_SAMPLES) if config.get('cutrun', {}).get('markdup', True) else expand(f"{config['outdir']}/cutrun/bam/{{sample}}.sorted.bam", sample=CUTRUN_SAMPLES),
        expand(f"{config['outdir']}/cutrun/bigwig/{{sample}}.bw", sample=CUTRUN_SAMPLES),
        expand(f"{config['outdir']}/cutrun/peaks/{{sample}}_peaks.narrowPeak", sample=get_cutrun_peak_samples())
    output:
        f"{GENOME_OUTDIR}/cutrun/tracks_complete.txt"
    params:
    shell:
        "echo 'CUT&RUN pipeline completed at $(date)' > {output}"

# Resolve mapping rule ambiguity based on configuration
if config.get('cutrun', {}).get('mapper', 'bowtie2') == "bowtie2":
    ruleorder: cutrun_map_bowtie2 > cutrun_map_bwa_mem2
else:
    ruleorder: cutrun_map_bwa_mem2 > cutrun_map_bowtie2
