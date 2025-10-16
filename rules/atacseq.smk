# ============================================================
# ATAC-seq analysis rules
# ATAC-seq: Assay for Transposase-Accessible Chromatin
# ============================================================

# Map ATAC-seq reads with bowtie2
rule atacseq_map_bowtie2:
    input:
        unpack(get_input_reads),
        idx = f"{GENOME_OUTDIR}/genome/bowtie2/genome.1.bt2"
    output:
        f"{config['outdir']}/atacseq/bam/{{sample}}.raw.bam"
    params:
        bowtie2_opts = config.get('atacseq', {}).get('bowtie2_opts', '--very-sensitive'),
        genome_prefix = f"{GENOME_OUTDIR}/genome/bowtie2/genome"
    threads: config.get('atacseq', {}).get('threads', config['resources']['mapping_threads'])
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

# Map ATAC-seq reads with bwa-mem2 
rule atacseq_map_bwa_mem2:
    input:
        unpack(get_input_reads),
        idx = f"{GENOME_OUTDIR}/genome/bwa_mem2/genome.fa.bwt.2bit.64"
    output:
        f"{config['outdir']}/atacseq/bam/{{sample}}.raw.bam"
    params:
        bwa_mem2_opts = config.get('atacseq', {}).get('bwa_mem2_opts', '-M'),
        genome_prefix = f"{GENOME_OUTDIR}/genome/bwa_mem2/genome.fa"
    threads: config.get('atacseq', {}).get('threads', config['resources']['mapping_threads'])
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

# Sort, mark duplicates, and index ATAC-seq BAM files
rule atacseq_sort_and_index:
    input:
        f"{config['outdir']}/atacseq/bam/{{sample}}.raw.bam"
    output:
        bam = f"{config['outdir']}/atacseq/bam/{{sample}}.sorted.bam",
        bai = f"{config['outdir']}/atacseq/bam/{{sample}}.sorted.bam.bai"
    params:
    threads: config.get('atacseq', {}).get('threads', config['resources']['mapping_threads'])
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

# Mark duplicates in ATAC-seq BAM files
rule atacseq_markdup:
    input:
        f"{config['outdir']}/atacseq/bam/{{sample}}.sorted.bam"
    output:
        bam = f"{config['outdir']}/atacseq/bam/{{sample}}.dedup.bam",
        metrics = f"{config['outdir']}/atacseq/qc/{{sample}}.markdup.metrics"
    threads: config.get('atacseq', {}).get('threads', config['resources']['mapping_threads'])
    shell:
        """
        mkdir -p $(dirname {output.metrics})
        # Mark duplicates but don't remove them (add -s flag for ATAC-seq)
        samtools markdup -@ {threads} -f {output.metrics} -s {input} {output.bam}
        samtools index {output.bam}
        """

# Call peaks with MACS3 for ATAC-seq (using ATAC-seq specific parameters)
rule atacseq_call_peaks:
    input:
        f"{config['outdir']}/atacseq/bam/{{sample}}.dedup.bam" if config.get('atacseq', {}).get('markdup', True) else f"{config['outdir']}/atacseq/bam/{{sample}}.sorted.bam"
    output:
        f"{config['outdir']}/atacseq/peaks/{{sample}}_peaks.narrowPeak"
    params:
        shift = config.get('atacseq', {}).get('shift', -75),
        extsize = config.get('atacseq', {}).get('extsize', 150),
        macs3_opts = config.get('atacseq', {}).get('macs3_opts', '--qval 0.05'),
        outdir = f"{config['outdir']}/atacseq/peaks",
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

# Parameter sweep: Call peaks with multiple q-values for ATAC-seq
rule atacseq_call_peaks_sweep:
    input:
        f"{config['outdir']}/atacseq/bam/{{sample}}.dedup.bam" if config.get('atacseq', {}).get('markdup', True) else f"{config['outdir']}/atacseq/bam/{{sample}}.sorted.bam"
    output:
        f"{config['outdir']}/atacseq/peaks_sweep/{{sample}}_q{{qval}}_peaks.narrowPeak"
    params:
        shift = config.get('atacseq', {}).get('shift', -75),
        extsize = config.get('atacseq', {}).get('extsize', 150),
        outdir = f"{config['outdir']}/atacseq/peaks_sweep",
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

# Convert ATAC-seq BAM to BigWig for visualization
rule atacseq_bam_to_bigwig:
    input:
        bam = f"{config['outdir']}/atacseq/bam/{{sample}}.dedup.bam" if config.get('atacseq', {}).get('markdup', True) else f"{config['outdir']}/atacseq/bam/{{sample}}.sorted.bam",
        fai = f"{GENOME_OUTDIR}/genome/genome.fa.fai"
    output:
        f"{config['outdir']}/atacseq/bigwig/{{sample}}.bw"
    params:
        normalize = config.get('atacseq', {}).get('bw_norm', 'CPM')
    threads: config.get('atacseq', {}).get('threads', config['resources']['mapping_threads'])
    shell:
        """
        mkdir -p $(dirname {output})
        
        # Check if bamCoverage is available
        if ! command -v bamCoverage &> /dev/null; then
            echo "ERROR: bamCoverage command not found (deepTools)"
            echo "Checking available tools:"
            which samtools || echo "samtools not found"
            exit 1
        fi
        
        echo "Using bamCoverage from: $(which bamCoverage)"
        echo "BAM file size: $(ls -lh {input.bam})"
        
        bamCoverage -b {input.bam} -o {output} \
                    --normalizeUsing {params.normalize} \
                    --numberOfProcessors {threads}
        """

# Complete ATAC-seq pipeline marker
rule atacseq_complete:
    input:
        expand(f"{config['outdir']}/atacseq/bam/{{sample}}.dedup.bam", sample=ATACSEQ_SAMPLES) if config.get('atacseq', {}).get('markdup', True) else expand(f"{config['outdir']}/atacseq/bam/{{sample}}.sorted.bam", sample=ATACSEQ_SAMPLES),
        expand(f"{config['outdir']}/atacseq/bigwig/{{sample}}.bw", sample=ATACSEQ_SAMPLES),
        expand(f"{config['outdir']}/atacseq/peaks/{{sample}}_peaks.narrowPeak", sample=ATACSEQ_SAMPLES)
    output:
        f"{GENOME_OUTDIR}/atacseq/tracks_complete.txt"
    params:
    shell:
        "echo 'ATAC-seq pipeline completed at $(date)' > {output}"

# Resolve mapping rule ambiguity based on configuration
if config.get('atacseq', {}).get('mapper', 'bowtie2') == "bowtie2":
    ruleorder: atacseq_map_bowtie2 > atacseq_map_bwa_mem2
else:
    ruleorder: atacseq_map_bwa_mem2 > atacseq_map_bowtie2
