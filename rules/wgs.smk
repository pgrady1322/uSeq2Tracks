# ============================================================
# WGS (Whole Genome Sequencing) analysis rules
# ============================================================

# Map WGS reads with bowtie2
rule wgs_map_bowtie2:
    input:
        unpack(get_input_reads),
        idx = f"{GENOME_OUTDIR}/genome/bowtie2/genome.1.bt2"
    output:
        f"{config['outdir']}/wgs/bam/{{sample}}.raw.bam"
    params:
        bowtie2_opts = config.get('wgs', {}).get('bowtie2_opts', '--very-sensitive'),
        genome_prefix = f"{GENOME_OUTDIR}/genome/bowtie2/genome"
    threads: config.get('wgs', {}).get('threads', config['resources']['mapping_threads'])
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

# Map WGS reads with bwa-mem2 (default for WGS)
rule wgs_map_bwa_mem2:
    input:
        unpack(get_input_reads),
        idx = f"{GENOME_OUTDIR}/genome/bwa_mem2/genome.fa.bwt.2bit.64"
    output:
        f"{config['outdir']}/wgs/bam/{{sample}}.raw.bam"
    params:
        bwa_mem2_opts = config.get('wgs', {}).get('bwa_mem2_opts', '-M'),
        genome_prefix = f"{GENOME_OUTDIR}/genome/bwa_mem2/genome.fa"
    threads: config.get('wgs', {}).get('threads', config['resources']['mapping_threads'])
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

# Sort and index WGS BAM files
rule wgs_sort_and_index:
    input:
        f"{config['outdir']}/wgs/bam/{{sample}}.raw.bam"
    output:
        bam = f"{config['outdir']}/wgs/bam/{{sample}}.sorted.bam",
        bai = f"{config['outdir']}/wgs/bam/{{sample}}.sorted.bam.bai"
    params:
    threads: config.get('wgs', {}).get('threads', config['resources']['mapping_threads'])
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input}
        samtools index -@ {threads} {output.bam}
        """

# Mark duplicates in WGS BAM files
rule wgs_markdup:
    input:
        f"{config['outdir']}/wgs/bam/{{sample}}.sorted.bam"
    output:
        bam = f"{config['outdir']}/wgs/bam/{{sample}}.dedup.bam",
        metrics = f"{config['outdir']}/wgs/qc/{{sample}}.markdup.metrics"
    threads: config.get('wgs', {}).get('threads', config['resources']['mapping_threads'])
    shell:
        """
        mkdir -p $(dirname {output.metrics})
        samtools markdup -@ {threads} -f {output.metrics} {input} {output.bam}
        samtools index {output.bam}
        """

# Convert WGS BAM to CRAM for space efficiency
rule wgs_bam_to_cram:
    input:
        bam = f"{config['outdir']}/wgs/bam/{{sample}}.dedup.bam" if config.get('wgs', {}).get('markdup', True) else f"{config['outdir']}/wgs/bam/{{sample}}.sorted.bam",
        ref = f"{GENOME_OUTDIR}/genome/genome.fa"
    output:
        f"{config['outdir']}/wgs/cram/{{sample}}.cram"
    params:
    threads: config.get('wgs', {}).get('threads', config['resources']['mapping_threads'])
    shell:
        """
        mkdir -p $(dirname {output})
        samtools view -C -T {input.ref} -@ {threads} -o {output} {input.bam}
        samtools index {output}
        """

# Convert WGS BAM to BigWig for visualization
rule wgs_bam_to_bigwig:
    input:
        bam = f"{config['outdir']}/wgs/bam/{{sample}}.dedup.bam" if config.get('wgs', {}).get('markdup', True) else f"{config['outdir']}/wgs/bam/{{sample}}.sorted.bam",
        fai = f"{GENOME_OUTDIR}/genome/genome.fa.fai"
    output:
        f"{config['outdir']}/wgs/bigwig/{{sample}}.bw"
    params:
        normalize = config.get('wgs', {}).get('bw_norm', 'CPM')
    threads: config.get('wgs', {}).get('threads', config['resources']['mapping_threads'])
    shell:
        """
        mkdir -p $(dirname {output})
        bamCoverage -b {input.bam} -o {output} \
                    --normalizeUsing {params.normalize} \
                    --numberOfProcessors {threads}
        """

# Run variant calling with bcftools (optional for WGS)
rule wgs_variant_calling:
    input:
        bam = f"{config['outdir']}/wgs/bam/{{sample}}.dedup.bam" if config.get('wgs', {}).get('markdup', True) else f"{config['outdir']}/wgs/bam/{{sample}}.sorted.bam",
        ref = f"{GENOME_OUTDIR}/genome/genome.fa"
    output:
        f"{config['outdir']}/wgs/variants/{{sample}}.vcf.gz"
    params:
    threads: config.get('wgs', {}).get('threads', config['resources']['mapping_threads'])
    shell:
        """
        mkdir -p $(dirname {output})
        bcftools mpileup -f {input.ref} {input.bam} | \
        bcftools call -mv -Oz -o {output}
        bcftools index {output}
        """

# Complete WGS pipeline marker
rule wgs_complete:
    input:
        expand(f"{config['outdir']}/wgs/bam/{{sample}}.dedup.bam", sample=WGS_SAMPLES) if config.get('wgs', {}).get('markdup', True) else expand(f"{config['outdir']}/wgs/bam/{{sample}}.sorted.bam", sample=WGS_SAMPLES),
        expand(f"{config['outdir']}/wgs/bigwig/{{sample}}.bw", sample=WGS_SAMPLES),
        expand(f"{config['outdir']}/wgs/cram/{{sample}}.cram", sample=WGS_SAMPLES)
    output:
        f"{GENOME_OUTDIR}/wgs/tracks_complete.txt"
    params:
    shell:
        "echo 'WGS pipeline completed at $(date)' > {output}"

# Resolve mapping rule ambiguity based on configuration
if config.get('wgs', {}).get('mapper', 'bwa_mem2') == "bowtie2":
    ruleorder: wgs_map_bowtie2 > wgs_map_bwa_mem2
else:
    ruleorder: wgs_map_bwa_mem2 > wgs_map_bowtie2
