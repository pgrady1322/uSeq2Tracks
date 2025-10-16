# ============================================================
# Ancient DNA analysis rules
# ============================================================

# Map ancient DNA reads with bwa aln (traditional approach for degraded DNA)
rule ancientdna_map_bwa_aln:
    input:
        unpack(get_input_reads),
        ref = f"{GENOME_OUTDIR}/genome/genome.fa"
    output:
        f"{config['outdir']}/ancientdna/bam/{{sample}}.raw.bam"
    params:
        bwa_aln_opts = config.get('ancientdna', {}).get('bwa_aln_opts', '-l 16500')
    threads: config.get('ancientdna', {}).get('threads', config['resources']['mapping_threads'])
    run:
        # Check if we have paired-end or single-end data
        has_r2 = hasattr(input, 'r2') and input.r2
        
        shell("mkdir -p $(dirname {output})")
        
        if has_r2:
            # Paired-end BWA aln
            shell("""
            bwa aln {params.bwa_aln_opts} -t {threads} {input.ref} {input.r1} > temp_{wildcards.sample}_1.sai
            bwa aln {params.bwa_aln_opts} -t {threads} {input.ref} {input.r2} > temp_{wildcards.sample}_2.sai
            bwa sampe {input.ref} temp_{wildcards.sample}_1.sai temp_{wildcards.sample}_2.sai {input.r1} {input.r2} \
              | samtools view -bS -o {output} -
            rm temp_{wildcards.sample}_*.sai
            """)
        else:
            # Single-end BWA aln
            shell("""
            bwa aln {params.bwa_aln_opts} -t {threads} {input.ref} {input.r1} > temp_{wildcards.sample}_1.sai
            bwa samse {input.ref} temp_{wildcards.sample}_1.sai {input.r1} \
              | samtools view -bS -o {output} -
            rm temp_{wildcards.sample}_*.sai
            """)

# Map ancient DNA reads with bwa mem (alternative, faster approach)
rule ancientdna_map_bwa_mem:
    input:
        unpack(get_input_reads),
        ref = f"{GENOME_OUTDIR}/genome/genome.fa"
    output:
        f"{config['outdir']}/ancientdna/bam/{{sample}}.raw.bam"
    params:
        bwa_mem_opts = config.get('ancientdna', {}).get('bwa_mem_opts', '-M')
    threads: config.get('ancientdna', {}).get('threads', config['resources']['mapping_threads'])
    run:
        # Check if we have paired-end or single-end data
        has_r2 = hasattr(input, 'r2') and input.r2
        
        shell("mkdir -p $(dirname {output})")
        
        if has_r2:
            # Paired-end
            shell("""
            bwa mem {params.bwa_mem_opts} -t {threads} {input.ref} {input.r1} {input.r2} \
              | samtools view -bS -o {output} -
            """)
        else:
            # Single-end
            shell("""
            bwa mem {params.bwa_mem_opts} -t {threads} {input.ref} {input.r1} \
              | samtools view -bS -o {output} -
            """)

# Sort and index ancient DNA BAM files
rule ancientdna_sort_and_index:
    input:
        f"{config['outdir']}/ancientdna/bam/{{sample}}.raw.bam"
    output:
        bam = f"{config['outdir']}/ancientdna/bam/{{sample}}.sorted.bam",
        bai = f"{config['outdir']}/ancientdna/bam/{{sample}}.sorted.bam.bai"
    params:
    threads: config.get('ancientdna', {}).get('threads', config['resources']['mapping_threads'])
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input}
        samtools index -@ {threads} {output.bam}
        """

# Remove duplicates (important for ancient DNA)
rule ancientdna_remove_duplicates:
    input:
        f"{config['outdir']}/ancientdna/bam/{{sample}}.sorted.bam"
    output:
        bam = f"{config['outdir']}/ancientdna/bam/{{sample}}.dedup.bam",
        metrics = f"{config['outdir']}/ancientdna/qc/{{sample}}.dedup.metrics"
    threads: config.get('ancientdna', {}).get('threads', config['resources']['mapping_threads'])
    shell:
        """
        mkdir -p $(dirname {output.metrics})
        samtools markdup -r -@ {threads} -f {output.metrics} {input} {output.bam}
        samtools index {output.bam}
        """

# Ancient DNA damage pattern analysis with mapDamage2
rule ancientdna_damage_analysis:
    input:
        bam = f"{config['outdir']}/ancientdna/bam/{{sample}}.dedup.bam" if config.get('ancientdna', {}).get('markdup', True) else f"{config['outdir']}/ancientdna/bam/{{sample}}.sorted.bam",
        ref = f"{GENOME_OUTDIR}/genome/genome.fa"
    output:
        stats = f"{config['outdir']}/ancientdna/damage/{{sample}}/3pGtoA_freq.txt",
        plot = f"{config['outdir']}/ancientdna/damage/{{sample}}/Fragmisincorporation_plot.pdf"
    params:
        outdir = f"{config['outdir']}/ancientdna/damage/{{sample}}"
    shell:
        """
        mkdir -p {params.outdir}
        mapDamage -i {input.bam} -r {input.ref} -d {params.outdir}
        """

# Variant calling optimized for ancient DNA with low coverage
rule ancientdna_variant_calling:
    input:
        bam = f"{config['outdir']}/ancientdna/bam/{{sample}}.dedup.bam" if config.get('ancientdna', {}).get('markdup', True) else f"{config['outdir']}/ancientdna/bam/{{sample}}.sorted.bam",
        ref = f"{GENOME_OUTDIR}/genome/genome.fa"
    output:
        f"{config['outdir']}/ancientdna/variants/{{sample}}.vcf.gz"
    params:
        min_mapq = config.get('ancientdna', {}).get('min_mapq', 20),
        min_baseq = config.get('ancientdna', {}).get('min_baseq', 20)
    threads: config.get('ancientdna', {}).get('threads', config['resources']['mapping_threads'])
    shell:
        """
        mkdir -p $(dirname {output})
        bcftools mpileup -f {input.ref} -q {params.min_mapq} -Q {params.min_baseq} {input.bam} | \
        bcftools call -c -v -Oz -o {output}
        bcftools index {output}
        """

# Ancient DNA authentication with ANGSD
rule ancientdna_authentication:
    input:
        bam = f"{config['outdir']}/ancientdna/bam/{{sample}}.dedup.bam" if config.get('ancientdna', {}).get('markdup', True) else f"{config['outdir']}/ancientdna/bam/{{sample}}.sorted.bam",
        ref = f"{GENOME_OUTDIR}/genome/genome.fa"
    output:
        depth = f"{config['outdir']}/ancientdna/auth/{{sample}}.depthSample",
        stats = f"{config['outdir']}/ancientdna/auth/{{sample}}.stats"
    params:
        outdir = f"{config['outdir']}/ancientdna/auth",
        prefix = "{sample}"
    threads: config.get('ancientdna', {}).get('threads', config['resources']['mapping_threads'])
    shell:
        """
        mkdir -p {params.outdir}
        angsd -i {input.bam} -ref {input.ref} \
              -doDepth 1 -doCounts 1 \
              -out {params.outdir}/{params.prefix} \
              -P {threads}
        """

# Convert ancient DNA BAM to BigWig for visualization
rule ancientdna_bam_to_bigwig:
    input:
        bam = f"{config['outdir']}/ancientdna/bam/{{sample}}.dedup.bam" if config.get('ancientdna', {}).get('markdup', True) else f"{config['outdir']}/ancientdna/bam/{{sample}}.sorted.bam",
        fai = f"{GENOME_OUTDIR}/genome/genome.fa.fai"
    output:
        f"{config['outdir']}/ancientdna/bigwig/{{sample}}.bw"
    params:
        normalize = config.get('ancientdna', {}).get('bw_norm', 'CPM')
    threads: config.get('ancientdna', {}).get('threads', config['resources']['mapping_threads'])
    shell:
        """
        mkdir -p $(dirname {output})
        bamCoverage -b {input.bam} -o {output} \
                    --normalizeUsing {params.normalize} \
                    --numberOfProcessors {threads}
        """

# Complete ancient DNA pipeline marker
rule ancientdna_complete:
    input:
        expand(f"{config['outdir']}/ancientdna/bam/{{sample}}.dedup.bam", sample=ADNA_SAMPLES) if config.get('ancientdna', {}).get('markdup', True) else expand(f"{config['outdir']}/ancientdna/bam/{{sample}}.sorted.bam", sample=ADNA_SAMPLES),
        expand(f"{config['outdir']}/ancientdna/bigwig/{{sample}}.bw", sample=ADNA_SAMPLES),
        expand(f"{config['outdir']}/ancientdna/damage/{{sample}}/3pGtoA_freq.txt", sample=ADNA_SAMPLES)
    output:
        f"{GENOME_OUTDIR}/ancientdna/tracks_complete.txt"
    params:
    shell:
        "echo 'Ancient DNA pipeline completed at $(date)' > {output}"
