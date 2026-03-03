# ============================================================
# PacBio long-read analysis rules
# ============================================================

# Map PacBio reads with minimap2
rule pacbio_map_minimap2:
    input:
        unpack(get_input_reads),
        ref = f"{GENOME_OUTDIR}/genome/genome.fa"
    output:
        f"{GENOME_OUTDIR}/pacbio/bam/{{sample}}.raw.bam"
    params:
        minimap2_opts = config.get('pacbio', {}).get('minimap2_opts', '-x map-pb')
    threads: config.get('pacbio', {}).get('threads', config['resources']['mapping_threads'])
    shell:
        """
        mkdir -p $(dirname {output})
        minimap2 -ax map-pb {params.minimap2_opts} -t {threads} \
                 {input.ref} {input.r1} \
          | samtools view -bS -o {output} -
        """

# Sort and index PacBio BAM files
rule pacbio_sort_and_index:
    input:
        f"{GENOME_OUTDIR}/pacbio/bam/{{sample}}.raw.bam"
    output:
        bam = f"{GENOME_OUTDIR}/pacbio/bam/{{sample}}.sorted.bam",
        bai = f"{GENOME_OUTDIR}/pacbio/bam/{{sample}}.sorted.bam.bai"
    params:
    threads: config.get('pacbio', {}).get('threads', config['resources']['mapping_threads'])
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input}
        samtools index -@ {threads} {output.bam}
        """

# Call variants with pbsv (PacBio structural variant caller)
rule pacbio_structural_variants:
    input:
        bam = f"{GENOME_OUTDIR}/pacbio/bam/{{sample}}.sorted.bam",
        ref = f"{GENOME_OUTDIR}/genome/genome.fa"
    output:
        svsig = f"{GENOME_OUTDIR}/pacbio/variants/{{sample}}.svsig.gz",
        vcf = f"{GENOME_OUTDIR}/pacbio/variants/{{sample}}.sv.vcf"
    params:
    threads: config.get('pacbio', {}).get('threads', config['resources']['mapping_threads'])
    shell:
        """
        mkdir -p $(dirname {output.svsig})
        pbsv discover {input.bam} {output.svsig}
        pbsv call {input.ref} {output.svsig} {output.vcf}
        """

# Convert PacBio BAM to BigWig for visualization
rule pacbio_bam_to_bigwig:
    input:
        bam = f"{GENOME_OUTDIR}/pacbio/bam/{{sample}}.sorted.bam",
        fai = f"{GENOME_OUTDIR}/genome/genome.fa.fai"
    output:
        f"{GENOME_OUTDIR}/pacbio/bigwig/{{sample}}.bw"
    params:
        normalize = config.get('pacbio', {}).get('bw_norm', 'CPM')
    threads: config.get('pacbio', {}).get('threads', config['resources']['mapping_threads'])
    shell:
        """
        mkdir -p $(dirname {output})
        bamCoverage -b {input.bam} -o {output} \
                    --normalizeUsing {params.normalize} \
                    --numberOfProcessors {threads}
        """

# Complete PacBio pipeline marker
rule pacbio_complete:
    input:
        expand(f"{GENOME_OUTDIR}/pacbio/bam/{{sample}}.sorted.bam", sample=PACBIO_SAMPLES),
        expand(f"{GENOME_OUTDIR}/pacbio/bigwig/{{sample}}.bw", sample=PACBIO_SAMPLES),
        expand(f"{GENOME_OUTDIR}/pacbio/variants/{{sample}}.sv.vcf", sample=PACBIO_SAMPLES)
    output:
        f"{GENOME_OUTDIR}/pacbio/tracks_complete.txt"
    params:
    shell:
        "echo 'PacBio pipeline completed at $(date)' > {output}"
