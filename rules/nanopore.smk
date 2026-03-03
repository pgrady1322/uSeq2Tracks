# ============================================================
# Nanopore long-read analysis rules
# ============================================================

# Map Nanopore reads with minimap2
rule nanopore_map_minimap2:
    input:
        unpack(get_input_reads),
        ref = f"{GENOME_OUTDIR}/genome/genome.fa"
    output:
        f"{GENOME_OUTDIR}/nanopore/bam/{{sample}}.raw.bam"
    params:
        minimap2_opts = config.get('nanopore', {}).get('minimap2_opts', '-x map-ont')
    threads: config.get('nanopore', {}).get('threads', config['resources']['mapping_threads'])
    shell:
        """
        mkdir -p $(dirname {output})
        minimap2 -ax map-ont {params.minimap2_opts} -t {threads} \
                 {input.ref} {input.r1} \
          | samtools view -bS -o {output} -
        """

# Sort and index Nanopore BAM files
rule nanopore_sort_and_index:
    input:
        f"{GENOME_OUTDIR}/nanopore/bam/{{sample}}.raw.bam"
    output:
        bam = f"{GENOME_OUTDIR}/nanopore/bam/{{sample}}.sorted.bam",
        bai = f"{GENOME_OUTDIR}/nanopore/bam/{{sample}}.sorted.bam.bai"
    params:
    threads: config.get('nanopore', {}).get('threads', config['resources']['mapping_threads'])
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input}
        samtools index -@ {threads} {output.bam}
        """

# Call variants with Clair3 (optimized for Nanopore) - only if enabled
rule nanopore_variant_calling:
    input:
        bam = f"{GENOME_OUTDIR}/nanopore/bam/{{sample}}.sorted.bam",
        ref = f"{GENOME_OUTDIR}/genome/genome.fa"
    output:
        f"{GENOME_OUTDIR}/nanopore/variants/{{sample}}.vcf.gz"
    params:
        model = config.get('nanopore', {}).get('clair3_model', '/path/to/model'),
        outdir = f"{GENOME_OUTDIR}/nanopore/variants/{{sample}}"
    threads: config.get('nanopore', {}).get('threads', config['resources']['mapping_threads'])
    shell:
        """
        mkdir -p {params.outdir}
        run_clair3.sh --bam_fn={input.bam} \
                      --ref_fn={input.ref} \
                      --threads={threads} \
                      --platform=ont \
                      --model_path={params.model} \
                      --output={params.outdir}
        mv {params.outdir}/merge_output.vcf.gz {output}
        """

# Convert Nanopore BAM to BigWig for visualization
rule nanopore_bam_to_bigwig:
    input:
        bam = f"{GENOME_OUTDIR}/nanopore/bam/{{sample}}.sorted.bam",
        fai = f"{GENOME_OUTDIR}/genome/genome.fa.fai"
    output:
        f"{GENOME_OUTDIR}/nanopore/bigwig/{{sample}}.bw"
    params:
        normalize = config.get('nanopore', {}).get('bw_norm', 'CPM')
    threads: config.get('nanopore', {}).get('threads', config['resources']['mapping_threads'])
    shell:
        """
        mkdir -p $(dirname {output})
        bamCoverage -b {input.bam} -o {output} \
                    --normalizeUsing {params.normalize} \
                    --numberOfProcessors {threads}
        """

# Nanopore read statistics and quality assessment
rule nanopore_read_stats:
    input:
        unpack(get_input_reads)
    output:
        stats = f"{GENOME_OUTDIR}/nanopore/qc/{{sample}}_read_stats.txt",
        plot = f"{GENOME_OUTDIR}/nanopore/qc/{{sample}}_read_length_dist.png"
    params:
    shell:
        """
        mkdir -p $(dirname {output.stats})
        NanoStat --fastq {input.r1} --name {output.stats}
        NanoPlot --fastq {input.r1} --plots dot --outdir $(dirname {output.plot}) --prefix {wildcards.sample}_
        """

# Complete Nanopore pipeline marker
rule nanopore_complete:
    input:
        expand(f"{GENOME_OUTDIR}/nanopore/bam/{{sample}}.sorted.bam", sample=NANOPORE_SAMPLES),
        expand(f"{GENOME_OUTDIR}/nanopore/bigwig/{{sample}}.bw", sample=NANOPORE_SAMPLES),
        # Only include variants if variant calling is enabled
        (expand(f"{GENOME_OUTDIR}/nanopore/variants/{{sample}}.vcf.gz", sample=NANOPORE_SAMPLES) 
         if config['nanopore']['variant_calling'] else [])
    output:
        f"{GENOME_OUTDIR}/nanopore/tracks_complete.txt"
    params:
    shell:
        "echo 'Nanopore pipeline completed at $(date)' > {output}"
