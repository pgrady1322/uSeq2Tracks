# ============================================================
# Nanopore long-read analysis rules
# ============================================================

# Map Nanopore reads with minimap2
rule nanopore_map_minimap2:
    input:
        reads = f"{config['outdir']}/sra/{{sample}}.fastq.gz",
        ref = f"{GENOME_OUTDIR}/genome/genome.fa"
    output:
        f"{config['outdir']}/nanopore/bam/{{sample}}.raw.bam"
    params:
        minimap2_opts = config.get('nanopore', {}).get('minimap2_opts', '-x map-ont')
    threads: config.get('nanopore', {}).get('threads', config['resources']['mapping_threads'])
    shell:
        """
        mkdir -p $(dirname {output})
        minimap2 -ax map-ont {params.minimap2_opts} -t {threads} \
                 {input.ref} {input.reads} \
          | samtools view -bS -o {output} -
        """

# Sort and index Nanopore BAM files
rule nanopore_sort_and_index:
    input:
        f"{config['outdir']}/nanopore/bam/{{sample}}.raw.bam"
    output:
        bam = f"{config['outdir']}/nanopore/bam/{{sample}}.sorted.bam",
        bai = f"{config['outdir']}/nanopore/bam/{{sample}}.sorted.bam.bai"
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
        bam = f"{config['outdir']}/nanopore/bam/{{sample}}.sorted.bam",
        ref = f"{GENOME_OUTDIR}/genome/genome.fa"
    output:
        f"{config['outdir']}/nanopore/variants/{{sample}}.vcf.gz"
    params:
        model = config.get('nanopore', {}).get('clair3_model', '/path/to/model'),
        outdir = f"{config['outdir']}/nanopore/variants/{{sample}}"
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
        bam = f"{config['outdir']}/nanopore/bam/{{sample}}.sorted.bam",
        fai = f"{GENOME_OUTDIR}/genome/genome.fa.fai"
    output:
        f"{config['outdir']}/nanopore/bigwig/{{sample}}.bw"
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
        reads = f"{config['outdir']}/sra/{{sample}}.fastq.gz"
    output:
        stats = f"{config['outdir']}/nanopore/qc/{{sample}}_read_stats.txt",
        plot = f"{config['outdir']}/nanopore/qc/{{sample}}_read_length_dist.png"
    params:
    shell:
        """
        mkdir -p $(dirname {output.stats})
        NanoStat --fastq {input.reads} --name {output.stats}
        NanoPlot --fastq {input.reads} --plots dot --outdir $(dirname {output.plot}) --prefix {wildcards.sample}_
        """

# Complete Nanopore pipeline marker
rule nanopore_complete:
    input:
        expand(f"{config['outdir']}/nanopore/bam/{{sample}}.sorted.bam", sample=NANOPORE_SAMPLES),
        expand(f"{config['outdir']}/nanopore/bigwig/{{sample}}.bw", sample=NANOPORE_SAMPLES),
        # Only include variants if variant calling is enabled
        (expand(f"{config['outdir']}/nanopore/variants/{{sample}}.vcf.gz", sample=NANOPORE_SAMPLES) 
         if config['nanopore']['variant_calling'] else [])
    output:
        f"{GENOME_OUTDIR}/nanopore/tracks_complete.txt"
    params:
    shell:
        "echo 'Nanopore pipeline completed at $(date)' > {output}"
