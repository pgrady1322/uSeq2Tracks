# ============================================================
# RNA-seq analysis rules
# ============================================================

# Map RNA-seq reads with STAR (gold standard for RNA-seq)
rule rnaseq_map_star:
    input:
        unpack(get_input_reads),
        idx = f"{GENOME_OUTDIR}/genome/star/Genome"
    output:
        bam = f"{GENOME_OUTDIR}/rnaseq/bam/{{sample}}.Aligned.sortedByCoord.out.bam",
        counts = f"{GENOME_OUTDIR}/rnaseq/counts/{{sample}}.ReadsPerGene.out.tab"
    params:
        star_opts = config['rnaseq']['star_opts'],
        genome_dir = f"{GENOME_OUTDIR}/genome/star",
        outdir = f"{GENOME_OUTDIR}/rnaseq/bam",
        prefix = "{sample}.",
        has_gtf = bool(config.get('gtf') and config['gtf'])
    threads: config['rnaseq']['threads']
    run:
        # Check if we have paired-end or single-end data
        has_r2 = hasattr(input, 'r2') and input.r2
        
        # Create output directories
        shell("mkdir -p {params.outdir}")
        shell("mkdir -p $(dirname {output.counts})")
        
        # Build STAR command with conditional gene counting
        quantmode = "--quantMode GeneCounts" if params.has_gtf else ""
        
        if has_r2:
            # Paired-end
            shell(f"""
            STAR --runThreadN {{threads}} \
                 --genomeDir {{params.genome_dir}} \
                 --readFilesIn {{input.r1}} {{input.r2}} \
                 --readFilesCommand zcat \
                 --outFileNamePrefix {{params.outdir}}/{{params.prefix}} \
                 --outSAMtype BAM SortedByCoordinate \
                 {quantmode} \
                 {{params.star_opts}}
            """)
        else:
            # Single-end
            shell(f"""
            STAR --runThreadN {{threads}} \
                 --genomeDir {{params.genome_dir}} \
                 --readFilesIn {{input.r1}} \
                 --readFilesCommand zcat \
                 --outFileNamePrefix {{params.outdir}}/{{params.prefix}} \
                 --outSAMtype BAM SortedByCoordinate \
                 {quantmode} \
                 {{params.star_opts}}
            """)
        
        # Handle counts file - create dummy if no GTF provided
        if params.has_gtf:
            shell("mv {params.outdir}/{params.prefix}ReadsPerGene.out.tab {output.counts}")
        else:
            # Create empty counts file as placeholder
            shell("touch {output.counts}")
            shell("echo '# No gene counts available - no GTF file provided' > {output.counts}")

# Index RNA-seq BAM files
rule rnaseq_index_bam:
    input:
        f"{GENOME_OUTDIR}/rnaseq/bam/{{sample}}.Aligned.sortedByCoord.out.bam"
    output:
        f"{GENOME_OUTDIR}/rnaseq/bam/{{sample}}.Aligned.sortedByCoord.out.bam.bai"
    params:
    threads: config['rnaseq']['threads']
    shell:
        """
        samtools index -@ {threads} {input}
        """

# Convert RNA-seq BAM to BigWig for visualization
rule rnaseq_bam_to_bigwig:
    input:
        bam = f"{GENOME_OUTDIR}/rnaseq/bam/{{sample}}.Aligned.sortedByCoord.out.bam",
        bai = f"{GENOME_OUTDIR}/rnaseq/bam/{{sample}}.Aligned.sortedByCoord.out.bam.bai",
        fai = f"{GENOME_OUTDIR}/genome/genome.fa.fai"
    output:
        f"{GENOME_OUTDIR}/rnaseq/bigwig/{{sample}}.bw"
    params:
        normalize = config['rnaseq']['bw_norm']
    threads: config['rnaseq']['threads']
    shell:
        """
        mkdir -p $(dirname {output})
        bamCoverage -b {input.bam} -o {output} \
                    --normalizeUsing {params.normalize} \
                    --numberOfProcessors {threads}
        """

# RNA-seq quality control with RSeQC (optional)
rule rnaseq_rseqc_qc:
    input:
        bam = f"{GENOME_OUTDIR}/rnaseq/bam/{{sample}}.Aligned.sortedByCoord.out.bam",
        bed = config['rnaseq']['gene_bed']
    output:
        txt = f"{GENOME_OUTDIR}/rnaseq/qc/{{sample}}.geneBodyCoverage.txt",
        pdf = f"{GENOME_OUTDIR}/rnaseq/qc/{{sample}}.geneBodyCoverage.curves.pdf"
    params:
    shell:
        """
        mkdir -p $(dirname {output.txt})
        geneBody_coverage.py -i {input.bam} -r {input.bed} -o $(dirname {output.txt})/{wildcards.sample}
        """

# Aggregate gene counts across all RNA-seq samples
rule rnaseq_aggregate_counts:
    input:
        expand(f"{GENOME_OUTDIR}/rnaseq/counts/{{sample}}.ReadsPerGene.out.tab", sample=RNASEQ_SAMPLES)
    output:
        f"{GENOME_OUTDIR}/rnaseq/counts/all_samples_counts.txt"
    params:
    shell:
        """
        mkdir -p $(dirname {output})
        python3 -c "
import pandas as pd
import sys

# Read sample files
count_files = sys.argv[1].split()
dfs = []

for f in count_files:
    sample = f.split('/')[-1].replace('.ReadsPerGene.out.tab', '')
    df = pd.read_csv(f, sep='\t', header=None, names=['gene_id', 'unstranded', 'strand1', 'strand2'])
    df = df[4:]  # Skip header lines
    df[sample] = df['unstranded']  # Use unstranded counts
    dfs.append(df[['gene_id', sample]])

# Merge all samples
result = dfs[0]
for df in dfs[1:]:
    result = result.merge(df, on='gene_id', how='outer')

result.to_csv(sys.argv[2], sep='\t', index=False)
        " "$(echo {input})" {output}
        """

# Complete RNA-seq pipeline marker
rule rnaseq_complete:
    input:
        expand(f"{GENOME_OUTDIR}/rnaseq/bam/{{sample}}.Aligned.sortedByCoord.out.bam", sample=RNASEQ_SAMPLES),
        expand(f"{GENOME_OUTDIR}/rnaseq/bigwig/{{sample}}.bw", sample=RNASEQ_SAMPLES),
        f"{GENOME_OUTDIR}/rnaseq/counts/all_samples_counts.txt"
    output:
        f"{GENOME_OUTDIR}/rnaseq/tracks_complete.txt"
    params:
    shell:
        "echo 'RNA-seq pipeline completed at $(date)' > {output}"
