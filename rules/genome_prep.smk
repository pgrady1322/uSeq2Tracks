# ============================================================
# Genome preparation and indexing rules
# ============================================================

rule copy_genome:
    input:
        config["genome"]
    output:
        f"{GENOME_OUTDIR}/genome/genome.fa"
    params:
    shell:
        "cp {input} {output}"

rule index_genome_fai:
    input:
        f"{GENOME_OUTDIR}/genome/genome.fa"
    output:
        f"{GENOME_OUTDIR}/genome/genome.fa.fai"
    params:
    shell:
        "samtools faidx {input}"

# Create UCSC chromosome sizes file from FAI
rule genome_sizes:
    input:
        f"{GENOME_OUTDIR}/genome/genome.fa.fai"
    output:
        f"{GENOME_OUTDIR}/genome/genome.sizes"
    params:
    shell:
        "cut -f1,2 {input} > {output}"

rule index_genome_dict:
    input:
        f"{GENOME_OUTDIR}/genome/genome.fa"
    output:
        f"{GENOME_OUTDIR}/genome/genome.dict"
    params:
    shell:
        "samtools dict {input} > {output}"

# BWA-MEM2 index (for WGS, ChIP-seq, CUT&RUN, ATAC-seq)
rule index_bwa_mem2:
    input:
        f"{GENOME_OUTDIR}/genome/genome.fa"
    output:
        f"{GENOME_OUTDIR}/genome/bwa_mem2/genome.fa.0123",
        f"{GENOME_OUTDIR}/genome/bwa_mem2/genome.fa.amb",
        f"{GENOME_OUTDIR}/genome/bwa_mem2/genome.fa.ann",
        f"{GENOME_OUTDIR}/genome/bwa_mem2/genome.fa.bwt.2bit.64",
        f"{GENOME_OUTDIR}/genome/bwa_mem2/genome.fa.pac"
    params:
        prefix=f"{GENOME_OUTDIR}/genome/bwa_mem2/genome.fa"
    threads: config["resources"]["index_threads"]
    shell:
        """
        mkdir -p {GENOME_OUTDIR}/genome/bwa_mem2
        cp {input} {params.prefix}
        bwa-mem2 index {params.prefix}
        """

# Bowtie2 index (for ChIP-seq, CUT&RUN, ATAC-seq)
rule index_bowtie2:
    input:
        f"{GENOME_OUTDIR}/genome/genome.fa"
    output:
        expand(f"{GENOME_OUTDIR}/genome/bowtie2/genome.{{ext}}", 
               ext=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])
    params:
        prefix=f"{GENOME_OUTDIR}/genome/bowtie2/genome",
        outdir=GENOME_OUTDIR
    threads: config["resources"]["index_threads"]
    shell:
        """
        mkdir -p {params.outdir}/genome/bowtie2
        bowtie2-build --threads {threads} {input} {params.prefix}
        """

# STAR index (for RNA-seq)
rule index_star:
    input:
        f"{GENOME_OUTDIR}/genome/genome.fa"
    output:
        f"{GENOME_OUTDIR}/genome/star/Genome",
        f"{GENOME_OUTDIR}/genome/star/SA",
        f"{GENOME_OUTDIR}/genome/star/SAindex"
    params:
        genomeDir=f"{GENOME_OUTDIR}/genome/star",
        gtf_opt=f"--sjdbGTFfile {config['gtf']}" if config.get('gtf') and config['gtf'] else ""
    threads: config["resources"]["index_threads"]
    shell:
        """
        mkdir -p {params.genomeDir}
        STAR --runMode genomeGenerate \
             --genomeDir {params.genomeDir} \
             --genomeFastaFiles {input} \
             --runThreadN {threads} \
             {params.gtf_opt}
        """

# HISAT2 index (alternative for RNA-seq)
rule index_hisat2:
    input:
        f"{GENOME_OUTDIR}/genome/genome.fa"
    output:
        expand(f"{GENOME_OUTDIR}/genome/hisat2/genome.{{ext}}", 
               ext=["1.ht2", "2.ht2", "3.ht2", "4.ht2", "5.ht2", "6.ht2", "7.ht2", "8.ht2"])
    params:
        prefix=f"{GENOME_OUTDIR}/genome/hisat2/genome",
        outdir=GENOME_OUTDIR
    threads: config["resources"]["index_threads"]
    shell:
        """
        mkdir -p {params.outdir}/genome/hisat2
        hisat2-build -p {threads} {input} {params.prefix}
        """

# Minimap2 index (for long reads)
rule index_minimap2:
    input:
        f"{GENOME_OUTDIR}/genome/genome.fa"
    output:
        f"{GENOME_OUTDIR}/genome/minimap2/genome.mmi"
    params:
        outdir=GENOME_OUTDIR
    threads: config["resources"]["index_threads"]
    shell:
        """
        mkdir -p {params.outdir}/genome/minimap2
        minimap2 -t {threads} -d {output} {input}
        """

# BWA index (for ancient DNA)
rule index_bwa_aln:
    input:
        f"{GENOME_OUTDIR}/genome/genome.fa"
    output:
        f"{GENOME_OUTDIR}/genome/bwa_aln/genome.fa.amb",
        f"{GENOME_OUTDIR}/genome/bwa_aln/genome.fa.ann",
        f"{GENOME_OUTDIR}/genome/bwa_aln/genome.fa.bwt",
        f"{GENOME_OUTDIR}/genome/bwa_aln/genome.fa.pac",
        f"{GENOME_OUTDIR}/genome/bwa_aln/genome.fa.sa"
    params:
        prefix=f"{GENOME_OUTDIR}/genome/bwa_aln/genome.fa",
        outdir=GENOME_OUTDIR
    shell:
        """
        mkdir -p {params.outdir}/genome/bwa_aln
        cp {input} {params.prefix}
        bwa index {params.prefix}
        """
