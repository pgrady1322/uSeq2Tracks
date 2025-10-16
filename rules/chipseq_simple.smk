# ============================================================
# ChIP-seq analysis rules - SIMPLIFIED
# ============================================================

rule chipseq_simple_map_bowtie2:
    input:
        unpack(get_input_reads),
        index=f"{GENOME_OUTDIR}/genome/bowtie2/genome.1.bt2"
    output:
        bam=f"{config['outdir']}/chipseq/bam/{{sample}}.raw.bam"
    params:
        index_prefix=f"{GENOME_OUTDIR}/genome/bowtie2/genome",
        bowtie2_opts=config["chipseq"]["bowtie2_opts"],
        outdir=config['outdir']
    threads: config["resources"]["mapping_threads"]
    run:
        # Check if we have paired-end or single-end data
        has_r2 = hasattr(input, 'r2') and input.r2
        
        shell("mkdir -p {params.outdir}/chipseq/bam")
        
        if has_r2:
            # Paired-end
            shell("""
            bowtie2 {params.bowtie2_opts} \
                    --threads {threads} \
                    -x {params.index_prefix} \
                    -1 {input.r1} -2 {input.r2} \
            | samtools view -bS - > {output.bam}
            """)
        else:
            # Single-end
            shell("""
            bowtie2 {params.bowtie2_opts} \
                    --threads {threads} \
                    -x {params.index_prefix} \
                    -U {input.r1} \
            | samtools view -bS - > {output.bam}
            """)

rule chipseq_simple_map_bwa_mem2:
    input:
        unpack(get_input_reads),
        index=f"{GENOME_OUTDIR}/genome/bwa_mem2/genome.fa.0123"
    output:
        bam=f"{config['outdir']}/chipseq/bam/{{sample}}.raw.bam"
    params:
        index_prefix=f"{GENOME_OUTDIR}/genome/bwa_mem2/genome.fa",
        bwa_opts=config["chipseq"]["bwa_mem2_opts"],
        outdir=config['outdir']
    threads: config["resources"]["mapping_threads"]
    run:
        # Check if we have paired-end or single-end data
        has_r2 = hasattr(input, 'r2') and input.r2
        
        shell("mkdir -p {params.outdir}/chipseq/bam")
        
        if has_r2:
            # Paired-end
            shell("""
            bwa-mem2 mem {params.bwa_opts} \
                     -t {threads} \
                     {params.index_prefix} \
                     {input.r1} {input.r2} \
            | samtools view -bS - > {output.bam}
            """)
        else:
            # Single-end  
            shell("""
            bwa-mem2 mem {params.bwa_opts} \
                     -t {threads} \
                     {params.index_prefix} \
                     {input.r1} \
            | samtools view -bS - > {output.bam}
            """)

# Use mapper selection based on config
if config["chipseq"]["mapper"] == "bowtie2":
    ruleorder: chipseq_map_bowtie2 > chipseq_map_bwa_mem2
else:
    ruleorder: chipseq_map_bwa_mem2 > chipseq_map_bowtie2

rule chipseq_sort_and_index:
    input:
        f"{config['outdir']}/chipseq/bam/{{sample}}.raw.bam"
    output:
        bam=f"{config['outdir']}/chipseq/bam/{{sample}}.sorted.bam",
        bai=f"{config['outdir']}/chipseq/bam/{{sample}}.sorted.bam.bai"
    threads: config["resources"]["sort_threads"]
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input}
        samtools index {output.bam}
        """

# Optional duplicate marking for ChIP-seq
rule chipseq_markdup:
    input:
        f"{config['outdir']}/chipseq/bam/{{sample}}.sorted.bam"
    output:
        bam=f"{config['outdir']}/chipseq/bam/{{sample}}.dedup.bam",
        bai=f"{config['outdir']}/chipseq/bam/{{sample}}.dedup.bam.bai",
        metrics=f"{config['outdir']}/chipseq/qc/{{sample}}.markdup.metrics"
    params:
        outdir=config['outdir']
    threads: config["resources"]["sort_threads"]
    shell:
        """
        mkdir -p {params.outdir}/chipseq/qc
        samtools markdup -@ {threads} -f {output.metrics} -s {input} {output.bam}
        samtools index -@ {threads} {output.bam}
        """

rule chipseq_bam_to_bigwig:
    input:
        bam=f"{config['outdir']}/chipseq/bam/{{sample}}.dedup.bam" if config["chipseq"]["markdup"] else f"{config['outdir']}/chipseq/bam/{{sample}}.sorted.bam",
        bai=f"{config['outdir']}/chipseq/bam/{{sample}}.dedup.bam.bai" if config["chipseq"]["markdup"] else f"{config['outdir']}/chipseq/bam/{{sample}}.sorted.bam.bai",
        fai=f"{GENOME_OUTDIR}/genome/genome.fa.fai",
        sizes=f"{GENOME_OUTDIR}/genome/genome.sizes"
    output:
        f"{config['outdir']}/chipseq/bigwig/{{sample}}.bw"
    params:
        outdir=config['outdir']
    threads: config["resources"]["sort_threads"]
    shell:
        """
        mkdir -p {params.outdir}/chipseq/bigwig
        bamCoverage --bam {input.bam} \
                    --outFileName {output} \
                    --outFileFormat bigwig \
                    --numberOfProcessors {threads} \
                    --normalizeUsing RPKM
        """

rule chipseq_call_peaks:
    input:
        treatment=f"{config['outdir']}/chipseq/bam/{{sample}}.dedup.bam" if config["chipseq"]["markdup"] else f"{config['outdir']}/chipseq/bam/{{sample}}.sorted.bam",
        control=get_control_for_sample
    output:
        peaks=f"{config['outdir']}/chipseq/peaks/{{sample}}_peaks.narrowPeak",
        summits=f"{config['outdir']}/chipseq/peaks/{{sample}}_summits.bed",
        xls=f"{config['outdir']}/chipseq/peaks/{{sample}}_peaks.xls"
    params:
        name="{sample}",
        outdir=f"{config['outdir']}/chipseq/peaks",
        macs_opts=config["chipseq"]["macs3_opts"]
    shell:
        """
        mkdir -p {params.outdir}
        
        # Skip peak calling for input/control samples
        if [[ "{wildcards.sample}" == *"input"* ]] || [[ "{wildcards.sample}" == *"control"* ]]; then
            echo "# Control/Input sample - no peaks called" > {output.peaks}
            touch {output.summits}
            echo "# Control/Input sample" > {output.xls}
            echo "Skipped peak calling for control sample: {wildcards.sample}"
            exit 0
        fi
        
        # Basic fragment count check
        fragment_count=$(samtools view -c {input.treatment})
        if [[ "$fragment_count" -lt 1000 ]]; then
            echo "# Too few fragments ($fragment_count)" > {output.peaks}
            touch {output.summits}
            echo "# Too few fragments" > {output.xls}
            echo "Warning: Too few fragments for {wildcards.sample}"
            exit 0
        fi
        
        # Auto-detect paired vs single-end
        paired_count=$(samtools view -f 1 -c {input.treatment})
        if [[ $paired_count -gt 0 ]]; then
            format_opt="--format BAMPE"
        else
            format_opt="--format BAM"
        fi
        
        # Build MACS3 command
        macs_cmd="macs3 callpeak -t {input.treatment} -n {params.name} --outdir {params.outdir} {params.macs_opts} --keep-dup all $format_opt"
        
        # Add control if available
        control_file="{input.control}"
        if [[ -n "$control_file" && -f "$control_file" ]]; then
            macs_cmd="$macs_cmd -c $control_file"
            echo "Using control: $control_file"
        else
            echo "No control file - running without control"
        fi
        
        echo "Running: $macs_cmd"
        eval $macs_cmd
        
        echo "Peak calling completed for {wildcards.sample}"
        """

rule chipseq_complete:
    input:
        bam=expand(f"{config['outdir']}/chipseq/bam/{{sample}}.dedup.bam", sample=CHIPSEQ_SAMPLES) if config["chipseq"]["markdup"] else expand(f"{config['outdir']}/chipseq/bam/{{sample}}.sorted.bam", sample=CHIPSEQ_SAMPLES),
        bw=expand(f"{config['outdir']}/chipseq/bigwig/{{sample}}.bw", sample=CHIPSEQ_SAMPLES),
        peaks=expand(f"{config['outdir']}/chipseq/peaks/{{sample}}_peaks.narrowPeak", sample=CHIPSEQ_SAMPLES)
    output:
        f"{GENOME_OUTDIR}/chipseq/tracks_complete.txt"
    shell:
        """
        echo "ChIP-seq analysis complete" > {output}
        echo "Samples processed: {CHIPSEQ_SAMPLES}" >> {output}
        echo "Generated BAM files: {input.bam}" >> {output}
        echo "Generated BigWig files: {input.bw}" >> {output}
        echo "Generated peak files: {input.peaks}" >> {output}
        """
