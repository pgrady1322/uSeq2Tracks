# ============================================================
# ChIP-seq analysis rules
# ============================================================

rule chipseq_map_bowtie2:
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

rule chipseq_map_bwa_mem2:
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
        
        samtools index {output.bam}
        
        # Clean up temporary directory
        rm -rf {resources.tmpdir}/samtools_sort_{wildcards.sample}
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
        # Mark duplicates but don't remove them for ChIP-seq (use -s flag)
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
        set +e  # Don't exit on errors in subshells
        
        # Check if this is a control/input sample
        sample_condition="TREATMENT"
        
        # Check sample name for input/control keywords
        if [[ "{wildcards.sample}" == *"input"* ]] || [[ "{wildcards.sample}" == *"control"* ]] || [[ "{wildcards.sample}" == *"Input"* ]] || [[ "{wildcards.sample}" == *"Control"* ]]; then
            sample_condition="CONTROL"
            echo "Detected control sample from name: {wildcards.sample}"
        fi
        
        # Also check CSV if available (but don't fail if it doesn't work)
        if [[ "$sample_condition" == "TREATMENT" ]]; then
            csv_condition=$(python3 -c "
import pandas as pd
try:
    samples = pd.read_csv('{config[samplesheet]}')
    sample_info = samples[samples['sample'] == '{wildcards.sample}']
    if not sample_info.empty:
        condition = str(sample_info.iloc[0].get('condition', '')).lower()
        if 'input' in condition or 'control' in condition:
            print('CONTROL')
        else:
            print('TREATMENT')
    else:
        print('TREATMENT')
except:
    print('TREATMENT')
" 2>/dev/null)
            
            if [[ "$csv_condition" == "CONTROL" ]]; then
                sample_condition="CONTROL"
                echo "Detected control sample from CSV: {wildcards.sample}"
            fi
        fi
        
        set -e  # Re-enable exit on error
        
        echo "Sample: {wildcards.sample}"
        echo "Final condition: $sample_condition"
        
        # Create output directory
        mkdir -p {params.outdir}
        
        if [[ "$sample_condition" == "CONTROL" ]]; then
            # This is a control sample - create empty peak files
            echo "# Control/Input sample - no peaks called" > {output.peaks}
            echo "# Sample: {wildcards.sample}" >> {output.peaks}
            echo "# Control samples are used as background for other ChIP samples" >> {output.peaks}
            
            touch {output.summits}
            
            echo "# Control/Input sample" > {output.xls}
            echo "# Sample: {wildcards.sample}" >> {output.xls}
            echo "# No peaks called for control samples" >> {output.xls}
            
            echo "Successfully created empty peak files for control sample: {wildcards.sample}"
            exit 0
        fi
        
        # This is a treatment sample - proceed with peak calling
        echo "Processing treatment sample: {wildcards.sample}"
        
        # Check fragment count
        fragment_count=$(samtools view -c {input.treatment})
        echo "Fragment count: $fragment_count"
        
        if [[ "$fragment_count" -lt 1000 ]]; then
            echo "Warning: Too few fragments ($fragment_count) for reliable peak calling"
            echo "# Insufficient fragments for peak calling" > {output.peaks}
            touch {output.summits}
            echo "# Insufficient fragments ($fragment_count)" > {output.xls}
            exit 0
        fi
        
        # Check if BAM has paired-end data
        paired_count=$(samtools view -f 1 -c {input.treatment})
        echo "Paired-end reads: $paired_count"
        
        # Prepare MACS3 options
        base_opts="{params.macs_opts}"
        
        # Clean up any conflicting options
        base_opts=$(echo "$base_opts" | sed 's/--keep-dup[[:space:]]\+[^[:space:]]*//g')
        base_opts=$(echo "$base_opts" | sed 's/--format[[:space:]]\+[^[:space:]]*//g')
        base_opts=$(echo "$base_opts" | sed 's/[[:space:]]\+/ /g' | sed 's/^[[:space:]]*//' | sed 's/[[:space:]]*$//')
        
        # Add our options
        base_opts="$base_opts --keep-dup all"
        
        if [[ $paired_count -gt 0 ]]; then
            echo "Detected paired-end data, using BAMPE format"
            base_opts="$base_opts --format BAMPE"
        else
            echo "Detected single-end data, using BAM format"
            base_opts="$base_opts --format BAM"
        fi
        
        echo "MACS3 options: $base_opts"
        
        # Run MACS3
        control_file="{input.control}"
        
        if [[ -n "$control_file" && "$control_file" != "" && -f "$control_file" ]]; then
            echo "Using control file: $control_file"
            macs3 callpeak -t {input.treatment} -c "$control_file" -n {params.name} --outdir {params.outdir} $base_opts
        else
            echo "No control file provided, running without control"
            macs3 callpeak -t {input.treatment} -n {params.name} --outdir {params.outdir} $base_opts
        fi
        
        echo "Peak calling completed successfully for {wildcards.sample}"
        """

# Parameter sweep: Call peaks with multiple q-values for ChIP-seq
rule chipseq_call_peaks_sweep:
    input:
        treatment=f"{config['outdir']}/chipseq/bam/{{sample}}.dedup.bam" if config["chipseq"]["markdup"] else f"{config['outdir']}/chipseq/bam/{{sample}}.sorted.bam",
        control=get_control_for_sample
    output:
        peaks=f"{config['outdir']}/chipseq/peaks_sweep/{{sample}}_q{{qval}}_peaks.narrowPeak",
        summits=f"{config['outdir']}/chipseq/peaks_sweep/{{sample}}_q{{qval}}_summits.bed",
        xls=f"{config['outdir']}/chipseq/peaks_sweep/{{sample}}_q{{qval}}_peaks.xls"
    params:
        name="{sample}_q{qval}",
        outdir=f"{config['outdir']}/chipseq/peaks_sweep",
        qval="{qval}"
    shell:
        """
        set +e  # Don't exit on errors in subshells
        
        # Check if this is a control/input sample
        sample_condition="TREATMENT"
        
        # Check sample name for input/control keywords
        if [[ "{wildcards.sample}" == *"input"* ]] || [[ "{wildcards.sample}" == *"control"* ]] || [[ "{wildcards.sample}" == *"Input"* ]] || [[ "{wildcards.sample}" == *"Control"* ]]; then
            sample_condition="CONTROL"
            echo "Detected control sample from name: {wildcards.sample}"
        fi
        
        # Also check CSV if available (but don't fail if it doesn't work)
        if [[ "$sample_condition" == "TREATMENT" ]]; then
            csv_condition=$(python3 -c "
import pandas as pd
try:
    samples = pd.read_csv('{config[samplesheet]}')
    sample_info = samples[samples['sample'] == '{wildcards.sample}']
    if not sample_info.empty:
        condition = str(sample_info.iloc[0].get('condition', '')).lower()
        if 'input' in condition or 'control' in condition:
            print('CONTROL')
        else:
            print('TREATMENT')
    else:
        print('TREATMENT')
except:
    print('TREATMENT')
" 2>/dev/null)
            
            if [[ "$csv_condition" == "CONTROL" ]]; then
                sample_condition="CONTROL"
                echo "Detected control sample from CSV: {wildcards.sample}"
            fi
        fi
        
        set -e  # Re-enable exit on error
        
        echo "Sample: {wildcards.sample}"
        echo "Final condition: $sample_condition"
        echo "Parameter sweep: Running with q-value = {params.qval}"
        
        # Create output directory
        mkdir -p {params.outdir}
        
        if [[ "$sample_condition" == "CONTROL" ]]; then
            # This is a control sample - create empty peak files
            echo "# Control/Input sample - no peaks called" > {output.peaks}
            echo "# Sample: {wildcards.sample}" >> {output.peaks}
            echo "# Q-value: {params.qval}" >> {output.peaks}
            echo "# Control samples are used as background for other ChIP samples" >> {output.peaks}
            
            touch {output.summits}
            
            echo "# Control/Input sample" > {output.xls}
            echo "# Sample: {wildcards.sample}" >> {output.xls}
            echo "# Q-value: {params.qval}" >> {output.xls}
            echo "# No peaks called for control samples" >> {output.xls}
            
            echo "Successfully created empty peak files for control sample: {wildcards.sample}"
            exit 0
        fi
        
        # This is a treatment sample - proceed with peak calling
        echo "Processing treatment sample: {wildcards.sample}"
        
        # Check fragment count
        fragment_count=$(samtools view -c {input.treatment})
        echo "Fragment count: $fragment_count"
        
        if [[ "$fragment_count" -lt 1000 ]]; then
            echo "Warning: Too few fragments ($fragment_count) for reliable peak calling"
            echo "# Insufficient fragments for peak calling" > {output.peaks}
            echo "# Q-value: {params.qval}" >> {output.peaks}
            touch {output.summits}
            echo "# Insufficient fragments ($fragment_count)" > {output.xls}
            exit 0
        fi
        
        # Check if BAM has paired-end data
        paired_count=$(samtools view -f 1 -c {input.treatment})
        echo "Paired-end reads: $paired_count"
        
        # Prepare MACS3 options with specific q-value
        base_opts="--qval {params.qval} --keep-dup all"
        
        if [[ $paired_count -gt 0 ]]; then
            echo "Detected paired-end data, using BAMPE format"
            base_opts="$base_opts --format BAMPE"
        else
            echo "Detected single-end data, using BAM format"
            base_opts="$base_opts --format BAM"
        fi
        
        echo "MACS3 options: $base_opts"
        
        # Run MACS3
        control_file="{input.control}"
        
        if [[ -n "$control_file" && "$control_file" != "" && -f "$control_file" ]]; then
            echo "Using control file: $control_file"
            macs3 callpeak -t {input.treatment} -c "$control_file" -n {params.name} --outdir {params.outdir} $base_opts
        else
            echo "No control file provided, running without control"
            macs3 callpeak -t {input.treatment} -n {params.name} --outdir {params.outdir} $base_opts
        fi
        
        echo "Parameter sweep peak calling completed successfully for {wildcards.sample} with q-value {params.qval}"
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