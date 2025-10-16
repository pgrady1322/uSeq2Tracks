# ============================================================
# Genrich Peak Calling Rules
# ============================================================

def get_genrich_samples_by_replicate_group(assay):
    """Get samples grouped by replicate groups for a specific assay"""
    if assay not in REPLICATE_GROUPS:
        return {}
    return REPLICATE_GROUPS[assay]

def get_genrich_treatment_bams_atacseq(wildcards):
    """Get ATAC-seq treatment BAM files for Genrich"""
    group = wildcards.group
    
    if 'atacseq' not in REPLICATE_GROUPS or group not in REPLICATE_GROUPS['atacseq']:
        return []
    
    samples = REPLICATE_GROUPS['atacseq'][group]
    
    # Determine BAM file path based on markdup setting
    bam_files = []
    for sample in samples:
        if config.get('atacseq', {}).get('markdup', False):
            bam_files.append(f"{config['outdir']}/atacseq/bam/{sample}.dedup.bam")
        else:
            bam_files.append(f"{config['outdir']}/atacseq/bam/{sample}.sorted.bam")
    
    return bam_files

def get_genrich_control_bams_atacseq(wildcards):
    """Get ATAC-seq control BAM files for Genrich (usually none for ATAC-seq)"""
    return []

def get_genrich_treatment_bams_cutrun(wildcards):
    """Get CUT&RUN treatment BAM files for Genrich"""
    group = wildcards.group
    
    if 'cutrun' not in REPLICATE_GROUPS or group not in REPLICATE_GROUPS['cutrun']:
        return []
    
    samples = REPLICATE_GROUPS['cutrun'][group]
    
    # Determine BAM file path based on markdup setting
    bam_files = []
    for sample in samples:
        if config.get('cutrun', {}).get('markdup', False):
            bam_files.append(f"{config['outdir']}/cutrun/bam/{sample}.dedup.bam")
        else:
            bam_files.append(f"{config['outdir']}/cutrun/bam/{sample}.sorted.bam")
    
    return bam_files

def get_genrich_control_bams_cutrun(wildcards):
    """Get CUT&RUN control BAM files for Genrich"""
    group = wildcards.group
    control_samples = []
    
    # Look for IgG or Input controls in CUT&RUN
    if 'cutrun' in REPLICATE_GROUPS:
        for control_group, samples in REPLICATE_GROUPS['cutrun'].items():
            if any(keyword in control_group.lower() for keyword in ['igg', 'input']):
                control_samples.extend(samples)
    
    if not control_samples:
        return []
    
    # Determine BAM file path based on markdup setting
    bam_files = []
    for sample in control_samples:
        if config.get('cutrun', {}).get('markdup', False):
            bam_files.append(f"{config['outdir']}/cutrun/bam/{sample}.dedup.bam")
        else:
            bam_files.append(f"{config['outdir']}/cutrun/bam/{sample}.sorted.bam")
    
    return bam_files

def get_genrich_treatment_bams_chipseq(wildcards):
    """Get ChIP-seq treatment BAM files for Genrich"""
    group = wildcards.group
    
    if 'chipseq' not in REPLICATE_GROUPS or group not in REPLICATE_GROUPS['chipseq']:
        return []
    
    samples = REPLICATE_GROUPS['chipseq'][group]
    
    # Determine BAM file path based on markdup setting
    bam_files = []
    for sample in samples:
        if config.get('chipseq', {}).get('markdup', False):
            bam_files.append(f"{config['outdir']}/chipseq/bam/{sample}.dedup.bam")
        else:
            bam_files.append(f"{config['outdir']}/chipseq/bam/{sample}.sorted.bam")
    
    return bam_files

def get_genrich_control_bams_chipseq(wildcards):
    """Get ChIP-seq control BAM files for Genrich"""
    group = wildcards.group
    control_samples = []
    
    # Look for Input or Control samples in ChIP-seq
    if 'chipseq' in REPLICATE_GROUPS:
        for control_group, samples in REPLICATE_GROUPS['chipseq'].items():
            if any(keyword in control_group.lower() for keyword in ['input', 'control']):
                control_samples.extend(samples)
    
    if not control_samples:
        return []
    
    # Determine BAM file path based on markdup setting
    bam_files = []
    for sample in control_samples:
        if config.get('chipseq', {}).get('markdup', False):
            bam_files.append(f"{config['outdir']}/chipseq/bam/{sample}.dedup.bam")
        else:
            bam_files.append(f"{config['outdir']}/chipseq/bam/{sample}.sorted.bam")
    
    return bam_files

# ATAC-seq Genrich peak calling with replicate merging
rule atacseq_genrich_peaks:
    input:
        treatment = get_genrich_treatment_bams_atacseq,
        control = get_genrich_control_bams_atacseq
    output:
        peaks = f"{config['outdir']}/atacseq/genrich_peaks/{{group}}_genrich.narrowPeak"
    params:
        atac_flag = "-j" if config.get('genrich', {}).get('atac_mode', True) else "",
        genrich_opts = config.get('genrich', {}).get('options', '-q 0.05'),
        outdir = f"{config['outdir']}/atacseq/genrich_peaks"
    shell:
        """
        mkdir -p {params.outdir}
        
        # Check if Genrich is available
        if ! command -v Genrich &> /dev/null; then
            echo "ERROR: Genrich command not found"
            echo "Please install with: conda install bioconda::genrich"
            exit 1
        fi
        
        echo "Using Genrich for ATAC-seq peak calling on replicate group: {wildcards.group}"
        echo "Treatment files: {input.treatment}"
        
        # Build Genrich command
        cmd="Genrich -t {input.treatment}"
        
        # Add controls if available (usually none for ATAC-seq)
        if [ ! -z "{input.control}" ]; then
            echo "Control files: {input.control}"
            cmd="$cmd -c {input.control}"
        fi
        
        # Add ATAC-seq flag and other options
        cmd="$cmd {params.atac_flag} {params.genrich_opts} -o {output.peaks}"
        
        echo "Running: $cmd"
        eval $cmd
        
        echo "Genrich ATAC-seq peak calling completed for {wildcards.group}"
        """

# CUT&RUN Genrich peak calling with replicate merging
rule cutrun_genrich_peaks:
    input:
        treatment = get_genrich_treatment_bams_cutrun,
        control = get_genrich_control_bams_cutrun
    output:
        peaks = f"{config['outdir']}/cutrun/genrich_peaks/{{group}}_genrich.narrowPeak"
    params:
        atac_flag = "-j" if config.get('genrich', {}).get('atac_mode', True) else "",
        genrich_opts = config.get('genrich', {}).get('options', '-q 0.05'),
        outdir = f"{config['outdir']}/cutrun/genrich_peaks"
    shell:
        """
        mkdir -p {params.outdir}
        
        # Check if Genrich is available
        if ! command -v Genrich &> /dev/null; then
            echo "ERROR: Genrich command not found"
            echo "Please install with: conda install bioconda::genrich"
            exit 1
        fi
        
        echo "Using Genrich for CUT&RUN peak calling on replicate group: {wildcards.group}"
        echo "Treatment files: {input.treatment}"
        
        # Build Genrich command
        cmd="Genrich -t {input.treatment}"
        
        # Add controls if available
        if [ ! -z "{input.control}" ]; then
            echo "Control files: {input.control}"
            cmd="$cmd -c {input.control}"
        fi
        
        # Add ATAC-seq flag (also good for CUT&RUN) and other options
        cmd="$cmd {params.atac_flag} {params.genrich_opts} -o {output.peaks}"
        
        echo "Running: $cmd"
        eval $cmd
        
        echo "Genrich CUT&RUN peak calling completed for {wildcards.group}"
        """

# ChIP-seq Genrich peak calling with replicate merging
rule chipseq_genrich_peaks:
    input:
        treatment = get_genrich_treatment_bams_chipseq,
        control = get_genrich_control_bams_chipseq
    output:
        peaks = f"{config['outdir']}/chipseq/genrich_peaks/{{group}}_genrich.narrowPeak"
    params:
        genrich_opts = config.get('genrich', {}).get('options', '-q 0.05'),
        outdir = f"{config['outdir']}/chipseq/genrich_peaks"
    shell:
        """
        mkdir -p {params.outdir}
        
        # Check if Genrich is available
        if ! command -v Genrich &> /dev/null; then
            echo "ERROR: Genrich command not found"
            echo "Please install with: conda install bioconda::genrich"
            exit 1
        fi
        
        echo "Using Genrich for ChIP-seq peak calling on replicate group: {wildcards.group}"
        echo "Treatment files: {input.treatment}"
        
        # Build Genrich command
        cmd="Genrich -t {input.treatment}"
        
        # Add controls if available
        if [ ! -z "{input.control}" ]; then
            echo "Control files: {input.control}"
            cmd="$cmd -c {input.control}"
        fi
        
        # ChIP-seq doesn't use -j flag
        cmd="$cmd {params.genrich_opts} -o {output.peaks}"
        
        echo "Running: $cmd"
        eval $cmd
        
        echo "Genrich ChIP-seq peak calling completed for {wildcards.group}"
        """

# Generate all Genrich peaks for all assays
rule genrich_all_peaks:
    input:
        atacseq = expand(f"{config['outdir']}/atacseq/genrich_peaks/{{group}}_genrich.narrowPeak", 
                        group=[g for g in get_genrich_samples_by_replicate_group('atacseq').keys()]) if ATACSEQ_SAMPLES and config.get('genrich', {}).get('enabled', False) else [],
        cutrun = expand(f"{config['outdir']}/cutrun/genrich_peaks/{{group}}_genrich.narrowPeak", 
                       group=[g for g in get_genrich_samples_by_replicate_group('cutrun').keys()]) if CUTRUN_SAMPLES and config.get('genrich', {}).get('enabled', False) else [],
        chipseq = expand(f"{config['outdir']}/chipseq/genrich_peaks/{{group}}_genrich.narrowPeak", 
                        group=[g for g in get_genrich_samples_by_replicate_group('chipseq').keys()]) if CHIPSEQ_SAMPLES and config.get('genrich', {}).get('enabled', False) else []
    output:
        f"{config['outdir']}/genrich/genrich_complete.txt"
    shell:
        """
        mkdir -p $(dirname {output})
        echo "Genrich peak calling complete" > {output}
        echo "ATAC-seq peaks: {input.atacseq}" >> {output}
        echo "CUT&RUN peaks: {input.cutrun}" >> {output}
        echo "ChIP-seq peaks: {input.chipseq}" >> {output}
        echo "Generated on: $(date)" >> {output}
        """

# Rule to run only Genrich analysis
rule genrich:
    input:
        f"{config['outdir']}/genrich/genrich_complete.txt"
    message:
        "Genrich peak calling analysis complete"

# ============================================================
# Genrich Parameter Sweep Rules
# ============================================================

def get_genrich_sweep_outputs(wildcards=None):
    """Generate output files for Genrich parameter sweep"""
    outputs = []
    if not config.get('parameter_sweep', {}).get('genrich_enabled', False):
        return outputs
    
    qvalues = config.get('parameter_sweep', {}).get('genrich_qvalues', [0.05, 0.01, 0.005, 0.001])
    
    for assay in ['atacseq', 'cutrun', 'chipseq']:
        if assay not in REPLICATE_GROUPS:
            continue
            
        assay_samples = globals().get(f"{assay.upper()}_SAMPLES", [])
        if not assay_samples:
            continue
            
        for group in REPLICATE_GROUPS[assay].keys():
            for qval in qvalues:
                # Convert qvalue to filename format (0.05 -> q005)
                qval_str = f"q{str(qval).replace('0.', '').zfill(3)}"
                outputs.append(f"{config['outdir']}/{assay}/genrich_sweep/{group}_genrich_{qval_str}.narrowPeak")
    
    return outputs

# ATAC-seq Genrich parameter sweep
rule atacseq_genrich_sweep:
    input:
        treatment = get_genrich_treatment_bams_atacseq,
        control = get_genrich_control_bams_atacseq
    output:
        peaks = f"{config['outdir']}/atacseq/genrich_sweep/{{group}}_genrich_{{qval}}.narrowPeak"
    params:
        atac_flag = "-j" if config.get('genrich', {}).get('atac_mode', True) else "",
        outdir = f"{config['outdir']}/atacseq/genrich_sweep"
    shell:
        """
        mkdir -p {params.outdir}
        
        # Convert qvalue format (q005 -> 0.05)
        qvalue=$(echo "{wildcards.qval}" | sed 's/q/0./' | sed 's/^0\./0.0/' | sed 's/00$//')
        
        echo "Using Genrich for ATAC-seq parameter sweep on group: {wildcards.group}, q-value: $qvalue"
        
        # Build Genrich command
        cmd="Genrich -t {input.treatment}"
        
        if [ ! -z "{input.control}" ]; then
            cmd="$cmd -c {input.control}"
        fi
        
        cmd="$cmd {params.atac_flag} -q $qvalue -o {output.peaks}"
        
        echo "Running: $cmd"
        eval $cmd
        """

# CUT&RUN Genrich parameter sweep
rule cutrun_genrich_sweep:
    input:
        treatment = get_genrich_treatment_bams_cutrun,
        control = get_genrich_control_bams_cutrun
    output:
        peaks = f"{config['outdir']}/cutrun/genrich_sweep/{{group}}_genrich_{{qval}}.narrowPeak"
    params:
        atac_flag = "-j" if config.get('genrich', {}).get('atac_mode', True) else "",
        outdir = f"{config['outdir']}/cutrun/genrich_sweep"
    shell:
        """
        mkdir -p {params.outdir}
        
        # Convert qvalue format (q005 -> 0.05)
        qvalue=$(echo "{wildcards.qval}" | sed 's/q/0./' | sed 's/^0\./0.0/' | sed 's/00$//')
        
        echo "Using Genrich for CUT&RUN parameter sweep on group: {wildcards.group}, q-value: $qvalue"
        
        # Build Genrich command
        cmd="Genrich -t {input.treatment}"
        
        if [ ! -z "{input.control}" ]; then
            cmd="$cmd -c {input.control}"
        fi
        
        cmd="$cmd {params.atac_flag} -q $qvalue -o {output.peaks}"
        
        echo "Running: $cmd"
        eval $cmd
        """

# ChIP-seq Genrich parameter sweep
rule chipseq_genrich_sweep:
    input:
        treatment = get_genrich_treatment_bams_chipseq,
        control = get_genrich_control_bams_chipseq
    output:
        peaks = f"{config['outdir']}/chipseq/genrich_sweep/{{group}}_genrich_{{qval}}.narrowPeak"
    params:
        outdir = f"{config['outdir']}/chipseq/genrich_sweep"
    shell:
        """
        mkdir -p {params.outdir}
        
        # Convert qvalue format (q005 -> 0.05)
        qvalue=$(echo "{wildcards.qval}" | sed 's/q/0./' | sed 's/^0\./0.0/' | sed 's/00$//')
        
        echo "Using Genrich for ChIP-seq parameter sweep on group: {wildcards.group}, q-value: $qvalue"
        
        # Build Genrich command
        cmd="Genrich -t {input.treatment}"
        
        if [ ! -z "{input.control}" ]; then
            cmd="$cmd -c {input.control}"
        fi
        
        # ChIP-seq doesn't use -j flag
        cmd="$cmd -q $qvalue -o {output.peaks}"
        
        echo "Running: $cmd"
        eval $cmd
        """

# Summary rule for Genrich parameter sweep
rule genrich_parameter_sweep_summary:
    input:
        get_genrich_sweep_outputs
    output:
        f"{config['outdir']}/genrich_sweep/genrich_sweep_summary.txt"
    shell:
        """
        mkdir -p $(dirname {output})
        echo "Genrich Parameter Sweep Summary" > {output}
        echo "Generated on: $(date)" >> {output}
        echo "" >> {output}
        
        echo "Peak files generated:" >> {output}
        for file in {input}; do
            if [ -f "$file" ]; then
                peaks=$(wc -l < "$file")
                echo "  $file: $peaks peaks" >> {output}
            fi
        done
        """
