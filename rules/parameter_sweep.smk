# ============================================================
# Parameter Sweep Rules
# ============================================================

def get_parameter_sweep_outputs():
    """Get all parameter sweep outputs if enabled"""
    if not config.get('parameter_sweep', {}).get('enabled', False):
        return []
    
    qvalues = config.get('parameter_sweep', {}).get('qvalues', [0.05, 0.01, 0.005, 0.001])
    outputs = []
    
    # Convert q-values to strings for file naming (replace dots with underscores)
    qval_strings = [str(q).replace('.', '') for q in qvalues]
    
    # ATAC-seq parameter sweep outputs
    if ATACSEQ_SAMPLES:
        for sample in ATACSEQ_SAMPLES:
            for qval in qval_strings:
                outputs.append(f"{config['outdir']}/atacseq/peaks_sweep/{sample}_q{qval}_peaks.narrowPeak")
    
    # CUT&RUN parameter sweep outputs (excluding control samples)
    if CUTRUN_SAMPLES:
        cutrun_peak_samples = [s for s in CUTRUN_SAMPLES if not any(keyword in s.lower() for keyword in ['igg', 'input', 'control'])]
        for sample in cutrun_peak_samples:
            for qval in qval_strings:
                outputs.append(f"{config['outdir']}/cutrun/peaks_sweep/{sample}_q{qval}_peaks.narrowPeak")
    
    # ChIP-seq parameter sweep outputs
    if CHIPSEQ_SAMPLES:
        for sample in CHIPSEQ_SAMPLES:
            for qval in qval_strings:
                outputs.append(f"{config['outdir']}/chipseq/peaks_sweep/{sample}_q{qval}_peaks.narrowPeak")
    
    return outputs

def get_parameter_sweep_samples():
    """Get samples that can have peaks called (excludes controls)"""
    peak_samples = []
    
    # ATAC-seq samples (all can have peaks called)
    peak_samples.extend(ATACSEQ_SAMPLES)
    
    # CUT&RUN samples (exclude controls)
    if CUTRUN_SAMPLES:
        cutrun_peak_samples = [s for s in CUTRUN_SAMPLES 
                              if not any(keyword in s.lower() for keyword in ['igg', 'input', 'control'])]
        peak_samples.extend(cutrun_peak_samples)
    
    # ChIP-seq samples (all - controls will be handled by the rule logic)
    peak_samples.extend(CHIPSEQ_SAMPLES)
    
    return peak_samples

# Create a summary report of parameter sweep results
rule parameter_sweep_summary:
    input:
        get_parameter_sweep_outputs()
    output:
        f"{config['outdir']}/parameter_sweep/summary_report.txt"
    params:
        qvalues = config.get('parameter_sweep', {}).get('qvalues', [0.05, 0.01, 0.005, 0.001])
    run:
        import os
        import glob
        
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        
        with open(output[0], 'w') as f:
            f.write("Parameter Sweep Summary Report\n")
            f.write("="*50 + "\n\n")
            
            f.write(f"Q-values tested: {params.qvalues}\n")
            f.write(f"Total samples analyzed: {len(get_parameter_sweep_samples())}\n")
            f.write(f"Total peak files generated: {len(input)}\n\n")
            
            # Analyze by assay type
            for assay in ['atacseq', 'cutrun', 'chipseq']:
                assay_files = [f for f in input if f"/{assay}/peaks_sweep/" in f]
                if assay_files:
                    f.write(f"\n{assay.upper()} Results:\n")
                    f.write("-" * 20 + "\n")
                    
                    samples = set()
                    for file_path in assay_files:
                        sample_name = os.path.basename(file_path).split('_q')[0]
                        samples.add(sample_name)
                    
                    f.write(f"Samples: {len(samples)}\n")
                    f.write(f"Peak files: {len(assay_files)}\n")
                    
                    # Peak count summary for each q-value
                    for qval in params.qvalues:
                        qval_str = str(qval).replace('.', '')
                        qval_files = [f for f in assay_files if f"_q{qval_str}_" in f]
                        
                        total_peaks = 0
                        for peak_file in qval_files:
                            if os.path.exists(peak_file) and os.path.getsize(peak_file) > 0:
                                try:
                                    with open(peak_file, 'r') as pf:
                                        peak_count = sum(1 for line in pf if not line.startswith('#'))
                                        total_peaks += peak_count
                                except:
                                    pass
                        
                        f.write(f"  Q-value {qval}: {total_peaks} peaks across {len(qval_files)} samples\n")
            
            f.write(f"\n\nGenerated on: {__import__('datetime').datetime.now()}\n")

# Rule to run only parameter sweep analysis
rule parameter_sweep:
    input:
        get_parameter_sweep_outputs() + [f"{config['outdir']}/parameter_sweep/summary_report.txt"]
    message:
        "Parameter sweep analysis complete"
