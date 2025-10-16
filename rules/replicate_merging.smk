# ============================================================
# Replicate merging and composite track generation rules
# ============================================================

def get_merge_inputs(wildcards):
    """Get input files for replicate merging with better error handling"""
    assay = wildcards.assay
    experiment = wildcards.experiment
    group = wildcards.group
    
    # Debug information
    print(f"DEBUG: Looking for assay={assay}, experiment={experiment}, group={group}")
    print(f"DEBUG: Available assays in EXPERIMENT_GROUPS: {list(EXPERIMENT_GROUPS.keys())}")
    
    if assay not in EXPERIMENT_GROUPS:
        raise ValueError(f"Assay '{assay}' not found in EXPERIMENT_GROUPS. Available: {list(EXPERIMENT_GROUPS.keys())}")
    
    print(f"DEBUG: Available experiments for {assay}: {list(EXPERIMENT_GROUPS[assay].keys())}")
    
    if experiment not in EXPERIMENT_GROUPS[assay]:
        raise ValueError(f"Experiment '{experiment}' not found for assay '{assay}'. Available: {list(EXPERIMENT_GROUPS[assay].keys())}")
    
    print(f"DEBUG: Available groups for {assay}/{experiment}: {list(EXPERIMENT_GROUPS[assay][experiment].keys())}")
    
    if group not in EXPERIMENT_GROUPS[assay][experiment]:
        raise ValueError(f"Group '{group}' not found for assay '{assay}' experiment '{experiment}'. Available: {list(EXPERIMENT_GROUPS[assay][experiment].keys())}")
    
    samples = EXPERIMENT_GROUPS[assay][experiment][group]
    print(f"DEBUG: Found samples for {assay}/{experiment}/{group}: {samples}")
    
    return expand(f"{config['outdir']}/{assay}/bigwig/{{sample}}.bw", sample=samples)

# Debug rule to inspect experiment groups structure
rule debug_experiment_groups:
    output:
        f"{config['outdir']}/debug/experiment_groups.txt"
    params:
    run:
        import os
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        
        with open(output[0], 'w') as f:
            f.write("EXPERIMENT_GROUPS structure:\n")
            f.write("="*50 + "\n")
            
            for assay, exp_dict in EXPERIMENT_GROUPS.items():
                f.write(f"\nAssay: {assay}\n")
                f.write(f"  Experiments: {list(exp_dict.keys())}\n")
                
                for exp, group_dict in exp_dict.items():
                    f.write(f"    Experiment: '{exp}'\n")
                    f.write(f"      Groups: {list(group_dict.keys())}\n")
                    
                    for group, samples in group_dict.items():
                        f.write(f"        Group: '{group}' -> Samples: {samples}\n")
            
            f.write("\n" + "="*50 + "\n")
            f.write("Generated merged file paths:\n")
            
            if 'atacseq' in EXPERIMENT_GROUPS:
                for exp in EXPERIMENT_GROUPS['atacseq'].keys():
                    for group in EXPERIMENT_GROUPS['atacseq'][exp].keys():
                        merged_path = f"{config['outdir']}/atacseq/bigwig_merged/{exp}.{group}_merged.bw"
                        f.write(f"  {merged_path}\n")

# Merge BigWig files for replicate groups using deepTools (within experiment groups)
rule merge_replicate_bigwigs:
    input:
        lambda wildcards: get_merge_inputs(wildcards)
    output:
        f"{config['outdir']}/{{assay}}/bigwig_merged/{{experiment}}.{{group}}_merged.bw"
    params:
        input_list=lambda wildcards, input: ' '.join(input),
        genome_sizes=f"{GENOME_OUTDIR}/genome/genome.sizes"
    shell:
        """
        mkdir -p $(dirname {output})
        
        # Use deeptools multiBigwigSummary and bigwigCompare for merging
        # Build bash array safely
        IFS=' ' read -r -a input_files <<< "{params.input_list}"
        
        if [ ${{#input_files[@]}} -eq 1 ]; then
            # Only one file, just copy
            cp "${{input_files[0]}}" {output}
        elif [ ${{#input_files[@]}} -eq 2 ]; then
            # Two files, use bigwigCompare
            bigwigCompare -b1 "${{input_files[0]}}" -b2 "${{input_files[1]}}" \
                         --operation mean \
                         -o {output} \
                         --binSize 10
        else
            # Multiple files, use multiBigwigSummary + custom averaging
            temp_matrix=$(mktemp).npz
            temp_bed=$(mktemp).bed
            
            # Create summary matrix
            multiBigwigSummary bins -b {params.input_list} \
                                  -o $temp_matrix \
                                  --outRawCounts $temp_bed \
                                  --binSize 10
            
            # Average the signals and create new BigWig
            python3 -c "
import numpy as np
import pandas as pd
import subprocess
import tempfile
import os

# Read the raw counts
df = pd.read_csv('$temp_bed', sep='\t')
score_cols = [col for col in df.columns if col not in ['chr', 'start', 'end']]

# Calculate mean across replicates
df['mean_score'] = df[score_cols].mean(axis=1)

# Remove zero/NaN regions
df = df[df['mean_score'] > 0]
df = df.dropna()

# Create bedGraph
bedgraph = df[['chr', 'start', 'end', 'mean_score']]
temp_bg = tempfile.NamedTemporaryFile(mode='w', suffix='.bedGraph', delete=False)
bedgraph.to_csv(temp_bg.name, sep='\t', header=False, index=False)
temp_bg.close()

# Convert to BigWig
subprocess.run([
    'bedGraphToBigWig',
    temp_bg.name,
    '{params.genome_sizes}',
    '{output}'
], check=True)

# Cleanup
os.unlink(temp_bg.name)
"
            
            # Cleanup temporary files
            rm -f $temp_matrix $temp_bed
        fi
        """

# Create merged BAM files for replicates (useful for peak calling, within experiment groups)
rule merge_replicate_bams:
    input:
        lambda wildcards: expand(f"{config['outdir']}/{{assay}}/bam/{{sample}}.sorted.bam", 
                                sample=EXPERIMENT_GROUPS[wildcards.assay][wildcards.experiment][wildcards.group])
    output:
        bam = f"{config['outdir']}/{{assay}}/bam_merged/{{experiment}}_{{group}}_merged.bam",
        bai = f"{config['outdir']}/{{assay}}/bam_merged/{{experiment}}_{{group}}_merged.bam.bai"
    params:
    threads: 4
    shell:
        """
        mkdir -p $(dirname {output.bam})
        samtools merge -@ {threads} {output.bam} {input}
        samtools index -@ {threads} {output.bam}
        """

# Enhanced peak calling on merged replicates for ATAC-seq (experiment-aware)
rule atacseq_merged_peaks:
    input:
        treatment = f"{config['outdir']}/atacseq/bam_merged/{{experiment}}_treatment_merged.bam",
        control = f"{config['outdir']}/atacseq/bam_merged/{{experiment}}_input_merged.bam" 
    output:
        narrowpeak = f"{config['outdir']}/atacseq/peaks_merged/{{experiment}}_treatment_vs_input_peaks.narrowPeak",
        summits = f"{config['outdir']}/atacseq/peaks_merged/{{experiment}}_treatment_vs_input_summits.bed"
    params:
        name = "{experiment}_treatment_vs_input",
        outdir = f"{config['outdir']}/atacseq/peaks_merged",
        shift = config['atacseq']['shift'],
        extsize = config['atacseq']['extsize'],
        macs3_opts = config['atacseq']['macs3_opts']
    shell:
        """
        mkdir -p {params.outdir}
        macs3 callpeak -t {input.treatment} -c {input.control} \
                       -f BAMPE -g hs -n {params.name} \
                       --outdir {params.outdir} \
                       --shift {params.shift} --extsize {params.extsize} \
                       --keep-dup all {params.macs3_opts}
        """

# Generate composite track configuration for UCSC (experiment-aware)
rule create_composite_tracks:
    input:
        # Ensure all merged BigWigs exist for all experiments - but only if we have valid experiment groups
        atacseq_merged = [f"{config['outdir']}/atacseq/bigwig_merged/{exp}.{group}_merged.bw" 
                         for exp in EXPERIMENT_GROUPS.get('atacseq', {}).keys()
                         for group in EXPERIMENT_GROUPS.get('atacseq', {}).get(exp, {}).keys() 
                         if exp and group] if 'atacseq' in EXPERIMENT_GROUPS and EXPERIMENT_GROUPS['atacseq'] else []
    output:
        f"{config['outdir']}/ucsc/composite_trackDb.txt"
    params:
    run:
        with open(output[0], 'w') as f:
            # ATAC-seq composite track
            if 'atacseq' in EXPERIMENT_GROUPS and EXPERIMENT_GROUPS['atacseq']:
                f.write("track atacseq_composite\n")
                f.write("compositeTrack on\n")
                f.write("shortLabel ATAC-seq\n")
                f.write("longLabel ATAC-seq Accessibility Tracks\n")
                f.write("type bigWig\n")
                f.write("visibility full\n")
                f.write("autoScale on\n")
                f.write("maxHeightPixels 100:50:20\n")
                f.write("subGroup1 experiment Experiment")
                for exp in EXPERIMENT_GROUPS['atacseq'].keys():
                    f.write(f" {exp}={exp}")
                f.write("\n")
                f.write("subGroup2 condition Condition")
                all_conditions = set()
                for exp_data in EXPERIMENT_GROUPS['atacseq'].values():
                    all_conditions.update(exp_data.keys())
                for condition in sorted(all_conditions):
                    f.write(f" {condition}={condition}")
                f.write("\n")
                f.write("subGroup3 replicate Replicate individual=Individual merged=Merged")
                f.write("\n\n")
                
                # Individual replicate tracks
                color_idx = 0
                colors = ["200,0,0", "0,150,0", "0,0,200", "150,100,0", "100,0,150"]
                
                for exp_name, exp_groups in EXPERIMENT_GROUPS['atacseq'].items():
                    for group_name, sample_list in exp_groups.items():
                        # Individual replicates as subtracks
                        for sample in sample_list:
                            sample_info = SAMPLES[sample]
                            f.write(f"    track atacseq_{sample}\n")
                            f.write(f"    parent atacseq_composite\n")
                            f.write(f"    bigDataUrl ../atacseq/bigwig/{sample}.bw\n")
                            f.write(f"    shortLabel {sample}\n")
                            f.write(f"    longLabel ATAC-seq: {sample} ({exp_name} - {group_name})\n")
                            f.write(f"    type bigWig\n")
                            f.write(f"    color {colors[color_idx % len(colors)]}\n")
                            f.write(f"    subGroups experiment={exp_name} condition={group_name} replicate=individual\n")
                            f.write("\n")
                        
                        # Merged track for the group within this experiment
                        f.write(f"    track atacseq_{exp_name}_{group_name}_merged\n")
                        f.write(f"    parent atacseq_composite\n")
                        f.write(f"    bigDataUrl ../atacseq/bigwig_merged/{exp_name}.{group_name}_merged.bw\n")
                        f.write(f"    shortLabel {exp_name}_{group_name} (merged)\n")
                        f.write(f"    longLabel ATAC-seq: {exp_name} {group_name} (merged replicates)\n")
                        f.write(f"    type bigWig\n")
                        f.write(f"    color 0,0,0\n")
                        f.write(f"    subGroups experiment={exp_name} condition={group_name} replicate=merged\n")
                        f.write("\n")
                        
                        color_idx += 1
