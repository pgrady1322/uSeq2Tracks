# ============================================================
# UCSC track hub generation rules
# ============================================================

rule create_ucsc_hub:
    input:
        # Depend on pipeline completion files instead of individual files
        chipseq_complete=f"{GENOME_OUTDIR}/chipseq/tracks_complete.txt" if CHIPSEQ_SAMPLES else [],
        rnaseq_complete=f"{GENOME_OUTDIR}/rnaseq/tracks_complete.txt" if RNASEQ_SAMPLES else [],
        atacseq_complete=f"{GENOME_OUTDIR}/atacseq/tracks_complete.txt" if ATACSEQ_SAMPLES else [],
        cutrun_complete=f"{GENOME_OUTDIR}/cutrun/tracks_complete.txt" if CUTRUN_SAMPLES else [],
        wgs_complete=f"{GENOME_OUTDIR}/wgs/tracks_complete.txt" if WGS_SAMPLES else [],
        ancientdna_complete=f"{GENOME_OUTDIR}/ancientdna/tracks_complete.txt" if ADNA_SAMPLES else [],
        nanopore_complete=f"{GENOME_OUTDIR}/nanopore/tracks_complete.txt" if NANOPORE_SAMPLES else [],
        pacbio_complete=f"{GENOME_OUTDIR}/pacbio/tracks_complete.txt" if PACBIO_SAMPLES else []
    output:
        hub=f"{GENOME_OUTDIR}/ucsc/hub.txt",
        genomes=f"{GENOME_OUTDIR}/ucsc/genomes.txt",
        trackdb=f"{GENOME_OUTDIR}/ucsc/trackDb.txt"
    params:
        ucsc_config=config["ucsc"],
        outdir=f"{GENOME_OUTDIR}/ucsc"
    run:
        import os
        
        # Create output directory
        os.makedirs(params.outdir, exist_ok=True)
        
        # Note: This rule depends on completion files rather than individual files
        # to avoid circular dependencies. The trackDb.txt references files by path
        # and assumes they exist after pipeline completion.
        
        # Write hub.txt
        with open(output.hub, 'w') as f:
            f.write(f"hub {params.ucsc_config['hub_name']}\n")
            f.write(f"shortLabel {params.ucsc_config['hub_short_label']}\n")
            f.write(f"longLabel {params.ucsc_config['hub_long_label']}\n")
            f.write(f"genomesFile genomes.txt\n")
            f.write(f"email {params.ucsc_config['hub_email']}\n")
        
        # Write genomes.txt
        with open(output.genomes, 'w') as f:
            f.write(f"genome {params.ucsc_config['genome_name']}\n")
            f.write(f"trackDb trackDb.txt\n")
        
        # Write trackDb.txt
        with open(output.trackdb, 'w') as f:
            track_num = 1
            
            # ChIP-seq tracks - hierarchical structure with supertracks
            if CHIPSEQ_SAMPLES:
                # Main supertrack for all ChIP-seq data
                f.write("track chipseq_supertrack\n")
                f.write("superTrack on show\n")
                f.write("shortLabel ChIP-seq\n")
                f.write("longLabel ChIP-seq Signal and Peak Tracks\n")
                f.write("\n")
                
                # ChIP-seq Peak tracks composite
                f.write("track chipseq_peaks\n")
                f.write("parent chipseq_supertrack\n")
                f.write("compositeTrack on\n")
                f.write("shortLabel ChIP Peaks\n")
                f.write("longLabel ChIP-seq Peak Calls\n")
                f.write("type narrowPeak\n")
                f.write("visibility dense\n")
                f.write("itemRgb on\n")
                f.write("\n")
                
                # Individual peak tracks
                peak_colors = ["200,0,0", "0,150,0", "0,0,200", "150,100,0", "100,0,150", "200,100,0", "0,200,100"]
                for i, sample in enumerate(CHIPSEQ_SAMPLES):
                    sample_info = SAMPLES[sample]
                    condition = sample_info.get('condition', sample)
                    
                    # Skip control samples for peak tracks
                    if any(control_term in condition.lower() for control_term in ['input', 'control']):
                        continue
                        
                    f.write(f"    track chipseq_peaks_{sample}\n")
                    f.write(f"    parent chipseq_peaks\n")
                    f.write(f"    bigDataUrl ../chipseq/peaks/{sample}_peaks.narrowPeak\n")
                    f.write(f"    shortLabel {condition} peaks\n")
                    f.write(f"    longLabel ChIP-seq peaks for {condition}\n")
                    f.write(f"    type narrowPeak\n")
                    f.write(f"    color {peak_colors[i % len(peak_colors)]}\n")
                    f.write("\n")
                
                # ChIP-seq Signal tracks - organized by experiment groups if available
                if 'chipseq' in EXPERIMENT_GROUPS and EXPERIMENT_GROUPS['chipseq']:
                    # Hierarchical organization by experiment and replicate groups
                    for exp_name, exp_groups in EXPERIMENT_GROUPS['chipseq'].items():
                        # Experiment-level composite track
                        f.write(f"track chipseq_signal_{exp_name}\n")
                        f.write(f"parent chipseq_supertrack\n")
                        f.write("compositeTrack on\n")
                        f.write(f"shortLabel {exp_name} Signal\n")
                        f.write(f"longLabel ChIP-seq Signal Tracks - {exp_name}\n")
                        f.write("type bigWig\n")
                        f.write("visibility full\n")
                        f.write("autoScale on\n")
                        f.write("maxHeightPixels 100:50:20\n")
                        
                        # Create subgroups for this experiment
                        f.write("subGroup1 replicate_group Replicate_Group")
                        for rep_group in exp_groups.keys():
                            f.write(f" {rep_group}={rep_group}")
                        f.write("\n")
                        f.write("subGroup2 sample Sample")
                        all_samples = [sample for sample_list in exp_groups.values() for sample in sample_list]
                        for sample in all_samples:
                            f.write(f" {sample}={sample}")
                        f.write("\n\n")
                        
                        # Individual sample tracks within this experiment
                        signal_colors = ["150,0,0", "0,100,0", "0,0,150", "100,75,0", "75,0,100", "150,75,0", "0,150,75"]
                        color_idx = 0
                        
                        for rep_group, sample_list in exp_groups.items():
                            for sample in sample_list:
                                sample_info = SAMPLES[sample]
                                condition = sample_info.get('condition', sample)
                                
                                f.write(f"        track chipseq_signal_{sample}\n")
                                f.write(f"        parent chipseq_signal_{exp_name}\n")
                                f.write(f"        bigDataUrl ../chipseq/bigwig/{sample}.bw\n")
                                f.write(f"        shortLabel {condition}\n")
                                f.write(f"        longLabel ChIP-seq signal for {condition} ({rep_group})\n")
                                f.write(f"        type bigWig\n")
                                f.write(f"        color {signal_colors[color_idx % len(signal_colors)]}\n")
                                f.write(f"        subGroups replicate_group={rep_group} sample={sample}\n")
                                f.write("\n")
                                color_idx += 1
                else:
                    # Fallback: simple flat structure if no experiment groups
                    f.write("track chipseq_signal\n")
                    f.write("parent chipseq_supertrack\n")
                    f.write("compositeTrack on\n")
                    f.write("shortLabel ChIP Signal\n")
                    f.write("longLabel ChIP-seq Signal Tracks\n")
                    f.write("type bigWig\n")
                    f.write("visibility full\n")
                    f.write("\n")
                    
                    for sample in CHIPSEQ_SAMPLES:
                        sample_info = SAMPLES[sample]
                        condition = sample_info.get('condition', sample)
                        f.write(f"    track chipseq_signal_{sample}\n")
                        f.write(f"    parent chipseq_signal\n")
                        f.write(f"    bigDataUrl ../chipseq/bigwig/{sample}.bw\n")
                        f.write(f"    shortLabel {condition}\n")
                        f.write(f"    longLabel ChIP-seq signal for {condition}\n")
                        f.write(f"    type bigWig\n")
                        f.write(f"    color {track_num*50 % 255},{track_num*80 % 255},{track_num*120 % 255}\n")
                        f.write("\n")
                        track_num += 1
            
            # RNA-seq tracks - hierarchical structure with experiment groups if available  
            if RNASEQ_SAMPLES:
                if 'rnaseq' in EXPERIMENT_GROUPS and EXPERIMENT_GROUPS['rnaseq']:
                    # Main supertrack for all RNA-seq data
                    f.write("track rnaseq_supertrack\n")
                    f.write("superTrack on show\n")
                    f.write("shortLabel RNA-seq\n")
                    f.write("longLabel RNA-seq Expression Tracks\n")
                    f.write("\n")
                    
                    # Hierarchical organization by experiment and replicate groups
                    for exp_name, exp_groups in EXPERIMENT_GROUPS['rnaseq'].items():
                        # Experiment-level composite track
                        f.write(f"track rnaseq_{exp_name}\n")
                        f.write(f"parent rnaseq_supertrack\n")
                        f.write("compositeTrack on\n")
                        f.write(f"shortLabel {exp_name}\n")
                        f.write(f"longLabel RNA-seq Expression Tracks - {exp_name}\n")
                        f.write("type bigWig\n")
                        f.write("visibility full\n")
                        f.write("autoScale on\n")
                        f.write("maxHeightPixels 100:50:20\n")
                        
                        # Create subgroups for this experiment
                        f.write("subGroup1 replicate_group Replicate_Group")
                        for rep_group in exp_groups.keys():
                            f.write(f" {rep_group}={rep_group}")
                        f.write("\n")
                        f.write("subGroup2 sample Sample")
                        all_samples = [sample for sample_list in exp_groups.values() for sample in sample_list]
                        for sample in all_samples:
                            f.write(f" {sample}={sample}")
                        f.write("\n\n")
                        
                        # Individual sample tracks within this experiment
                        rnaseq_colors = ["0,100,200", "200,100,0", "100,200,0", "200,0,100", "0,200,100", "100,0,200", "200,150,0"]
                        color_idx = 0
                        
                        for rep_group, sample_list in exp_groups.items():
                            for sample in sample_list:
                                sample_info = SAMPLES[sample]
                                condition = sample_info.get('condition', sample)
                                
                                f.write(f"        track rnaseq_{sample}\n")
                                f.write(f"        parent rnaseq_{exp_name}\n")
                                f.write(f"        bigDataUrl ../rnaseq/bigwig/{sample}.bw\n")
                                f.write(f"        shortLabel {condition}\n")
                                f.write(f"        longLabel RNA-seq expression for {condition} ({rep_group})\n")
                                f.write(f"        type bigWig\n")
                                f.write(f"        color {rnaseq_colors[color_idx % len(rnaseq_colors)]}\n")
                                f.write(f"        subGroups replicate_group={rep_group} sample={sample}\n")
                                f.write("\n")
                                color_idx += 1
                else:
                    # Fallback: simple flat structure if no experiment groups
                    f.write("track rnaseq\n")
                    f.write("compositeTrack on\n")
                    f.write("shortLabel RNA-seq\n")
                    f.write("longLabel RNA-seq Expression Tracks\n")
                    f.write("type bigWig\n")
                    f.write("visibility full\n")
                    f.write("\n")
                    
                    for sample in RNASEQ_SAMPLES:
                        sample_info = SAMPLES[sample]
                        condition = sample_info.get('condition', sample)
                        f.write(f"    track rnaseq_{sample}\n")
                        f.write(f"    parent rnaseq\n")
                        f.write(f"    bigDataUrl ../rnaseq/bigwig/{sample}.bw\n")
                        f.write(f"    shortLabel {condition}\n")
                        f.write(f"    longLabel RNA-seq expression for {condition}\n")
                        f.write(f"    type bigWig\n")
                        f.write(f"    color {track_num*50 % 255},{track_num*80 % 255},{track_num*120 % 255}\n")
                        f.write("\n")
                        track_num += 1

            # ATAC-seq tracks - hierarchical structure with experiment groups if available
            if ATACSEQ_SAMPLES:
                if 'atacseq' in EXPERIMENT_GROUPS and EXPERIMENT_GROUPS['atacseq']:
                    # Main supertrack for all ATAC-seq data
                    f.write("track atacseq_supertrack\n")
                    f.write("superTrack on show\n")
                    f.write("shortLabel ATAC-seq\n")
                    f.write("longLabel ATAC-seq Accessibility Tracks\n")
                    f.write("\n")
                    
                    # Hierarchical organization by experiment and replicate groups
                    for exp_name, exp_groups in EXPERIMENT_GROUPS['atacseq'].items():
                        # Experiment-level composite track
                        f.write(f"track atacseq_{exp_name}\n")
                        f.write(f"parent atacseq_supertrack\n")
                        f.write("compositeTrack on\n")
                        f.write(f"shortLabel {exp_name}\n")
                        f.write(f"longLabel ATAC-seq Accessibility Tracks - {exp_name}\n")
                        f.write("type bigWig\n")
                        f.write("visibility full\n")
                        f.write("autoScale on\n")
                        f.write("maxHeightPixels 100:50:20\n")
                        
                        # Create subgroups for this experiment
                        f.write("subGroup1 replicate_group Replicate_Group")
                        for rep_group in exp_groups.keys():
                            f.write(f" {rep_group}={rep_group}")
                        f.write("\n")
                        f.write("subGroup2 sample Sample")
                        all_samples = [sample for sample_list in exp_groups.values() for sample in sample_list]
                        for sample in all_samples:
                            f.write(f" {sample}={sample}")
                        f.write("\n\n")
                        
                        # Individual sample tracks within this experiment
                        atacseq_colors = ["0,150,0", "150,0,0", "0,0,150", "150,100,0", "100,0,150", "0,150,100", "150,0,100"]
                        color_idx = 0
                        
                        for rep_group, sample_list in exp_groups.items():
                            for sample in sample_list:
                                sample_info = SAMPLES[sample]
                                condition = sample_info.get('condition', sample)
                                
                                f.write(f"        track atacseq_{sample}\n")
                                f.write(f"        parent atacseq_{exp_name}\n")
                                f.write(f"        bigDataUrl ../atacseq/bigwig/{sample}.bw\n")
                                f.write(f"        shortLabel {condition}\n")
                                f.write(f"        longLabel ATAC-seq accessibility for {condition} ({rep_group})\n")
                                f.write(f"        type bigWig\n")
                                f.write(f"        color {atacseq_colors[color_idx % len(atacseq_colors)]}\n")
                                f.write(f"        subGroups replicate_group={rep_group} sample={sample}\n")
                                f.write("\n")
                                color_idx += 1
                else:
                    # Fallback: simple flat structure if no experiment groups
                    f.write("track atacseq\n")
                    f.write("compositeTrack on\n")
                    f.write("shortLabel ATAC-seq\n")
                    f.write("longLabel ATAC-seq Accessibility Tracks\n")
                    f.write("type bigWig\n")
                    f.write("visibility full\n")
                    f.write("\n")
                    
                    for sample in ATACSEQ_SAMPLES:
                        sample_info = SAMPLES[sample]
                        condition = sample_info.get('condition', sample)
                        f.write(f"    track atacseq_{sample}\n")
                        f.write(f"    parent atacseq\n")
                        f.write(f"    bigDataUrl ../atacseq/bigwig/{sample}.bw\n")
                        f.write(f"    shortLabel {condition}\n")
                        f.write(f"    longLabel ATAC-seq accessibility for {condition}\n")
                        f.write(f"    type bigWig\n")
                        f.write(f"    color {track_num*50 % 255},{track_num*80 % 255},{track_num*120 % 255}\n")
                        f.write("\n")
                        track_num += 1

            # CUT&RUN tracks - hierarchical structure with experiment groups if available
            if CUTRUN_SAMPLES:
                if 'cutrun' in EXPERIMENT_GROUPS and EXPERIMENT_GROUPS['cutrun']:
                    # Main supertrack for all CUT&RUN data
                    f.write("track cutrun_supertrack\n")
                    f.write("superTrack on show\n")
                    f.write("shortLabel CUT&RUN\n")
                    f.write("longLabel CUT&RUN Signal Tracks\n")
                    f.write("\n")
                    
                    # Hierarchical organization by experiment and replicate groups
                    for exp_name, exp_groups in EXPERIMENT_GROUPS['cutrun'].items():
                        # Experiment-level composite track
                        f.write(f"track cutrun_{exp_name}\n")
                        f.write(f"parent cutrun_supertrack\n")
                        f.write("compositeTrack on\n")
                        f.write(f"shortLabel {exp_name}\n")
                        f.write(f"longLabel CUT&RUN Signal Tracks - {exp_name}\n")
                        f.write("type bigWig\n")
                        f.write("visibility full\n")
                        f.write("autoScale on\n")
                        f.write("maxHeightPixels 100:50:20\n")
                        
                        # Create subgroups for this experiment
                        f.write("subGroup1 replicate_group Replicate_Group")
                        for rep_group in exp_groups.keys():
                            f.write(f" {rep_group}={rep_group}")
                        f.write("\n")
                        f.write("subGroup2 sample Sample")
                        all_samples = [sample for sample_list in exp_groups.values() for sample in sample_list]
                        for sample in all_samples:
                            f.write(f" {sample}={sample}")
                        f.write("\n\n")
                        
                        # Individual sample tracks within this experiment
                        cutrun_colors = ["100,50,150", "150,100,0", "0,150,50", "150,0,100", "50,150,0", "0,100,150", "150,50,0"]
                        color_idx = 0
                        
                        for rep_group, sample_list in exp_groups.items():
                            for sample in sample_list:
                                sample_info = SAMPLES[sample]
                                condition = sample_info.get('condition', sample)
                                
                                f.write(f"        track cutrun_{sample}\n")
                                f.write(f"        parent cutrun_{exp_name}\n")
                                f.write(f"        bigDataUrl ../cutrun/bigwig/{sample}.bw\n")
                                f.write(f"        shortLabel {condition}\n")
                                f.write(f"        longLabel CUT&RUN signal for {condition} ({rep_group})\n")
                                f.write(f"        type bigWig\n")
                                f.write(f"        color {cutrun_colors[color_idx % len(cutrun_colors)]}\n")
                                f.write(f"        subGroups replicate_group={rep_group} sample={sample}\n")
                                f.write("\n")
                                color_idx += 1
                else:
                    # Fallback: simple flat structure if no experiment groups
                    f.write("track cutrun\n")
                    f.write("compositeTrack on\n")
                    f.write("shortLabel CUT&RUN\n")
                    f.write("longLabel CUT&RUN Signal Tracks\n")
                    f.write("type bigWig\n")
                    f.write("visibility full\n")
                    f.write("\n")
                    
                    for sample in CUTRUN_SAMPLES:
                        sample_info = SAMPLES[sample]
                        condition = sample_info.get('condition', sample)
                        f.write(f"    track cutrun_{sample}\n")
                        f.write(f"    parent cutrun\n")
                        f.write(f"    bigDataUrl ../cutrun/bigwig/{sample}.bw\n")
                        f.write(f"    shortLabel {condition}\n")
                        f.write(f"    longLabel CUT&RUN signal for {condition}\n")
                        f.write(f"    type bigWig\n")
                        f.write(f"    color {track_num*50 % 255},{track_num*80 % 255},{track_num*120 % 255}\n")
                        f.write("\n")
                        track_num += 1

            # WGS tracks
            if WGS_SAMPLES:
                f.write("track wgs\n")
                f.write("compositeTrack on\n")
                f.write("shortLabel WGS\n")
                f.write("longLabel Whole Genome Sequencing Tracks\n")
                f.write("type bigWig\n")
                f.write("visibility dense\n")
                f.write("\n")
                
                for sample in WGS_SAMPLES:
                    sample_info = SAMPLES[sample]
                    condition = sample_info.get('condition', sample)
                    f.write(f"    track wgs_{sample}\n")
                    f.write(f"    parent wgs\n")
                    f.write(f"    bigDataUrl ../wgs/bigwig/{sample}.bw\n")
                    f.write(f"    shortLabel {condition}\n")
                    f.write(f"    longLabel WGS coverage for {condition}\n")
                    f.write(f"    type bigWig\n")
                    f.write(f"    color {track_num*50 % 255},{track_num*80 % 255},{track_num*120 % 255}\n")
                    f.write("\n")
                    track_num += 1

            # Ancient DNA tracks
            if ADNA_SAMPLES:
                f.write("track ancientdna\n")
                f.write("compositeTrack on\n")
                f.write("shortLabel Ancient DNA\n")
                f.write("longLabel Ancient DNA Coverage Tracks\n")
                f.write("type bigWig\n")
                f.write("visibility dense\n")
                f.write("\n")
                
                for sample in ADNA_SAMPLES:
                    sample_info = SAMPLES[sample]
                    condition = sample_info.get('condition', sample)
                    f.write(f"    track ancientdna_{sample}\n")
                    f.write(f"    parent ancientdna\n")
                    f.write(f"    bigDataUrl ../ancientdna/bigwig/{sample}.bw\n")
                    f.write(f"    shortLabel {condition}\n")
                    f.write(f"    longLabel Ancient DNA coverage for {condition}\n")
                    f.write(f"    type bigWig\n")
                    f.write(f"    color {track_num*50 % 255},{track_num*80 % 255},{track_num*120 % 255}\n")
                    f.write("\n")
                    track_num += 1

            # Long-read tracks (Nanopore and PacBio)
            if NANOPORE_SAMPLES:
                f.write("track nanopore\n")
                f.write("compositeTrack on\n")
                f.write("shortLabel Nanopore\n")
                f.write("longLabel Nanopore Long-read Coverage Tracks\n")
                f.write("type bigWig\n")
                f.write("visibility dense\n")
                f.write("\n")
                
                for sample in NANOPORE_SAMPLES:
                    sample_info = SAMPLES[sample]
                    condition = sample_info.get('condition', sample)
                    f.write(f"    track nanopore_{sample}\n")
                    f.write(f"    parent nanopore\n")
                    f.write(f"    bigDataUrl ../nanopore/bigwig/{sample}.bw\n")
                    f.write(f"    shortLabel {condition}\n")
                    f.write(f"    longLabel Nanopore coverage for {condition}\n")
                    f.write(f"    type bigWig\n")
                    f.write(f"    color {track_num*50 % 255},{track_num*80 % 255},{track_num*120 % 255}\n")
                    f.write("\n")
                    track_num += 1

            if PACBIO_SAMPLES:
                f.write("track pacbio\n")
                f.write("compositeTrack on\n")
                f.write("shortLabel PacBio\n")
                f.write("longLabel PacBio Long-read Coverage Tracks\n")
                f.write("type bigWig\n")
                f.write("visibility dense\n")
                f.write("\n")
                
                for sample in PACBIO_SAMPLES:
                    sample_info = SAMPLES[sample]
                    condition = sample_info.get('condition', sample)
                    f.write(f"    track pacbio_{sample}\n")
                    f.write(f"    parent pacbio\n")
                    f.write(f"    bigDataUrl ../pacbio/bigwig/{sample}.bw\n")
                    f.write(f"    shortLabel {condition}\n")
                    f.write(f"    longLabel PacBio coverage for {condition}\n")
                    f.write(f"    type bigWig\n")
                    f.write(f"    color {track_num*50 % 255},{track_num*80 % 255},{track_num*120 % 255}\n")
                    f.write("\n")
                    track_num += 1
