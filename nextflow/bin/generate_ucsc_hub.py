#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
uSeq2Tracks v1.0.0

UCSC Genome Browser track hub generator (Nextflow).

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: MIT License - See LICENSE
"""

import sys
import argparse
from pathlib import Path
import json


def parse_args(args=None):
    parser = argparse.ArgumentParser(
        description="Generate UCSC track hub files"
    )
    parser.add_argument("--genome-id", required=True, help="Genome identifier")
    parser.add_argument("--hub-name", required=True, help="Hub name")
    parser.add_argument("--hub-short-label", required=True, help="Hub short label")
    parser.add_argument("--hub-long-label", required=True, help="Hub long label")
    parser.add_argument("--genome-name", required=True, help="Genome name")
    parser.add_argument("--hub-email", required=True, help="Hub email")
    parser.add_argument("--bigwigs", nargs='+', help="BigWig files")
    parser.add_argument("--peaks", nargs='*', help="Peak files")
    parser.add_argument("--output-dir", default=".", help="Output directory")
    return parser.parse_args(args)


def write_hub_txt(args, output_dir):
    """Write hub.txt file"""
    hub_file = output_dir / "hub.txt"
    with open(hub_file, 'w') as f:
        f.write(f"hub {args.hub_name}\\n")
        f.write(f"shortLabel {args.hub_short_label}\\n")
        f.write(f"longLabel {args.hub_long_label}\\n")
        f.write(f"genomesFile genomes.txt\\n")
        f.write(f"email {args.hub_email}\\n")
    print(f"✓ Created {hub_file}")


def write_genomes_txt(args, output_dir):
    """Write genomes.txt file"""
    genomes_file = output_dir / "genomes.txt"
    with open(genomes_file, 'w') as f:
        f.write(f"genome {args.genome_name}\\n")
        f.write(f"trackDb trackDb.txt\\n")
    print(f"✓ Created {genomes_file}")


def write_trackdb_txt(args, output_dir):
    """Write trackDb.txt file"""
    trackdb_file = output_dir / "trackDb.txt"
    
    # Organize tracks by assay type
    bigwig_paths = [Path(bw) for bw in args.bigwigs] if args.bigwigs else []
    peak_paths = [Path(p) for p in args.peaks] if args.peaks else []
    
    # Extract sample info from filenames
    samples = {}
    for bw in bigwig_paths:
        # Assuming format: sample_id.bw
        sample_id = bw.stem
        samples[sample_id] = {'bigwig': bw}
    
    for peak in peak_paths:
        # Assuming format: sample_id_peaks.narrowPeak
        sample_id = peak.stem.replace('_peaks', '')
        if sample_id in samples:
            samples[sample_id]['peaks'] = peak
    
    with open(trackdb_file, 'w') as f:
        # Create a composite track for all samples
        f.write("track all_samples\\n")
        f.write("compositeTrack on\\n")
        f.write("shortLabel All Tracks\\n")
        f.write("longLabel All uSeq2Tracks Signal and Peak Tracks\\n")
        f.write("type bigWig\\n")
        f.write("visibility full\\n")
        f.write("\\n")
        
        # Individual BigWig tracks
        colors = ["200,0,0", "0,150,0", "0,0,200", "150,100,0", "100,0,150", "200,100,0", "0,200,100"]
        for i, (sample_id, files) in enumerate(samples.items()):
            color = colors[i % len(colors)]
            
            # BigWig track
            if 'bigwig' in files:
                rel_path = f"../bigwig/{files['bigwig'].name}"
                f.write(f"    track {sample_id}_signal\\n")
                f.write(f"    parent all_samples\\n")
                f.write(f"    bigDataUrl {rel_path}\\n")
                f.write(f"    shortLabel {sample_id}\\n")
                f.write(f"    longLabel {sample_id} Signal Track\\n")
                f.write(f"    type bigWig\\n")
                f.write(f"    color {color}\\n")
                f.write(f"    autoScale on\\n")
                f.write("\\n")
            
            # Peak track (if available)
            if 'peaks' in files:
                rel_path = f"../peaks/{files['peaks'].name}"
                f.write(f"    track {sample_id}_peaks\\n")
                f.write(f"    parent all_samples\\n")
                f.write(f"    bigDataUrl {rel_path}\\n")
                f.write(f"    shortLabel {sample_id} Peaks\\n")
                f.write(f"    longLabel {sample_id} Peak Calls\\n")
                f.write(f"    type narrowPeak\\n")
                f.write(f"    color {color}\\n")
                f.write("\\n")
    
    print(f"✓ Created {trackdb_file} with {len(samples)} samples")


def main(args=None):
    args = parse_args(args)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    write_hub_txt(args, output_dir)
    write_genomes_txt(args, output_dir)
    write_trackdb_txt(args, output_dir)
    
    print("\\n✓ UCSC track hub generation complete!")


if __name__ == '__main__':
    main()

# uSeq2Tracks v1.0.0
# Any usage is subject to this software's license.
