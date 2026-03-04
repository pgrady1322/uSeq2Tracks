#!/usr/bin/env python3
"""
uSeq2Tracks – UCSC Genome Browser Track Hub Generator

Generates the three files required for a UCSC track hub:
  hub.txt, genomes.txt, trackDb.txt

Author: Patrick Grady
License: MIT License – See LICENSE
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

# ── Constants ─────────────────────────────────────────────────────
TRACK_COLORS: list[str] = [
    "200,0,0",  # red
    "0,150,0",  # green
    "0,0,200",  # blue
    "150,100,0",  # olive
    "100,0,150",  # purple
    "200,100,0",  # orange
    "0,200,100",  # teal
]


# ── CLI ───────────────────────────────────────────────────────────
def parse_args(args: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments for hub generation.

    Args:
        args: Explicit argument list (defaults to ``sys.argv[1:]``).

    Returns:
        Parsed :class:`argparse.Namespace`.
    """
    parser = argparse.ArgumentParser(
        description="Generate UCSC track hub files",
    )
    parser.add_argument("--genome-id", required=True, help="Genome identifier")
    parser.add_argument("--hub-name", required=True, help="Hub name (no spaces)")
    parser.add_argument("--hub-short-label", required=True, help="Hub short label")
    parser.add_argument("--hub-long-label", required=True, help="Hub long label")
    parser.add_argument("--genome-name", required=True, help="Genome name")
    parser.add_argument("--hub-email", required=True, help="Contact email for the hub")
    parser.add_argument("--bigwigs", nargs="+", help="BigWig files to include")
    parser.add_argument("--peaks", nargs="*", help="Peak files to include (optional)")
    parser.add_argument("--output-dir", default=".", help="Output directory (default: .)")
    return parser.parse_args(args)


# ── Writers ───────────────────────────────────────────────────────
def write_hub_txt(args: argparse.Namespace, output_dir: Path) -> None:
    """Write the ``hub.txt`` descriptor file.

    Args:
        args: Parsed CLI arguments containing hub metadata.
        output_dir: Directory in which to create the file.
    """
    hub_file = output_dir / "hub.txt"
    with open(hub_file, "w", encoding="utf-8") as fh:
        fh.write(f"hub {args.hub_name}\n")
        fh.write(f"shortLabel {args.hub_short_label}\n")
        fh.write(f"longLabel {args.hub_long_label}\n")
        fh.write("genomesFile genomes.txt\n")
        fh.write(f"email {args.hub_email}\n")
    logger.info("Created %s", hub_file)


def write_genomes_txt(args: argparse.Namespace, output_dir: Path) -> None:
    """Write the ``genomes.txt`` file linking genome to trackDb.

    Args:
        args: Parsed CLI arguments containing genome metadata.
        output_dir: Directory in which to create the file.
    """
    genomes_file = output_dir / "genomes.txt"
    with open(genomes_file, "w", encoding="utf-8") as fh:
        fh.write(f"genome {args.genome_name}\n")
        fh.write("trackDb trackDb.txt\n")
    logger.info("Created %s", genomes_file)


def write_trackdb_txt(args: argparse.Namespace, output_dir: Path) -> None:
    """Write the ``trackDb.txt`` with signal and peak tracks.

    Organises BigWig and peak files into a composite track.  Peak files
    are matched to BigWig files by stripping the ``_peaks`` suffix from
    the peak file stem.

    Args:
        args: Parsed CLI arguments containing bigwig/peak file lists.
        output_dir: Directory in which to create the file.
    """
    bigwig_paths = [Path(bw) for bw in args.bigwigs] if args.bigwigs else []
    peak_paths = [Path(p) for p in args.peaks] if args.peaks else []

    if not bigwig_paths:
        logger.warning("No bigwig files provided — trackDb.txt will be empty")

    # Build sample dict: {sample_id: {bigwig: Path, peaks?: Path}}
    samples: dict[str, dict[str, Path]] = {}
    for bw in bigwig_paths:
        samples[bw.stem] = {"bigwig": bw}

    for peak in peak_paths:
        sample_id = peak.stem.replace("_peaks", "")
        if sample_id in samples:
            samples[sample_id]["peaks"] = peak
        else:
            logger.warning("Peak file '%s' has no matching bigwig — skipping", peak.name)

    trackdb_file = output_dir / "trackDb.txt"
    with open(trackdb_file, "w", encoding="utf-8") as fh:
        # Composite track header
        fh.write("track all_samples\n")
        fh.write("compositeTrack on\n")
        fh.write("shortLabel All Tracks\n")
        fh.write("longLabel All uSeq2Tracks Signal and Peak Tracks\n")
        fh.write("type bigWig\n")
        fh.write("visibility full\n\n")

        for idx, (sample_id, files) in enumerate(samples.items()):
            color = TRACK_COLORS[idx % len(TRACK_COLORS)]

            if "bigwig" in files:
                rel_path = f"../bigwig/{files['bigwig'].name}"
                fh.write(f"    track {sample_id}_signal\n")
                fh.write("    parent all_samples\n")
                fh.write(f"    bigDataUrl {rel_path}\n")
                fh.write(f"    shortLabel {sample_id}\n")
                fh.write(f"    longLabel {sample_id} Signal Track\n")
                fh.write("    type bigWig\n")
                fh.write(f"    color {color}\n")
                fh.write("    autoScale on\n\n")

            if "peaks" in files:
                rel_path = f"../peaks/{files['peaks'].name}"
                fh.write(f"    track {sample_id}_peaks\n")
                fh.write("    parent all_samples\n")
                fh.write(f"    bigDataUrl {rel_path}\n")
                fh.write(f"    shortLabel {sample_id} Peaks\n")
                fh.write(f"    longLabel {sample_id} Peak Calls\n")
                fh.write("    type narrowPeak\n")
                fh.write(f"    color {color}\n\n")

    logger.info("Created %s with %d samples", trackdb_file, len(samples))


# ── Main ──────────────────────────────────────────────────────────
def main(args: list[str] | None = None) -> None:
    """Entry point for the UCSC hub generator CLI.

    Args:
        args: Explicit argument list (defaults to ``sys.argv[1:]``).
    """
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    parsed = parse_args(args)
    output_dir = Path(parsed.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    write_hub_txt(parsed, output_dir)
    write_genomes_txt(parsed, output_dir)
    write_trackdb_txt(parsed, output_dir)

    print(f"\n✓ UCSC track hub generated in {output_dir}")


if __name__ == "__main__":
    main()
