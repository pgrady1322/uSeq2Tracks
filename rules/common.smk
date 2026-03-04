# ============================================================
# Common utility functions and rules
# ============================================================

import gzip
import os

import pandas as pd
from pathlib import Path

# ── Shared sentinel values ────────────────────────────────────────
# Canonical set of strings that should be treated as "missing / NA".
# Used by is_valid_sra_id(), has_local_files(), and get_trimmed_reads().
INVALID_SENTINELS: set[str] = frozenset(
    {"", "nan", "NaN", "none", "None", "null", "NULL", "<NA>", "N/A", "n/a"}
)


def is_valid_sra_id(sra_id: object) -> bool:
    """Check whether *sra_id* is a meaningful SRA accession.

    Handles ``None``, pandas ``NaN``, empty strings, whitespace, and
    common textual representations of missing data.

    Args:
        sra_id: Value from the samplesheet's ``sra_id`` column.

    Returns:
        ``True`` if *sra_id* looks like a real accession (e.g. SRR…).
    """
    if sra_id is None:
        return False
    if pd.isna(sra_id):
        return False
    return str(sra_id).strip() not in INVALID_SENTINELS


def has_local_files(sample_info: dict) -> bool:
    """Check whether *sample_info* contains a valid local ``read1`` path.

    Args:
        sample_info: Row dict from the parsed samplesheet.

    Returns:
        ``True`` when a non-sentinel ``read1`` value is present.
    """
    read1 = sample_info.get("read1", "")
    if not read1 or pd.isna(read1):
        return False
    return str(read1).strip() not in INVALID_SENTINELS


def get_sample_reads(wildcards) -> dict[str, str]:
    """Return input read paths for *wildcards.sample*.

    Priority: local files first, then SRA downloads.  For SRA data the
    function auto-detects single-end vs paired-end by inspecting whether
    the R2 file is a placeholder.

    Args:
        wildcards: Snakemake wildcards (must include ``sample``).

    Returns:
        Dict with ``"r1"`` (and optionally ``"r2"``) mapped to file paths.

    Raises:
        KeyError: If ``wildcards.sample`` is not found in ``SAMPLES``.
        ValueError: If the sample lacks both local files and an SRA ID.
    """
    if wildcards.sample not in SAMPLES:
        raise KeyError(
            f"Sample '{wildcards.sample}' not found in the samplesheet. "
            f"Known samples: {list(SAMPLES.keys())[:10]}"
        )
    sample_info = SAMPLES[wildcards.sample]

    # ── Local files (highest priority) ──
    if has_local_files(sample_info):
        reads: dict[str, str] = {"r1": sample_info["read1"]}
        read2 = sample_info.get("read2", "")
        if read2 and not pd.isna(read2) and str(read2).strip() not in INVALID_SENTINELS:
            reads["r2"] = read2
        return reads

    # ── SRA data ──
    if is_valid_sra_id(sample_info.get("sra_id", "")):
        r2_file = f"data/raw/{wildcards.sample}_R2.fastq.gz"

        if os.path.exists(r2_file):
            try:
                with gzip.open(r2_file, "rt") as fh:
                    first_line = fh.readline().strip()
                    if first_line == "# SINGLE_END_PLACEHOLDER":
                        return {"r1": f"data/raw/{wildcards.sample}_R1.fastq.gz"}
            except (OSError, gzip.BadGzipFile):
                # If R2 is tiny it's likely a placeholder
                if os.path.getsize(r2_file) < 100:
                    return {"r1": f"data/raw/{wildcards.sample}_R1.fastq.gz"}

        # Default: assume paired-end
        return {
            "r1": f"data/raw/{wildcards.sample}_R1.fastq.gz",
            "r2": r2_file,
        }

    raise ValueError(
        f"Sample '{wildcards.sample}' has neither valid local files nor a valid SRA ID. "
        f"read1={sample_info.get('read1', 'missing')}, "
        f"sra_id={sample_info.get('sra_id', 'missing')}"
    )


def get_genome_indices() -> dict[str, str]:
    """Return a dict of genome index files needed for the current sample set.

    Keys are descriptive labels (``fai``, ``bwa_mem2``, ``star``, …);
    values are paths under ``GENOME_OUTDIR``.

    Returns:
        Mapping of index label → expected file path.
    """
    genome_base = Path(config["genome"]).stem
    indices: dict[str, str] = {
        "fai": f"{GENOME_OUTDIR}/genome/{genome_base}.fa.fai",
        "dict": f"{GENOME_OUTDIR}/genome/{genome_base}.dict",
    }

    if WGS_SAMPLES or CHIPSEQ_SAMPLES or CUTRUN_SAMPLES or ATACSEQ_SAMPLES:
        indices["bwa_mem2"] = f"{GENOME_OUTDIR}/genome/bwa_mem2/{genome_base}.fa.0123"
        indices["bowtie2"] = f"{GENOME_OUTDIR}/genome/bowtie2/{genome_base}.1.bt2"

    if RNASEQ_SAMPLES:
        indices["star"] = f"{GENOME_OUTDIR}/genome/star/SA"
        indices["hisat2"] = f"{GENOME_OUTDIR}/genome/hisat2/{genome_base}.1.ht2"

    if NANOPORE_SAMPLES or PACBIO_SAMPLES:
        indices["minimap2"] = f"{GENOME_OUTDIR}/genome/minimap2/{genome_base}.mmi"

    if ADNA_SAMPLES:
        indices["bwa_aln"] = f"{GENOME_OUTDIR}/genome/bwa_aln/{genome_base}.fa.bwt"

    return indices


def get_samples_by_type(assay_type: str) -> list[str]:
    """Return sample IDs for a given assay type.

    Args:
        assay_type: One of the recognised assay strings (e.g. ``"atacseq"``).

    Returns:
        List of matching sample IDs (may be empty).
    """
    return [s for s, info in SAMPLES.items() if info.get("type") == assay_type]


def get_control_for_sample(wildcards) -> str | list:
    """Find the control BAM for a ChIP-seq / CUT&RUN sample.

    The control is identified by matching the ``condition`` field with the
    ``control_tag`` suffix from the config.

    Args:
        wildcards: Snakemake wildcards (must include ``sample``).

    Returns:
        Path to the control BAM file, or ``[]`` if no control sample
        was found (MACS3 will run without a control).
    """
    if wildcards.sample not in SAMPLES:
        raise KeyError(f"Sample '{wildcards.sample}' not in samplesheet")
    sample_info = SAMPLES[wildcards.sample]
    condition = sample_info.get("condition", "")
    control_tag = config.get("chipseq", {}).get("control_tag", "input")
    control_condition = f"{condition}_{control_tag}"

    for sample_id, info in SAMPLES.items():
        if info.get("condition") == control_condition:
            suffix = "dedup.bam" if config.get("chipseq", {}).get("markdup", False) else "sorted.bam"
            return f"{GENOME_OUTDIR}/chipseq/bam/{sample_id}.{suffix}"

    # No control found — warn but don't fail (MACS3 can run without)
    print(f"  ⚠ No control sample found for '{wildcards.sample}' "
          f"(expected condition '{control_condition}')")
    return []


def get_fastq_for_fastqc(wildcards) -> list[str]:
    """Return FASTQ file paths suitable for FastQC.

    Args:
        wildcards: Snakemake wildcards (must include ``sample``).

    Returns:
        List containing one (single-end) or two (paired-end) FASTQ paths.
    """
    reads = get_sample_reads(wildcards)
    return [reads["r1"]] + ([reads["r2"]] if "r2" in reads else [])


def get_trimmed_reads(wildcards) -> dict[str, str]:
    """Return expected trimmed-read paths for *wildcards.sample*.

    For local paired-end data, both R1 and R2 paths are returned.  For
    SRA data the R2 path is always included (actual layout is
    auto-detected at download time).

    Args:
        wildcards: Snakemake wildcards (must include ``sample``).

    Returns:
        Dict with ``"r1"`` (and optionally ``"r2"``) trimmed FASTQ paths.
    """
    sample_info = SAMPLES[wildcards.sample]
    reads: dict[str, str] = {
        "r1": f"{GENOME_OUTDIR}/trimmed/{wildcards.sample}_R1_trimmed.fastq.gz"
    }

    if has_local_files(sample_info):
        read2 = sample_info.get("read2", "")
        if read2 and not pd.isna(read2) and str(read2).strip() not in INVALID_SENTINELS:
            reads["r2"] = f"{GENOME_OUTDIR}/trimmed/{wildcards.sample}_R2_trimmed.fastq.gz"
    elif is_valid_sra_id(sample_info.get("sra_id", "")):
        reads["r2"] = f"{GENOME_OUTDIR}/trimmed/{wildcards.sample}_R2_trimmed.fastq.gz"

    return reads


def get_input_reads(wildcards) -> dict[str, str]:
    """Return trimmed reads if trimming is enabled, otherwise raw reads.

    Args:
        wildcards: Snakemake wildcards (must include ``sample``).

    Returns:
        Dict with ``"r1"`` (and optionally ``"r2"``) file paths.
    """
    if config.get("adapter_trimming", {}).get("enabled", False):
        return get_trimmed_reads(wildcards)
    return get_sample_reads(wildcards)
