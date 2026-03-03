#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
uSeq2Tracks – Samplesheet Validator

Validates sample sheet format for the uSeq2Tracks pipeline, checking
required columns, assay types, data sources (SRA / local), and
duplicate sample IDs.

Author: Patrick Grady
License: MIT License – See LICENSE
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

# ── Shared constants ──────────────────────────────────────────────
VALID_ASSAY_TYPES: list[str] = [
    "atacseq", "chipseq", "cutrun", "rnaseq",
    "wgs", "nanopore", "pacbio", "ancientdna",
]

INVALID_SENTINELS: set[str] = frozenset(
    {"", "nan", "NaN", "none", "None", "null", "NULL", "<NA>", "N/A", "n/a"}
)


class SamplesheetError(Exception):
    """Raised when samplesheet validation fails."""


# ── Helpers ───────────────────────────────────────────────────────
def is_valid_value(value: Any) -> bool:
    """Return *True* if *value* is a meaningful, non-sentinel string.

    Args:
        value: Any cell value read from the CSV.

    Returns:
        ``True`` when the value is neither ``None`` nor a common
        representation of missing data.
    """
    if value is None:
        return False
    return str(value).strip() not in INVALID_SENTINELS


# ── Core validation ───────────────────────────────────────────────
def check_samplesheet(samplesheet_path: Path, output_path: Path) -> int:
    """Validate samplesheet format, contents, and data sources.

    Args:
        samplesheet_path: Path to the input CSV samplesheet.
        output_path: Path where the validated samplesheet will be written.

    Returns:
        Number of validated samples.

    Raises:
        FileNotFoundError: If *samplesheet_path* does not exist.
        SamplesheetError: On any validation failure (missing columns,
            invalid types, duplicate IDs, missing data sources, etc.).
    """
    if not samplesheet_path.exists():
        raise FileNotFoundError(f"Samplesheet not found: {samplesheet_path}")

    required_cols = ["sample_id", "type"]
    seen_ids: set[str] = set()
    samples: list[dict[str, str]] = []

    with open(samplesheet_path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)

        if reader.fieldnames is None:
            raise SamplesheetError("Samplesheet is empty or has no header row")

        # Check required columns
        missing = [c for c in required_cols if c not in reader.fieldnames]
        if missing:
            raise SamplesheetError(
                f"Missing required columns: {', '.join(missing)}"
            )

        for lineno, row in enumerate(reader, start=2):
            # -- sample_id --
            if not is_valid_value(row.get("sample_id")):
                raise SamplesheetError(
                    f"Line {lineno}: sample_id is empty or invalid"
                )
            sample_id = row["sample_id"].strip()

            if sample_id in seen_ids:
                raise SamplesheetError(
                    f"Line {lineno}: Duplicate sample_id '{sample_id}'"
                )
            seen_ids.add(sample_id)

            # -- assay type --
            assay_type = row.get("type", "").strip().lower()
            if assay_type not in VALID_ASSAY_TYPES:
                raise SamplesheetError(
                    f"Line {lineno}: Invalid type '{assay_type}'. "
                    f"Must be one of: {', '.join(VALID_ASSAY_TYPES)}"
                )

            # -- data source --
            has_sra = is_valid_value(row.get("sra_id"))
            has_local = is_valid_value(row.get("read1"))

            if not has_sra and not has_local:
                raise SamplesheetError(
                    f"Line {lineno}: Sample '{sample_id}' has neither "
                    "SRA ID nor local files (read1)"
                )

            # -- optional file-existence warnings --
            if has_local:
                read1 = Path(row["read1"])
                if not read1.exists():
                    logger.warning("Line %d: read1 does not exist: %s", lineno, read1)
                if is_valid_value(row.get("read2")):
                    read2 = Path(row["read2"])
                    if not read2.exists():
                        logger.warning("Line %d: read2 does not exist: %s", lineno, read2)

            samples.append(row)

    if not samples:
        raise SamplesheetError("No valid samples found in samplesheet")

    # Write validated output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=reader.fieldnames)
        writer.writeheader()
        writer.writerows(samples)

    assay_types = sorted({s["type"] for s in samples})
    logger.info("Validated %d samples (%s)", len(samples), ", ".join(assay_types))
    return len(samples)


# ── CLI ───────────────────────────────────────────────────────────
def parse_args(args: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments.

    Args:
        args: Explicit argument list (defaults to ``sys.argv[1:]``).

    Returns:
        Parsed :class:`argparse.Namespace`.
    """
    parser = argparse.ArgumentParser(
        description="Validate uSeq2Tracks samplesheet and check for required fields",
    )
    parser.add_argument("samplesheet", type=Path, help="Input samplesheet CSV")
    parser.add_argument("output", type=Path, help="Output validated samplesheet")
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable debug logging"
    )
    return parser.parse_args(args)


def main(args: list[str] | None = None) -> None:
    """Entry point for the samplesheet validator CLI.

    Args:
        args: Explicit argument list (defaults to ``sys.argv[1:]``).
    """
    parsed = parse_args(args)
    logging.basicConfig(
        level=logging.DEBUG if parsed.verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )
    try:
        n = check_samplesheet(parsed.samplesheet, parsed.output)
        print(f"✓ Validated {n} samples")
    except (SamplesheetError, FileNotFoundError) as exc:
        logger.error(str(exc))
        sys.exit(1)


if __name__ == "__main__":
    main()
