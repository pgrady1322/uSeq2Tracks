"""Tests for bin/check_samplesheet.py"""

from __future__ import annotations

import csv

# The script lives in bin/ — import by manipulating sys.path lazily.
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "bin"))

from check_samplesheet import (
    INVALID_SENTINELS,
    VALID_ASSAY_TYPES,
    SamplesheetError,
    check_samplesheet,
    is_valid_value,
    main,
)


# ── Fixtures ──────────────────────────────────────────────────────
@pytest.fixture()
def tmp_dir(tmp_path: Path) -> Path:
    """Return a clean temp directory."""
    return tmp_path


def _write_csv(path: Path, rows: list[dict[str, str]]) -> Path:
    """Helper to write a CSV samplesheet."""
    fieldnames = list(rows[0].keys()) if rows else ["sample_id", "type"]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return path


# ── is_valid_value ────────────────────────────────────────────────
class TestIsValidValue:
    def test_none_is_invalid(self):
        assert is_valid_value(None) is False

    @pytest.mark.parametrize("sentinel", sorted(INVALID_SENTINELS))
    def test_sentinels_are_invalid(self, sentinel: str):
        assert is_valid_value(sentinel) is False

    def test_real_value_is_valid(self):
        assert is_valid_value("SRR1234567") is True

    def test_whitespace_stripped(self):
        assert is_valid_value("  ") is False

    def test_numeric_string_is_valid(self):
        assert is_valid_value("42") is True


# ── check_samplesheet — happy path ───────────────────────────────
class TestCheckSamplesheetHappy:
    def test_valid_sra_samples(self, tmp_dir: Path):
        ss = _write_csv(
            tmp_dir / "input.csv",
            [
                {"sample_id": "s1", "type": "atacseq", "sra_id": "SRR111", "read1": ""},
                {"sample_id": "s2", "type": "chipseq", "sra_id": "SRR222", "read1": ""},
            ],
        )
        out = tmp_dir / "validated.csv"
        n = check_samplesheet(ss, out)
        assert n == 2
        assert out.exists()

    def test_valid_local_samples(self, tmp_dir: Path):
        # Create dummy FASTQ files so warnings are not logged
        r1 = tmp_dir / "sample_R1.fq.gz"
        r1.write_text("@read")
        ss = _write_csv(
            tmp_dir / "input.csv",
            [{"sample_id": "s1", "type": "wgs", "sra_id": "", "read1": str(r1)}],
        )
        out = tmp_dir / "validated.csv"
        n = check_samplesheet(ss, out)
        assert n == 1

    def test_output_preserves_columns(self, tmp_dir: Path):
        ss = _write_csv(
            tmp_dir / "input.csv",
            [
                {
                    "sample_id": "s1",
                    "type": "rnaseq",
                    "sra_id": "SRR999",
                    "read1": "",
                    "extra": "val",
                }
            ],
        )
        out = tmp_dir / "validated.csv"
        check_samplesheet(ss, out)
        with open(out, encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            row = next(reader)
            assert row["extra"] == "val"

    def test_output_dir_created_automatically(self, tmp_dir: Path):
        ss = _write_csv(
            tmp_dir / "input.csv",
            [{"sample_id": "s1", "type": "pacbio", "sra_id": "SRR555", "read1": ""}],
        )
        out = tmp_dir / "sub" / "dir" / "validated.csv"
        check_samplesheet(ss, out)
        assert out.exists()

    @pytest.mark.parametrize("assay", VALID_ASSAY_TYPES)
    def test_all_assay_types_accepted(self, tmp_dir: Path, assay: str):
        ss = _write_csv(
            tmp_dir / "input.csv",
            [{"sample_id": "s1", "type": assay, "sra_id": "SRR000", "read1": ""}],
        )
        out = tmp_dir / "validated.csv"
        assert check_samplesheet(ss, out) == 1


# ── check_samplesheet — error paths ──────────────────────────────
class TestCheckSamplesheetErrors:
    def test_file_not_found(self, tmp_dir: Path):
        with pytest.raises(FileNotFoundError, match="not found"):
            check_samplesheet(tmp_dir / "nope.csv", tmp_dir / "out.csv")

    def test_empty_file(self, tmp_dir: Path):
        ss = tmp_dir / "empty.csv"
        ss.write_text("")
        with pytest.raises(SamplesheetError, match="empty"):
            check_samplesheet(ss, tmp_dir / "out.csv")

    def test_missing_required_column(self, tmp_dir: Path):
        ss = tmp_dir / "bad.csv"
        ss.write_text("sample_id,foo\ns1,bar\n")
        with pytest.raises(SamplesheetError, match="Missing required columns.*type"):
            check_samplesheet(ss, tmp_dir / "out.csv")

    def test_empty_sample_id(self, tmp_dir: Path):
        ss = _write_csv(
            tmp_dir / "input.csv",
            [{"sample_id": "", "type": "wgs", "sra_id": "SRR111", "read1": ""}],
        )
        with pytest.raises(SamplesheetError, match="sample_id is empty"):
            check_samplesheet(ss, tmp_dir / "out.csv")

    def test_duplicate_sample_id(self, tmp_dir: Path):
        ss = _write_csv(
            tmp_dir / "input.csv",
            [
                {"sample_id": "dup", "type": "wgs", "sra_id": "SRR111", "read1": ""},
                {"sample_id": "dup", "type": "wgs", "sra_id": "SRR222", "read1": ""},
            ],
        )
        with pytest.raises(SamplesheetError, match="Duplicate sample_id 'dup'"):
            check_samplesheet(ss, tmp_dir / "out.csv")

    def test_invalid_assay_type(self, tmp_dir: Path):
        ss = _write_csv(
            tmp_dir / "input.csv",
            [{"sample_id": "s1", "type": "metagenomics", "sra_id": "SRR111", "read1": ""}],
        )
        with pytest.raises(SamplesheetError, match="Invalid type 'metagenomics'"):
            check_samplesheet(ss, tmp_dir / "out.csv")

    def test_no_data_source(self, tmp_dir: Path):
        ss = _write_csv(
            tmp_dir / "input.csv",
            [{"sample_id": "s1", "type": "wgs", "sra_id": "", "read1": ""}],
        )
        with pytest.raises(SamplesheetError, match="neither SRA ID nor local"):
            check_samplesheet(ss, tmp_dir / "out.csv")

    def test_nan_sentinels_rejected(self, tmp_dir: Path):
        ss = _write_csv(
            tmp_dir / "input.csv",
            [{"sample_id": "s1", "type": "wgs", "sra_id": "NaN", "read1": "nan"}],
        )
        with pytest.raises(SamplesheetError, match="neither SRA ID nor local"):
            check_samplesheet(ss, tmp_dir / "out.csv")

    def test_header_only_no_samples(self, tmp_dir: Path):
        ss = tmp_dir / "header_only.csv"
        ss.write_text("sample_id,type,sra_id,read1\n")
        with pytest.raises(SamplesheetError, match="No valid samples"):
            check_samplesheet(ss, tmp_dir / "out.csv")


# ── CLI (main) ────────────────────────────────────────────────────
class TestCLI:
    def test_main_success(self, tmp_dir: Path):
        ss = _write_csv(
            tmp_dir / "input.csv",
            [{"sample_id": "s1", "type": "atacseq", "sra_id": "SRR111", "read1": ""}],
        )
        out = tmp_dir / "out.csv"
        main([str(ss), str(out)])
        assert out.exists()

    def test_main_exits_on_error(self, tmp_dir: Path):
        with pytest.raises(SystemExit) as exc_info:
            main([str(tmp_dir / "missing.csv"), str(tmp_dir / "out.csv")])
        assert exc_info.value.code == 1
