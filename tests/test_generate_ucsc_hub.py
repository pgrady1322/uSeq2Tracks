"""Tests for useq2tracks.generate_ucsc_hub."""

from __future__ import annotations

from pathlib import Path

import pytest

from useq2tracks.generate_ucsc_hub import (
    TRACK_COLORS,
    main,
    parse_args,
    write_genomes_txt,
    write_hub_txt,
    write_trackdb_txt,
)


# ── Helpers ───────────────────────────────────────────────────────
def _base_args(**overrides) -> list[str]:
    """Build a minimal valid CLI arg list, with optional overrides."""
    defaults = {
        "--genome-id": "galGal6",
        "--hub-name": "TestHub",
        "--hub-short-label": "Test",
        "--hub-long-label": "Test Hub Long Label",
        "--genome-name": "galGal6",
        "--hub-email": "test@example.com",
    }
    defaults.update(overrides)
    args: list[str] = []
    for k, v in defaults.items():
        args.extend([k, v])
    return args


# ── parse_args ────────────────────────────────────────────────────
class TestParseArgs:
    def test_required_args(self):
        ns = parse_args(_base_args())
        assert ns.genome_id == "galGal6"
        assert ns.hub_name == "TestHub"
        assert ns.hub_email == "test@example.com"

    def test_defaults(self):
        ns = parse_args(_base_args())
        assert ns.output_dir == "."
        assert ns.bigwigs is None
        assert ns.peaks is None

    def test_bigwigs_and_peaks(self):
        extra = _base_args() + ["--bigwigs", "a.bw", "b.bw", "--peaks", "a_peaks.narrowPeak"]
        ns = parse_args(extra)
        assert ns.bigwigs == ["a.bw", "b.bw"]
        assert ns.peaks == ["a_peaks.narrowPeak"]

    def test_missing_required_exits(self):
        with pytest.raises(SystemExit):
            parse_args(["--genome-id", "x"])  # missing other required args


# ── write_hub_txt ─────────────────────────────────────────────────
class TestWriteHubTxt:
    def test_creates_hub_file(self, tmp_path: Path):
        ns = parse_args(_base_args())
        write_hub_txt(ns, tmp_path)
        hub = tmp_path / "hub.txt"
        assert hub.exists()
        content = hub.read_text()
        assert "hub TestHub" in content
        assert "email test@example.com" in content
        assert "genomesFile genomes.txt" in content

    def test_labels_written(self, tmp_path: Path):
        ns = parse_args(_base_args(**{"--hub-short-label": "Short", "--hub-long-label": "Long"}))
        write_hub_txt(ns, tmp_path)
        content = (tmp_path / "hub.txt").read_text()
        assert "shortLabel Short" in content
        assert "longLabel Long" in content


# ── write_genomes_txt ─────────────────────────────────────────────
class TestWriteGenomesTxt:
    def test_creates_genomes_file(self, tmp_path: Path):
        ns = parse_args(_base_args())
        write_genomes_txt(ns, tmp_path)
        content = (tmp_path / "genomes.txt").read_text()
        assert "genome galGal6" in content
        assert "trackDb trackDb.txt" in content


# ── write_trackdb_txt ─────────────────────────────────────────────
class TestWriteTrackDbTxt:
    def test_empty_bigwigs(self, tmp_path: Path):
        ns = parse_args(_base_args())
        write_trackdb_txt(ns, tmp_path)
        content = (tmp_path / "trackDb.txt").read_text()
        assert "compositeTrack on" in content
        # No sample tracks
        assert "track " in content  # header exists
        assert "_signal" not in content

    def test_bigwig_tracks_written(self, tmp_path: Path):
        args = _base_args() + ["--bigwigs", "sample1.bw", "sample2.bw"]
        ns = parse_args(args)
        write_trackdb_txt(ns, tmp_path)
        content = (tmp_path / "trackDb.txt").read_text()
        assert "track sample1_signal" in content
        assert "track sample2_signal" in content
        assert "bigDataUrl ../bigwig/sample1.bw" in content
        assert "autoScale on" in content

    def test_peak_tracks_matched(self, tmp_path: Path):
        args = _base_args() + [
            "--bigwigs",
            "s1.bw",
            "--peaks",
            "s1_peaks.narrowPeak",
        ]
        ns = parse_args(args)
        write_trackdb_txt(ns, tmp_path)
        content = (tmp_path / "trackDb.txt").read_text()
        assert "track s1_signal" in content
        assert "track s1_peaks" in content
        assert "type narrowPeak" in content

    def test_unmatched_peak_skipped(self, tmp_path: Path):
        args = _base_args() + [
            "--bigwigs",
            "s1.bw",
            "--peaks",
            "orphan_peaks.narrowPeak",
        ]
        ns = parse_args(args)
        write_trackdb_txt(ns, tmp_path)
        content = (tmp_path / "trackDb.txt").read_text()
        assert "track orphan_peaks" not in content

    def test_color_cycling(self, tmp_path: Path):
        bws = [f"s{i}.bw" for i in range(len(TRACK_COLORS) + 2)]
        args = _base_args() + ["--bigwigs"] + bws
        ns = parse_args(args)
        write_trackdb_txt(ns, tmp_path)
        content = (tmp_path / "trackDb.txt").read_text()
        # First and (N+1)th sample should share a color
        assert content.count(TRACK_COLORS[0]) >= 2


# ── main (integration) ───────────────────────────────────────────
class TestMain:
    def test_full_run(self, tmp_path: Path):
        args = _base_args(**{"--output-dir": str(tmp_path)}) + [
            "--bigwigs",
            "a.bw",
            "b.bw",
            "--peaks",
            "a_peaks.narrowPeak",
        ]
        main(args)
        assert (tmp_path / "hub.txt").exists()
        assert (tmp_path / "genomes.txt").exists()
        assert (tmp_path / "trackDb.txt").exists()

    def test_output_dir_created(self, tmp_path: Path):
        out = tmp_path / "nested" / "hub"
        args = _base_args(**{"--output-dir": str(out)}) + ["--bigwigs", "x.bw"]
        main(args)
        assert (out / "hub.txt").exists()
