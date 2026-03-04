"""Tests for useq2tracks.cli – the unified CLI entry point."""

from __future__ import annotations

import csv
import textwrap
from pathlib import Path
from unittest.mock import patch

import pytest

from useq2tracks.cli import (
    _build_nextflow_cmd,
    _build_snakemake_cmd,
    _count_samples,
    _validate_config,
    build_parser,
    main,
)

# ── Helpers ───────────────────────────────────────────────────────


def _write_csv(path: Path, rows: list[dict[str, str]]) -> Path:
    fieldnames = list(rows[0].keys()) if rows else ["sample_id", "type"]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return path


def _write_config(
    path: Path,
    *,
    ss: str = "samples.csv",
    genome_id: str = "galGal6",
    genome: str = "genome.fa",
) -> Path:
    path.write_text(
        textwrap.dedent(f"""\
        samplesheet: "{ss}"
        genome: "{genome}"
        genome_id: "{genome_id}"
        outdir: "./results"
        rapid_mode: false
    """),
        encoding="utf-8",
    )
    return path


def _write_minimal_repo(tmp_path: Path) -> Path:
    """Create a minimal repo layout with Snakefile + config + samplesheet."""
    (tmp_path / "Snakefile").write_text("# stub", encoding="utf-8")
    ss = _write_csv(
        tmp_path / "samples.csv",
        [
            {"sample_id": "s1", "type": "atacseq", "sra_id": "SRR111", "read1": ""},
            {"sample_id": "s2", "type": "chipseq", "sra_id": "SRR222", "read1": ""},
        ],
    )
    _write_config(tmp_path / "config.yaml", ss=str(ss))
    (tmp_path / "main.nf").write_text("// stub", encoding="utf-8")
    return tmp_path


# ── build_parser ──────────────────────────────────────────────────


class TestBuildParser:
    def test_no_args_shows_help(self, capsys: pytest.CaptureFixture[str]):
        with pytest.raises(SystemExit, match="0"):
            main([])
        assert "commands" in capsys.readouterr().out.lower()

    def test_version(self, capsys: pytest.CaptureFixture[str]):
        with pytest.raises(SystemExit, match="0"):
            main(["--version"])
        assert "useq2tracks" in capsys.readouterr().out

    def test_subcommands_exist(self):
        parser = build_parser()
        # Subcommands that need no positional args
        for cmd in ("run", "info"):
            ns = parser.parse_known_args([cmd])[0]
            assert ns.command == cmd
        # Subcommands with required positional — just check the
        # parser.choices dict directly
        assert "validate" in parser._subparsers._group_actions[0].choices
        assert "hub" in parser._subparsers._group_actions[0].choices


# ── _count_samples ────────────────────────────────────────────────


class TestCountSamples:
    def test_basic(self, tmp_path: Path):
        ss = _write_csv(
            tmp_path / "s.csv",
            [
                {"sample_id": "s1", "type": "atacseq"},
                {"sample_id": "s2", "type": "atacseq"},
                {"sample_id": "s3", "type": "rnaseq"},
            ],
        )
        result = _count_samples(ss)
        assert result == {"atacseq": ["s1", "s2"], "rnaseq": ["s3"]}

    def test_empty(self, tmp_path: Path):
        ss = tmp_path / "empty.csv"
        ss.write_text("sample_id,type\n", encoding="utf-8")
        assert _count_samples(ss) == {}


# ── _validate_config ─────────────────────────────────────────────


class TestValidateConfig:
    def test_valid(self, tmp_path: Path):
        cfg = _write_config(tmp_path / "c.yaml")
        issues = _validate_config(cfg)
        # genome doesn't exist on disk, that's 1 issue
        assert any("genome path" in i for i in issues)
        # but required keys are present
        assert not any("Missing required" in i for i in issues)

    def test_missing_keys(self, tmp_path: Path):
        cfg = tmp_path / "bad.yaml"
        cfg.write_text("outdir: results\n", encoding="utf-8")
        issues = _validate_config(cfg)
        assert any("genome_id" in i for i in issues)
        assert any("samplesheet" in i for i in issues)

    def test_bad_genome_id(self, tmp_path: Path):
        cfg = _write_config(tmp_path / "c.yaml", genome_id="bad id!")
        issues = _validate_config(cfg)
        assert any("invalid characters" in i for i in issues)

    def test_valid_genome_id_chars(self, tmp_path: Path):
        cfg = _write_config(tmp_path / "c.yaml", genome_id="hg38.p14_v2+patch")
        issues = _validate_config(cfg)
        assert not any("invalid characters" in i for i in issues)


# ── _build_snakemake_cmd ─────────────────────────────────────────


class TestBuildSnakemakeCmd:
    @patch("useq2tracks.cli._require_tool", return_value="/usr/bin/snakemake")
    def test_default(self, _mock, tmp_path: Path):
        parser = build_parser()
        args = parser.parse_args(["run", "--configfile", "config.yaml"])
        cmd = _build_snakemake_cmd(args, tmp_path)
        assert cmd[0] == "/usr/bin/snakemake"
        assert "--cores" in cmd
        assert "all" in cmd
        assert "--snakefile" in cmd

    @patch("useq2tracks.cli._require_tool", return_value="/usr/bin/snakemake")
    def test_dryrun(self, _mock, tmp_path: Path):
        parser = build_parser()
        args = parser.parse_args(["run", "--dryrun"])
        cmd = _build_snakemake_cmd(args, tmp_path)
        assert "--dryrun" in cmd

    @patch("useq2tracks.cli._require_tool", return_value="/usr/bin/snakemake")
    def test_cores(self, _mock, tmp_path: Path):
        parser = build_parser()
        args = parser.parse_args(["run", "--cores", "16"])
        cmd = _build_snakemake_cmd(args, tmp_path)
        assert "16" in cmd

    @patch("useq2tracks.cli._require_tool", return_value="/usr/bin/snakemake")
    def test_rapid(self, _mock, tmp_path: Path):
        parser = build_parser()
        args = parser.parse_args(["run", "--rapid"])
        cmd = _build_snakemake_cmd(args, tmp_path)
        assert "rapid_mode=true" in cmd

    @patch("useq2tracks.cli._require_tool", return_value="/usr/bin/snakemake")
    def test_resume(self, _mock, tmp_path: Path):
        parser = build_parser()
        args = parser.parse_args(["run", "--resume"])
        cmd = _build_snakemake_cmd(args, tmp_path)
        assert "--rerun-incomplete" in cmd

    @patch("useq2tracks.cli._require_tool", return_value="/usr/bin/snakemake")
    def test_profile(self, _mock, tmp_path: Path):
        parser = build_parser()
        args = parser.parse_args(["run", "--profile", "slurm"])
        cmd = _build_snakemake_cmd(args, tmp_path)
        assert "--profile" in cmd
        assert "slurm" in cmd

    @patch("useq2tracks.cli._require_tool", return_value="/usr/bin/snakemake")
    def test_extra_args(self, _mock, tmp_path: Path):
        parser = build_parser()
        args = parser.parse_args(["run", "--", "--printshellcmds"])
        cmd = _build_snakemake_cmd(args, tmp_path)
        assert "--printshellcmds" in cmd


# ── _build_nextflow_cmd ──────────────────────────────────────────


class TestBuildNextflowCmd:
    @patch("useq2tracks.cli._require_tool", return_value="/usr/bin/nextflow")
    def test_default(self, _mock, tmp_path: Path):
        cfg = _write_config(tmp_path / "config.yaml")
        parser = build_parser()
        args = parser.parse_args(["run", "--engine", "nextflow", "--configfile", str(cfg)])
        cmd = _build_nextflow_cmd(args, tmp_path)
        assert cmd[0] == "/usr/bin/nextflow"
        assert "run" in cmd
        assert "--genome_id" in cmd

    @patch("useq2tracks.cli._require_tool", return_value="/usr/bin/nextflow")
    def test_resume(self, _mock, tmp_path: Path):
        cfg = _write_config(tmp_path / "config.yaml")
        parser = build_parser()
        args = parser.parse_args(
            ["run", "--engine", "nextflow", "--resume", "--configfile", str(cfg)]
        )
        cmd = _build_nextflow_cmd(args, tmp_path)
        assert "-resume" in cmd

    @patch("useq2tracks.cli._require_tool", return_value="/usr/bin/nextflow")
    def test_dryrun(self, _mock, tmp_path: Path):
        cfg = _write_config(tmp_path / "config.yaml")
        parser = build_parser()
        args = parser.parse_args(
            ["run", "--engine", "nextflow", "--dryrun", "--configfile", str(cfg)]
        )
        cmd = _build_nextflow_cmd(args, tmp_path)
        assert "-preview" in cmd

    @patch("useq2tracks.cli._require_tool", return_value="/usr/bin/nextflow")
    def test_profile(self, _mock, tmp_path: Path):
        cfg = _write_config(tmp_path / "config.yaml")
        parser = build_parser()
        args = parser.parse_args(
            ["run", "--engine", "nextflow", "--profile", "docker", "--configfile", str(cfg)]
        )
        cmd = _build_nextflow_cmd(args, tmp_path)
        assert "-profile" in cmd
        assert "docker" in cmd


# ── validate subcommand ──────────────────────────────────────────


class TestValidateSubcommand:
    def test_valid_samplesheet(self, tmp_path: Path, capsys: pytest.CaptureFixture[str]):
        ss = _write_csv(
            tmp_path / "s.csv",
            [{"sample_id": "s1", "type": "atacseq", "sra_id": "SRR111", "read1": ""}],
        )
        main(["validate", str(ss)])
        out = capsys.readouterr().out
        assert "OK" in out
        assert "1 samples" in out

    def test_invalid_samplesheet(self, tmp_path: Path):
        ss = tmp_path / "bad.csv"
        ss.write_text("sample_id,type\n", encoding="utf-8")
        with pytest.raises(SystemExit, match="1"):
            main(["validate", str(ss)])

    def test_with_config(self, tmp_path: Path, capsys: pytest.CaptureFixture[str]):
        ss = _write_csv(
            tmp_path / "s.csv",
            [{"sample_id": "s1", "type": "wgs", "sra_id": "SRR111", "read1": ""}],
        )
        # Create a dummy genome file so config validation doesn't flag it
        genome = tmp_path / "genome.fa"
        genome.write_text(">chr1\nACGT\n", encoding="utf-8")
        cfg = _write_config(tmp_path / "c.yaml", ss=str(ss), genome=str(genome))
        main(["validate", str(ss), "--configfile", str(cfg)])
        out = capsys.readouterr().out
        assert "Samplesheet OK" in out

    def test_missing_samplesheet(self, tmp_path: Path):
        with pytest.raises(SystemExit, match="1"):
            main(["validate", str(tmp_path / "nope.csv")])


# ── info subcommand ──────────────────────────────────────────────


class TestInfoSubcommand:
    def test_missing_config(self, tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch):
        monkeypatch.chdir(tmp_path)
        main(["info", "--configfile", str(tmp_path / "nope.yaml")])
        assert "not found" in capsys.readouterr().out

    def test_with_config_and_samples(self, tmp_path: Path, capsys: pytest.CaptureFixture[str]):
        ss = _write_csv(
            tmp_path / "samples.csv",
            [
                {"sample_id": "s1", "type": "atacseq"},
                {"sample_id": "s2", "type": "atacseq"},
                {"sample_id": "s3", "type": "rnaseq"},
            ],
        )
        cfg = _write_config(tmp_path / "c.yaml", ss=str(ss))
        main(["info", "--configfile", str(cfg)])
        out = capsys.readouterr().out
        assert "galGal6" in out
        assert "atacseq" in out
        assert "rnaseq" in out
        assert "3 samples" in out


# ── hub subcommand ───────────────────────────────────────────────


class TestHubSubcommand:
    def test_generates_files(self, tmp_path: Path):
        outdir = tmp_path / "hub_out"
        main(
            [
                "hub",
                "--genome-id",
                "galGal6",
                "--hub-name",
                "TestHub",
                "--hub-short-label",
                "Test",
                "--hub-long-label",
                "Test Hub",
                "--genome-name",
                "galGal6",
                "--hub-email",
                "test@example.com",
                "--output-dir",
                str(outdir),
            ]
        )
        assert (outdir / "hub.txt").exists()
        assert (outdir / "genomes.txt").exists()
        assert (outdir / "trackDb.txt").exists()
