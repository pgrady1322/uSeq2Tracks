"""useq2tracks CLI – unified entry point for the pipeline.

Subcommands
-----------
run       Execute the pipeline (Snakemake *or* Nextflow).
validate  Check samplesheet + config without running anything.
hub       Generate a UCSC track-hub from existing bigWig/peak files.
info      Print version, detected samples and assay types.

Author: Patrick Grady
License: MIT License – See LICENSE
"""

from __future__ import annotations

import argparse
import csv
import logging
import shutil
import subprocess
from pathlib import Path
from typing import NoReturn

import useq2tracks
from useq2tracks.check_samplesheet import (
    VALID_ASSAY_TYPES,
    SamplesheetError,
    check_samplesheet,
)
from useq2tracks.generate_ucsc_hub import main as hub_main

logger = logging.getLogger("useq2tracks")

# ── Helpers ───────────────────────────────────────────────────────

_SNAKEMAKE_CMD = "snakemake"
_NEXTFLOW_CMD = "nextflow"


def _find_pipeline_root() -> Path:
    """Walk up from *cwd* to locate the repo root (contains ``Snakefile``).

    Returns:
        Absolute path to the pipeline root directory.

    Raises:
        SystemExit: If no ``Snakefile`` is found in *cwd* or any parent.
    """
    cwd = Path.cwd().resolve()
    for parent in [cwd, *cwd.parents]:
        if (parent / "Snakefile").exists():
            return parent
    _die(
        "Cannot locate pipeline root (no Snakefile found).\n"
        "Run this command from inside the uSeq2Tracks repository."
    )


def _die(msg: str, code: int = 1) -> NoReturn:
    logger.error(msg)
    raise SystemExit(code)


def _require_tool(name: str) -> str:
    """Return the resolved path to *name* or exit with a helpful message."""
    path = shutil.which(name)
    if path is None:
        _die(f"'{name}' not found on PATH.  Install it first.")
    return path


def _load_yaml(path: Path) -> dict:
    """Load a YAML file using the stdlib-only approach (simple key: value).

    Falls back to PyYAML if available; otherwise does a minimal parse
    that handles the flat top-level keys the pipeline cares about.
    """
    try:
        import yaml  # type: ignore[import-untyped]

        with open(path, encoding="utf-8") as fh:
            return yaml.safe_load(fh) or {}
    except ImportError:
        pass

    # Minimal stdlib fallback — handles flat scalars only
    data: dict[str, str] = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            stripped = line.split("#", 1)[0].strip()
            if ":" in stripped:
                key, _, val = stripped.partition(":")
                key = key.strip()
                val = val.strip().strip('"').strip("'")
                if val in ("~", "null", ""):
                    continue
                data[key] = val
    return data


def _count_samples(samplesheet: Path) -> dict[str, list[str]]:
    """Parse the CSV samplesheet and return ``{assay_type: [sample_ids]}``."""
    result: dict[str, list[str]] = {}
    with open(samplesheet, encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            sid = (row.get("sample_id") or "").strip()
            assay = (row.get("type") or "").strip().lower()
            if sid and assay:
                result.setdefault(assay, []).append(sid)
    return result


# ── Subcommand: run ───────────────────────────────────────────────


def _add_run_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "run",
        help="Execute the pipeline via Snakemake or Nextflow",
        description=(
            "Wraps 'snakemake' or 'nextflow run' with sensible defaults.  "
            "Extra arguments after '--' are forwarded verbatim to the engine."
        ),
    )
    p.add_argument(
        "--engine",
        choices=["snakemake", "nextflow"],
        default="snakemake",
        help="Workflow engine (default: snakemake)",
    )
    p.add_argument(
        "-c",
        "--cores",
        type=int,
        default=None,
        help="CPU cores (Snakemake --cores / Nextflow -qs). Default: all",
    )
    p.add_argument(
        "--configfile",
        type=Path,
        default=Path("config.yaml"),
        help="Config YAML (default: config.yaml)",
    )
    p.add_argument(
        "--profile",
        default=None,
        help="Snakemake --profile or Nextflow -profile (e.g. slurm, docker)",
    )
    p.add_argument(
        "-n",
        "--dryrun",
        action="store_true",
        help="Dry-run / preview (Snakemake -n / Nextflow -preview)",
    )
    p.add_argument(
        "--rapid",
        action="store_true",
        help="Enable rapid mode (skip QC, essential tracks only)",
    )
    p.add_argument(
        "--resume",
        action="store_true",
        help="Resume a previous run (Snakemake --rerun-incomplete / Nextflow -resume)",
    )
    p.add_argument(
        "extra",
        nargs=argparse.REMAINDER,
        help="Extra args forwarded to the engine (place after '--')",
    )
    p.set_defaults(func=_cmd_run)


def _build_snakemake_cmd(args: argparse.Namespace, root: Path) -> list[str]:
    cmd = [_require_tool(_SNAKEMAKE_CMD)]
    cmd += ["--snakefile", str(root / "Snakefile")]
    cmd += ["--configfile", str(args.configfile)]

    if args.cores is not None:
        cmd += ["--cores", str(args.cores)]
    else:
        cmd += ["--cores", "all"]

    if args.dryrun:
        cmd.append("--dryrun")
    if args.resume:
        cmd.append("--rerun-incomplete")
    if args.rapid:
        cmd += ["--config", "rapid_mode=true"]
    if args.profile:
        cmd += ["--profile", args.profile]

    # Forward extra arguments (strip leading '--' separator)
    extras = args.extra or []
    if extras and extras[0] == "--":
        extras = extras[1:]
    cmd.extend(extras)
    return cmd


def _build_nextflow_cmd(args: argparse.Namespace, root: Path) -> list[str]:
    cmd = [_require_tool(_NEXTFLOW_CMD), "run", str(root / "main.nf")]
    cmd += ["--samplesheet", str(args.configfile)]

    # Load the config to pass top-level params Nextflow expects
    cfg = _load_yaml(args.configfile)
    for key in ("genome", "genome_id", "outdir", "gtf"):
        val = cfg.get(key)
        if val:
            cmd += [f"--{key}", str(val)]

    if args.cores is not None:
        cmd += ["-qs", str(args.cores)]
    if args.dryrun:
        cmd.append("-preview")
    if args.resume:
        cmd.append("-resume")
    if args.rapid:
        cmd += ["-profile", "rapid"]
    if args.profile:
        cmd += ["-profile", args.profile]

    extras = args.extra or []
    if extras and extras[0] == "--":
        extras = extras[1:]
    cmd.extend(extras)
    return cmd


def _cmd_run(args: argparse.Namespace) -> None:
    root = _find_pipeline_root()

    if not args.configfile.exists():
        _die(f"Config file not found: {args.configfile}")

    if args.engine == "snakemake":
        cmd = _build_snakemake_cmd(args, root)
    else:
        cmd = _build_nextflow_cmd(args, root)

    logger.info("Engine : %s", args.engine)
    logger.info("Command: %s", " ".join(cmd))
    raise SystemExit(subprocess.call(cmd))


# ── Subcommand: validate ─────────────────────────────────────────


def _add_validate_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "validate",
        help="Validate samplesheet and config without running the pipeline",
    )
    p.add_argument(
        "samplesheet",
        type=Path,
        help="Input samplesheet CSV",
    )
    p.add_argument(
        "--configfile",
        type=Path,
        default=None,
        help="Config YAML to check (optional — validates schema keys)",
    )
    p.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Write validated samplesheet to this path (default: stdout summary only)",
    )
    p.set_defaults(func=_cmd_validate)


def _validate_config(cfg_path: Path) -> list[str]:
    """Return a list of issues found in the config YAML (empty = OK)."""
    issues: list[str] = []
    cfg = _load_yaml(cfg_path)

    for key in ("samplesheet", "genome", "genome_id", "outdir"):
        if not cfg.get(key):
            issues.append(f"Missing required config key: '{key}'")

    genome = cfg.get("genome")
    if genome and not Path(genome).exists():
        issues.append(f"genome path does not exist: {genome}")

    gid = cfg.get("genome_id", "")
    if gid:
        import re

        if not re.fullmatch(r"[A-Za-z0-9_.+-]+", gid):
            issues.append(
                f"genome_id '{gid}' contains invalid characters "
                "(allowed: alphanumerics, dots, underscores, hyphens, plus)"
            )
    return issues


def _cmd_validate(args: argparse.Namespace) -> None:
    ok = True

    # ── Config validation ──
    if args.configfile:
        if not args.configfile.exists():
            _die(f"Config file not found: {args.configfile}")
        issues = _validate_config(args.configfile)
        if issues:
            ok = False
            print(f"Config issues ({args.configfile}):")
            for issue in issues:
                print(f"  - {issue}")
        else:
            print(f"Config OK: {args.configfile}")

    # ── Samplesheet validation ──
    out_path = args.output or Path("/dev/null")
    try:
        n = check_samplesheet(args.samplesheet, out_path)
        print(f"Samplesheet OK: {n} samples validated")
        by_type = _count_samples(args.samplesheet)
        for assay, ids in sorted(by_type.items()):
            print(f"  {assay}: {len(ids)} samples")
    except (SamplesheetError, FileNotFoundError) as exc:
        ok = False
        print(f"Samplesheet ERROR: {exc}")

    if not ok:
        raise SystemExit(1)


# ── Subcommand: hub ───────────────────────────────────────────────


def _add_hub_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "hub",
        help="Generate a UCSC track hub from bigWig / peak files",
        description="Delegates to the generate_ucsc_hub module.",
    )
    p.add_argument("--genome-id", required=True, help="Genome identifier")
    p.add_argument("--hub-name", required=True, help="Hub name (no spaces)")
    p.add_argument("--hub-short-label", required=True, help="Hub short label")
    p.add_argument("--hub-long-label", required=True, help="Hub long label")
    p.add_argument("--genome-name", required=True, help="Genome name")
    p.add_argument("--hub-email", required=True, help="Contact email")
    p.add_argument("--bigwigs", nargs="+", help="BigWig files")
    p.add_argument("--peaks", nargs="*", help="Peak files (optional)")
    p.add_argument("--output-dir", default=".", help="Output directory (default: .)")
    p.set_defaults(func=_cmd_hub)


def _cmd_hub(args: argparse.Namespace) -> None:
    # Re-pack into the flat argv list that generate_ucsc_hub.main() expects
    argv: list[str] = []
    argv += ["--genome-id", args.genome_id]
    argv += ["--hub-name", args.hub_name]
    argv += ["--hub-short-label", args.hub_short_label]
    argv += ["--hub-long-label", args.hub_long_label]
    argv += ["--genome-name", args.genome_name]
    argv += ["--hub-email", args.hub_email]
    argv += ["--output-dir", args.output_dir]
    if args.bigwigs:
        argv += ["--bigwigs"] + args.bigwigs
    if args.peaks:
        argv += ["--peaks"] + args.peaks
    hub_main(argv)


# ── Subcommand: info ──────────────────────────────────────────────


def _add_info_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "info",
        help="Show version, detected samples and assay types",
    )
    p.add_argument(
        "--configfile",
        type=Path,
        default=Path("config.yaml"),
        help="Config YAML (default: config.yaml)",
    )
    p.set_defaults(func=_cmd_info)


def _cmd_info(args: argparse.Namespace) -> None:
    print(f"useq2tracks  v{useq2tracks.__version__}")
    print()

    # ── Engine availability ──
    for tool in (_SNAKEMAKE_CMD, _NEXTFLOW_CMD):
        loc = shutil.which(tool)
        status = loc if loc else "not found"
        print(f"  {tool:12s}  {status}")
    print()

    # ── Config summary ──
    cfg_path = args.configfile
    if not cfg_path.exists():
        print(f"Config file not found: {cfg_path}")
        return

    cfg = _load_yaml(cfg_path)
    genome_id = cfg.get("genome_id", "<not set>")
    genome = cfg.get("genome", "<not set>")
    outdir = cfg.get("outdir", "./results")
    rapid = cfg.get("rapid_mode", False)

    print(f"  genome_id    {genome_id}")
    print(f"  genome       {genome}")
    print(f"  outdir       {outdir}")
    print(f"  rapid_mode   {rapid}")
    print()

    # ── Samplesheet summary ──
    ss_path = cfg.get("samplesheet")
    if not ss_path:
        print("  samplesheet  <not set in config>")
        return

    ss = Path(ss_path)
    if not ss.exists():
        print(f"  samplesheet  {ss_path} (file not found)")
        return

    by_type = _count_samples(ss)
    total = sum(len(ids) for ids in by_type.values())
    print(f"  samplesheet  {ss_path}  ({total} samples)")
    for assay in VALID_ASSAY_TYPES:
        ids = by_type.get(assay, [])
        if ids:
            print(
                f"    {assay:14s}  {len(ids):>3d}  {', '.join(ids[:5])}"
                + (" ..." if len(ids) > 5 else "")
            )


# ── Top-level CLI ─────────────────────────────────────────────────


def build_parser() -> argparse.ArgumentParser:
    """Build and return the top-level argument parser."""
    parser = argparse.ArgumentParser(
        prog="useq2tracks",
        description="uSeq2Tracks – Universal Sequencing-to-Genome-Browser-Tracks Pipeline",
    )
    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version=f"%(prog)s {useq2tracks.__version__}",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable debug logging",
    )

    sub = parser.add_subparsers(dest="command", title="commands")
    _add_run_parser(sub)
    _add_validate_parser(sub)
    _add_hub_parser(sub)
    _add_info_parser(sub)

    return parser


def main(argv: list[str] | None = None) -> None:
    """CLI entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    if args.command is None:
        parser.print_help()
        raise SystemExit(0)

    args.func(args)


if __name__ == "__main__":
    main()
