# uSeq2Tracks Pipeline - Clean Setup

## Overview
This pipeline processes sequencing data (RNA-seq, ChIP-seq, WGS) for UCSC track generation, optimized for SLURM cluster execution without conda environment issues.

## Quick Start

### 1. Initial Setup (run once)
```bash
chmod +x manage_rules.sh create_no_conda_rules.sh run_simple_no_conda.sh
./manage_rules.sh setup
```

### 2. Run Pipeline
```bash
./run_simple_no_conda.sh
```

## File Structure

### Core Files
- `Snakefile` - Main workflow definition
- `config.yaml` - Pipeline configuration
- `samples.csv` - Sample metadata with layout information
- `validate_pipeline.sh` - Pipeline validation and setup guide

### Execution Scripts
- `run_simple_no_conda.sh` - **RECOMMENDED**: Uses existing conda environment
- `run_slurm.sh` - Standard execution (may have conda issues on some clusters)
- `manage_rules.sh` - Switch between conda/no-conda rule versions

### Rule Management
- `rules/` - Active rule files (currently: no-conda version)
- `rules_with_conda_archived/` - Original rules with conda directives
- `rules_no_conda/` - Clean rules without conda directives
- `create_no_conda_rules.sh` - Creates the archive structure

### Configuration
- `profiles/` - SLURM execution profiles
- `envs/` - Conda environment files (for reference)

## Rule Version Management

```bash
# Check current version
./manage_rules.sh status

# Switch to no-conda rules (for cluster)
./manage_rules.sh use-noconda

# Switch to conda rules (for development)
./manage_rules.sh use-conda
```

## Troubleshooting

### Conda Environment Issues
If you encounter conda file locking issues:
1. Ensure you're using the no-conda rules: `./manage_rules.sh use-noconda`
2. Run with existing environment: `./run_simple_no_conda.sh`

### SLURM Issues
- Check QoS settings in `profiles/slurm_general_general/config.yaml`
- Verify SLURM partition availability
- Use `--dry-run` to test without execution

## Data Types Supported
- Single-end and paired-end RNA-seq
- Single-end and paired-end ChIP-seq  
- Whole genome sequencing (WGS)
- SRA data download and processing

## Output
- Processed BAM files
- BigWig tracks for UCSC
- Quality control reports
- Peak calling results (ChIP-seq)
