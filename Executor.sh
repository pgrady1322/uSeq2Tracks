#!/bin/bash
#SBATCH --job-name=uSeq2Tracks
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=10g
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

# uSeq2Tracks Pipeline Executor
# Runs the Snakemake workflow with appropriate cluster configuration

echo "🧬 Starting uSeq2Tracks Pipeline..."
echo "📅 $(date)"
echo "📁 Working directory: $(pwd)"

# Check for required files
if [[ ! -f "config.yaml" ]]; then
    echo "❌ ERROR: config.yaml not found in current directory"
    exit 1
fi

if [[ ! -f "Snakefile" ]]; then
    echo "❌ ERROR: Snakefile not found in current directory"
    exit 1
fi

# Check if genome_id is set in config
GENOME_ID=$(grep "^genome_id:" config.yaml | cut -d'"' -f2 | tr -d ' ')
if [[ -z "$GENOME_ID" || "$GENOME_ID" == "" ]]; then
    echo "❌ ERROR: genome_id must be set in config.yaml"
    echo "Example: genome_id: \"galGal6\""
    exit 1
fi

echo "🔬 Genome ID: $GENOME_ID"

# Load required modules (adjust for your cluster)
# Note: Comment out module loading if running locally
# module load miniconda/22.11.1-1
# module load snakemake/7.32.4

# Check if snakemake is available
if ! command -v snakemake &> /dev/null; then
    echo "❌ ERROR: snakemake command not found"
    echo "Please install snakemake or activate a conda environment with snakemake"
    echo "Example: conda activate snakemake_env"
    exit 1
fi

echo "🐍 Using snakemake: $(which snakemake)"

# Run Snakemake with your global cluster configuration
snakemake \
    --use-conda \
    --conda-frontend mamba \
    --jobs 500 \
    --latency-wait 60 \
    --rerun-incomplete \
    --printshellcmds \
    "$@"

echo "✅ Pipeline execution completed!"
echo "📅 $(date)"