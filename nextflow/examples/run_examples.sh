#!/bin/bash
#
# Example execution script for uSeq2Tracks Nextflow pipeline
# This demonstrates various ways to run the pipeline

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}uSeq2Tracks Nextflow Pipeline - Execution Examples${NC}\n"

# Function to print section headers
print_section() {
    echo -e "\n${YELLOW}========================================${NC}"
    echo -e "${YELLOW}$1${NC}"
    echo -e "${YELLOW}========================================${NC}\n"
}

# Check if nextflow is installed
if ! command -v nextflow &> /dev/null; then
    echo -e "${RED}ERROR: Nextflow is not installed!${NC}"
    echo "Install with: curl -s https://get.nextflow.io | bash"
    exit 1
fi

# Example 1: Basic local execution with Docker
print_section "Example 1: Local Execution with Docker"
cat << 'EOF'
nextflow run main.nf \
    -profile docker \
    -params-file params.yaml
EOF

# Example 2: SLURM cluster execution with Singularity
print_section "Example 2: SLURM Cluster with Singularity"
cat << 'EOF'
nextflow run main.nf \
    -profile slurm,singularity \
    -params-file params.yaml \
    -qs 100
EOF

# Example 3: Rapid mode (skip QC)
print_section "Example 3: Rapid Mode (Skip QC)"
cat << 'EOF'
nextflow run main.nf \
    -profile rapid,docker \
    -params-file params.yaml
EOF

# Example 4: Resume failed run
print_section "Example 4: Resume Failed Run"
cat << 'EOF'
nextflow run main.nf \
    -profile docker \
    -params-file params.yaml \
    -resume
EOF

# Example 5: Custom resource limits
print_section "Example 5: Custom Resource Limits"
cat << 'EOF'
nextflow run main.nf \
    -profile docker \
    -params-file params.yaml \
    --max_cpus 32 \
    --max_memory 128.GB
EOF

# Example 6: With custom configuration
print_section "Example 6: With Custom Config"
cat << 'EOF'
nextflow run main.nf \
    -profile docker \
    -params-file params.yaml \
    -c custom.config
EOF

# Example 7: Test run
print_section "Example 7: Test Run (Minimal Dataset)"
cat << 'EOF'
nextflow run main.nf \
    -profile test,docker
EOF

# Example 8: Generate execution report
print_section "Example 8: With Detailed Reports"
cat << 'EOF'
nextflow run main.nf \
    -profile docker \
    -params-file params.yaml \
    -with-report report.html \
    -with-timeline timeline.html \
    -with-dag dag.svg
EOF

# Example 9: Specify output directory
print_section "Example 9: Custom Output Directory"
cat << 'EOF'
nextflow run main.nf \
    -profile docker \
    -params-file params.yaml \
    --outdir /path/to/custom/output
EOF

# Example 10: Use specific genome and samplesheet
print_section "Example 10: Command-line Override of Parameters"
cat << 'EOF'
nextflow run main.nf \
    -profile docker \
    --samplesheet my_samples.csv \
    --genome my_genome.fa \
    --genome_id myGenome \
    --outdir results_myGenome
EOF

# Actual execution section
print_section "Ready to Run?"
echo "Edit the examples above with your specific parameters and execute."
echo ""
echo "Quick start:"
echo -e "${GREEN}nextflow run main.nf -profile docker -params-file params.yaml${NC}"
echo ""
echo "For help:"
echo -e "${GREEN}nextflow run main.nf --help${NC}"
