# uSeq2Tracks Nextflow - Quick Reference

## Quick Start Commands

### Basic Execution
```bash
# With Docker
nextflow run main.nf -profile docker -params-file examples/params.yaml

# With Singularity
nextflow run main.nf -profile singularity -params-file examples/params.yaml

# With SLURM + Singularity
nextflow run main.nf -profile slurm,singularity -params-file examples/params.yaml

# Rapid mode (skip QC)
nextflow run main.nf -profile rapid,docker -params-file examples/params.yaml

# Resume failed run
nextflow run main.nf -profile docker -resume
```

## File Organization

### Input Files Required
```
samples.csv          # Samplesheet with sample metadata
genome.fa            # Reference genome FASTA
annotations.gtf      # Gene annotations (optional, for RNA-seq)
```

### Configuration Files
```
params.yaml          # Pipeline parameters (recommended)
nextflow.config      # Main configuration (edit for advanced usage)
custom.config        # Optional custom configuration
```

## Samplesheet Format

### Required Columns
- `sample_id`: Unique sample identifier
- `type`: Assay type (atacseq, chipseq, cutrun, rnaseq, wgs, etc.)

### Data Source (choose one)
- `sra_id`: SRA accession number OR
- `read1`, `read2`: Local FASTQ file paths

### Optional Columns
- `experiment_group`: Group samples by experiment
- `replicate_group`: Group replicates for merging
- `condition`: Sample condition/treatment

### Example
```csv
sample_id,type,sra_id,read1,read2,experiment_group,replicate_group,condition
sample1,atacseq,,data/s1_R1.fq.gz,data/s1_R2.fq.gz,exp1,atac,treatment
sample2,chipseq,SRR123456,,,exp1,H3K27ac,K562
sample3,chipseq,SRR123457,,,exp1,input,input
```

## Key Parameters

### Core Settings
```yaml
genome_id: "galGal6"           # REQUIRED - genome identifier
rapid_mode: false              # Skip QC for faster processing
outdir: "./results"            # Output directory
```

### ATAC-seq
```yaml
atacseq:
  mapper: "bowtie2"            # bowtie2 or bwa_mem2
  markdup: false               # Mark/remove duplicates
  shift: -75                   # MACS3 shift parameter
  extsize: 150                 # MACS3 extension size
  bw_norm: "CPM"               # BigWig normalization
```

### ChIP-seq
```yaml
chipseq:
  mapper: "bowtie2"
  control_tag: "input"         # Identifier for control samples
  markdup: false
  macs3_opts: "--qval 0.05"
  bw_norm: "CPM"
```

### CUT&RUN
```yaml
cutrun:
  mapper: "bowtie2"
  markdup: false
  shift: 0
  extsize: 160
  bw_norm: "CPM"
```

## Profile Reference

### Container Profiles
- `docker`: Use Docker (recommended for local)
- `singularity`: Use Singularity (recommended for HPC)
- `conda`: Use Conda environments

### Executor Profiles
- `local`: Run on local machine (default)
- `slurm`: Submit to SLURM scheduler
- `sge`: Submit to SGE scheduler
- `pbs`: Submit to PBS scheduler

### Special Profiles
- `test`: Use test dataset
- `rapid`: Skip QC steps
- `debug`: Enable debugging output

### Combining Profiles
Use comma to combine profiles:
```bash
-profile slurm,singularity      # SLURM with Singularity
-profile docker,rapid            # Docker with rapid mode
```

## Output Structure

```
results/
└── <genome_id>/
    ├── atacseq/
    │   ├── bam/                # Aligned reads
    │   ├── peaks/              # Peak calls
    │   └── bigwig/             # Signal tracks
    ├── chipseq/
    │   ├── bam/
    │   ├── peaks/
    │   └── bigwig/
    ├── qc/
    │   ├── fastqc/             # FastQC reports
    │   └── multiqc_report.html # MultiQC summary
    ├── ucsc/
    │   ├── hub.txt             # UCSC hub file
    │   ├── genomes.txt
    │   └── trackDb.txt         # Track definitions
    └── pipeline_info/          # Execution reports
```

## Resource Management

### Set Maximum Resources
```bash
nextflow run main.nf --max_cpus 32 --max_memory 128.GB --max_time 48.h
```

### In params.yaml
```yaml
max_cpus: 64
max_memory: "256.GB"
max_time: "240.h"
```

## Common Use Cases

### Process Local FASTQ Files
```csv
sample_id,type,read1,read2
sample1,atacseq,data/sample1_R1.fastq.gz,data/sample1_R2.fastq.gz
```

### Download from SRA
```csv
sample_id,type,sra_id
ENCODE_sample,chipseq,SRR1536404
```

### ChIP-seq with Input Control
```csv
sample_id,type,sra_id,experiment_group,replicate_group,condition
ChIP_H3K27ac,chipseq,SRR1,exp1,H3K27ac,K562
ChIP_Input,chipseq,SRR2,exp1,input,input
```

## Troubleshooting

### Pipeline Won't Start
```bash
# Check configuration
nextflow config main.nf

# Validate samplesheet
head -n 5 samples.csv

# Check Nextflow version
nextflow -version  # Should be >=21.10.3
```

### Out of Memory
```bash
# Increase memory limit
nextflow run main.nf --max_memory 512.GB

# Or retry with more resources
nextflow run main.nf -resume
```

### Resume Not Working
```bash
# Clean work directory and restart
rm -rf work/
nextflow run main.nf
```

### Container Issues
```bash
# Pull containers beforehand
nextflow pull main.nf -profile docker

# Or use Singularity
nextflow run main.nf -profile singularity
```

## Advanced Features

### Custom Configuration
```bash
# Use custom config file
nextflow run main.nf -c custom.config
```

### Override Parameters
```bash
# Command-line parameter override
nextflow run main.nf --genome my_genome.fa --genome_id myGenome
```

### Generate Reports
```bash
nextflow run main.nf \\
    -with-report report.html \\
    -with-timeline timeline.html \\
    -with-dag dag.svg
```

### Limit Concurrent Jobs
```bash
# Limit to 50 simultaneous jobs
nextflow run main.nf -profile slurm -qs 50
```

## Getting Help

### View Help Message
```bash
nextflow run main.nf --help
```

### Check Configuration
```bash
nextflow config main.nf -show-profiles
```

### Validate Pipeline
```bash
nextflow run main.nf -profile test
```

## Next Steps After Installation

1. **Install nf-core modules** (required):
   ```bash
   # Install nf-core tools
   pip install nf-core
   
   # Install required modules
   cd nextflow/
   nf-core modules install fastqc
   nf-core modules install bowtie2/align
   ...
   ```

2. **Prepare your data**:
   - Create samplesheet
   - Prepare reference genome
   - Configure parameters

3. **Test run**:
   ```bash
   nextflow run main.nf -profile test,docker
   ```

4. **Production run**:
   ```bash
   nextflow run main.nf -profile docker -params-file params.yaml
   ```

## Support

- Documentation: See [README.md](README.md)
- Examples: See [examples/](examples/) directory

---

**Quick Tip**: Start with the test profile to verify your installation before running production data.
