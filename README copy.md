# uSeq2Tracks: Universal Sequencing to Browser Tracks Pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.0.0-brightgreen.svg)](https://snakemake.github.io)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

A comprehensive Snakemake pipeline for processing diverse sequencing datasets and generating standardized genomic tracks for UCSC Genome Browser visualization. uSeq2Tracks handles everything from raw sequencing data to publication-ready browser tracks with minimal user intervention.

## 🌟 Overview

uSeq2Tracks is designed to standardize the processing of heterogeneous sequencing datasets, transforming raw sequencing data into organized, genome browser-ready tracks. Whether you're working with public datasets from ENCODE or your own experimental data, uSeq2Tracks provides a unified workflow that handles quality control, alignment, track generation, and browser hub creation.

### Key Features

- **Universal Input Support**: Handles both local FASTQ files and SRA accessions
- **Multiple Assay Types**: ChIP-seq, ATAC-seq, CUT&RUN, RNA-seq, WGS, Ancient DNA, Long-reads
- **Two Processing Modes**: Standard (full QC) and Rapid (streamlined for public data)
- **Genome-ID Organization**: All outputs tagged with unique genome identifiers
- **UCSC Integration**: Automatic track hub generation for browser visualization
- **Quality Control**: Comprehensive QC reports with FastQC, MultiQC, and assay-specific metrics
- **Flexible Configuration**: Extensive parameter customization for each assay type

## 🧬 Supported Assay Types

| Assay Type | Purpose | Key Outputs | Peak Calling |
|------------|---------|-------------|--------------|
| **ChIP-seq** | Histone modifications, TF binding | BigWig tracks, narrowPeak files | MACS3 |
| **ATAC-seq** | Chromatin accessibility | BigWig tracks, narrowPeak files | MACS3 |
| **CUT&RUN** | Low-input chromatin profiling | BigWig tracks, narrowPeak files | MACS3 |
| **RNA-seq** | Gene expression | BigWig tracks, count matrices | N/A |
| **WGS** | Genome-wide sequencing | BigWig coverage, variant calls | N/A |
| **Ancient DNA** | Historical/archaeological samples | BigWig tracks, damage analysis | N/A |
| **Nanopore** | Long-read sequencing | BigWig tracks, structural variants | N/A |
| **PacBio** | Long-read sequencing | BigWig tracks, high-accuracy variants | N/A |

## 🚀 Quick Start

### 1. Installation

```bash
# Clone the repository
git clone https://github.com/your-org/uSeq2Tracks.git
cd uSeq2Tracks

# Install dependencies with conda/mamba
conda env create -f envs/useq2tracks.yml
conda activate useq2tracks
```

### 2. Configuration

Edit the main configuration file:

```yaml
# config.yaml
samplesheet: "samples.csv"           # Your sample metadata
genome: "/path/to/genome.fa"         # Reference genome
genome_id: "galGal6"                 # REQUIRED: Unique genome identifier
gtf: "/path/to/annotations.gtf"      # Gene annotations (optional)
outdir: "./results"                  # Output directory
rapid_mode: false                    # Set to true for streamlined processing
```

### 3. Sample Sheet Setup

Create a sample sheet describing your data:

```csv
sample_id,type,sra_id,read1,read2,experiment_group,replicate_group,condition
sample1,chipseq,SRR123456,,,H3K27ac,replicate1,treatment
sample2,atacseq,,data/sample2_R1.fastq.gz,data/sample2_R2.fastq.gz,accessibility,replicate1,control
sample3,rnaseq,,data/sample3_R1.fastq.gz,data/sample3_R2.fastq.gz,expression,timepoint1,control
```

### 4. Run the Pipeline

```bash
# Execute with the provided script
./Executor.sh

# Or run directly with Snakemake
snakemake --use-conda --jobs 100
```

## 📊 Pipeline Modes

### Standard Mode (Default)
**Purpose**: Comprehensive analysis with full quality control

**Features**:
- Complete FastQC and MultiQC reports
- Adapter trimming with detailed summaries
- Comprehensive track generation
- Parameter sweep analysis (optional)
- Replicate-merged composite tracks
- Detailed QC metrics and visualizations

**Use Cases**:
- Novel experimental datasets
- Publication-ready analysis
- Unknown data quality
- Parameter optimization studies

### Rapid Mode
**Purpose**: Streamlined processing for high-confidence datasets

**Features**:
- Essential track generation only
- Skipped QC reports (FastQC, MultiQC)
- No adapter trimming summaries
- Core outputs: BigWig tracks, peak calls, UCSC hubs
- Faster processing time
- Simplified output structure

**Use Cases**:
- Public datasets (ENCODE, TCGA, GEO)
- Known high-quality data
- Quick browser track generation
- Time-sensitive analyses

**Activation**:
```yaml
rapid_mode: true
```

## 📁 Output Structure

The pipeline generates genome-tagged outputs with the following structure:

```
results/
└── {genome_id}/                    # e.g., galGal6/
    ├── genome/                     # Genome indices
    │   ├── genome.fa              # Reference genome
    │   ├── star/                  # STAR index
    │   ├── bowtie2/               # Bowtie2 index
    │   └── bwa_mem2/              # BWA-MEM2 index
    ├── qc/                        # Quality control
    │   ├── fastqc/                # FastQC reports
    │   └── multiqc_report.html    # Aggregated QC report
    ├── {assay_type}/              # Per-assay outputs
    │   ├── bam/                   # Aligned reads
    │   ├── bigwig/                # Coverage tracks
    │   └── peaks/                 # Peak calls (when applicable)
    ├── ucsc/                      # UCSC track hubs
    │   ├── hub.txt                # Main hub file
    │   ├── genomes.txt            # Genome specification
    │   └── trackDb.txt            # Track definitions
    └── rapid/                     # Rapid mode summary (if enabled)
        └── rapid_tracks_complete.txt
```

## 🔧 Configuration Details

### Required Parameters

```yaml
# Essential settings - must be configured
samplesheet: "samples.csv"          # Sample metadata file
genome: "/path/to/genome.fa"        # Reference genome FASTA
genome_id: "your_genome_id"         # REQUIRED: Unique identifier
outdir: "./results"                 # Output directory
```

### Sample Sheet Format

| Column | Description | Example | Required |
|--------|-------------|---------|----------|
| `sample_id` | Unique sample identifier | `H3K27ac_rep1` | Yes |
| `type` | Assay type | `chipseq`, `atacseq`, `rnaseq` | Yes |
| `sra_id` | SRA accession (if downloading) | `SRR123456` | If no local files |
| `read1` | Path to R1 FASTQ | `data/sample_R1.fastq.gz` | If no SRA |
| `read2` | Path to R2 FASTQ | `data/sample_R2.fastq.gz` | For paired-end |
| `experiment_group` | Experimental grouping | `H3K27ac`, `timepoint1` | No |
| `replicate_group` | Replicate grouping | `replicate1` | No |
| `condition` | Sample condition | `treatment`, `control` | No |

### Assay-Specific Parameters

#### ChIP-seq Configuration
```yaml
chipseq:
  mapper: "bowtie2"                 # bowtie2 or bwa_mem2
  bowtie2_opts: "--very-sensitive"
  markdup: false                    # Mark duplicates
  macs3_opts: "--qval 0.05 --keep-dup all"
  control_tag: "input"              # Identify control samples
  bw_norm: "CPM"                    # BigWig normalization
```

#### ATAC-seq Configuration
```yaml
atacseq:
  mapper: "bowtie2"
  markdup: false                    # Recommended: false for accessibility
  shift: -75                        # MACS3 shift for ATAC-seq
  extsize: 150                      # MACS3 extension size
  macs3_opts: "--qval 0.05"
```

#### RNA-seq Configuration
```yaml
rnaseq:
  mapper: "star"                    # star or hisat2
  star_opts: "--outFilterMultimapNmax 20"
  strand_specific: false
  gene_bed: "genes.bed"             # For QC analysis
```

### Advanced Features

#### Parameter Sweep
Test multiple peak-calling thresholds:
```yaml
parameter_sweep:
  enabled: true
  qvalues: [0.05, 0.01, 0.005, 0.001]
```

#### Adapter Trimming
Enable quality-based trimming:
```yaml
adapter_trimming:
  enabled: true
  min_length: 20
  quality_cutoff: 20
```

#### Alternative Peak Callers
Enable Genrich for ATAC-seq:
```yaml
## 🧮 Computational Requirements

### Resource Recommendations

| Dataset Size | CPU Cores | Memory | Storage | Time Estimate |
|--------------|-----------|---------|---------|---------------|
| Small (< 10 samples) | 8-16 | 32 GB | 100 GB | 2-6 hours |
| Medium (10-50 samples) | 16-32 | 64 GB | 500 GB | 6-24 hours |
| Large (50+ samples) | 32-64 | 128 GB | 1 TB+ | 1-3 days |

### Cluster Configuration

The pipeline includes SLURM integration via `Executor.sh`. Customize for your cluster:

```bash
#SBATCH --job-name=uSeq2Tracks
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=10g
#SBATCH --cpus-per-task=4
```

## 📖 Detailed Workflows

### ChIP-seq/ATAC-seq/CUT&RUN Workflow
1. **Quality Control**: FastQC analysis of raw reads
2. **Adapter Trimming**: Optional quality-based trimming with fastp
3. **Genome Indexing**: Build mapper-specific indices (Bowtie2/BWA-MEM2)
4. **Read Mapping**: Align reads to reference genome
5. **Post-processing**: Sort, index, optional duplicate marking
6. **Coverage Tracks**: Generate normalized BigWig files
7. **Peak Calling**: Identify enriched regions with MACS3
8. **Quality Metrics**: Generate assay-specific QC reports

### RNA-seq Workflow
1. **Quality Control**: FastQC analysis of raw reads
2. **Adapter Trimming**: Optional preprocessing with fastp
3. **Genome Indexing**: Build STAR or HISAT2 indices
4. **Read Mapping**: Splice-aware alignment to reference
5. **Quantification**: Generate gene-level count matrices
6. **Coverage Tracks**: Create strand-specific BigWig files
7. **Quality Assessment**: RNA-seq specific QC with RSeQC

### Long-read Workflow (Nanopore/PacBio)
1. **Quality Assessment**: Basic statistics and length distributions
2. **Genome Indexing**: Build Minimap2 index
3. **Read Mapping**: Long-read aware alignment
4. **Coverage Analysis**: Generate coverage tracks
5. **Variant Calling**: Optional structural variant detection

## 🎨 UCSC Browser Integration

### Automatic Track Hub Generation

uSeq2Tracks automatically creates UCSC-compatible track hubs:

**Hub Structure**:
```
ucsc/
├── hub.txt                 # Hub metadata
├── genomes.txt             # Genome definitions
└── trackDb.txt             # Track configurations
```

**Track Organization**:
- **Composite Tracks**: Group related samples (e.g., same experiment)
- **Subgroups**: Organize by condition, replicate, timepoint
- **Color Coding**: Consistent color schemes per assay type
- **Metadata Integration**: Sample information in track descriptions

### Loading in UCSC Browser

1. Upload track hub files to web-accessible location
2. In UCSC Browser: `My Data` → `Track Hubs` → `My Hubs`
3. Enter hub URL: `https://your-server.com/path/to/ucsc/hub.txt`
4. Browse your data with full metadata integration

## 🔍 Quality Control Features

### Standard QC Reports
- **FastQC**: Per-sample quality metrics
- **MultiQC**: Aggregated quality dashboard
- **Mapping Statistics**: Alignment rates and quality scores
- **Library Complexity**: Duplication rates and insert sizes

### Assay-Specific QC
- **ChIP-seq**: Fragment length distributions, enrichment metrics
- **ATAC-seq**: TSS enrichment, fragment size profiles
- **RNA-seq**: Gene body coverage, junction analysis
- **WGS**: Coverage uniformity, variant quality metrics

### Quality Thresholds
The pipeline includes built-in quality checks:
- Minimum mapping rates
- Fragment count requirements for peak calling
- Insert size validation
- Strand specificity assessment

## 🛠 Troubleshooting

### Common Issues

#### Genome ID Not Set
```
ERROR: genome_id must be set in config.yaml
```
**Solution**: Add `genome_id: "your_genome"` to config.yaml

#### Sample Sheet Formatting
```
ERROR: Missing required columns in sample sheet
```
**Solution**: Ensure sample sheet includes `sample_id`, `type`, and either `sra_id` or local file paths

#### Memory Issues
```
ERROR: Job exceeded memory limit
```
**Solution**: Increase memory allocation in config.yaml:
```yaml
memory:
  large: 128000    # Increase for memory-intensive jobs
```

#### Disk Space
```
ERROR: No space left on device
```
**Solution**: 
- Clean up intermediate files: `snakemake --delete-temp-output`
- Use scratch storage for temporary files
- Monitor disk usage during execution

### Performance Optimization

#### Speed Up Processing
1. **Enable Rapid Mode** for public datasets
2. **Increase Parallelization**: More `--jobs` in Snakemake
3. **Use SSDs** for scratch space
4. **Optimize Resource Allocation**: Match CPU/memory to job requirements

#### Reduce Storage
1. **Delete Intermediate Files**: Use `--delete-temp-output`
2. **Compress Outputs**: Enable compression for BAM files
3. **Archive Unused Data**: Move completed analyses to long-term storage

## 📚 Examples

### Example 1: ENCODE ChIP-seq Analysis
```yaml
# config.yaml
genome_id: "hg38_ENCODE_H3K27ac"
genome: "/data/genomes/hg38.fa"
samplesheet: "encode_chipseq.csv"
rapid_mode: true                     # Skip QC for public data

chipseq:
  mapper: "bowtie2"
  macs3_opts: "--qval 0.01"
```

```csv
# encode_chipseq.csv
sample_id,type,sra_id,read1,read2,experiment_group,replicate_group,condition
ENCSR000EWQ_rep1,chipseq,SRR1536404,,,H3K27ac,ENCSR000EWQ,K562
ENCSR000EWQ_rep2,chipseq,SRR1536405,,,H3K27ac,ENCSR000EWQ,K562
ENCSR000EWQ_input,chipseq,SRR1536406,,,H3K27ac,input,input
```

### Example 2: Multi-assay Developmental Study
```yaml
# config.yaml
genome_id: "mm10_development"
genome: "/data/genomes/mm10.fa"
gtf: "/data/annotations/mm10.gtf"
samplesheet: "development_study.csv"
rapid_mode: false                    # Full QC for novel data

parameter_sweep:
  enabled: true
  qvalues: [0.1, 0.05, 0.01, 0.001]
```

```csv
# development_study.csv
sample_id,type,sra_id,read1,read2,experiment_group,replicate_group,condition
E10_ATAC_rep1,atacseq,,data/E10_ATAC_1_R1.fq.gz,data/E10_ATAC_1_R2.fq.gz,E10,ATAC_rep1,E10
E10_ATAC_rep2,atacseq,,data/E10_ATAC_2_R1.fq.gz,data/E10_ATAC_2_R2.fq.gz,E10,ATAC_rep2,E10
E10_RNA_rep1,rnaseq,,data/E10_RNA_1_R1.fq.gz,data/E10_RNA_1_R2.fq.gz,E10,RNA_rep1,E10
E12_ATAC_rep1,atacseq,,data/E12_ATAC_1_R1.fq.gz,data/E12_ATAC_1_R2.fq.gz,E12,ATAC_rep1,E12
E12_RNA_rep1,rnaseq,,data/E12_RNA_1_R1.fq.gz,data/E12_RNA_1_R2.fq.gz,E12,RNA_rep1,E12
```

### Example 3: Ancient DNA Analysis
```yaml
# config.yaml
genome_id: "ancientDNA_sample"
genome: "/data/reference/ancient_ref.fa"
samplesheet: "ancient_samples.csv"

ancientdna:
  mapper: "bwa_aln"
  markdup: true
  damage_analysis: true
  min_mapq: 30
```

## 🤝 Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Development Setup
```bash
git clone https://github.com/your-org/uSeq2Tracks.git
cd uSeq2Tracks
conda env create -f envs/development.yml
conda activate useq2tracks-dev
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📞 Support

- **Documentation**: [GitHub Wiki](https://github.com/your-org/uSeq2Tracks/wiki)
- **Issues**: [GitHub Issues](https://github.com/your-org/uSeq2Tracks/issues)
- **Discussions**: [GitHub Discussions](https://github.com/your-org/uSeq2Tracks/discussions)

## 🏆 Citation

If you use uSeq2Tracks in your research, please cite:

```
Your Name et al. (2025). uSeq2Tracks: A universal pipeline for sequencing data to genome browser tracks. 
Journal Name, Volume(Issue), pages. DOI: 10.xxxx/xxxxx
```

## 🙏 Acknowledgments

- **Snakemake Community**: For the excellent workflow management system
- **Bioconda**: For streamlined software distribution
- **UCSC Genome Browser**: For track hub specifications
- **Tool Developers**: FastQC, MultiQC, STAR, Bowtie2, BWA, MACS3, and all other integrated tools

---

**uSeq2Tracks**: From raw sequencing data to publication-ready browser tracks in one streamlined workflow.
