/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Local to the pipeline
//
include { SAMPLESHEET_CHECK       } from '../modules/local/samplesheet_check'
include { GENOME_PREP             } from '../modules/local/genome_prep'
include { UCSC_HUB                } from '../modules/local/ucsc_hub'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { RNASEQ                  } from '../workflows/subworkflows/rnaseq'
// TODO: Uncomment when nf-core modules are installed
// include { ATACSEQ                 } from '../workflows/subworkflows/atacseq'
// include { CHIPSEQ                 } from '../workflows/subworkflows/chipseq'
// include { CUTRUN                  } from '../workflows/subworkflows/cutrun'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
// TODO: Uncomment when nf-core modules are installed
// include { FASTQC                  } from '../modules/nf-core/fastqc/main'
// include { MULTIQC                 } from '../modules/nf-core/multiqc/main'
// include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow USEQ2TRACKS {

    // Validate input parameters
    WorkflowUseq2tracks.initialise(params, log)

    // Check mandatory parameters
    if (!params.samplesheet) { 
        error 'Input samplesheet not specified! Use --samplesheet <path/to/samplesheet.csv>'
    }

    if (!params.genome) { 
        error 'Reference genome not specified! Use --genome <path/to/genome.fa>'
    }

    if (!params.genome_id) {
        error 'Genome ID not specified! Please provide a genome identifier with --genome_id <ID>'
    }

    def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)
    ch_samplesheet = file(params.samplesheet, checkIfExists: true)
    ch_genome = file(params.genome, checkIfExists: true)

    // Config files
    ch_multiqc_config        = Channel.fromPath("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

    ch_versions = Channel.empty()

    //
    // MODULE: Check and validate samplesheet
    //
    SAMPLESHEET_CHECK (
        ch_samplesheet
    )
    ch_versions = ch_versions.mix(SAMPLESHEET_CHECK.out.versions)
    
    //
    // Create channels from validated samplesheet
    //
    SAMPLESHEET_CHECK.out.csv
        .splitCsv(header:true, sep:',')
        .map { row ->
            def meta = [:]
            meta.id               = row.sample_id
            meta.type             = row.type
            meta.single_end       = false  // Default to paired-end
            meta.experiment_group = row.experiment_group ?: 'exp1'
            meta.replicate_group  = row.replicate_group ?: (row.condition ?: row.type)
            meta.condition        = row.condition ?: row.sample_id
            
            // Determine data source
            if (row.sra_id && row.sra_id != '' && row.sra_id != 'NA') {
                meta.source = 'sra'
                return [meta, row.sra_id]
            } else if (row.read1 && row.read1 != '' && row.read1 != 'NA') {
                meta.source = 'local'
                def reads = []
                reads << file(row.read1, checkIfExists: true)
                if (row.read2 && row.read2 != '' && row.read2 != 'NA') {
                    reads << file(row.read2, checkIfExists: true)
                } else {
                    meta.single_end = true
                }
                return [meta, reads]
            } else {
                error "Sample ${row.sample_id} has neither SRA ID nor local files!"
            }
        }
        .branch { meta, data ->
            sra: meta.source == 'sra'
            local: meta.source == 'local'
        }
        .set { ch_input }
    
    //
    // Handle SRA downloads (TODO: implement SRA download module)
    //
    ch_sra_reads = ch_input.sra
        .map { meta, sra_id ->
            // For now, assume SRA download is handled externally
            // In full implementation, would call SRA_DOWNLOAD module here
            log.warn "SRA download not yet implemented. Please provide local FASTQ files."
            return [meta, []]
        }
    
    //
    // Combine local and SRA reads
    //
    ch_reads = ch_input.local.mix(ch_sra_reads)
    
    //
    // MODULE: Prepare reference genome indices
    //
    GENOME_PREP (
        ch_genome,
        params.genome_id
    )
    ch_versions = ch_versions.mix(GENOME_PREP.out.versions)
    
    //
    // MODULE: Run FastQC (unless in rapid mode)
    //
    // TODO: Uncomment when FASTQC module is installed
    ch_fastqc_raw_multiqc = Channel.empty()
    // if (!params.rapid_mode && !params.skip_fastqc) {
    //     FASTQC (
    //         ch_reads
    //     )
    //     ch_fastqc_raw_multiqc = FASTQC.out.zip
    //     ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    // }
    
    //
    // Branch samples by assay type
    //
    ch_reads
        .branch { meta, reads ->
            atacseq:    meta.type == 'atacseq'
            chipseq:    meta.type == 'chipseq'
            cutrun:     meta.type == 'cutrun'
            rnaseq:     meta.type == 'rnaseq'
            wgs:        meta.type == 'wgs'
            nanopore:   meta.type == 'nanopore'
            pacbio:     meta.type == 'pacbio'
            ancientdna: meta.type == 'ancientdna'
        }
        .set { ch_by_assay }
    
    //
    // Initialize output channels
    //
    ch_all_bams     = Channel.empty()
    ch_all_peaks    = Channel.empty()
    ch_all_bigwigs  = Channel.empty()
    
    //
    // SUBWORKFLOW: Process ATAC-seq samples
    //
    // TODO: Uncomment when ATACSEQ subworkflow modules are installed
    // ch_atacseq_results = Channel.empty()
    // ch_by_assay.atacseq.ifEmpty([])
    //     .map { meta, reads -> 
    //         log.info "Processing ATAC-seq sample: ${meta.id}"
    //         [meta, reads]
    //     }
    //     .set { ch_atacseq_input }
    // 
    // if (ch_atacseq_input) {
    //     ATACSEQ (
    //         ch_atacseq_input,
    //         GENOME_PREP.out.bowtie2_index.collect().ifEmpty([]),
    //         GENOME_PREP.out.bwa_index.collect().ifEmpty([]),
    //         GENOME_PREP.out.fasta,
    //         GENOME_PREP.out.fai,
    //         params.atacseq.mapper
    //     )
    //     ch_all_bams     = ch_all_bams.mix(ATACSEQ.out.bam)
    //     ch_all_peaks    = ch_all_peaks.mix(ATACSEQ.out.peaks)
    //     ch_all_bigwigs  = ch_all_bigwigs.mix(ATACSEQ.out.bigwig)
    //     ch_versions     = ch_versions.mix(ATACSEQ.out.versions)
    // }
    
    //
    // SUBWORKFLOW: Process ChIP-seq samples
    //
    // TODO: Uncomment when CHIPSEQ subworkflow modules are installed
    // ch_by_assay.chipseq.ifEmpty([])
    //     .map { meta, reads -> 
    //         log.info "Processing ChIP-seq sample: ${meta.id}"
    //         [meta, reads]
    //     }
    //     .set { ch_chipseq_input }
    // 
    // if (ch_chipseq_input) {
    //     CHIPSEQ (
    //         ch_chipseq_input,
    //         GENOME_PREP.out.bowtie2_index.collect().ifEmpty([]),
    //         GENOME_PREP.out.bwa_index.collect().ifEmpty([]),
    //         GENOME_PREP.out.fasta,
    //         GENOME_PREP.out.fai,
    //         params.chipseq.mapper
    //     )
    //     ch_all_bams     = ch_all_bams.mix(CHIPSEQ.out.bam)
    //     ch_all_peaks    = ch_all_peaks.mix(CHIPSEQ.out.peaks)
    //     ch_all_bigwigs  = ch_all_bigwigs.mix(CHIPSEQ.out.bigwig)
    //     ch_versions     = ch_versions.mix(CHIPSEQ.out.versions)
    // }
    
    //
    // SUBWORKFLOW: Process CUT&RUN samples
    //
    // TODO: Uncomment when CUTRUN subworkflow modules are installed
    // ch_by_assay.cutrun.ifEmpty([])
    //     .map { meta, reads -> 
    //         log.info "Processing CUT&RUN sample: ${meta.id}"
    //         [meta, reads]
    //     }
    //     .set { ch_cutrun_input }
    // 
    // if (ch_cutrun_input) {
    //     CUTRUN (
    //         ch_cutrun_input,
    //         GENOME_PREP.out.bowtie2_index.collect().ifEmpty([]),
    //         GENOME_PREP.out.bwa_index.collect().ifEmpty([]),
    //         GENOME_PREP.out.fasta,
    //         GENOME_PREP.out.fai,
    //         params.cutrun.mapper
    //     )
    //     ch_all_bams     = ch_all_bams.mix(CUTRUN.out.bam)
    //     ch_all_peaks    = ch_all_peaks.mix(CUTRUN.out.peaks)
    //     ch_all_bigwigs  = ch_all_bigwigs.mix(CUTRUN.out.bigwig)
    //     ch_versions     = ch_versions.mix(CUTRUN.out.versions)
    // }
    
    //
    // SUBWORKFLOW: Process RNA-seq samples
    //
    ch_by_assay.rnaseq.ifEmpty([])
        .map { meta, reads -> 
            log.info "Processing RNA-seq sample: ${meta.id}"
            [meta, reads]
        }
        .set { ch_rnaseq_input }
    
    if (ch_rnaseq_input) {
        RNASEQ (
            ch_rnaseq_input,
            GENOME_PREP.out.hisat2_index.collect().ifEmpty([]),
            GENOME_PREP.out.fasta,
            GENOME_PREP.out.fai,
            GENOME_PREP.out.sizes
        )
        ch_all_bams     = ch_all_bams.mix(RNASEQ.out.bam)
        ch_all_bigwigs  = ch_all_bigwigs.mix(RNASEQ.out.bigwig)
        ch_versions     = ch_versions.mix(RNASEQ.out.versions)
    }
    
    // TODO: Add other assay type subworkflows (WGS, Nanopore, PacBio, Ancient DNA, etc.)
    
    //
    // MODULE: Generate UCSC track hub
    //
    ch_all_bigwigs
        .map { meta, bigwig -> bigwig }
        .collect()
        .set { ch_bigwigs_for_hub }
    
    ch_all_peaks
        .map { meta, peaks -> peaks }
        .collect()
        .ifEmpty([])
        .set { ch_peaks_for_hub }
    
    // Collect all sample metadata
    ch_reads
        .map { meta, reads -> meta }
        .collect()
        .set { ch_sample_metadata }
    
    UCSC_HUB (
        params.genome_id,
        ch_bigwigs_for_hub,
        ch_peaks_for_hub,
        ch_sample_metadata,
        GENOME_PREP.out.sizes
    )
    ch_versions = ch_versions.mix(UCSC_HUB.out.versions)
    
    //
    // MODULE: Pipeline reporting
    //
    // TODO: Uncomment when modules are installed
    // CUSTOM_DUMPSOFTWAREVERSIONS (
    //     ch_versions.unique().collectFile(name: 'collated_versions.yml')
    // )

    //
    // MODULE: MultiQC
    //
    // TODO: Uncomment when MULTIQC module is installed
    // if (!params.rapid_mode && !params.skip_multiqc) {
    //     workflow_summary    = WorkflowUseq2tracks.paramsSummaryMultiqc(workflow, summary_params)
    //     ch_workflow_summary = Channel.value(workflow_summary)

    //     ch_multiqc_files = Channel.empty()
    //     ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    //     ch_multiqc_files = ch_multiqc_files.mix(ch_fastqc_raw_multiqc.collect{it[1]}.ifEmpty([]))
    //     ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    //     MULTIQC (
    //         ch_multiqc_files.collect()
    //     )
    //     multiqc_report = MULTIQC.out.report.toList()
    //     ch_versions    = ch_versions.mix(MULTIQC.out.versions)
    // }
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
