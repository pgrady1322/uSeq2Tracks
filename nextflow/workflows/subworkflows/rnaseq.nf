/*
========================================================================================
    RNA-seq Analysis Subworkflow
========================================================================================
    Processes RNA-seq samples: alignment and BigWig generation
*/

// Local modules
include { HISAT2_ALIGN           } from '../../modules/local/hisat2_align'
include { SAMTOOLS_SORT_INDEX    } from '../../modules/local/samtools_sort_index'
include { BIGWIG_GENERATE        } from '../../modules/local/bigwig_generate'

workflow RNASEQ {
    take:
    ch_reads            // channel: [meta, [reads]]
    ch_hisat2_index     // channel: /path/to/hisat2/index/*
    ch_fasta            // channel: /path/to/genome.fa
    ch_fai              // channel: /path/to/genome.fa.fai
    ch_sizes            // channel: /path/to/genome.sizes

    main:
    ch_versions = Channel.empty()

    //
    // Align reads with HISAT2
    //
    HISAT2_ALIGN (
        ch_reads,
        ch_hisat2_index
    )
    ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions)

    //
    // Sort and index BAM files
    //
    SAMTOOLS_SORT_INDEX (
        HISAT2_ALIGN.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_INDEX.out.versions)

    //
    // Generate BigWig coverage tracks
    //
    BIGWIG_GENERATE (
        SAMTOOLS_SORT_INDEX.out.bam,
        SAMTOOLS_SORT_INDEX.out.bai,
        ch_sizes
    )
    ch_versions = ch_versions.mix(BIGWIG_GENERATE.out.versions)

    emit:
    bam      = SAMTOOLS_SORT_INDEX.out.bam      // channel: [ meta, bam ]
    bai      = SAMTOOLS_SORT_INDEX.out.bai      // channel: [ meta, bai ]
    bigwig   = BIGWIG_GENERATE.out.bigwig       // channel: [ meta, bigwig ]
    versions = ch_versions                       // channel: [ versions.yml ]
}
