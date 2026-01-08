/*
========================================================================================
    ChIP-seq Analysis Subworkflow
========================================================================================
    Processes ChIP-seq samples: alignment, peak calling with controls, BigWig generation
*/

include { BOWTIE2_ALIGN           } from '../../modules/nf-core/bowtie2/align/main'
include { BWA_MEM                 } from '../../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_SORT           } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX          } from '../../modules/nf-core/samtools/index/main'
include { PICARD_MARKDUPLICATES   } from '../../modules/nf-core/picard/markduplicates/main'
include { MACS3_CALLPEAK          } from '../../modules/nf-core/macs3/callpeak/main'
include { DEEPTOOLS_BAMCOVERAGE   } from '../../modules/nf-core/deeptools/bamcoverage/main'

workflow CHIPSEQ {
    take:
    ch_reads            // channel: [meta, [reads]]
    ch_bowtie2_index    // channel: /path/to/bowtie2/index/*
    ch_bwa_index        // channel: /path/to/bwa/index/*
    ch_genome           // channel: /path/to/genome.fa
    ch_genome_fai       // channel: /path/to/genome.fa.fai
    mapper              // val: 'bowtie2' or 'bwa_mem2'
    
    main:
    ch_versions = Channel.empty()
    
    //
    // Separate controls from experimental samples
    //
    def control_tag = params.chipseq.control_tag.toLowerCase()
    ch_reads
        .branch { meta, reads ->
            control: meta.condition.toLowerCase().contains(control_tag) ||
                     meta.replicate_group.toLowerCase().contains(control_tag)
            experimental: true
        }
        .set { ch_branched }
    
    //
    // Align reads with selected mapper
    //
    if (mapper == 'bowtie2') {
        ch_index = ch_bowtie2_index
            .map { files -> 
                def index_file = files[0]
                def index_dir = index_file.getParent()
                def index_base = index_file.getName().toString().replaceAll(/\\.\\d+\\.bt2$/, '')
                [index_dir, index_base]
            }
            .first()
        
        BOWTIE2_ALIGN (
            ch_reads,
            ch_index,
            false,
            false
        )
        ch_bam = BOWTIE2_ALIGN.out.bam
        ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())
    } else {
        BWA_MEM (
            ch_reads,
            ch_bwa_index.first(),
            false
        )
        ch_bam = BWA_MEM.out.bam
        ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())
    }
    
    //
    // Sort BAM files
    //
    SAMTOOLS_SORT (
        ch_bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())
    
    //
    // Index sorted BAM files
    //
    SAMTOOLS_INDEX (
        SAMTOOLS_SORT.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())
    
    //
    // Optional: Mark duplicates
    //
    if (params.chipseq.markdup) {
        PICARD_MARKDUPLICATES (
            SAMTOOLS_SORT.out.bam,
            ch_genome,
            ch_genome_fai
        )
        ch_final_bam = PICARD_MARKDUPLICATES.out.bam
        ch_final_bai = PICARD_MARKDUPLICATES.out.bai
        ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.first())
    } else {
        ch_final_bam = SAMTOOLS_SORT.out.bam
        ch_final_bai = SAMTOOLS_INDEX.out.bai
    }
    
    //
    // Separate control and experimental BAMs
    //
    ch_final_bam
        .branch { meta, bam ->
            control: meta.condition.toLowerCase().contains(control_tag) ||
                     meta.replicate_group.toLowerCase().contains(control_tag)
            experimental: true
        }
        .set { ch_bam_branched }
    
    //
    // Match experimental samples with controls for peak calling
    //
    ch_bam_branched.experimental
        .map { meta, bam ->
            // Match with control based on experiment group
            def control_key = meta.experiment_group
            [control_key, meta, bam]
        }
        .set { ch_exp_with_key }
    
    ch_bam_branched.control
        .map { meta, bam ->
            def control_key = meta.experiment_group
            [control_key, bam]
        }
        .set { ch_ctrl_with_key }
    
    ch_exp_with_key
        .combine(ch_ctrl_with_key, by: 0)
        .map { key, meta, exp_bam, ctrl_bam ->
            [meta, exp_bam, ctrl_bam]
        }
        .set { ch_for_macs }
    
    //
    // Call peaks with MACS3 (with control)
    //
    MACS3_CALLPEAK (
        ch_for_macs.map { meta, exp, ctrl -> [meta, exp] },
        ch_for_macs.map { meta, exp, ctrl -> ctrl }
    )
    ch_versions = ch_versions.mix(MACS3_CALLPEAK.out.versions.first())
    
    //
    // Generate BigWig tracks for all samples
    //
    ch_final_bam
        .join(ch_final_bai, by: 0)
        .set { ch_bam_bai }
    
    DEEPTOOLS_BAMCOVERAGE (
        ch_bam_bai,
        ch_genome_fai
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE.out.versions.first())
    
    emit:
    bam      = ch_final_bam
    bai      = ch_final_bai
    peaks    = MACS3_CALLPEAK.out.peak
    bigwig   = DEEPTOOLS_BAMCOVERAGE.out.bigwig
    versions = ch_versions
}
