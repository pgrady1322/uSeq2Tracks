/*
========================================================================================
    nf-core Schema functions
========================================================================================
*/

class NfcoreSchema {

    //
    // Generate a parameter summary map
    //
    public static LinkedHashMap paramsSummaryMap(workflow, params) {
        def summary_params = [:]
        
        // Input/Output
        def io_params = [:]
        io_params['samplesheet'] = params.samplesheet
        io_params['genome'] = params.genome
        io_params['genome_id'] = params.genome_id
        io_params['gtf'] = params.gtf
        io_params['outdir'] = params.outdir
        summary_params['Input/Output'] = io_params
        
        // Pipeline Mode
        def mode_params = [:]
        mode_params['rapid_mode'] = params.rapid_mode
        mode_params['skip_fastqc'] = params.skip_fastqc
        mode_params['skip_multiqc'] = params.skip_multiqc
        summary_params['Pipeline Mode'] = mode_params
        
        // ATAC-seq
        def atacseq_params = [:]
        atacseq_params['mapper'] = params.atacseq.mapper
        atacseq_params['markdup'] = params.atacseq.markdup
        atacseq_params['bw_norm'] = params.atacseq.bw_norm
        summary_params['ATAC-seq'] = atacseq_params
        
        // ChIP-seq
        def chipseq_params = [:]
        chipseq_params['mapper'] = params.chipseq.mapper
        chipseq_params['control_tag'] = params.chipseq.control_tag
        chipseq_params['markdup'] = params.chipseq.markdup
        summary_params['ChIP-seq'] = chipseq_params
        
        // Resources
        def resource_params = [:]
        resource_params['max_cpus'] = params.max_cpus
        resource_params['max_memory'] = params.max_memory
        resource_params['max_time'] = params.max_time
        summary_params['Max Resources'] = resource_params
        
        return summary_params
    }

    //
    // Generate help message
    //
    public static String paramsHelp(workflow, params, command) {
        def help_string = """
        Usage:
            ${command}

        Mandatory arguments:
            --samplesheet [file]      Path to comma-separated file containing information about the samples
            --genome [file]           Path to reference genome FASTA file
            --genome_id [string]      Genome identifier (e.g., 'galGal6')

        Optional arguments:
            --outdir [file]           Output directory (default: ./results)
            --gtf [file]              Path to GTF annotation file (required for RNA-seq)
            
        Pipeline modes:
            --rapid_mode [bool]       Skip QC and generate only essential tracks
            --skip_fastqc [bool]      Skip FastQC
            --skip_multiqc [bool]     Skip MultiQC
            
        ATAC-seq options:
            --atacseq.mapper [string]     Alignment tool (bowtie2 | bwa_mem2)
            --atacseq.markdup [bool]      Mark duplicates
            --atacseq.bw_norm [string]    BigWig normalization (CPM | RPKM | BPM | RPGC)
            
        ChIP-seq options:
            --chipseq.mapper [string]     Alignment tool (bowtie2 | bwa_mem2)
            --chipseq.control_tag [string] Control sample identifier
            --chipseq.markdup [bool]      Mark duplicates
            
        Resource limits:
            --max_cpus [int]          Maximum number of CPUs
            --max_memory [string]     Maximum memory
            --max_time [string]       Maximum time

        Other options:
            --help                    Display this help message
        """.stripIndent()
        
        return help_string
    }
}
