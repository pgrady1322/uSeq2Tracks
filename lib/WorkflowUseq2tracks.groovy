/*
========================================================================================
    Workflow-specific functions for uSeq2Tracks
========================================================================================
*/

class WorkflowUseq2tracks {

    //
    // Check and validate pipeline parameters
    //
    public static void initialise(params, log) {
        
        // Check assay-specific parameters
        if (params.atacseq.mapper !in ['bowtie2', 'bwa_mem2']) {
            log.warn "Invalid ATAC-seq mapper '${params.atacseq.mapper}'. Using 'bowtie2'."
        }
        
        if (params.chipseq.mapper !in ['bowtie2', 'bwa_mem2']) {
            log.warn "Invalid ChIP-seq mapper '${params.chipseq.mapper}'. Using 'bowtie2'."
        }
        
        // Validate normalization methods
        def valid_norms = ['CPM', 'RPKM', 'BPM', 'RPGC', 'None']
        ['atacseq', 'chipseq', 'cutrun', 'rnaseq', 'wgs', 'nanopore', 'pacbio', 'ancientdna'].each { assay ->
            def norm = params."${assay}".bw_norm
            if (norm && !(norm in valid_norms)) {
                log.warn "Invalid BigWig normalization '${norm}' for ${assay}. Valid options: ${valid_norms.join(', ')}"
            }
        }
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)
            if (group_params) {
                summary_section += "    <p style=\\"font-size:110%\\"><b>$group</b></p>\\n"
                summary_section += "    <dl class=\\"dl-horizontal\\">\\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\\"color:#999999;\\">N/A</a>'}</samp></dd>\\n"
                }
                summary_section += "    </dl>\\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\\n"
        yaml_file_text        += "plot_type: 'html'\\n"
        yaml_file_text        += "data: |\\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }
}
