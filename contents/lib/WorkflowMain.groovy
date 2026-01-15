/*
========================================================================================
    Workflow Initialization and Validation
========================================================================================
*/

class WorkflowMain {

    //
    // Print help message
    //
    public static String help(workflow, params, log) {
        def command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity> --samplesheet <file> --genome <file> --genome_id <string>"
        def help_string = ''
        help_string += NfcoreTemplate.logo(workflow, params.monochrome_logs)
        help_string += NfcoreSchema.paramsHelp(workflow, params, command)
        help_string += '\n' + WorkflowMain.citation(workflow) + '\n'
        help_string += NfcoreTemplate.dashedLine(params.monochrome_logs)
        return help_string
    }

    //
    // Citation string
    //
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
            "* The pipeline\n" +
            "  https://github.com/${workflow.manifest.name}\n\n" +
            "* The nf-core framework\n" +
            "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
            "* Software dependencies\n" +
            "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
    }

    //
    // Validate parameters and print summary
    //
    public static void initialise(workflow, params, log) {
        // Print help to screen if required
        if (params.help) {
            log.info help(workflow, params, log)
            System.exit(0)
        }

        // Print workflow version and exit on --version
        if (params.version) {
            String version_string = ""
            version_string += "${workflow.manifest.name} ${workflow.manifest.version}\n"
            log.info version_string
            System.exit(0)
        }

        // Check that required parameters are provided
        if (!params.samplesheet) {
            log.error "Please provide an input samplesheet with --samplesheet"
            System.exit(1)
        }
        if (!params.genome) {
            log.error "Please provide a reference genome with --genome"
            System.exit(1)
        }
        if (!params.genome_id) {
            log.error "Please provide a genome identifier with --genome_id (e.g., 'galGal6')"
            System.exit(1)
        }
    }
}
