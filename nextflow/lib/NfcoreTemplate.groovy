/*
========================================================================================
    nf-core template functions
========================================================================================
*/

class NfcoreTemplate {

    //
    // Generate a workflow summary for the terminal
    //
    public static void summary(workflow, params, log) {
        Map colors = logColours(params.monochrome_logs)
        log.info "${colors.cyan}${workflow.manifest.name} ${colors.green}v${workflow.manifest.version}${colors.reset}"
        log.info "${colors.purple}=======================================${colors.reset}"
        
        def summary_params = [:]
        
        // Core parameters
        if (params.samplesheet)   summary_params['Samplesheet'] = params.samplesheet
        if (params.genome)        summary_params['Genome'] = params.genome
        if (params.genome_id)     summary_params['Genome ID'] = params.genome_id
        if (params.gtf)           summary_params['GTF'] = params.gtf
        if (params.outdir)        summary_params['Output dir'] = params.outdir
        
        // Pipeline mode
        summary_params['Rapid mode'] = params.rapid_mode ? 'Yes' : 'No'
        
        // Resources
        summary_params['Max CPUs'] = params.max_cpus
        summary_params['Max Memory'] = params.max_memory
        summary_params['Max Time'] = params.max_time
        
        log.info summary_params.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
        log.info "${colors.purple}=======================================${colors.reset}"
    }

    //
    // Generate a completion email
    //
    public static void email(workflow, params, summary_params, projectDir, log, multiqc_report=[]) {
        // TODO: Implement email functionality
        log.info "Email notification not yet implemented"
    }

    //
    // ASCII art logo
    //
    public static String logo(workflow, monochrome_logs) {
        Map colors = logColours(monochrome_logs)
        String.format(
            """\n
            ${colors.cyan}
            ┬ ┬┌─┐┌─┐┌─┐┌─┐┌─┐╔╦╗┬─┐┌─┐┌─┐┬┌─┌─┐
            │ │└─┐├┤ │─┼┤ ┴┐ ║ ├┬┘├─┤│  ├┴┐└─┐
            └─┘└─┘└─┘└─┘└─┘└─┘ ╩ ┴└─┴ ┴└─┘┴ ┴└─┘
            ${colors.reset}
            """.stripIndent()
        )
    }

    //
    // Dashed line
    //
    public static String dashedLine(monochrome_logs) {
        Map colors = logColours(monochrome_logs)
        return "-${colors.dim}----------------------------------------------------${colors.reset}-"
    }

    //
    // ANSI colour codes for terminal output
    //
    private static Map logColours(Boolean monochrome_logs) {
        Map colorcodes = [:]
        
        // Reset / Meta
        colorcodes['reset']      = monochrome_logs ? '' : "\033[0m"
        colorcodes['bold']       = monochrome_logs ? '' : "\033[1m"
        colorcodes['dim']        = monochrome_logs ? '' : "\033[2m"
        colorcodes['underlined'] = monochrome_logs ? '' : "\033[4m"
        colorcodes['blink']      = monochrome_logs ? '' : "\033[5m"
        colorcodes['reverse']    = monochrome_logs ? '' : "\033[7m"
        colorcodes['hidden']     = monochrome_logs ? '' : "\033[8m"
        
        // Regular Colors
        colorcodes['black']      = monochrome_logs ? '' : "\033[0;30m"
        colorcodes['red']        = monochrome_logs ? '' : "\033[0;31m"
        colorcodes['green']      = monochrome_logs ? '' : "\033[0;32m"
        colorcodes['yellow']     = monochrome_logs ? '' : "\033[0;33m"
        colorcodes['blue']       = monochrome_logs ? '' : "\033[0;34m"
        colorcodes['purple']     = monochrome_logs ? '' : "\033[0;35m"
        colorcodes['cyan']       = monochrome_logs ? '' : "\033[0;36m"
        colorcodes['white']      = monochrome_logs ? '' : "\033[0;37m"
        
        return colorcodes
    }
}
