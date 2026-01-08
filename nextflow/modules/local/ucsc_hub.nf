process UCSC_HUB {
    tag "$genome_id"
    label 'process_low'
    publishDir "${params.outdir}/${genome_id}/ucsc", mode: params.publish_dir_mode

    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    val genome_id
    path bigwigs
    path peaks
    val sample_metadata
    path genome_sizes

    output:
    path "hub.txt"        , emit: hub
    path "genomes.txt"    , emit: genomes
    path "trackDb.txt"    , emit: trackdb
    path "versions.yml"   , emit: versions

    script:
    def hub_name        = params.ucsc.hub_name
    def hub_short_label = params.ucsc.hub_short_label
    def hub_long_label  = params.ucsc.hub_long_label
    def genome_name     = params.ucsc.genome_name
    def hub_email       = params.ucsc.hub_email
    
    """
    generate_ucsc_hub.py \\
        --genome-id ${genome_id} \\
        --hub-name "${hub_name}" \\
        --hub-short-label "${hub_short_label}" \\
        --hub-long-label "${hub_long_label}" \\
        --genome-name "${genome_name}" \\
        --hub-email "${hub_email}" \\
        --bigwigs ${bigwigs} \\
        --peaks ${peaks} \\
        --output-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
