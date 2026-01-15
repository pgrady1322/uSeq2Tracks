process BIGWIG_GENERATE {
    tag "$meta.id"
    label 'process_medium'
    
    publishDir "${params.outdir}/${params.genome_id}/rnaseq/bigwig", mode: params.publish_dir_mode, pattern: "*.bw"

    conda (params.enable_conda ? "bioconda::deeptools=3.5.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.1--py_0' :
        'biocontainers/deeptools:3.5.1--py_0' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta), path(bai)
    path sizes

    output:
    tuple val(meta), path("*.bigWig"), emit: bigwig
    path "versions.yml"              , emit: versions

    script:
    def prefix = "${meta.id}"
    def norm_method = params.rnaseq?.bw_norm ?: 'CPM'
    
    """
    bamCoverage \\
        -b $bam \\
        -o ${prefix}.bigWig \\
        --binSize 10 \\
        --normalizeUsing $norm_method \\
        --numberOfProcessors ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(bamCoverage --version | sed -e "s/bamCoverage //g")
    END_VERSIONS
    """
}
