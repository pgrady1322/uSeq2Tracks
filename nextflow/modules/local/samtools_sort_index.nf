process SAMTOOLS_SORT_INDEX {
    tag "$meta.id"
    label 'process_medium'
    
    publishDir "${params.outdir}/${params.genome_id}/rnaseq/bam", mode: params.publish_dir_mode, pattern: "*.{bam,bai}"

    conda (params.enable_conda ? "bioconda::samtools=1.16.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
        'biocontainers/samtools:1.16.1--h6899075_1' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.sorted.bam"), emit: bam
    tuple val(meta), path("*.sorted.bam.bai"), emit: bai
    path "versions.yml"                  , emit: versions

    script:
    def prefix = "${meta.id}"
    
    """
    samtools sort \\
        -@ ${task.cpus} \\
        -o ${prefix}.sorted.bam \\
        $bam

    samtools index \\
        -@ ${task.cpus} \\
        ${prefix}.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
