process HISAT2_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::hisat2=2.2.1 bioconda::samtools=1.16.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' :
        'biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' }"

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    script:
    def prefix = "${meta.id}"
    def read_input = meta.single_end ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    
    """
    INDEX=`find -L ./ -name "*.1.ht2" | sed 's/.1.ht2//'`
    
    hisat2 \\
        -x \$INDEX \\
        $read_input \\
        --threads ${task.cpus} \\
        --summary-file ${prefix}.hisat2_summary.txt \\
        | samtools view -@ ${task.cpus} -bS -o ${prefix}.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: \$(hisat2 --version 2>&1 | grep -o 'version [0-9.]*' | cut -d' ' -f2)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
