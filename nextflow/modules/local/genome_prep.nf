process GENOME_PREP {
    tag "$genome_id"
    label 'process_high'
    
    conda (params.enable_conda ? "bioconda::bowtie2=2.4.5 bioconda::bwa=0.7.17 bioconda::samtools=1.16.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:4d8b9f5f5c8b5f8e6c8c8d8e6c8c8c8c8c8c8c8c' :
        'quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:4d8b9f5f5c8b5f8e6c8c8d8e6c8c8c8c8c8c8c8c' }"

    input:
    path genome
    val genome_id

    output:
    path "genome.fa"               , emit: fasta
    path "genome.fa.fai"           , emit: fai
    path "genome.sizes"            , emit: sizes
    path "bowtie2/*"               , emit: bowtie2_index, optional: true
    path "bwa_mem2/*"              , emit: bwa_index, optional: true
    path "versions.yml"            , emit: versions

    script:
    def build_bowtie2 = params.atacseq.mapper == 'bowtie2' || 
                        params.chipseq.mapper == 'bowtie2' || 
                        params.cutrun.mapper == 'bowtie2' ||
                        params.wgs.mapper == 'bowtie2'
    def build_bwa = params.atacseq.mapper == 'bwa_mem2' || 
                    params.chipseq.mapper == 'bwa_mem2' ||
                    params.cutrun.mapper == 'bwa_mem2' ||
                    params.wgs.mapper == 'bwa_mem2' ||
                    params.ancientdna.mapper == 'bwa_aln' ||
                    params.ancientdna.mapper == 'bwa_mem'
    
    """
    # Copy and prepare genome
    cp $genome genome.fa
    
    # Create FASTA index
    samtools faidx genome.fa
    
    # Create chromosome sizes file
    cut -f1,2 genome.fa.fai > genome.sizes
    
    # Build Bowtie2 index if needed
    if [ "$build_bowtie2" = "true" ]; then
        mkdir -p bowtie2
        bowtie2-build --threads ${task.cpus} genome.fa bowtie2/genome
    fi
    
    # Build BWA index if needed
    if [ "$build_bwa" = "true" ]; then
        mkdir -p bwa_mem2
        cp genome.fa bwa_mem2/
        bwa index bwa_mem2/genome.fa
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*version //; s/ .*\$//')
        bwa: \$(echo \$(bwa 2>&1) | grep -i version | sed 's/Version: //g')
    END_VERSIONS
    """
}
