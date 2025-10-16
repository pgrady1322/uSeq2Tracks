# ============================================================
# Quality control rules using fastp (replacing FastQC)
# ============================================================

rule fastqc:
    input:
        get_fastq_for_fastqc
    output:
        html=f"{GENOME_OUTDIR}/qc/fastqc/{{sample}}_fastqc.html",
        zip=f"{GENOME_OUTDIR}/qc/fastqc/{{sample}}_fastqc.zip"
    params:
        outdir=f"{GENOME_OUTDIR}/qc/fastqc",
        json=f"{GENOME_OUTDIR}/qc/fastqc/{{sample}}_fastp.json"
    threads: 2
    shell:
        """
        mkdir -p {params.outdir}
        
        # Determine if input is single-end or paired-end
        input_files=({input})
        if [ "${{#input_files[@]}}" -eq 2 ]; then
            # Paired-end mode
            fastp -i ${{input_files[0]}} -I ${{input_files[1]}} \
                  --html {output.html} \
                  --json {params.json} \
                  --thread {threads} \
                  --dont_eval_duplication \
                  --disable_adapter_trimming \
                  --disable_trim_poly_g \
                  --disable_quality_filtering \
                  --disable_length_filtering
        else
            # Single-end mode
            fastp -i ${{input_files[0]}} \
                  --html {output.html} \
                  --json {params.json} \
                  --thread {threads} \
                  --dont_eval_duplication \
                  --disable_adapter_trimming \
                  --disable_trim_poly_g \
                  --disable_quality_filtering \
                  --disable_length_filtering
        fi
        
        # Create a dummy zip file to match FastQC output expectations
        cd {params.outdir}
        echo "{{" > {wildcards.sample}_fastqc_data.txt
        echo '  "summary": "fastp quality analysis",' >> {wildcards.sample}_fastqc_data.txt
        echo '  "sample": "{wildcards.sample}",' >> {wildcards.sample}_fastqc_data.txt
        echo '  "tool": "fastp"' >> {wildcards.sample}_fastqc_data.txt
        echo "}}" >> {wildcards.sample}_fastqc_data.txt
        zip {wildcards.sample}_fastqc.zip {wildcards.sample}_fastqc_data.txt
        rm {wildcards.sample}_fastqc_data.txt
        """
