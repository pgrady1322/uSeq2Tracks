# ============================================================
# MultiQC aggregation rules
# ============================================================

rule multiqc:
    input:
        fastqc=expand(f"{GENOME_OUTDIR}/qc/fastqc/{{sample}}_fastqc.zip", sample=SAMPLES.keys()),
        # Add markdup metrics only if markdup is enabled for each assay type
        chipseq_markdup=expand(f"{GENOME_OUTDIR}/chipseq/qc/{{sample}}.markdup.metrics", sample=CHIPSEQ_SAMPLES) if CHIPSEQ_SAMPLES and config.get("chipseq", {}).get("markdup", False) else [],
        atacseq_markdup=expand(f"{GENOME_OUTDIR}/atacseq/qc/{{sample}}.markdup.metrics", sample=ATACSEQ_SAMPLES) if ATACSEQ_SAMPLES and config.get("atacseq", {}).get("markdup", False) else [],
        cutrun_markdup=expand(f"{GENOME_OUTDIR}/cutrun/qc/{{sample}}.markdup.metrics", sample=CUTRUN_SAMPLES) if CUTRUN_SAMPLES and config.get("cutrun", {}).get("markdup", False) else [],
        wgs_markdup=expand(f"{GENOME_OUTDIR}/wgs/qc/{{sample}}.markdup.metrics", sample=WGS_SAMPLES) if WGS_SAMPLES and config.get("wgs", {}).get("markdup", True) else [],
        adna_markdup=expand(f"{GENOME_OUTDIR}/ancientdna/qc/{{sample}}.dedup.metrics", sample=ADNA_SAMPLES) if ADNA_SAMPLES and config.get("ancientdna", {}).get("markdup", True) else []
    output:
        f"{GENOME_OUTDIR}/qc/multiqc_report.html"
    params:
        outdir=f"{GENOME_OUTDIR}/qc",
        search_dirs=f"{GENOME_OUTDIR}"
    shell:
        """
        multiqc --outdir {params.outdir} \
                --filename multiqc_report.html \
                --force \
                {params.search_dirs}
        """
