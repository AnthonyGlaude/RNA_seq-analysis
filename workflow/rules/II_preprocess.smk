# RNA-Seq quality control followed by alignment with STAR 
rule fastqc:
    """Assess the FASTQ quality using FastQC BEFORE TRIMMING"""
    input:
        fq1 = os.path.join(config["path"]["RNAseq_dir"], "{id}", "{id}_1.fastq.gz"),
        fq2 = os.path.join(config["path"]["RNAseq_dir"], "{id}", "{id}_2.fastq.gz")

    output:
        qc_fq1_out = "data/qc/{id}/{id}_1_fastqc.html",
        qc_fq2_out = "data/qc/{id}/{id}_2_fastqc.html"

    params:
        out_dir = "data/qc/{id}"
    log:
        "logs/{id}/fastqc.log"
    threads: 8
    conda:
        "../envs/fastqc.yml" 
    shell:
        "mkdir -p {params.out_dir} && "
        "fastqc --outdir {params.out_dir} --format fastq --threads {threads} {input.fq1} {input.fq2} &> {log} "


# Cutadapt trimming (calling by trim_galore)
rule trim_reads:
    input:
        fq1 = os.path.join(config["path"]["RNAseq_dir"], "{id}",  "{id}_1.fastq.gz"),
        fq2 = os.path.join(config["path"]["RNAseq_dir"], "{id}",  "{id}_2.fastq.gz")
    output:
        gal_trim1 = "data/trim_galore/{id}/{id}_1_val_1.fq.gz", 
        gal_trim2 = "data/trim_galore/{id}/{id}_2_val_2.fq.gz",
        unpaired1 = "data/trim_galore/{id}/{id}_1_unpaired_1.fq.gz", 
        unpaired2 = "data/trim_galore/{id}/{id}_2_unpaired_2.fq.gz"
    params:
        out_dir = "data/trim_galore/{id}"
    threads:
        6
    conda:
        "../envs/trim_galore.yml"
    log:
        "logs/{id}/trim.log"
    shell:
        """
        mkdir -p {params.out_dir} &&\
        trim_galore --paired \
        --retain_unpaired \
        --cores {threads} \
        --gzip \
        --output_dir {params.out_dir} \
        {input.fq1} {input.fq2} \
        &> {log}
        """

rule qc_fastq:
    """ Assess the FASTQ quality using FastQC AFTER TRIMMING"""
    input:
        trimm_fq1 = rules.trim_reads.output.gal_trim1,
        trimm_fq2 = rules.trim_reads.output.gal_trim2,
        unpaired1 = rules.trim_reads.output.unpaired1,
        unpaired2 = rules.trim_reads.output.unpaired2
    output:
        qc_trimm_fq1_out = "data/qc_after_trim/{id}/{id}_1_val_1_fastqc.html",
        qc_trimm_fq2_out = "data/qc_after_trim/{id}/{id}_2_val_2_fastqc.html"
    params:
        out_dir = "data/qc_after_trim/{id}"
    log:
        "logs/{id}/FASTQC2.log"
    threads:
        8
    conda:
        "../envs/fastqc.yml"
    shell:
        """
        mkdir -p {params.out_dir} &&
        fastqc \
            --outdir {params.out_dir} \
            --format fastq \
            --threads {threads} \
            {input.trimm_fq1} {input.trimm_fq2} {input.unpaired1} {input.unpaired2}\
            &> {log}
        """
