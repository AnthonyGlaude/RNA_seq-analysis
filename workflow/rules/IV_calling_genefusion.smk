rule call_fusion_gene:
    input:
        bam = rules.star_chimeric.output.bam,
    output:
        fusions   = "results/fusion_gene/{id}/fusions.tsv",
        discarded = "results/fusion_gene/{id}/discarded.tsv"
    conda:
        "../envs/arriba.yml"
    params:
        genome    = rules.download_human_genome.output.genome,
        gtf       = rules.download_human_gtf.output.gtf,
        blacklist = rules.prepare_arriba_database.output.blacklist,
        known     = rules.prepare_arriba_database.output.known,
        domains   = rules.prepare_arriba_database.output.domains
    log:
        "logs/arriba/{id}.log"
    threads: 8
    shell:
        r"""
        mkdir -p $(dirname {output.fusions})

        arriba \
          -x {input.bam} \
          -a {params.genome} \
          -g {params.gtf} \
          -b {params.blacklist} \
          -k {params.known} \
          -p {params.domains} \
          -o {output.fusions} \
          -O {output.discarded} \
          -@ {threads} \
          2> {log}
        """
rule index_star_chimeric_bam:
    input:
        bam = rules.star_chimeric.output.bam
    output:
        bai = "results/STAR_fusion/{id}/Aligned.sortedByCoord.out.bam.bai"
    shell:
        """
        module load samtools
        samtools index {input.bam}
        """

rule draw_fusions:
    input:
        fusions = rules.call_fusion_gene.output.fusions,
        bam     = rules.star_chimeric.output.bam,
        bai     = rules.index_star_chimeric_bam.output.bai,
    output:
        pdf = "results/fusion_gene/{id}/fusions.pdf"
    params:
        script     = "data/references/arriba/arriba_v2.5.0/draw_fusions.R",
        annotation = "data/references/gtf/Homo_sapiens.GRCh38.113.gtf",
        cytobands  = "data/references/arriba/arriba_v2.5.0/database/cytobands_hg38_GRCh38_v2.5.0.tsv",
        domains    = "data/references/arriba/protein_domains_hg38_GRCh38_v2.5.0.gff3"
    log:
        "logs/arriba_draw/{id}.log"
    conda:
        "../envs/arriba.yml"
    shell:
        r"""
        mkdir -p $(dirname {output.pdf})

        Rscript {params.script} \
          --fusions={input.fusions} \
          --alignments={input.bam} \
          --output={output.pdf} \
          --annotation={params.annotation} \
          --cytobands={params.cytobands} \
          --proteinDomains={params.domains} \
          2> {log}
        """

#rule gene_fusion_transcripts:
#    input:
#        fusion_gene = rules.call_fusion_gene.output.fusions
#    output:
#        fusion_transcriptome = "results/fusion_gene/{id}/fusion_transcriptome.fa"
#    conda:
#        "../envs/python.yml"
#    script: "../scripts/IV_creation_fusion_transcriptome.py"

#ANNOFUSE ?
#### STAR-FUSION ###

rule star_fusion:
    input:
        fq1 = rules.trim_reads.output.gal_trim1,
        fq2 = rules.trim_reads.output.gal_trim2
    output:
        predictions = "results/star_fusion/{id}/star-fusion.fusion_predictions.tsv",
        abridged    = "results/star_fusion/{id}/star-fusion.fusion_predictions.abridged.tsv"
    params:
        genome_lib_dir = "data/references/ctat/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play/ctat_genome_lib_build_dir",
        outdir = "results/star_fusion/{id}"
    threads: 32
    conda:
        "../envs/starfusion.yml"
    log:
        "logs/star_fusion/{id}.log"
    shell:
        r"""
        set -euo pipefail

        mkdir -p {params.outdir} logs/star_fusion

        STAR-Fusion \
            --genome_lib_dir {params.genome_lib_dir} \
            --left_fq {input.fq1} \
            --right_fq {input.fq2} \
            --output_dir {params.outdir} \
            --CPU {threads} \
            > {log} 2>&1
        """

rule apply_transcriptome_arriba_star_fusion:
    """Create fusion transcripts from Arriba and STAR-Fusion predictions."""
    input:
        arriba_fusion=rules.call_fusion_gene.output.fusions,
        star_fusion=rules.star_fusion.output.predictions,
        gtf=config["reference"]["gtf"],
        transcriptome_fasta=rules.build_transcriptome.output.transcriptome,
        kallisto_abundance=rules.kallisto_quant.output.abundance_tsv
    output:
        fusions_transcripts="results/final/{id}/fusion_transcripts.fa",
        summary="results/final/{id}/fusion_transcripts_summary.tsv"
    params:
        min_tpm=0.0,
        allow_intronic_snap=True
    conda:
        "../envs/python.yml"
    script:
        "../scripts/IV_creation_fusion_transcript.py"


