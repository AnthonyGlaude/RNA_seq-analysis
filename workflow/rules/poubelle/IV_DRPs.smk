rule extract_discordant_bam:
    input:
        #bam = rules.star_alignreads.output.bam
        bam = "/home/glaudea/scratch/glaudea/AGlaude_2024/workflow/results/STAR/171992_SIMG0590_T_totalRNA_sarcoma_43378_S9_L002/Aligned.sortedByCoord.out.bam"
    output:
        disc_bam = "results/DRP/{id}/drps_data.disc.bam",
        disc_bai = "results/DRP/{id}/drps_data.disc.bam.bai"
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        samtools view -b -F 1294 {input.bam} > {output.disc_bam}
        samtools index {output.disc_bam}
        """

rule remapping:
    input:
        bam = rules.star_alignreads.output.bam,
        disc_bam = rules.extract_discordant_bam.output.disc_bam,
        genome_dict = rules.build_exon_dataframe.output.parquet,
        genome = rules.download_human_genome.output.genome
    output:
        same_gene_remap = "results/DRP/{id}/pipeline_output/remapped_same_gene.pkl", 
        diff_gene_remap = "results/DRP/{id}/pipeline_output/remapped_diff_gene.pkl"
    conda:
        "../envs/python.yml"
    script:
        "../scripts/drp_maps.py"
 
rule interpretation:
    input:
        remapped_drps = rules.remapping.output.diff_gene_remap
    output:
        compilation = "results/DRP/{id}/figures_WES/strong_evidences_compilation_T.pkl"
    conda:
        "../envs/python.yml"##
    script:
        "../scripts/DRP_interpretation.py"
