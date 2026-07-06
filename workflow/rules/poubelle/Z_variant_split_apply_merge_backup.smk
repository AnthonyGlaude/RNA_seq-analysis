CHR_10_TO_13 = [str(c) for c in range(10, 14)]
CHR_14_TO_18 = [str(c) for c in range(14, 19)]
CHR_1_TO_18 = [str(c) for c in range(1, 19)]

rule build_exon_dataframe:
    """Extracts exons from the GTF file and saves them in a Parquet file."""
    input:
        gtf=rules.download_human_gtf.output.gtf  
    output:
        parquet="data/references/exon_data.parquet"  
    conda:
        "../envs/python.yml"  
    script:
        "../scripts/build_exon_dataframe.py"

rule compress_vcf_filtered:
    input:
        vcf=rules.filter_variants.output.vcf_filtered
    output:
        vcf_gz="results/variants/{id}/20QC_variant.vcf.gz"
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        bgzip -c {input.vcf} > {output.vcf_gz}
        tabix -p vcf {output.vcf_gz}
        """

rule split_transcripts_all:
    input:
        transcripts_fasta=rules.build_filtered_transcriptome.output.transcriptome_final_custom,
        exon_parquet=rules.build_exon_dataframe.output.parquet
    output:
        expand("results/variants/{{id}}/split_transcripts/{group}_transcripts.fa", group=[
            "1_to_2", "3_to_4", "5_to_7", "7_to_10", "10_to_13", "14_to_18", "rest"
        ])
    conda:
        "../envs/python.yml"
    shell:
        """
        mkdir -p results/variants/{wildcards.id}/split_transcripts
        python scripts/split_transcripts_all_groups.py \
            {input.transcripts_fasta} {input.exon_parquet} results/variants/{wildcards.id}/split_transcripts
        """

rule split_vcf_1_to_2:
    input:
        vcf=rules.compress_vcf_filtered.output.vcf_gz
    output:
        vcf_split="results/variants/{id}/split_vcf/1_to_2.vcf"
    params:
        chromosomes=["1", "2"]
    conda:
        "../envs/bcftools.yml"  
    script:
        "../scripts/split_vcf_by_group.py"

rule split_vcf_3_to_4:
    input:
        vcf=rules.compress_vcf_filtered.output.vcf_gz
    output:
        vcf_split="results/variants/{id}/split_vcf/3_to_4.vcf"
    params:
        chromosomes=["3", "4"]
    conda:
        "../envs/bcftools.yml"  
    script:
        "../scripts/split_vcf_by_group.py"

rule split_vcf_5_to_7:
    input:
        vcf=rules.compress_vcf_filtered.output.vcf_gz
    output:
        vcf_split="results/variants/{id}/split_vcf/5_to_7.vcf"
    params:
        chromosomes=["5", "6", "7"]
    conda:
        "../envs/bcftools.yml"
    script:
        "../scripts/split_vcf_by_group.py"

rule split_vcf_7_to_10:
    input:
        vcf=rules.compress_vcf_filtered.output.vcf_gz
    output:
        vcf_split="results/variants/{id}/split_vcf/7_to_10.vcf"
    params:
        chromosomes=["7", "8", "9", "10"]
    conda:
        "../envs/bcftools.yml"
    script:
        "../scripts/split_vcf_by_group.py"

rule split_vcf_10_to_13:
    input:
        vcf=rules.compress_vcf_filtered.output.vcf_gz
    output:
        vcf_split="results/variants/{id}/split_vcf/10_to_13.vcf"
    params:
        chromosomes=CHR_10_TO_13
    conda:
        "../envs/bcftools.yml"
    script:
        "../scripts/split_vcf_by_group.py"

rule split_vcf_14_to_18:
    input:
        vcf=rules.compress_vcf_filtered.output.vcf_gz
    output:
        vcf_split="results/variants/{id}/split_vcf/14_to_18.vcf"
    params:
        chromosomes=CHR_14_TO_18
    conda:
        "../envs/bcftools.yml"
    script:
        "../scripts/split_vcf_by_group.py"
        
rule split_vcf_rest:
    input:
        vcf=rules.compress_vcf_filtered.output.vcf_gz
    output:
        vcf_split="results/variants/{id}/split_vcf/rest.vcf"
    params:
        chromosomes=CHR_1_TO_18,
        include_rest=True
    conda:
        "../envs/bcftools.yml"
    script:
        "../scripts/split_vcf_by_group.py"

rule apply_variants_1_to_2:
    input:
        exon_parquet=rules.build_exon_dataframe.output.parquet,
        vcf=rules.split_vcf_1_to_2.output.vcf_split,
        transcripts_fasta="results/variants/{id}/split_transcripts/1_to_2_transcripts.fa"
    output:
        mutated_output_fasta="results/final/{id}/1_to_4_mutated_transcripts.fa"
    conda:
        "../envs/python.yml"
    script:
        "../scripts/add_variants.py"

rule apply_variants_3_to_4:
    input:
        exon_parquet=rules.build_exon_dataframe.output.parquet,
        vcf=rules.split_vcf_3_to_4.output.vcf_split,
        transcripts_fasta="results/variants/{id}/split_transcripts/1_to_2_transcripts.fa"
    output:
        mutated_output_fasta="results/final/{id}/1_to_4_mutated_transcripts.fa"
    conda:
        "../envs/python.yml"
    script:
        "../scripts/add_variants.py"

rule apply_variants_5_to_7:
    input:
        exon_parquet=rules.build_exon_dataframe.output.parquet,
        vcf=rules.split_vcf_5_to_7.output.vcf_split,
        transcripts_fasta="results/variants/{id}/split_transcripts/5_to_7_transcripts.fa"
    output:
        mutated_output_fasta="results/final/{id}/5_to_7_mutated_transcripts.fa"
    conda:
        "../envs/python.yml"
    script:
        "../scripts/add_variants.py"

rule apply_variants_7_to_10:
    input:
        exon_parquet=rules.build_exon_dataframe.output.parquet,
        vcf=rules.split_vcf_7_to_10.output.vcf_split,
        transcripts_fasta="results/variants/{id}/split_transcripts/7_to_10_transcripts.fa"
    output:
        mutated_output_fasta="results/final/{id}/7_to_10_mutated_transcripts.fa"
    conda:
        "../envs/python.yml"
    script:
        "../scripts/add_variants.py"

rule apply_variants_10_to_13:
    input:
        exon_parquet=rules.build_exon_dataframe.output.parquet,
        vcf=rules.split_vcf_10_to_13.output.vcf_split,
        transcripts_fasta="results/variants/{id}/split_transcripts/10_to_13_transcripts.fa"
    output:
        mutated_output_fasta="results/final/{id}/10_to_13_mutated_transcripts.fa"
    conda:
        "../envs/python.yml"
    script:
        "../scripts/add_variants.py"

rule apply_variants_14_to_18:
    input:
        exon_parquet=rules.build_exon_dataframe.output.parquet,
        vcf=rules.split_vcf_14_to_18.output.vcf_split,
        transcripts_fasta="results/variants/{id}/split_transcripts/14_to_18_transcripts.fa"
    output:
        mutated_output_fasta="results/final/{id}/14_to_18_mutated_transcripts.fa"
    conda:
        "../envs/python.yml"
    script:
        "../scripts/add_variants.py"

rule apply_variants_rest:
    input:
        exon_parquet=rules.build_exon_dataframe.output.parquet,
        vcf=rules.split_vcf_rest.output.vcf_split,
        transcripts_fasta="results/variants/{id}/split_transcripts/rest_transcripts.fa"
    output:
        mutated_output_fasta="results/final/{id}/rest_mutated_transcripts.fa"
    conda:
        "../envs/python.yml"
    script:
        "../scripts/add_variants.py"


rule fusion_transcripts:
    input:
        mutated_output_fasta1=rules.apply_variants_1_to_2.output.mutated_output_fasta,
        mutated_output_fasta2=rules.apply_variants_3_to_4.output.mutated_output_fasta,
        mutated_output_fasta3=rules.apply_variants_5_to_7.output.mutated_output_fasta,
        mutated_output_fasta4=rules.apply_variants_7_to_10.output.mutated_output_fasta,
        mutated_output_fasta5=rules.apply_variants_10_to_13.output.mutated_output_fasta,
        mutated_output_fasta6=rules.apply_variants_14_to_18.output.mutated_output_fasta,
        mutated_output_fasta7=rules.apply_variants_rest.output.mutated_output_fasta,
        transcripts_fasta=rules.build_filtered_transcriptome.output.transcriptome_final_custom
    output:
        combined_transcripts="results/final/{id}/combined_transcripts.fa"
    conda:
        "../envs/python.yml"
    script:
        "../scripts/combine_transcripts.py"

rule filter_transcripts:
    "Fichier de sortie pour la portion proteomique"
    input:
        combined_fasta=rules.fusion_transcripts.output.combined_transcripts
    output:
        filtered_fasta="results/final/{id}/filtered_transcripts.fa"
    conda:
        "../envs/python.yml"
    script:
        "../scripts/filter_transcripts.py"
