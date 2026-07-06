CHR_1_TO_4 = [str(c) for c in range(1, 5)]
CHR_5_TO_9 = [str(c) for c in range(5, 10)]
CHR_10_TO_16 = [str(c) for c in range(10, 17)]
CHR_1_TO_16 = [str(c) for c in range(1, 17)]

rule build_exon_dataframe:
    """
    Extrait les exons du GTF, les écrit en Parquet,
    et génère les fichiers pkl pour le génome et les arbres d'intervalles.
    """
    input:
        gtf=rules.download_human_gtf.output.gtf,
        transcripts_fasta=rules.build_transcriptome.output.transcriptome
    output:
        parquet="data/references/{id}/exon_data.parquet",
        genome_pkl="data/references/{id}/genome_dict.pkl", # pas utile
        tree_pkl="data/references/{id}/trees.pkl"
    conda:
        "../envs/python.yml"
    log:
        "logs/build_exon_{id}.log"
    script:
        "../scripts/V_build_exon_dataframe.py"

rule compress_vcf_filtered:
    """Compress and index VCF file."""
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
    """Split transcripts by chromosome groups."""
    input:
        transcripts_fasta=rules.build_transcriptome.output.transcriptome,
        exon_parquet=rules.build_exon_dataframe.output.parquet
    output:
        expand("results/variants/{{id}}/split_transcripts/{group}_transcripts.fa", group=["1_to_4","5_to_9", "10_to_16", "rest"])
    conda:
        "../envs/python.yml"
    shell:
        """
        mkdir -p results/variants/{wildcards.id}/split_transcripts
        python scripts/V_split_transcripts_all_groups.py \
            {input.transcripts_fasta} {input.exon_parquet} results/variants/{wildcards.id}/split_transcripts
        """

rule split_vcf_1_to_4:
    """Split VCF file for chromosomes 1 to 4."""
    input:
        vcf=rules.compress_vcf_filtered.output.vcf_gz
    output:
        vcf_split="results/variants/{id}/split_vcf/1_to_4.vcf"
    params:
        chromosomes=CHR_1_TO_4
    conda:
        "../envs/bcftools.yml"
    script:
        "../scripts/V_split_vcf_by_group.py"

rule split_vcf_5_to_9:
    """Split VCF file for chromosomes 1 to 10."""
    input:
        vcf=rules.compress_vcf_filtered.output.vcf_gz
    output:
        vcf_split="results/variants/{id}/split_vcf/5_to_9.vcf"
    params:
        chromosomes=CHR_5_TO_9
    conda:
        "../envs/bcftools.yml"
    script:
        "../scripts/V_split_vcf_by_group.py"

rule split_vcf_10_to_16:
    """Split VCF file for chromosomes 10 to 16."""
    input:
        vcf=rules.compress_vcf_filtered.output.vcf_gz
    output:
        vcf_split="results/variants/{id}/split_vcf/10_to_16.vcf"
    params:
        chromosomes=CHR_10_TO_16
    conda:
        "../envs/bcftools.yml"
    script:
        "../scripts/V_split_vcf_by_group.py"

rule split_vcf_rest:
    """Split VCF file for rest of the chromosomes (non 1-16)."""
    input:
        vcf=rules.compress_vcf_filtered.output.vcf_gz
    output:
        vcf_split="results/variants/{id}/split_vcf/rest.vcf"
    params:
        chromosomes=CHR_1_TO_16,
        include_rest=True
    conda:
        "../envs/bcftools.yml"
    script:
        "../scripts/V_split_vcf_by_group.py"

rule apply_variants_1_to_4:
    """Apply variants from chr 1 to 4."""
    input:
        exon_parquet=rules.build_exon_dataframe.output.parquet,
        vcf=rules.split_vcf_1_to_4.output.vcf_split,
#        vcf="results/variants/SRR4787043/20QC_variant_test_5000.vcf",
        transcripts_fasta="results/variants/{id}/split_transcripts/1_to_4_transcripts.fa",
        genome_pkl = "data/references/{id}/genome_dict.pkl",  
        tree_pkl = "data/references/{id}/trees.pkl"           
    output:
        mutated_output_fasta="results/final/{id}/1_to_4_mutated_transcripts.fa"
    log:
        "logs/variants/{id}/apply_variants_1_to_4.log"
    conda:
        "../envs/apply_variant.yml"
    script:
        "../scripts/V_combinatoire_mut_transcript.py"

rule apply_variants_5_to_9:
    """Apply variants from chr 5 to 9."""
    input:
        exon_parquet=rules.build_exon_dataframe.output.parquet,
        vcf=rules.split_vcf_5_to_9.output.vcf_split,
        transcripts_fasta="results/variants/{id}/split_transcripts/5_to_9_transcripts.fa",
        genome_pkl = "data/references/{id}/genome_dict.pkl",  
        tree_pkl = "data/references/{id}/trees.pkl"           
    output:
        mutated_output_fasta="results/final/{id}/5_to_9_mutated_transcripts.fa"
    log:
        "logs/variants/{id}/apply_variants_5_to_9.log"
    conda:
        "../envs/apply_variant.yml"
    script:
        "../scripts/V_combinatoire_mut_transcript.py"

rule apply_variants_10_to_16:
    """Apply variants from chr 10 to 16."""
    input:
        exon_parquet=rules.build_exon_dataframe.output.parquet,
        vcf=rules.split_vcf_10_to_16.output.vcf_split,
        transcripts_fasta="results/variants/{id}/split_transcripts/10_to_16_transcripts.fa",
        genome_pkl = "data/references/{id}/genome_dict.pkl",  
        tree_pkl = "data/references/{id}/trees.pkl"     
    output:
        mutated_output_fasta="results/final/{id}/10_to_16_mutated_transcripts.fa"
    log:
        "logs/variants/{id}/apply_variants_10_to_16.log"
    conda:
        "../envs/apply_variant.yml"
    script:
        "../scripts/V_combinatoire_mut_transcript.py"

rule apply_variants_rest:
    """Apply variants to the rest of the chromosomes."""
    input:
        exon_parquet=rules.build_exon_dataframe.output.parquet,
        vcf=rules.split_vcf_rest.output.vcf_split,
        transcripts_fasta="results/variants/{id}/split_transcripts/rest_transcripts.fa",
        genome_pkl = "data/references/{id}/genome_dict.pkl",  
        tree_pkl = "data/references/{id}/trees.pkl"     
    log:
        "logs/variants/{id}/apply_variants_rest.log"
    output:
        mutated_output_fasta="results/final/{id}/rest_mutated_transcripts.fa"
    conda:
        "../envs/apply_variant.yml"
    script:
        "../scripts/V_combinatoire_mut_transcript.py"

rule fusion_transcripts:
    """Merge all mutated transcript FASTAs with full transcriptome."""
    input:
        mutated_output_fasta1=rules.apply_variants_1_to_4.output.mutated_output_fasta,
        mutated_output_fasta2=rules.apply_variants_5_to_9.output.mutated_output_fasta,
        mutated_output_fasta3=rules.apply_variants_10_to_16.output.mutated_output_fasta,
        mutated_output_fasta4=rules.apply_variants_rest.output.mutated_output_fasta,
        transcripts_fasta=rules.build_transcriptome.output.transcriptome
    output:
        combined_transcripts="results/final/{id}/combined_transcripts.fa",
        report="results/final/{id}/combined_transcripts_dedup_report.tsv"
    conda:
        "../envs/python.yml"
    script:
        "../scripts/V_combine_transcripts.py"

#rule filter_transcripts:
#    """Filter transcripts for proteomics."""
#    input:
#        combined_fasta=rules.fusion_transcripts.output.combined_transcripts
#    output:
#        filtered_fasta="results/final/{id}/filtered_transcripts.fa"
#    conda:
#        "../envs/python.yml"
#    script:
#        "../scripts/V_filter_transcripts.py"

# split_filtered_transcript.py script utilisé pour découper 
